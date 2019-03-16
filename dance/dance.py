#!/usr/bin/env python3
# Author: Bryon Tjanaka
# See README for more information.

import argparse
import glob
import logging
import math
import os
import pickle
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from openeye import oechem
from dance import dancegenerator
from dance import danceprops
from dance import dancesaver
from dance import dancewibhist


def parse_commandline_flags() -> {str: "argument value"}:
    """Uses argparse to handle parsing of command line flags."""
    parser = argparse.ArgumentParser(
        description=(
            "Performs various functions for selecting molecules from a "
            "database. It will do the following based on the mode. "
            "|GENERATE| - Take in directories of mol2 files, generate the "
            "initial set of molecules with a single trivalent nitrogen, and "
            "write results to a directory with the following files: "
            "mols.smi - molecules stored in SMILES format, "
            "mols.oeb - the same molecules stored in OEB (Openeye Binary) "
            "format, "
            "tri_n_data.csv - data about the trivalent nitrogen in each "
            "molecule, "
            "tri_n_bonds.csv - data about the bonds around the trivalent "
            "nitrogen in each molecule, "
            "props.binary - binary storage of DanceProperties for the "
            "molecules. "
            "|PLOTHIST| - Take in data files from the previous step and use "
            "matplotlib to generate histograms of the Wiberg bond orders. "
            "|SELECT| - Make a final selection of molecules from the ones "
            "generated in the GENERATE step. Writes the smallest molecules "
            "of a given \"bin\" to a SMILES file. (See README for "
            "more info about bins.) "
            "|SELECT-ANALYZE| - Provide statistics and visualizations of the "
            "output from SELECT mode. Writes to the following files: "
            "statistics.txt - facts about the number of molecules in each bin, "
            "visualization.pdf - a bar graph of numbers of molecules in each "
            "bin."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    mode_agnostic = parser.add_argument_group(
        "Mode Agnostic args", "Arguments which apply to every mode of DANCE")
    mode_agnostic.add_argument(
        "--mode",
        default="GENERATE",
        metavar="MODE",
        help=("The mode in which to run DANCE - one of GENERATE, PLOTHIST, "
              "or SELECT. See README for more info"))
    mode_agnostic.add_argument(
        "--log",
        default="info",
        metavar="LEVEL",
        help=("logging level - one of DEBUG, INFO, WARNING, ERROR, and CRITICAL"
              " - See https://docs.python.org/3/howto/logging.html for more "
              "information"))

    generate_group = parser.add_argument_group("GENERATE args")
    generate_group.add_argument(
        "--mol2dirs",
        default="",
        metavar="DIR1,DIR2,...",
        help=("a comma-separated list of directories with mol2 files to be "
              "filtered and saved"))
    generate_group.add_argument(
        "--generate-output-dir",
        default="generate-output",
        metavar="DIRNAME",
        help="directory for saving the output - refer to beginning of this msg")

    plothist_group = parser.add_argument_group("PLOTHIST args")
    plothist_group.add_argument(
        "--wiberg-csvs",
        default="",
        metavar="CSV1,CSV2,...",
        help=("a comma-separated list of CSV files with a column containing "
              "wiberg bond orders - these files are likely generated "
              "in the GENERATE step"))
    plothist_group.add_argument(
        "--wiberg-csv-col",
        default=0,
        metavar="INT",
        type=int,
        help=("Column in the CSV files holding the Wiberg bond orders "
              "(0-indexed)"))
    plothist_group.add_argument(
        "--output-histograms",
        default="output-histograms.pdf",
        metavar="FILENAME.pdf",
        help="location of PDF file for histograms")
    plothist_group.add_argument(
        "--hist-min",
        default=2.0,
        metavar="FLOAT",
        type=float,
        help="Minimum bin for histogram")
    plothist_group.add_argument(
        "--hist-max",
        default=3.4,
        metavar="FLOAT",
        type=float,
        help="Maximum bin for histogram")
    plothist_group.add_argument(
        "--hist-step",
        default=0.1,
        metavar="FLOAT",
        type=float,
        help="Step/bin size for histogram")

    select_group = parser.add_argument_group("SELECT args")
    select_group.add_argument(
        "--input-binaries",
        default="",
        metavar="OEB,BINARY,OEB,BINARY,...",
        help=("a comma-separated list of pairs of OEB and DanceProperties "
              "binary files - each OEB should correspond to the binary file "
              "next to it"))
    select_group.add_argument(
        "--select-bin-size",
        default=0.02,
        metavar="FLOAT",
        type=float,
        help=("bin size for separating molecules by total Wiberg bond order "
              "around the trivalent nitrogen"))
    select_group.add_argument(
        "--wiberg-precision",
        default=0.05,
        metavar="FLOAT",
        type=float,
        help=("value to which to round the Wiberg bond orders in the "
              "fingerprints; e.g. round to the nearest 0.02"))
    select_group.add_argument(
        "--bin-select",
        default=5,
        metavar="INT",
        type=int,
        help=("specifies how many of the smallest molecules to select from "
              "each bin, e.g. select the 5 smallest"))
    select_group.add_argument(
        "--select-output-dir",
        default="select-output",
        metavar="DIRNAME",
        help=("directory for writing SMILES files with molecules of each "
              "fingerprint"))

    select_analyze_group = parser.add_argument_group("SELECT-ANALYZE args")
    select_analyze_group.add_argument(
        "--select-analyze-dir",
        default="select-output",
        metavar="DIR",
        help="directory containing output from SELECT mode")
    select_analyze_group.add_argument(
        "--select-analyze-output-dir",
        default="select-analyze-output",
        metavar="FILENAME.pdf",
        help="directory for saving analysis")

    # Check for no arguments
    if len(sys.argv) == 1:
        parser.print_usage(sys.stderr)
        sys.exit(1)

    args = vars(parser.parse_args())

    # Make mode case-insensitive
    args["mode"] = args["mode"].upper()

    # Parse comma-separated lists
    for comma_sep_list in ["mol2dirs", "wiberg_csvs", "input_binaries"]:
        comma_sep_list = comma_sep_list.strip()
        args[comma_sep_list] = [] if args[comma_sep_list] == "" else args[
            comma_sep_list].split(",")

    # Check file extensions and append them if necessary
    for arg, extension in (("output_histograms", "pdf"),):
        if not args[arg].endswith("." + extension):
            args[arg] += "." + extension

    return args


def configure_logging(loglevel: str):
    """Configures Python's logging library"""
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f"Invalid log level: {loglevel}")
    logging.basicConfig(
        format="%(levelname)s: %(message)s", level=numeric_level)


def print_welcome(mode: str):
    """Print a fun little message"""
    print()
    print(r"    |                         ")
    print(r"  __|   __,   _  _    __   _  ")
    print(r" /  |  /  |  / |/ |  /    |/  ")
    print(r" \_/|_/\_/|_/  |  |_/\___/|__/")
    print()
    print(f"Welcome! You are using DANCE in {mode} mode")
    print()


#
# GENERATE mode
#


def run_generator(args):
    """Generates and saves molecules from the database."""
    dgenerator = dancegenerator.DanceGenerator(args["mol2dirs"])
    dgenerator.run()
    mols, properties = dgenerator.get_data()
    dancesaver.mkdir_and_save(mols, properties, args["generate_output_dir"])


#
# PLOTHIST mode
#


def run_plothist(args):
    """Plots Wiberg bond order histograms."""
    dwibhist = dancewibhist.DanceWibHist(
        args["wiberg_csvs"], args["wiberg_csv_col"], args["output_histograms"],
        args["hist_min"], args["hist_max"], args["hist_step"])
    dwibhist.run()


#
# SELECT mode
#


def write_mol_to_fingerprint_file(
        mol: oechem.OEMol,
        properties: [danceprops.DanceProperties],
        select_output_dir: str,
        select_bin_size: float,
        wiberg_precision: float,
):
    """Writes a molecule to its appropriate SMILES fingerprint file"""

    #  Some of the molecules coming in may be invalid. DanceGenerator may find
    #  there was an error in charge calculations, in which case the charged
    #  copy was not assigned to the molecule. This function checks for that.
    is_valid_molecule = \
            lambda mol: mol.HasData(danceprops.DANCE_CHARGED_COPY_KEY)

    if not is_valid_molecule(mol):
        logging.debug(f"Ignored molecule {mol.GetTitle()}")
        return

    charged_copy = mol.GetData(danceprops.DANCE_CHARGED_COPY_KEY)
    for atom in charged_copy.GetAtoms(oechem.OEIsInvertibleNitrogen()):
        tri_n = atom
        break
    fingerprint = danceprops.DanceFingerprint(tri_n, wiberg_precision)

    # Retrieve the total bond order around the trivalent nitrogen
    bond_order = danceprops.get_dance_property(mol, properties).tri_n_bond_order

    # Round the total bond order down to the lowest multiple of bin_size. For
    # instance, if bin_size is 0.02, and the bond_order is 2.028, it becomes
    # 2.02. This works because (bond_order / self._bin_size) generates a
    # multiple of the bin_size. Then floor() finds the next integer less than
    # the multiple. Finally, multiplying back by bin_size obtains the nearest
    # actual value.
    bond_order = math.floor(bond_order / select_bin_size) * select_bin_size

    filename = f"{select_output_dir}/{bond_order},{fingerprint}.smi"
    with open(filename, "a") as f:
        f.write(f"{oechem.OEMolToSmiles(mol)} {mol.GetTitle()}\n")
    logging.debug(f"Wrote {mol.GetTitle()} to {filename}")


def run_select(args):
    """Makes final selections of molecules."""
    if os.path.exists(args["select_output_dir"]):
        logging.error("Output directory already exists")
        sys.exit(1)
    os.mkdir(args["select_output_dir"])

    logging.info("STARTING SELECTING")
    logging.info("Labeling molecules with fingerprints")
    for oebf, binf in zip(args["input_binaries"][0::2],
                          args["input_binaries"][1::2]):
        with open(binf, "rb") as f:
            properties = pickle.load(f)
            istream = oechem.oemolistream(oebf)
            mol = oechem.OEMol()
            while oechem.OEReadMolecule(istream, mol):
                write_mol_to_fingerprint_file(
                    mol, properties, args["select_output_dir"],
                    args["select_bin_size"], args["wiberg_precision"])
    logging.info("FINISHED SELECTING")


#
# SELECT-ANALYZE mode
#


def write_data_to_pdf(bin_count: {str: int}, dirname: str):
    """
    Writes statistics and visualizations about the molecule bins to a text and
    and PDF file in the given directory.
    """
    logging.info(f"Writing results to {dirname}")

    os.makedirs(dirname, exist_ok=True)
    text = open(dirname + "/statistics.txt", "w")
    pdf = PdfPages(dirname + "/visualization.pdf")

    max_count = max(bin_count.values())
    sorted_bins = sorted(bin_count, key=lambda b: float(b.split(',')[0]))

    text.write(f"Max mols in a bin: {max_count}\n")
    text.write(f"These bin(s) have the max mols:\n")
    for bin_name in sorted_bins:
        if bin_count[bin_name] == max_count:
            text.write(f"  {bin_name}\n")
    text.write("\n")
    name_width = max(len(b) for b in bin_count)
    count_width = len(str(max(c for c in bin_count.values())))
    text.write("Bins and number of molecules in each bin:\n")
    for bin_name in sorted_bins:
        text.write(
            f"{bin_name:<{name_width}}: {bin_count[bin_name]:>{count_width}}\n")

    plt.rcParams.update({'font.size': 6})
    plt.rcParams.update({'figure.figsize': (15, int(len(bin_count) * .15))})
    fig, ax = plt.subplots()
    y_pos = np.arange(len(sorted_bins))
    counts = np.array([bin_count[b] for b in sorted_bins], dtype=int)
    ax.barh(y_pos, counts)
    ax.set_title("Number of molecules in each bin")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(sorted_bins)
    ax.invert_yaxis()
    ax.set_xlabel("Count")
    pdf.savefig(fig)

    text.close()
    pdf.close()


def run_select_analyze(args):
    """Runs some analysis on output from SELECT mode and generates a PDF."""
    logging.info("STARTING SELECT-ANALYZE")

    bin_count = {}
    for smi_file in glob.glob(args["select_analyze_dir"] + "/*.smi"):
        bin_name = smi_file[len(args["select_analyze_dir"]) +
                            1:-4]  # exclude directory name and ".smi"
        istream = oechem.oemolistream(smi_file)
        istream.SetFormat(oechem.OEFormat_USM)
        bin_count[bin_name] = sum(1 for mol in istream.GetOEMols())
    write_data_to_pdf(bin_count, args["select_analyze_output_dir"])

    logging.info("FINISHED SELECT-ANALYZE")


#
# Main
#


def main():
    """Runs everything."""
    args = parse_commandline_flags()
    configure_logging(args["log"])
    run_mode = {
        "GENERATE": run_generator,
        "PLOTHIST": run_plothist,
        "SELECT": run_select,
        "SELECT-ANALYZE": run_select_analyze,
    }
    if args["mode"] in run_mode:
        print_welcome(args["mode"])
        run_mode[args["mode"]](args)
    else:
        raise RuntimeError(f"Invalid mode: {args['mode']}")


if __name__ == "__main__":
    main()
