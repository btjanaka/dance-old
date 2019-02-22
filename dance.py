#!/usr/bin/env python3
# Author: Bryon Tjanaka
# See README for more information.

import argparse
import logging
import pickle
from openeye import oechem
from dancelib import dancegenerator
from dancelib import danceprops
from dancelib import dancesaver
from dancelib import danceselector
from dancelib import dancewibhist


def parse_commandline_flags() -> {str: "argument value"}:
    """Uses argparse to handle all parsing of command line flags."""
    parser = argparse.ArgumentParser(
        description=(
            "Performs various functions for selecting molecules from a "
            "database. It will do the following based on the mode. "
            "GENERATE - Take in directories of mol2 files, generate the "
            "initial set of molecules with a single trivalent nitrogen, and "
            "write results to a directory. "
            "PLOTHIST - Take in data files from the previous step and use "
            "matplotlib to generate histograms of the Wiberg bond orders. "
            "SELECT - Make a final selection of molecules from the ones "
            "generated in the GENERATE step and write results to a directory. "
            "Note that when a part of this script \"writes results to a "
            "directory\", that means it generates a directory with the "
            "following files: mols.smi - molecules from that step stored "
            "in SMILES format, mols.oeb - the same molecules stored in OEB "
            "(Openeye Binary) format, tri_n_data.csv - data about the "
            "trivalent nitrogen in each molecule, tri_n_bonds.csv - data about "
            "the bonds around the trivalent nitrogen in each molecule, "
            "props.binary - binary storage of DanceProperties for the "
            "molecules"),
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
        help="bin size for separating molecules by Wiberg bond order")
    select_group.add_argument(
        "--wiberg-precision",
        default=0.02,
        metavar="FLOAT",
        type=float,
        help=("value to which to round the Wiberg bond orders in the "
              "fingerprints; e.g. round to the nearest 0.02"))
    select_group.add_argument(
        "--select-output-dir",
        default="select-output",
        metavar="DIRNAME",
        help="directory for saving the output - refer to beginning of this msg")

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


def read_binaries(oebfiles: [str], binfiles: [str]
                 ) -> ([oechem.OEMol], [danceprops.DanceProperties]):
    """
    Reads in molecules and DanceProperties from OEB files and binaries of
    DanceProperties generated by pickle.
    """
    logging.info(f"Reading oeb files {oebfiles}")
    logging.info(f"Reading binary files {binfiles}")

    mols = []
    properties = []

    for oebf, binf in zip(oebfiles, binfiles):
        mols2 = []
        istream = oechem.oemolistream(oebf)
        mol = oechem.OEMol()
        while oechem.OEReadMolecule(istream, mol):
            mols2.append(oechem.OEMol(mol))

        properties2 = None
        with open(binf, "rb") as f:
            properties2 = pickle.load(f)

        danceprops.append_properties_list(mols, properties, mols2, properties2)

    return mols, properties


def run_generator(args):
    """Generates and saves molecules from the database."""
    dgenerator = dancegenerator.DanceGenerator(args["mol2dirs"])
    dgenerator.run()
    mols, properties = dgenerator.get_data()
    dancesaver.mkdir_and_save(mols, properties, args["generate_output_dir"])


def run_plothist(args):
    """Plots Wiberg bond order histograms."""
    dwibhist = dancewibhist.DanceWibHist(
        args["wiberg_csvs"], args["wiberg_csv_col"], args["output_histograms"],
        args["hist_min"], args["hist_max"], args["hist_step"])
    dwibhist.run()


def run_select(args):
    """Makes final selections of molecules."""
    mols, properties = read_binaries(args["input_binaries"][0::2],
                                     args["input_binaries"][1::2])
    dselector = danceselector.DanceSelector(
        mols, properties, args["select_bin_size"], args["wiberg_precision"])
    dselector.run()
    mols, properties = dselector.get_data()
    dancesaver.mkdir_and_save(mols, properties, args["select_output_dir"])


def main():
    """Runs everything."""
    args = parse_commandline_flags()
    configure_logging(args["log"])
    run_mode = {
        "GENERATE": run_generator,
        "PLOTHIST": run_plothist,
        "SELECT": run_select,
    }
    if args["mode"] in run_mode:
        run_mode[args["mode"]](args)
    else:
        raise RuntimeError(f"Invalid mode: {args['mode']}")


if __name__ == "__main__":
    main()
