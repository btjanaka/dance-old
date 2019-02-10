#!/usr/bin/env python3
# Author: Bryon Tjanaka
# See README for more information.

import argparse
import logging
from dancelib import dancegenerator
from dancelib import dancesaver
from dancelib import dancewibhist


def parse_commandline_flags() -> {str: "argument value"}:
    """Uses argparse to handle all parsing of command line flags."""
    parser = argparse.ArgumentParser(
        description=(
            "Performs various functions for selecting molecules from a "
            "database. It will do the following based on the mode. "
            "GENERATE - Take in directories of mol2 files, generate the "
            "initial set of molecules with a single trivalent nitrogen, and "
            "write the molecules and accompanying data to various files. "
            "PLOTHIST - Take in data files from the previous step and use "
            "matplotlib to generate histograms of the Wiberg bond orders. "
            "SELECT - Make a final selection of molecules from the ones "
            "generated in the first step. "
            "See README for more info."),
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
        "--output-mols",
        default="output-mols.smi",
        metavar="FILENAME.smi",
        help="location of SMILES file holding final generated molecules")
    generate_group.add_argument(
        "--output-tri-n-data",
        default="output-tri-n-data.csv",
        metavar="FILENAME.csv",
        help=("location of CSV file holding data about trivalent nitrogens "
              "(with molecules in the same order as the SMILES file"))
    generate_group.add_argument(
        "--output-tri-n-bonds",
        default="output-tri-n-bonds.csv",
        metavar="FILENAME.csv",
        help=("location of CSV file holding data about individual bonds around "
              "trivalent nitrogens"))

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

    args = vars(parser.parse_args())

    # Make mode case-insensitive
    args["mode"] = args["mode"].upper()

    # Parse comma-separated lists
    for comma_sep_list in ["mol2dirs", "wiberg_csvs"]:
        comma_sep_list = comma_sep_list.strip()
        args[comma_sep_list] = [] if args[comma_sep_list] == "" else args[
            comma_sep_list].split(",")

    # Check file extensions and append them if necessary
    for arg, extension in (
        ("output_mols", "smi"),
        ("output_tri_n_data", "csv"),
        ("output_tri_n_bonds", "csv"),
        ("output_histograms", "pdf"),
    ):
        if not args[arg].endswith("." + extension):
            args[arg] += "." + extension

    return args


def configure_logging(loglevel: str):
    """Configure Python's logging library"""
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f"Invalid log level: {loglevel}")
    logging.basicConfig(
        format="%(levelname)s: %(message)s", level=numeric_level)


def run_generator(args):
    """Generate molecules from the database."""
    dgenerator = dancegenerator.DanceGenerator(args["mol2dirs"])
    dgenerator.run()
    mols, properties = dgenerator.get_data()
    dsaver = dancesaver.DanceSaver(mols, properties, args["output_mols"],
                                   args["output_tri_n_data"],
                                   args["output_tri_n_bonds"])
    dsaver.run()


def run_plothist(args):
    """Plots Wiberg bond order histograms."""
    dwibhist = dancewibhist.DanceWibHist(
        args["wiberg_csvs"], args["wiberg_csv_col"], args["output_histograms"],
        args["hist_min"], args["hist_max"], args["hist_step"])
    dwibhist.run()


def run_select(args):
    """Makes final selections of molecules."""
    pass


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
