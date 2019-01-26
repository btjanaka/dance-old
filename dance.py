#!/usr/bin/env python3
# Author: Bryon Tjanaka
# See README for more information.

import argparse
import logging
from dancelib import danceanalyzer
from dancelib import dancefilter


def parse_commandline_flags() -> {}:
    """Uses argparse to handle all parsing of command line flags.

    Returns:
        a dictionary holding the values of all the flags
    Raises:
        IndexError: not enough arguments were passed in
    """
    parser = argparse.ArgumentParser(
        description=(
            "Filters molecules stored in mol2 files from a list of directories "
            "and stores them as SMILES strings. Also generates several "
            "additional files with useful data. Note that for filenames, if "
            "the appropriate extension is not given, it will be added on."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        "--mol2dirs",
        default="",
        metavar="DIR1,DIR2,...",
        help=
        "a comma-separated list of directories with mol2 files to be filtered")
    parser.add_argument(
        "--output-mols",
        default="output-mols.smi",
        metavar="FILENAME.smi",
        help="location of SMILES file holding final filtered molecules")
    parser.add_argument(
        "--output-tri-n-data",
        default="output-tri-n-data.csv",
        metavar="FILENAME.csv",
        help=("location of CSV file holding data about trivalent nitrogens "
              "(with molecules in the same order as the SMILES file"))
    parser.add_argument(
        "--output-tri-n-bonds",
        default="output-tri-n-bonds.csv",
        metavar="FILENAME.csv",
        help=("location of CSV file holding data about individual bonds around "
              "trivalent nitrogens"))
    parser.add_argument(
        "--log",
        default="info",
        metavar="LEVEL",
        help=("logging level - one of DEBUG, INFO, WARNING, ERROR, and CRITICAL"
              " - See https://docs.python.org/3/howto/logging.html for more "
              "information"))

    args = vars(parser.parse_args())

    # Parse comma-separated lists
    for comma_sep_list in ["mol2dirs"]:
        comma_sep_list = comma_sep_list.strip()
        args[comma_sep_list] = [] if args[comma_sep_list] == "" else args[
            comma_sep_list].split(",")

    # Check file extensions
    for arg, extension in (
        ("output_mols", "smi"),
        ("output_tri_n_data", "csv"),
        ("output_tri_n_bonds", "csv"),
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


def main():
    """Run everything."""
    args = parse_commandline_flags()
    configure_logging(args["log"])
    dfilter = dancefilter.DanceFilter(args["mol2dirs"], args["output_mols"])
    dfilter.run()
    mols, properties = dfilter.get_data()
    danalyzer = danceanalyzer.DanceAnalyzer(
        mols, properties, args["output_tri_n_data"], args["output_tri_n_bonds"])
    danalyzer.run()


if __name__ == "__main__":
    main()
