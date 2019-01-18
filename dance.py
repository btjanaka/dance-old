#!/usr/bin/env python3
# Author: Bryon Tjanaka
# See README for more information.

from dancelib import dancefilter
import argparse
import logging


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
            "and stores them as SMILES strings"),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        "--mol2dirs",
        default="",
        help=
        "a comma-separated list of directories with mol2 files to be filtered")
    parser.add_argument(
        "--output-mols",
        default="output-mols.smi",
        metavar="FILENAME.smi",
        help="location of SMILES file holding final filtered molecules")
    parser.add_argument(
        "--log",
        default="info",
        metavar="LEVEL",
        help=("logging level - one of DEBUG, INFO, WARNING, ERROR, and CRITICAL"
              " - See https://docs.python.org/3/howto/logging.html for more "
              "information"))

    args = vars(parser.parse_args())
    for comma_sep_list in ["mol2dirs"]:
        comma_sep_list = comma_sep_list.strip()
        args[comma_sep_list] = [] if args[comma_sep_list] == "" else args[
            comma_sep_list].split(",")
    return args


def configure_logging(loglevel: str):
    """Configure Python's logging library"""
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f"Invalid log level: {loglevel}")
    logging.basicConfig(
        format="%(levelname)s: %(message)s", level=numeric_level)


if __name__ == "__main__":
    args = parse_commandline_flags()
    configure_logging(args["log"])
    dfilter = dancefilter.DanceFilter(args["mol2dirs"], args["output_mols"])
    dfilter.generate()
