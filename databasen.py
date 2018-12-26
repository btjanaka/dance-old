#!/usr/bin/env python3
# Author: Bryon Tjanaka
# Purpose: Filters mol2 files in a list of directories to find molecules with
#          trivalent nitrogens within them. Stores these filtered molecules as
#          SMILES strings.
# Usage: python databasen.py [-h] [-o OUTPUT] DIRS [DIRS ...]

import argparse
import glob
import sys
from openeye import oechem


def print_status(status: str):
    """Prints a short status message to stderr."""
    print(status, file=sys.stderr)


def parse_commandline_flags() -> (str, [str]):
    """Uses argparse to handle all parsing of command line flags.

    Returns:
        name of the output file to which to write the SMILES
        strings
    Raises:
        IndexError: not enough arguments were passed in
    """
    parser = argparse.ArgumentParser(
        description=
        "Finds molecules with trivalent nitrogens in a directory and stores them as SMILES strings."
    )
    parser.add_argument(
        "-o",
        "--output",
        help=
        "the name of an output file for the SMILES strings - defaults to output.smi"
    )
    parser.set_defaults(output="output.smi")

    req_args = parser.add_argument_group("required arguments")
    req_args.add_argument(
        "DIRS", nargs='+', help="directories with mol2 files to be filtered")

    args = parser.parse_args()

    return args.output, args.DIRS


def generate_mol2files(input_dirs: [str]) -> [str]:
    """Generator that yields mol2 files from the given input directories.

    A generator has to be used to store as few items as possible. We know that
    the dataset can have millions of files.
    """
    for d in input_dirs:
        for f in glob.iglob(d + "/*.mol2"):
            yield f


def check_one_molecule(mol2file: str) -> oechem.OEGraphMol or None:
    """Checks if the molecule in the given file has a trivalent nitrogen.

    Returns:
        The molecule itself if any of the atoms are trivalent nitrogens,
        otherwise None.
    """
    in_stream = oechem.oemolistream(mol2file)
    in_stream.SetFormat(oechem.OEFormat_MOL2)
    mol = oechem.OEGraphMol()
    oechem.OEReadMolecule(in_stream, mol)
    is_trivalent_n = oechem.OEIsInvertibleNitrogen()
    return mol if any(is_trivalent_n(atom) for atom in mol.GetAtoms()) else None


def generate_output_file(out_filename: str, input_dirs: [str]):
    """Performs the main workflow for filtering out the molecules."""
    out_stream = oechem.oemolostream(out_filename)
    out_stream.SetFormat(oechem.OEFormat_USM)

    for mol2file in generate_mol2files(input_dirs):
        checked_mol = check_one_molecule(mol2file)
        if checked_mol is not None:
            oechem.OEWriteMolecule(out_stream, checked_mol)

    out_stream.close()


if __name__ == "__main__":
    out_filename, input_dirs = parse_commandline_flags()
    print_status("Filtering molecules in directories " + str(input_dirs))
    print_status("Writing to " + out_filename)
    generate_output_file(out_filename, input_dirs)
    print_status("Done")
