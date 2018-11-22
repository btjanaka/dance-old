#!/usr/bin/env python3
# Script for filtering out trivalent nitrogen molecules from the database.
# See README for usage.

import glob
import sys
from openeye import oechem

# Constants
USAGE_MSG = "Usage: databasen.py OUTPUT_FILENAME [INPUT_DIRS...]"


def print_status(status: str):
    """Prints a short status message to stderr."""
    print(status, file=sys.stderr)


def parse_commandline_flags() -> (str, [str]):
    """Handles all parsing of command line flags.

    Returns:
        name of the output file to which to write the SMILES
        strings
    Raises:
        IndexError: not enough arguments were passed in
    """
    if len(sys.argv) == 1: raise IndexError(USAGE_MSG)
    return sys.argv[1], sys.argv[2:]


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
    print_status("Filtering molecules in directories " + " ".join(input_dirs))
    print_status("Writing to " + out_filename)
    generate_output_file(out_filename, input_dirs)
    print_status("Done")
