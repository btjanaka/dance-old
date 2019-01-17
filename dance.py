#!/usr/bin/env python3
# Author: Bryon Tjanaka
# See README for more information.

import argparse
import glob
import logging
import sys
from openeye import oechem
from openeye import oequacpac


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
        "--output",
        default="output.smi",
        metavar="FILENAME.smi",
        help="location of output file")
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


class DanceMoleculeFilter:
    """
    DanceMoleculeFilter is a class for performing all the various filtering
    actions. To use, instantiate it with the list of mol2dirs and output. Then
    call the generate() method to filter the molecules and write them to a file
    as SMILES strings.

    Attributes:
        _mol2dirs: list of names of directories with mol2 files
        _output: the name of the output file to which to write the molecules
        _mols: a list storing the molecules the class is currently handling
    """

    #
    # static constants
    #

    # Name of the label used for storing total Wiberg bond order in molecules
    WIBERG_MOL = "TotalWibergBondOrder"

    # Name of the label used for storing Wiberg bond order in inidividual bonds
    WIBERG_BOND = "WibergBondOrder"

    # Trivalent Nitrogen checker - instead of repeatedly constructing this class
    IS_TRI_N = oechem.OEIsInvertibleNitrogen()

    #
    # Public
    #

    def __init__(self, mol2dirs: [str], output: str):
        self._mol2dirs = mol2dirs
        self._output = output
        self._mols = []

    def generate(self):
        """Pushes all the molecules through the various filtering steps."""
        logging.info("STARTING FILTERING")
        self._filter_tri_n()
        self._apply_wiberg()
        self._write_to_smiles_file()
        logging.info("FINISHED FILTERING")

    #
    # Private
    #

    def _filter_tri_n(self):
        """
        Filters out the molecules with one trivalent nitrogen from the molecules
        stored as mol2 files.
        """
        logging.info(f"Filtering trivalent nitrogen molecules "
                     f"in directories {self._mol2dirs}")

        for mol2file in self._generate_mol2files():
            mol, ok = self._check_one_molecule(mol2file)
            logging.debug(f"tri-n-check: molecule {mol2file}: {ok}")
            if ok: self._mols.append(mol)

    def _apply_wiberg(self):
        """
        Adds Wiberg bond order annotations and sorts the molecules by them.
        """
        logging.info("Annotating molecules with Wiberg bond order")

        for mol in self._mols:
            mol.SetData(self.WIBERG_MOL, self._calc_wiberg(mol))
            logging.debug(f"wiberg-annotation: molecule {mol.GetTitle()}: "
                          f"{mol.GetData(self.WIBERG_MOL)}")

        self._mols.sort(key=lambda m: m.GetData(self.WIBERG_MOL))

    def _write_to_smiles_file(self):
        """Writes the molecules to the SMILES file"""
        logging.info(f"Writing molecules to SMILES file {self._output}")
        ostream = oechem.oemolostream(self._output)
        ostream.SetFormat(oechem.OEFormat_USM)
        for mol in self._mols:
            oechem.OEWriteMolecule(ostream, mol)
        ostream.close()

    #
    # Private (utility)
    #

    def _generate_mol2files(self) -> str:
        """Generator that yields mol2 files from the input directories.

        A generator is used to try to store as few items as possible. We know
        that the dataset can have millions of files.
        """
        for mol2dir in self._mol2dirs:
            for mol2file in glob.iglob(mol2dir + "/*.mol2"):
                yield mol2file

    def _check_one_molecule(self, mol2file: str) -> (oechem.OEMol, bool):
        """Checks if the molecule in the given file has only one trivalent
        nitrogen.

        Returns:
            The molecule itself, as well as a bool telling if there is only one
            trivalent nitrogen.
        """
        istream = oechem.oemolistream(mol2file)
        istream.SetFormat(oechem.OEFormat_MOL2)
        mol = oechem.OEMol()
        oechem.OEReadMolecule(istream, mol)

        tri_n_count = sum(
            1 if self.IS_TRI_N(atom) else 0 for atom in mol.GetAtoms())
        return mol, tri_n_count == 1

    def _calc_wiberg(self, mol: oechem.OEMol) -> float:
        """
        Calculates the sum of the Wiberg bond order among the bonds surrounding
        the trivalent nitrogen in the molecule (there should be only one at this
        point). Only considers the first conformation of the molecule.

        Based on Victoria Lim's am1wib.py - see
        https://github.com/vtlim/misc/blob/master/oechem/am1wib.py
        """
        total = 0.0
        for conf in mol.GetConfs():
            charged_copy = oechem.OEMol(mol)
            # TODO: Deprecation warning here, switch to OEAssignCharges
            status = oequacpac.OEAssignPartialCharges(
                charged_copy, oequacpac.OECharges_AM1BCCSym, False, False)
            if not status:
                raise RuntimeError(
                    f"OEAssignPartialCharges returned error code {status}")

            # Copy over bonds and charges from the charged copy to our copy
            for charged_atom, atom in zip(charged_copy.GetAtoms(),
                                          mol.GetAtoms()):
                atom.SetPartialCharge(charged_atom.GetPartialCharge())
            for charged_bond, bond in zip(charged_copy.GetBonds(),
                                          mol.GetBonds()):
                bond.SetData(self.WIBERG_BOND,
                             charged_bond.GetData(self.WIBERG_BOND))

            # Sum up bond orders around the trivalent nitrogen
            for atom in conf.GetAtoms(self.IS_TRI_N):
                for bond in atom.GetBonds():
                    total += bond.GetData(self.WIBERG_BOND)
                break  # terminate after one trivalent nitrogen
            break  # terminate after one conformation
        return total


if __name__ == "__main__":
    args = parse_commandline_flags()
    configure_logging(args["log"])
    dance_filter = DanceMoleculeFilter(args["mol2dirs"], args["output"])
    dance_filter.generate()
