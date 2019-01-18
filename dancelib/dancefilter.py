"""Provides a class which is used to perform molecule filtering."""

from dancelib import danceutil
import glob
import logging
from openeye import oechem
from openeye import oequacpac


class DanceFilter:
    """
    DanceFilter is a class for performing all the various filtering
    actions. To use, instantiate it with the list of mol2dirs and output. Then
    call the generate() method to filter the molecules and write them to a file
    as SMILES strings. After generate() has been called, the molecules will each
    have the WIBERG_MOL data (where WIBERG_MOL is defined in danceutil) set
    to the Wiberg bond order. mol() can then be called to access all the
    molecules.

    Attributes:
        _mol2dirs: list of names of directories with mol2 files
        _output: the name of the output file to which to write the molecules
        _mols: a list storing the molecules the class is currently handling
    """

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

    def mols(self) -> [oechem.OEMol]:
        """Accessor for the molecules stored in this class."""
        return self._mols

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
            mol.SetData(danceutil.WIBERG_MOL, self._calc_wiberg(mol))
            logging.debug(f"wiberg-annotation: molecule {mol.GetTitle()}: "
                          f"{mol.GetData(danceutil.WIBERG_MOL)}")

        self._mols.sort(key=lambda m: m.GetData(danceutil.WIBERG_MOL))

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
            1 if danceutil.IS_TRI_N(atom) else 0 for atom in mol.GetAtoms())
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
                logging.debug(
                    f"failed to assign partial charges to {mol.GetTitle()}")
                return -1

            # Copy over bonds and charges from the charged copy to our copy
            for charged_atom, atom in zip(charged_copy.GetAtoms(),
                                          mol.GetAtoms()):
                atom.SetPartialCharge(charged_atom.GetPartialCharge())
            for charged_bond, bond in zip(charged_copy.GetBonds(),
                                          mol.GetBonds()):
                bond.SetData(danceutil.WIBERG_BOND,
                             charged_bond.GetData(danceutil.WIBERG_BOND))

            # Sum up bond orders around the trivalent nitrogen
            for atom in conf.GetAtoms(danceutil.IS_TRI_N):
                for bond in atom.GetBonds():
                    total += bond.GetData(danceutil.WIBERG_BOND)
                break  # terminate after one trivalent nitrogen
            break  # terminate after one conformation
        return total
