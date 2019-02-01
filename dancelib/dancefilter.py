"""Provides a class for performing molecule filtering."""

import glob
import logging
import math
from openeye import oechem
from openeye import oequacpac
from dancelib import danceprops
from dancelib import dancerunbase

#
# Constants
#

# Trivalent Nitrogen checker - instead of repeatedly constructing this class
IS_TRI_N = oechem.OEIsInvertibleNitrogen()

# Object for calculating AM1 charges - instead of reconstructing it everywhere
AM1 = oequacpac.OEAM1()


class DanceFilter(dancerunbase.DanceRunBase):
    """Performs various filtering actions on molecules

    After initializing DanceFilter, run it by calling the run() method.
    After run() has been called:
        - each molecule will have a danceprops.DANCE_PROPS_KEY tag (see
          danceprops.py for more info)
        - the molecules and their properties will be available via the
          get_data() method - the molecules returned will be sorted by Wiberg
          bond order and only have one trivalent nitrogen
        - a SMILES file with the molecules will have been written at the
          location specified upon initialization
    Note that the class is only meant to be run once, and attempting to call
    run() again will result in a RuntimeError.

    Attributes:
        _mol2dirs: list of names of directories with mol2 files
        _output: the name of the output file to which to write the molecules
        _mols: a list storing the molecules the class is currently handling
        _properties: a list storing properties of the molecules (see
                     danceprops.py for more info)
    """

    #
    # Public
    #

    def __init__(self, mol2dirs: [str], output: str):
        super().__init__()
        self._mol2dirs = mol2dirs
        self._output = output
        self._mols = []
        self._properties = []

    def run(self):
        """Pushes all the molecules through the various filtering steps"""
        super().check_run_fatal()
        logging.info("STARTING FILTERING")
        self._filter_tri_n()
        self._apply_properties()
        self._sort_by_wiberg()
        self._write_to_smiles_file()
        logging.info("FINISHED FILTERING")

    def get_data(self) -> ([oechem.OEMol], [danceprops.DanceProperties]):
        """Return the molecules and properties associated with this class."""
        return self._mols, self._properties

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
            mol, valid = self._check_one_molecule(mol2file)
            logging.debug(f"tri-n-check: molecule {mol2file}: {valid}")
            if valid:
                self._mols.append(mol)

    def _apply_properties(self):
        """Adds various properties to the molecules"""
        logging.info("Calculating properties for molecules")

        for mol in self._mols:
            danceprops.add_dance_property(mol, self._calc_properties(mol),
                                          self._properties)
            logging.debug(f"adding properties to molecule {mol.GetTitle()}")

    def _sort_by_wiberg(self):
        """
        Sorts the molecules by total Wiberg bond order around the trivalent
        nitrogen.
        """
        logging.info("Sorting molecules by Wiberg bond order")

        self._mols.sort(
            key=
            lambda m: danceprops.get_dance_property(m, self._properties).tri_n_bond_order
        )

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

    @staticmethod
    def _check_one_molecule(mol2file: str) -> (oechem.OEMol, bool):
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

        tri_n_count = sum(1 if IS_TRI_N(atom) else 0 for atom in mol.GetAtoms())
        return mol, tri_n_count == 1

    @staticmethod
    def _calc_properties(mol: oechem.OEMol) -> danceprops.DanceProperties:
        """
        Calculates properties of the given molecule and returns a
        DanceProperties object holding them.

        Based on Victoria Lim's am1wib.py - see
        https://github.com/vtlim/misc/blob/master/oechem/am1wib.py
        """
        props = danceprops.DanceProperties()

        for conf in mol.GetConfs():
            charged_copy = oechem.OEMol(conf)
            results = oequacpac.OEAM1Results()
            if not AM1.CalcAM1(results, charged_copy):
                logging.debug(
                    f"failed to assign partial charges to {mol.GetTitle()}")
                return props

            # Sum bond orders, bond lengths, and bond angles
            for atom in charged_copy.GetAtoms(IS_TRI_N):
                nbors = list(atom.GetAtoms())  # (neighbors)
                ang1 = math.degrees(
                    oechem.OEGetAngle(charged_copy, nbors[0], atom, nbors[1]))
                ang2 = math.degrees(
                    oechem.OEGetAngle(charged_copy, nbors[1], atom, nbors[2]))
                ang3 = math.degrees(
                    oechem.OEGetAngle(charged_copy, nbors[2], atom, nbors[0]))
                props.tri_n_bond_angle = ang1 + ang2 + ang3

                for nbor in nbors:
                    bond_order = results.GetBondOrder(atom.GetIdx(),
                                                      nbor.GetIdx())
                    bond_length = oechem.OEGetDistance(charged_copy, atom, nbor)
                    element = nbor.GetAtomicNum()

                    props.tri_n_bonds.append(
                        danceprops.DanceTriNBond(bond_order, bond_length,
                                                 element))
                    props.tri_n_bond_order += bond_order
                    props.tri_n_bond_length += bond_length

                break  # terminate after one trivalent nitrogen
            break  # terminate after one conformation

        return props
