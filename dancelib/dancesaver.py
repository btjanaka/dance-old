"""Provides a class for saving data about molecules."""

import logging
from openeye import oechem
from dancelib import danceprops
from dancelib import dancerunbase


class DanceSaver(dancerunbase.DanceRunBase):
    """Saves molecule data from DanceGenerator as SMILES strings and in CSV's.

    Molecules passed in should have the danceprops.DANCE_PROPS_KEY data set;
    failure to have this results in undefined behavior.

    After initializing DanceSaver with a list of molecules, run it by calling
    the run() method. Note that the class is only meant to be run once, and
    attempting to call run() again will result in a RuntimeError.

    Attributes:
        _mols: the molecules which are being saved
        _properties: a list storing properties of the molecules (see
                     danceprops.py for more info)
        _output_mols: SMILES file for SMILES strings of molecules
        _output_tri_n_data: csv file for data about trivalent nitrogens
        _output_tri_n_bonds: csv file for data about trivalent nitrogen bonds
    """

    #
    # Public
    #

    def __init__(self, mols: [oechem.OEMol],
                 properties: [danceprops.DanceProperties], output_mols: str,
                 output_tri_n_data: str, output_tri_n_bonds: str):
        super().__init__()
        self._mols = mols
        self._properties = properties
        self._output_mols = output_mols
        self._output_tri_n_data = output_tri_n_data
        self._output_tri_n_bonds = output_tri_n_bonds

    def run(self):
        """Perform all saving actions"""
        super().check_run_fatal()
        logging.info("STARTING SAVE")
        self._write_to_smiles_file()
        self._write_to_csv()
        logging.info("FINISHED SAVE")

    #
    # Private
    #

    def _write_to_smiles_file(self):
        """Writes the molecules to the SMILES file"""
        logging.info(f"Writing molecules to SMILES file {self._output_mols}")
        ostream = oechem.oemolostream(self._output_mols)
        ostream.SetFormat(oechem.OEFormat_USM)
        for mol in self._mols:
            oechem.OEWriteMolecule(ostream, mol)
        ostream.close()

    def _write_to_csv(self):
        """Writes the data about trivalent nitrogens and their bonds to CSVs"""
        logging.info(f"Writing tri-n data to {self._output_tri_n_data}")
        logging.info(f"Writing bond data to {self._output_tri_n_bonds}")

        with open(self._output_tri_n_data, 'w') as datacsv:
            with open(self._output_tri_n_bonds, 'w') as bondcsv:
                datacsv.write(
                    "tri_n_bond_order,tri_n_bond_angle,tri_n_bond_length\n")
                bondcsv.write("bond_order,bond_length,element\n")

                for mol in self._mols:
                    prop = danceprops.get_dance_property(mol, self._properties)
                    datacsv.write(f"{prop.tri_n_bond_order},"
                                  f"{prop.tri_n_bond_angle},"
                                  f"{prop.tri_n_bond_length}\n")
                    for bond in prop.tri_n_bonds:
                        bondcsv.write(f"{bond.bond_order},"
                                      f"{bond.bond_length},"
                                      f"{bond.element}\n")
