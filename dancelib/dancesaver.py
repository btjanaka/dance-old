"""Provides a class and functions for saving data about molecules."""

import logging
import os
import pickle
from openeye import oechem
from dancelib import danceprops
from dancelib import dancerunbase


def create_filenames_for_dir(dirname: str) -> {str: str}:
    """
    Generates a dict of parameters that can be passed in to DanceSaver as a
    kwarg. Each value in the dict will be a file that ends up in the given
    directory
    """
    return {
        "output_mols_smi": f"{dirname}/mols.smi",
        "output_mols_oeb": f"{dirname}/mols.oeb",
        "output_tri_n_data": f"{dirname}/tri-n-data.csv",
        "output_tri_n_bonds": f"{dirname}/tri-n-bonds.csv",
        "output_props_binary": f"{dirname}/props.binary",
    }


def mkdir_and_save(mols: [oechem.OEMol],
                   properties: [danceprops.DanceProperties], dirname: str):
    """
    Creates a directory and saves files related to the given mols and
    properites to it.
    """
    os.makedirs(dirname, exist_ok=True)
    dsaver = DanceSaver(mols, properties, **create_filenames_for_dir(dirname))
    dsaver.run()


class DanceSaver(dancerunbase.DanceRunBase):
    """Saves molecule data from DanceGenerator several files.

    See the attributes starting with _output for a list of files saved. Note
    that if any of these are passed as None, they will not be saved.

    Molecules passed in should have the danceprops.DANCE_PROPS_KEY data set;
    failure to have this results in undefined behavior.

    After initializing DanceSaver with a list of molecules, run it by calling
    the run() method. Note that the class is only meant to be run once, and
    attempting to call run() again will result in a RuntimeError.

    Attributes:
        _mols: the molecules which are being saved
        _properties: a list storing properties of the molecules (see
                     danceprops.py for more info)
        _output_mols_smi: SMILES file for SMILES strings of molecules
        _output_mols_oeb: OEB (Openeye Binary) file for raw molecule data
        _output_tri_n_data: CSV file for data about trivalent nitrogens
        _output_tri_n_bonds: CSV file for data about trivalent nitrogen bonds
        _output_props_binary: binary file for storing list of DanceProperties
                              with data about the molecules
    """

    #
    # Public
    #

    # yapf: disable (otherwise has weird formatting for params)
    def __init__(self, mols: [oechem.OEMol],
                 properties: [danceprops.DanceProperties],
                 output_mols_smi=None, output_mols_oeb=None,
                 output_tri_n_data=None, output_tri_n_bonds=None,
                 output_props_binary=None):
        super().__init__()
        self._mols = mols
        self._properties = properties
        self._output_mols_smi = output_mols_smi
        self._output_mols_oeb = output_mols_oeb
        self._output_tri_n_data = output_tri_n_data
        self._output_tri_n_bonds = output_tri_n_bonds
        self._output_props_binary = output_props_binary
    # yapf: enable

    def run(self):
        """Perform all saving actions"""
        super().check_run_fatal()
        logging.info("STARTING SAVE")
        if self._output_mols_smi is not None: self._write_mols_to_smi()
        if self._output_mols_oeb is not None: self._write_mols_to_oeb()
        if self._output_tri_n_data is not None: self._write_data_to_csv()
        if self._output_tri_n_bonds is not None: self._write_bonds_to_csv()
        if self._output_props_binary is not None: self._write_props_to_binary()
        logging.info("FINISHED SAVE")

    #
    # Private
    #

    def _write_mols_with_oemolostream(self, filename: str,
                                      oeformat: "oechem.OEFormat"):
        """Writes molecules to the given file using an oemolostream"""
        ostream = oechem.oemolostream(filename)
        ostream.SetFormat(oeformat)
        for mol in self._mols:
            oechem.OEWriteMolecule(ostream, mol)
        ostream.close()

    def _write_mols_to_smi(self):
        """Writes the molecules to the SMILES file"""
        logging.info(
            f"Writing molecules to SMILES file {self._output_mols_smi}")
        self._write_mols_with_oemolostream(self._output_mols_smi,
                                           oechem.OEFormat_USM)

    def _write_mols_to_oeb(self):
        """Writes the molecules to an OEB file"""
        logging.info(f"Writing molecules to OEB file {self._output_mols_oeb}")
        self._write_mols_with_oemolostream(self._output_mols_oeb,
                                           oechem.OEFormat_OEB)

    def _write_data_to_csv(self):
        """Writes the data about trivalent nitrogens to a CSV"""
        logging.info(f"Writing tri-n data to {self._output_tri_n_data}")

        with open(self._output_tri_n_data, 'w') as datacsv:
            datacsv.write(
                "tri_n_bond_order,tri_n_bond_angle,tri_n_bond_length\n")
            for mol in self._mols:
                prop = danceprops.get_dance_property(mol, self._properties)
                datacsv.write(f"{prop.tri_n_bond_order},"
                              f"{prop.tri_n_bond_angle},"
                              f"{prop.tri_n_bond_length}\n")

    def _write_bonds_to_csv(self):
        """Writes data about bonds around trivalent nitrogens to a CSV"""
        logging.info(f"Writing bond data to {self._output_tri_n_bonds}")
        with open(self._output_tri_n_bonds, 'w') as bondcsv:
            bondcsv.write("bond_order,bond_length,element\n")

            for mol in self._mols:
                prop = danceprops.get_dance_property(mol, self._properties)
                for bond in prop.tri_n_bonds:
                    bondcsv.write(f"{bond.bond_order},"
                                  f"{bond.bond_length},"
                                  f"{bond.element}\n")

    def _write_props_to_binary(self):
        """Writes the properties objects to a binary file"""
        with open(self._output_props_binary, "wb") as bin_file:
            pickle.dump(self._properties, bin_file)
