"""Provides a class for analyzing and outputting data about molecules."""

import logging
from openeye import oechem
from dancelib import danceprops


class DanceAnalyzer:
    """Saves and analyzes molecule data from DanceFilter.

    Molecules passed in should have the danceprops.DANCE_PROPS_KEY data set;
    failure to have this results in undefined behavior.

    After initializing DanceAnalyzer with a list of molecules, run it by calling
    the run() method. Note that the class is only meant to be run once, and
    attempting to call run() again will result in a RuntimeError.

    Attributes:
        _mols: the molecules which are being analyzed
        _properties: a list storing properties of the molecules (see
                     danceprops.py for more info)
        _output_tri_n_data: csv file for data about trivalent nitrogens
        _output_tri_n_bonds: csv file for data about trivalent nitrogen bonds
        _run_yet: whether run() has been called yet
    """

    #
    # Public
    #

    def __init__(self, mols: [oechem.OEMol],
                 properties: [danceprops.DanceProperties],
                 output_tri_n_data: str, output_tri_n_bonds: str):
        self._mols = mols
        self._properties = properties
        self._output_tri_n_data = output_tri_n_data
        self._output_tri_n_bonds = output_tri_n_bonds
        self._run_yet = False

    def run(self):
        """Perform all saving and analysis actions on the molecules"""
        if self._run_yet:
            raise RuntimeError("This DanceAnalyzer has already been run")
        self._run_yet = True

        logging.info("STARTING ANALYSIS")
        self._write_to_csv()
        # TODO: future functions of the analyzer
        #  self._plot_order_vs_angle()
        #  self._plot_total_order_vs_length()
        #  self._plot_individual_order_vs_length()
        logging.info("FINISHED ANALYSIS")

    #
    # Private
    #

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
