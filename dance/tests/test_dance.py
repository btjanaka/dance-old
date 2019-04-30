"""Tests for DANCE"""

import os
from dance import dance
from openeye import oechem


class TestSelectFinal:
    """Integration tests of SELECT-FINAL mode"""

    SMILES = ["N", "N#N", "O=C=O"]

    @staticmethod
    def invoke_select_final(n, directory, output_file):
        """Wrapper for invoking SELECT-FINAL"""
        args = {
            "select_final_n": n,
            "select_final_dir": directory,
            "select_final_output_file": output_file,
        }
        dance.run_select_final(args)

    @staticmethod
    def create_smiles_in_tmpdir(tmpdir):
        """
        Creates 3 files, each with all the molecules defined in the SMILES
        variable, but arranged in different orders.
        """
        for i in range(3):
            f = (tmpdir / f"{i}.smi").open("w")
            counter = i
            for _ in range(3):
                f.write(TestSelectFinal.SMILES[counter] + "\n")
                counter = (counter + 1) % 3
            f.close()

    def test_selects_just_one_smallest_molecule(self, tmpdir):
        output_file = str(tmpdir / "select-final.smi")
        self.create_smiles_in_tmpdir(tmpdir)
        self.invoke_select_final(1, str(tmpdir), output_file)

        ifs = oechem.oemolistream(output_file)
        mol = oechem.OEMol()
        while oechem.OEReadMolecule(ifs, mol):
            assert oechem.OEMolToSmiles(mol) == "N"
