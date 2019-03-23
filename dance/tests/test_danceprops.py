"""Tests for functions and classes in danceprops"""

import pytest
from collections import defaultdict
from dance import danceprops
from openeye import oechem

# No tests for DanceProperties class, as it only has data

# No tests for DanceTriNBond class, as it only has data


@pytest.fixture(scope="function")
def ammonia_fake() -> (oechem.OEMol, oechem.OEAtomBase):
    """
    Creates an ammonia molecule with fake Wiberg bond orders. Also returns
    the trivalent nitrogen in the molecule
    """
    ammonia = oechem.OEMol()
    oechem.OESmilesToMol(ammonia, "N")
    oechem.OEAddExplicitHydrogens(ammonia)

    fake_wbo = [1.01, 1.03, 1.05]
    for i, bond in enumerate(ammonia.GetBonds()):
        bond.SetData(danceprops.DANCE_BOND_ORDER_KEY, fake_wbo[i])
    for atom in ammonia.GetAtoms(oechem.OEIsInvertibleNitrogen()):
        tri_n = atom
    return ammonia, tri_n


@pytest.fixture(scope="function")
def two_ammonia_fake_print(ammonia_fake) -> (oechem.OEMol, oechem.OEMol):
    """
    Returns two fingerprints for ammonia molecules with fake Wiberg bond
    orders
    """
    fingerprint1 = danceprops.DanceFingerprint(ammonia_fake[1], 0.05)
    fingerprint2 = danceprops.DanceFingerprint(ammonia_fake[1], 0.05)
    return (fingerprint1, fingerprint2)


# yapf: disable (otherwise has weird formatting for params)
@pytest.fixture(scope="function") # important since tests may modify the lists
def small_mols_and_props_list() -> ([oechem.OEMol],
                                    [danceprops.DanceProperties]):
    """
    Returns a list of 4 molecules, as well as another list with their 4
    corresponding properties. Note that the properties are actually empty.

    Molecules will be as follows - the index of their properties in the props
    list is given in parens:
    [ N (2), N#N (3), O=C=O (0), C#N (1) ]
    """
    # yapf: enable
    smiles = ["N", "N#N", "O=C=O", "C#N"]
    mols = [oechem.OEMol() for _ in range(4)]
    for i in range(4):
        oechem.OESmilesToMol(mols[i], smiles[i])

    mols[0].SetData(danceprops.DANCE_PROPS_KEY, 2)
    mols[1].SetData(danceprops.DANCE_PROPS_KEY, 3)
    mols[2].SetData(danceprops.DANCE_PROPS_KEY, 0)
    mols[3].SetData(danceprops.DANCE_PROPS_KEY, 1)

    props = [danceprops.DanceProperties() for _ in range(4)]

    return (mols, props)


class TestDanceFingerprint:
    """Tests the DanceFingerprint class"""

    def test_init_from_ammonia(self, ammonia_fake):
        fingerprint = danceprops.DanceFingerprint(ammonia_fake[1], 0.05)
        assert fingerprint.data == ((1, 1, False, 1.00), (1, 1, False, 1.05),
                                    (1, 1, False, 1.05))

    def test_equal_for_same_mols(self, two_ammonia_fake_print):
        assert two_ammonia_fake_print[0] == two_ammonia_fake_print[1]

    def test_same_hash_for_same_mols(self, two_ammonia_fake_print):
        assert hash(two_ammonia_fake_print[0]) == hash(
            two_ammonia_fake_print[1])

    def test_same_dict_bin_for_same_mols(self, two_ammonia_fake_print):
        prints = defaultdict(int)
        prints[two_ammonia_fake_print[0]] += 1
        prints[two_ammonia_fake_print[1]] += 1
        assert prints[two_ammonia_fake_print[0]] == 2

    def test_str_for_ammonia_is_correct(self, two_ammonia_fake_print):
        assert str(
            two_ammonia_fake_print[0]) == "#1x1(1.00)#1x1(1.05)#1x1(1.05)"


class TestDancePropertiesLists:
    """Tests functions for handling OEMol's and lists of DanceProperties"""

    def test_set_dance_property_persists(self):
        mol = oechem.OEMol()
        danceprops.set_dance_property(mol, 5)
        assert mol.GetData(danceprops.DANCE_PROPS_KEY) == 5

    def test_add_dance_property_appends_to_empty_list(self):
        mol, prop = oechem.OEMol(), danceprops.DanceProperties()
        props = []
        danceprops.add_dance_property(mol, prop, props)

        assert len(props) == 1
        assert mol.GetData(danceprops.DANCE_PROPS_KEY) == 0
        assert props[-1] == prop

    def test_add_dance_property_appends_to_small_list(
            self, small_mols_and_props_list):
        mols, props = small_mols_and_props_list
        mol, prop = oechem.OEMol(), danceprops.DanceProperties()
        mols.append(mol)
        danceprops.add_dance_property(mol, prop, props)

        assert len(props) == 5
        assert mol.GetData(danceprops.DANCE_PROPS_KEY) == 4
        assert props[-1] == prop

    def assert_props_equal(self, mols: [oechem.OEMol], mol_indices: [int],
                           props: [danceprops.DanceProperties], correct_props: [
                               danceprops.DanceProperties
                           ], correct_prop_indices: [int]):
        """
        Assert that, for each molecule index given, the property returned with
        danceprops.get_dance_property is equal to the property at the given
        index in correct_props.
        """
        if len(mol_indices) != len(correct_prop_indices):
            raise RuntimeError("There must be an equal number of molecule and "
                               "correct property indices")
        for i in range(len(mol_indices)):
            assert danceprops.get_dance_property(
                mols[mol_indices[i]],
                props) == correct_props[correct_prop_indices[i]]

    def test_get_dance_property_retrieves_correctly(self,
                                                    small_mols_and_props_list):
        mols, props = small_mols_and_props_list
        self.assert_props_equal(mols, [0, 1, 2, 3], props, props, [2, 3, 0, 1])

    def test_clean_properties_list_reorders(self, small_mols_and_props_list):
        old_mols, old_props = small_mols_and_props_list
        mols, props = old_mols.copy(), old_props.copy()
        danceprops.clean_properties_list(mols, props)

        for i in range(4):
            assert mols[i].GetData(danceprops.DANCE_PROPS_KEY) == i

        self.assert_props_equal(mols, [0, 1, 2, 3], props, old_props,
                                [2, 3, 0, 1])

    def test_clean_properties_list_removes_unused(self,
                                                  small_mols_and_props_list):
        old_mols, old_props = small_mols_and_props_list
        mols, props = old_mols.copy(), old_props.copy()
        del mols[2]
        danceprops.clean_properties_list(mols, props)

        for i in range(3):
            assert mols[i].GetData(danceprops.DANCE_PROPS_KEY) == i

        self.assert_props_equal(mols, [0, 1, 2], props, old_props, [2, 3, 1])

    def test_append_properties_list_does_nothing_for_empty_append(
            self, small_mols_and_props_list):
        old_mols, old_props = small_mols_and_props_list
        mols, props = old_mols.copy(), old_props.copy()
        mols2, props2 = [], []
        danceprops.append_properties_list(mols, props, mols2, props2)

        assert mols == old_mols
        assert props == old_props

    def test_append_properties_list_appends_small_list(
            self, small_mols_and_props_list):
        old_mols, old_props = small_mols_and_props_list
        mols, props = old_mols.copy(), old_props.copy()

        mols2 = [oechem.OEMol() for i in range(2)]
        oechem.OESmilesToMol(mols2[0], "CN=C=O")
        oechem.OESmilesToMol(mols2[1], "[C-]#[O+]")
        danceprops.set_dance_property(mols2[0], 1)
        danceprops.set_dance_property(mols2[1], 0)
        old_mols2 = [oechem.OEMol(m) for m in mols2]
        props2 = [danceprops.DanceProperties() for i in range(2)]

        danceprops.append_properties_list(mols, props, mols2, props2)

        # The properties obtained from the molecules put in the new list should
        # be the same as the properties in the old list.
        assert danceprops.get_dance_property(
            mols[4], props) == danceprops.get_dance_property(
                old_mols2[0], props2)
        assert danceprops.get_dance_property(
            mols[5], props) == danceprops.get_dance_property(
                old_mols2[1], props2)
