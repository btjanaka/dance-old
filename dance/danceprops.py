"""Provides classes and functions for handling molecule properties

Within a program, DanceProperties are intended to be used as follows. First, one
keeps a list where all the properties are stored. Then, when a molecule is
assigned a DanceProperties object, the molecule actually stores a data tag
called DANCE_PROPS_KEY telling where to find the properties in the list. This is
necessary because OEMol's are not allowed to store custom objects in their data.
The process of adding a molecule to and retrieving a molecule from this list can
be streamlined by using the add_dance_property and get_dance_property functions.

EXAMPLE
=======

Consider the following two lists of molecules and their properties.

Molecules (list of oechem.OEMol)

Index | Data
------|-------------------------
  0   | ..., DANCE_PROPS_KEY: 2
  1   | ..., DANCE_PROPS_KEY: 3
  2   | ..., DANCE_PROPS_KEY: 0
  3   | ..., DANCE_PROPS_KEY: 1

Properties (list of DanceProperties)

Index | Attributes
-----------------------------------------------------------------------
  0   | tri_n_bond_order: 3.3, tri_n_bond_angle: 360,
      | tri_n_bond_length: 3.245,
      | tri_n_bonds: [<DanceTriNBond>,<DanceTriNBond>,<DanceTriNBond>],
      | fingerprint: [<DanceFingerprint>]
  1   | ...
  2   | ...
  3   | ...

Here, each molecule stores a DANCE_PROPS_KEY that tells where in the properties
list we can find the matching DanceProperties object. For instance, if we wish
to obtain the properties for molecule 0, we can do:

    index = molecules[0].GetData(DANCE_PROPS_KEY)
    mol0_props = properties[index]

We can also use a shortcut with get_dance_property:

    mol0_props = get_dance_property(molecules[0], properties)

Then we can access attributes of the properties; for instance, the fingerprint
would be accessed with:

    mol0_props.fingerprint
"""

from openeye import oechem

#
# Constants
#

# Tag used for storing the key to DanceProperties data in a molecule
DANCE_PROPS_KEY = "DancePropertiesKey"

# Tag used for storing the charged copy (using AM1 charges) of a molecule's
# first conformation
DANCE_CHARGED_COPY_KEY = "DanceChargedCopy"

# Tag used for storing Wiberg Bond Orders in molecules - needed because we want
# to save molecules and their properties, but the default object holding the
# bond orders - OEAM1Results - is not save-able using either pickle or Openeye
# tools
DANCE_BOND_ORDER_KEY = "DanceWibergBondOrder"

#
# Classes
#


class DanceProperties:
    """
    DanceProperties stores properties for a given molecule. It assumes the
    molecule has only one trivalent nitrogen.

    All attributes are public.

    Attributes:
        tri_n_bond_order: total bond order around the trivalent nitrogen
                          (float)
        tri_n_bond_angle: total bond angle around the trivalent nitrogen
                          (float)
        tri_n_bond_length: total bond length around the trivalent nitrogen
                           (float)
        tri_n_bonds: data about the bonds around the trivalent nitrogen
                     ([DanceTriNBond])
        fingerprint: provides greater detail about the bonds around the
                     trivalent nitrogen (DanceFingerprint)
    """

    def __init__(self):
        self.tri_n_bond_order = 0.0
        self.tri_n_bond_angle = 0.0
        self.tri_n_bond_length = 0.0
        self.tri_n_bonds = []
        self.fingerprint = None


class DanceTriNBond:
    """
    DanceBond stores information about a bond which has a trivalent nitrogen on
    one end. All attributes are public.

    Attributes:
        bond_order: Wiberg bond order (float)
        bond_length: length of the bond (float)
        element: the atomic number of the atom at the other end of the bond
                 (int)
    """

    def __init__(self, bond_order: float, bond_length: float, element: int):
        self.bond_order = bond_order
        self.bond_length = bond_length
        self.element = element


class DanceFingerprint:
    """
    DanceFingerprint is used to more uniquely identify the environment around
    the trivalent nitrogen in a molecule.

    Attributes:
        data: a tuple of size 3 representing each atom connected to the
              trivalent nitrogen, each tuple has the form (atomic
              number, connectivity, aromaticity, Wiberg bond order (rounded to
              nearest self._precision))
        _precision: the value to which to round the Wiberg bond orders; e.g.
                    "round to the nearest 0.02"
    """

    def __init__(self, tri_n: oechem.OEAtomBase, precision):
        """
        Creates the data based on a single trivalent nitrogen atom and its
        connections
        """
        self.data = [None, None, None]
        self._precision = precision

        for i, bond in enumerate(tri_n.GetBonds()):
            atom = bond.GetEnd() if bond.GetBgn() == tri_n else bond.GetBgn()
            atomic_num = atom.GetAtomicNum()
            connectivity = atom.GetDegree()
            aromaticity = atom.IsAromatic()
            wbo = round(bond.GetData(DANCE_BOND_ORDER_KEY) /
                        self._precision) * self._precision

            self.data[i] = (atomic_num, connectivity, aromaticity, wbo)

        self.data.sort()
        self.data = tuple(self.data)

    def __hash__(self):
        return hash(self.data)

    def __eq__(self, rhs):
        return self.data == rhs.data

    def __lt__(self, rhs):
        return self.data < rhs.data

    def __str__(self):
        return ''.join(
            f"#{d[0]}x{d[1]}{'a' if d[2] else ''}({d[3]})" for d in self.data)


#
# Functions
#


def add_dance_property(mol: oechem.OEMol, prop: DanceProperties,
                       properties: [DanceProperties]):
    """
    Appends the given molecule's properties to the list of DanceProperties and
    assigns the molecule a DANCE_PROPS_KEY data tag telling the index of its
    properties in the array.
    """
    properties.append(prop)
    set_dance_property(mol, len(properties) - 1)


def get_dance_property(mol: oechem.OEMol,
                       properties: [DanceProperties]) -> DanceProperties:
    """
    Returns the DanceProperties associated with a given molecule from the array.
    """
    key = mol.GetData(DANCE_PROPS_KEY)
    return properties[key]


def set_dance_property(mol: oechem.OEMol, key: int):
    """
    Sets the DANCE_PROPS_KEY data of a molecule.
    """
    mol.SetData(DANCE_PROPS_KEY, key)


def clean_properties_list(mols: [oechem.OEMol], properties: [DanceProperties]):
    """
    Modifies the given molecules and properties in-place to remove any
    properties not being used anymore.
    """
    props_copy = properties.copy()
    properties.clear()
    for mol in mols:
        add_dance_property(mol, get_dance_property(mol, props_copy), properties)


def append_properties_list(mols: [oechem.OEMol], properties: [DanceProperties],
                           mols2: [oechem.OEMol],
                           properties2: [DanceProperties]):
    """
    Adds on the molecules and properties in the second set of molecules and
    properties to the first set. This is not necessarily a trivial task because
    the keys in the second list of molecules have to be modified to point to the
    correct properties.
    """
    for i in range(len(mols2)):
        set_dance_property(mols2[i], len(mols) + i)
    mols.extend(mols2)
    properties.extend(properties2)
