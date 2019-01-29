"""Provides classes and functions for handling molecule properties

Within a program, DanceProperties are intended to be used as follows. First, one
keeps a list where all the properties are stored. Then, when a molecule is
assigned a DanceProperties object, the molecule actually stores a data tag
called DANCE_PROPS_KEY telling where to find the properties in the list. This is
necessary because OEMol's are not allowed to store custom objects in their data.
The process of adding a molecule to and retrieving a molecule from this list can
be streamlined by using the add_dance_property and get_dance_property functions.
"""

from openeye import oechem

#
# Constants
#

# Tag used for storing the key to DanceProperties data in a molecule
DANCE_PROPS_KEY = "DancePropertiesKey"

#
# Classes
#


class DanceProperties:
    """
    DanceProperties stores properties for a given molecule. It assumes the
    molecule has only one trivalent nitrogen. DanceProperties is intended to be
    stored in a molecule's data under the DANCE_PROPS tag.

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
    """

    def __init__(self):
        self.tri_n_bond_order = 0.0
        self.tri_n_bond_angle = 0.0
        self.tri_n_bond_length = 0.0
        self.tri_n_bonds = []


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
    mol.SetData(DANCE_PROPS_KEY, len(properties) - 1)


def get_dance_property(mol: oechem.OEMol,
                       properties: [DanceProperties]) -> DanceProperties:
    """
    Returns the DanceProperties associated with a give molecule from the array.
    """
    key = mol.GetData(DANCE_PROPS_KEY)
    return properties[key]
