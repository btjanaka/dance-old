"""Common code for use by all other parts of DANCE"""

from openeye import oechem

#
# Constants
#

# Name of the label used for storing total Wiberg bond order in molecules
WIBERG_MOL = "TotalWibergBondOrder"

# Name of the label used for storing Wiberg bond order in inidividual bonds
WIBERG_BOND = "WibergBondOrder"

# Trivalent Nitrogen checker - instead of repeatedly constructing this class
IS_TRI_N = oechem.OEIsInvertibleNitrogen()
