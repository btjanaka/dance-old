# DANCE - DAtabase of Nitrogen CEnters

<!-- toc -->

- [Usage](#usage)
  * [Example](#example)
- [A Note on Logging](#a-note-on-logging)

<!-- tocstop -->

This tool filters molecules passed in as directories consisting of mol2 files.
The molecules are selected in several steps:
1. Only molecules with a single trivalent/invertible nitrogen are selected
2. Those molecules are sorted by Wiberg bond order and separated into bins
3. In each bin, the smallest molecules are chosen
Thus, the final set of molecules is "small molecules with a single trivalent
nitrogen, with a wide range of Wiberg bond orders". These molecules are stored
as SMILES strings. Ultimately, this tool is used for filtering molecules from
various databases.

*Note: So far only step 1 is fully implemented. Step 2 is partially implemented,
but we are still working on determining the bins. Step 3 has yet to be
implemented.*

Along the way, the tool also generates several additional files with interesting
data about the molecules. Currently, it generates these files:
- `output-tri-n-data.csv`: holds data about the trivalent nitrogen in each
  molecule - the total Wiberg bond order, total bond angle, and total bond
  length of the bonds surrounding the nitrogen
- `output-tri-n-bonds.csv`: holds data about the individual bonds connected to
  the trivalent nitrogen - the Wiberg bond order, bond length, and element of
  each bond


## Usage
Below is the help message for dance.py. Adding `python` before the
invocation of dance.py is optional.
```
usage: dance.py [-h] [--mol2dirs DIR1,DIR2,...] [--output-mols FILENAME.smi]
                [--output-tri-n-data FILENAME.csv]
                [--output-tri-n-bonds FILENAME.csv] [--log LEVEL]

Filters molecules stored in mol2 files from a list of directories and stores
them as SMILES strings. Also generates several additional files with useful
data. Note that for filenames, if the appropriate extension is not given, it
will be added on.

optional arguments:
  -h, --help            show this help message and exit
  --mol2dirs DIR1,DIR2,...
                        a comma-separated list of directories with mol2 files
                        to be filtered (default: )
  --output-mols FILENAME.smi
                        location of SMILES file holding final filtered
                        molecules (default: output-mols.smi)
  --output-tri-n-data FILENAME.csv
                        location of CSV file holding data about trivalent
                        nitrogens (with molecules in the same order as the
                        SMILES file (default: output-tri-n-data.csv)
  --output-tri-n-bonds FILENAME.csv
                        location of CSV file holding data about individual
                        bonds around trivalent nitrogens (default: output-tri-
                        n-bonds.csv)
  --log LEVEL           logging level - one of DEBUG, INFO, WARNING, ERROR,
                        and CRITICAL - See
                        https://docs.python.org/3/howto/logging.html for more
                        information (default: info)
```

### Example
```
dance.py --mol2dirs dir1,dir2,dir3 \
         --output-mols output-mols.smi \
         --output-tri-n-data output-tri-n-data.csv \
         --output-tri-n-bonds output-tri-n-bonds.csv \
         --log debug
```
Reads in molecules from dir1, dir2, and dir3, filters them, and writes the
resulting molecules to output-mols.smi. Additional data are stored in
output-tri-n-data.csv and output-tri-n-bonds.csv. Prints log messages as low as
DEBUG to stderr.

## A Note on Logging
Python's standard logging library is used to write log messages of varying
severity to stderr. The severity level required for a message to be printed can
be adjusted with the `--log` flag. To capture the messages in a file, you will
have to redirect stderr to a file. For example, the following command will
redirect stderr to a file called status.txt when running dance.py.
```
dance.py --mol2dirs dir1,dir2,dir3 2> status.txt
```
