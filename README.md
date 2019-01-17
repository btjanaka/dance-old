# DANCE - DAtabase of Nitrogen CEnters

<!-- toc -->

- [Usage](#usage)
  * [Example](#example)
  * [As a Library](#as-a-library)
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

## Usage
Below is the help message for dance.py. Adding `python` before the
invocation of dance.py is optional.
```
usage: dance.py [-h] [--mol2dirs MOL2DIRS] [--output FILENAME.smi]
                [--log LEVEL]

Filters molecules stored in mol2 files from a list of directories and stores
them as SMILES strings

optional arguments:
  -h, --help            show this help message and exit
  --mol2dirs MOL2DIRS   a comma-separated list of directories with mol2 files
                        to be filtered (default: )
  --output FILENAME.smi
                        location of output file (default: output.smi)
  --log LEVEL           logging level - one of DEBUG, INFO, WARNING, ERROR,
                        and CRITICAL - See
                        https://docs.python.org/3/howto/logging.html for more
                        information (default: info)
```

### Example
```
dance.py --mol2dirs dir1,dir2,dir3 --output ouput.smi --log debug
```
Reads in molecules from dir1, dir2, and dir3, filters them, and writes the
results to output.smi. Log messages as low as DEBUG are printed to stderr.

### As a Library
It is possible to use dance.py as a library. Simply import the module,
construct a DanceMoleculeFilter object, call the `.generate()` method, and you
will have filtered the appropriate molecules and stored them in a file.
```Python
import dance

mol2dirs = ["dir1", "dir2"]
output = "output.smi"

dance_filter = dance.DanceMoleculeFilter(mol2dirs, output)
dance_filter.generate()
```

## A Note on Logging
Python's standard logging library is used to write log messages of varying
severity to stderr. The severity level required for a message to be printed can
be adjusted with the `--log` flag. To capture the messages in a file, you will
have to redirect stderr to a file. For example, the following command will
redirect stderr to a file called status.txt when running dance.py.
```
dance.py --mol2dirs dir1,dir2,dir3 2> status.txt
```
