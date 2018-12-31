# DANCE - DAtabase of Nitrogen CEnters

This tool filters out trivalent nitrogen molecules from the eMolecules database.
The molecules are passed in as mol2 files, and a file with SMILES strings
representing the trivalent nitrogen molecules is outputted.

## Usage
Below is the help message for dance.py. Adding `python` before the
invocation of dance.py is optional.
```
usage: dance.py [-h] [-o OUTPUT] DIRS [DIRS ...]

Finds molecules with trivalent nitrogens in a directory and stores them as
SMILES strings.

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        the name of an output file for the SMILES strings -
                        defaults to output.smi

required arguments:
  DIRS                  directories with mol2 files to be filtered
```

### Example
```
dance.py -o output.smi ./database1 ./database2 ./database3
```
Writes SMILES representations of the molecules with trivalent nitrogens in
database1, database2, and database3 to output.smi.

## Other Notes
Status messages printed during the program are printed to stderr to avoid
interfering with any real output, so if you wish to capture them, you will have
to capture stderr. For example, the following command will redirect the status
messages to a file called status.txt.
```
dance.py dev 2> status.txt
```
