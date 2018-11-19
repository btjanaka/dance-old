# DatabaseN.py

This tool filters out trivalent nitrogen molecules from the eMolecules database.
The molecules are passed in as mol2 files, and a file with SMILES strings
representing the molecules with trivalent nitrogen's is outputted.

## Usage
```
databasen.py output_filename [input_directories...]
```

### Parameters
* `output_filename`: The file to which to output the resulting SMILES strings.
* `input_directories`: A (space-separated) list of directories containing the
  mol2 files. All the mol2 files in each directory are processed. This list may
  be empty, but the output file will definitely be as well.

### Example
```
databasen.py output.smi ./database1 ./database2 ./database3
```
Writes SMILES representations of the molecules with trivalent nitrogens in
database1, database2, and database3 to output.smi.
