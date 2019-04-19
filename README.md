# DANCE - DAtabase of Nitrogen CEnters

<!-- toc -->

- [Installation](#installation)
- [Output Directories](#output-directories)
- [Bins](#bins)
- [Usage](#usage)
  - [Example](#example)
    - [GENERATE mode](#generate-mode)
    - [PLOTHIST mode](#plothist-mode)
    - [SELECT mode](#select-mode)
    - [SELECT-ANALYZE mode](#select-analyze-mode)
- [A Note on Logging](#a-note-on-logging)

<!-- tocstop -->

This tool filters molecules from a database. These molecules are then used to
generate parameters for smirnoff99frosst.

The heart of the tool is `dance.py`. It should be used as follows:

1. Use GENERATE mode to generate an initial set of trivalent nitrogen molecules
   from the database. The database must be represented as directories consisting
   of mol2 files. GENERATE mode ultimately generates an output directory with
   several files (see [Output Directories](#output-directories))
2. Use PLOTHIST mode to visualize the Wiberg bond orders from the previous step.
   This requires either the `output-tri-n-data.csv` file or the
   `output-tri-n-bonds.csv` from the GENERATE step. Note that if you choose the
   `output-tri-n-bonds.csv` file, you will have to change some of the command
   line arguments, as the defaults are for `output-tri-n-data.csv`.
   Specifically, `hist-min`, `hist-max`, and `hist-step` should be adjusted,
   most likely to around 0.5, 1.5, and 0.1, respectively. This step ultimately
   outputs the following file:
   - `output-histogram.pdf`: a PDF file holding histograms of the bond orders in
     every output-tri-n-data.csv file you pass in, as well as a histogram of the
     bond orders in all files combined (put on one plot together)
3. Use SELECT mode to make a final selection of molecules. This mode creates a
   directory with SMILES files corresponding to certain molecule
   [bins](#bins). Each file contains the n (specified on the
   command line) smallest molecules with that bin.
4. Use SELECT-ANALYZE mode to provide some statistics about the output of SELECT
   mode, as well as a bar graph of the counts of molecules in each bin.

## Installation

After cloning this repo, run the following command to install DANCE. You may
want to set up a [virtualenv](https://virtualenv.pypa.io/en/stable/) first. You
will also need an Openeye license to be able to use the Openeye toolkits.

```
pip install --extra-index-url https://pypi.anaconda.org/openeye/simple -e .
```

## Output Directories

The following files are generated whenever DANCE generates an "output directory"
of files.

- `mols.smi`: SMILES strings representing the molecules
- `mols.oeb`: an OEB (Openeye Binary) file for raw molecule data
- `tri-n-data.csv`: holds data about the trivalent nitrogen in each
  molecule - the total Wiberg bond order, total bond angle, and total bond
  length of the bonds surrounding the nitrogen
- `tri-n-bonds.csv`: holds data about the individual bonds connected
  to the trivalent nitrogen - the Wiberg bond order, bond length, and element
  of each bond
- `props.binary`: binary file for storing list of DanceProperties with
  data about the molecules

## Bins

During SELECT mode, molecules are separated into bins by two criteria. The first
is the total of the Wiberg bond orders of the bonds around the trivalent
nitrogen in the molecule. The total bond orders are rounded down to the nearest
multiple of the `--select-bin-size` command line arg. The second is a
"fingerprint", which attempts to describe the environment around the trivalent
nitrogen via several characteristics of each bond:

- atomic number of the atom at the end of the bond - what element the atom is
- connectivity of the atom at the end of the bond - how many other atoms the
  atom is connected to (including 1 for the trivalent nitrogen)
- rounded Wiberg bond order of the bond

## Usage

```
usage: dance [-h] [--mode MODE] [--log LEVEL] [--mol2dirs DIR1,DIR2,...]
             [--generate-output-dir DIRNAME] [--wiberg-csvs CSV1,CSV2,...]
             [--wiberg-csv-col INT] [--output-histograms FILENAME.pdf]
             [--hist-min FLOAT] [--hist-max FLOAT] [--hist-step FLOAT]
             [--input-binaries OEB,BINARY,OEB,BINARY,...]
             [--select-bin-size FLOAT] [--wiberg-precision FLOAT]
             [--bin-select INT] [--select-output-dir DIRNAME]
             [--select-analyze-dir DIR]
             [--select-analyze-output-dir FILENAME.pdf]

Performs various functions for selecting molecules from a database. It will do
the following based on the mode. |GENERATE| - Take in directories of mol2
files, generate the initial set of molecules with a single trivalent nitrogen,
and write results to a directory with the following files: mols.smi -
molecules stored in SMILES format, mols.oeb - the same molecules stored in OEB
(Openeye Binary) format, tri_n_data.csv - data about the trivalent nitrogen in
each molecule, tri_n_bonds.csv - data about the bonds around the trivalent
nitrogen in each molecule, props.binary - binary storage of DanceProperties
for the molecules. |PLOTHIST| - Take in data files from the previous step and
use matplotlib to generate histograms of the Wiberg bond orders. |SELECT| -
Make a final selection of molecules from the ones generated in the GENERATE
step. Writes the smallest molecules of a given "bin" to a SMILES file. (See
README for more info about bins.) |SELECT-ANALYZE| - Provide statistics and
visualizations of the output from SELECT mode. Writes to the following files:
statistics.txt - facts about the number of molecules in each bin,
visualization.pdf - a bar graph of numbers of molecules in each bin.

optional arguments:
  -h, --help            show this help message and exit

Mode Agnostic args:
  Arguments which apply to every mode of DANCE

  --mode MODE           The mode in which to run DANCE - one of GENERATE,
                        PLOTHIST, or SELECT. See README for more info
                        (default: GENERATE)
  --log LEVEL           logging level - one of DEBUG, INFO, WARNING, ERROR,
                        and CRITICAL - See
                        https://docs.python.org/3/howto/logging.html for more
                        information (default: info)

GENERATE args:
  --mol2dirs DIR1,DIR2,...
                        a comma-separated list of directories with mol2 files
                        to be filtered and saved (default: )
  --generate-output-dir DIRNAME
                        directory for saving the output - refer to beginning
                        of this msg (default: generate-output)

PLOTHIST args:
  --wiberg-csvs CSV1,CSV2,...
                        a comma-separated list of CSV files with a column
                        containing wiberg bond orders - these files are likely
                        generated in the GENERATE step (default: )
  --wiberg-csv-col INT  Column in the CSV files holding the Wiberg bond orders
                        (0-indexed) (default: 0)
  --output-histograms FILENAME.pdf
                        location of PDF file for histograms (default: output-
                        histograms.pdf)
  --hist-min FLOAT      Minimum bin for histogram (default: 2.0)
  --hist-max FLOAT      Maximum bin for histogram (default: 3.4)
  --hist-step FLOAT     Step/bin size for histogram (default: 0.1)

SELECT args:
  --input-binaries OEB,BINARY,OEB,BINARY,...
                        a comma-separated list of pairs of OEB and
                        DanceProperties binary files - each OEB should
                        correspond to the binary file next to it (default: )
  --select-bin-size FLOAT
                        bin size for separating molecules by total Wiberg bond
                        order around the trivalent nitrogen (default: 0.02)
  --wiberg-precision FLOAT
                        value to which to round the Wiberg bond orders in the
                        fingerprints; e.g. round to the nearest 0.02 (default:
                        0.05)
  --bin-select INT      specifies how many of the smallest molecules to select
                        from each bin, e.g. select the 5 smallest (default: 5)
  --select-output-dir DIRNAME
                        directory for writing SMILES files with molecules of
                        each fingerprint (default: select-output)

SELECT-ANALYZE args:
  --select-analyze-dir DIR
                        directory containing output from SELECT mode (default:
                        select-output)
  --select-analyze-output-dir FILENAME.pdf
                        directory for saving analysis (default: select-
                        analyze-output)
```

### Example

#### GENERATE mode

```
dance.py --mode GENERATE \
         --mol2dirs dir1,dir2,dir3 \
         --generate-output-dir my-output \
         --log debug
```

Reads in molecules from dir1, dir2, and dir3, filters out the ones with a single
trivalent nitrogen atom, and writes the results to files in a directory called
my-output. Prints log messages as low as DEBUG to stderr.

#### PLOTHIST mode

```
dance.py --mode PLOTHIST \
         --wiberg-csvs data1.csv,data2.csv,data3.csv \
         --output-histograms output-histograms.pdf \
         --log debug
```

Reads in Wiberg bond orders from data1.csv, data2.csv, and data3.csv and
generates histograms of the bond orders in each file. Also generates a histogram
for the bond orders from all files combined together. Writes the histograms to
output-histograms.pdf. Prints log messages as low as DEBUG to stderr.

#### SELECT mode

```
dance.py --mode SELECT \
         --input-binaries mol1.oeb,prop1.binary,mol2.oeb,prop2.binary \
         --bin-select 3 \
         --select-output-dir my-output
```

Reads in molecules and their properties from mol1.oeb, prop1.binary, mol2.oeb,
and prop2.binary, selects the 3 smallest molecules from each bin, and writes
molecules of each bin to files in a directory called my-output.

#### SELECT-ANALYZE mode

```
dance.py --mode SELECT-ANALYZE \
         --select-analyze-dir select-output \
         --select-analyze-output-dir select-analyze-output
```

Looks at SMILES files in the directory select-output and writes analysis to the
directory select-analyze-output.

## A Note on Logging

Python's standard logging library is used to write log messages of varying
severity to stderr. The severity level required for a message to be printed can
be adjusted with the `--log` flag. To capture the messages in a file, you will
have to redirect stderr to a file. For example, the following command will
redirect stderr to a file called status.txt when running dance.py.

```
dance.py --mode GENERATE --mol2dirs dir1,dir2,dir3 2> status.txt
```
