RescoreBN
A Python tool for rescoring Bayesian Networks using user-provided datasets and DOT-formatted network structures.

This repository provides utilities for rescoring Bayesian Networks (BNs) based on input trajectories or a combined dataset. The main entry point is run.py, which demonstrates how to configure and execute the rescoring workflow using the BN_Rescore class from the rescore_BN module.

Features
Load Bayesian Network structures from DOT files

Rescore networks using either:

Multiple trajectory files, or

A single combined dataset

Support for discrete and continuous data

Flexible labeling options (union of labels, user-provided subsets, etc.)

Parallel execution across multiple CPU cores

Output results as a simple CSV table

Repository Structure
arduino
￼Copy code
your-repo/
├── run.py
├── rescore_BN.py
└── README.md
Installation
Python requirements
Typical dependencies include:

Python ≥ 3.8

numpy

pandas

networkx

pydot or graphviz

Install dependencies:

nginx
￼Copy code
pip install -r requirements.txt
Environment setup
Modify this line inside run.py to match the location of your module:

go
￼Copy code
sys.path.append('/path/to/rescore/')
Input Requirements
1. Bayesian Network DOT File
Provide the path to a .dot file:

ini
￼Copy code
dotfile = '/path/to/dotfile'
2. Data Input Options
You may provide data in one of two ways:

Option A — Multiple trajectory files
ini
￼Copy code
files = ['Rep1.csv', 'Rep2.csv', 'Rep3.csv']
union = False
filename = None
subsets = None
Use when your dataset is split across multiple replicates.
If union = True, labels will be merged across files.

Option B — Single combined file
ini
￼Copy code
filename = 'filename.csv'
files = None
subsets = None
Use when your data is already merged.
The file must contain a label column.
If not, provide subset indices manually through subsets.

3. Data Type
Specify whether your data is discrete or continuous:

ini
￼Copy code
discrete = True
4. Output File
ini
￼Copy code
outfile = 'out.csv'
5. Parallelization
Choose the number of CPU cores:

ini
￼Copy code
ncore = 3
Suggested: use one core per subset.

Running the Script
cpp
￼Copy code
Rescore = rBN.BN_Rescore(
    dotfile=dotfile,
    data=filename,
    files=files,
    union=union,
    subsets=subsets,
    discrete=discrete
)

Rescore.rescore(ncore=ncore)
Rescore.table_write(outfile)
Run from the command line:

arduino
￼Copy code
python run.py
Output
The script writes a CSV file (e.g., out.csv) containing rescoring results for each node in the Bayesian Network.

Example Minimal Configuration
ini
￼Copy code
ncore = 4
dotfile = 'net.dot'
filename = 'combined.csv'
files = None
union = False
subsets = None
discrete = True
outfile = 'rescore_output.csv'
Contact / Issues
If you encounter any issues or have feature requests, please open an Issue on GitHub.


