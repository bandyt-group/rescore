# RescoreBN

A Python tool for rescoring Bayesian Networks using user-provided datasets and DOT-formatted network structures.

This repository provides utilities for rescoring Bayesian Networks (BNs) based on input trajectories or a combined dataset. The main entry point is `run.py`, which demonstrates how to configure and execute the rescoring workflow using the `BN_Rescore` class from the `rescore_BN` module.

---

## Features

- Load Bayesian Network structures from DOT files  
- Rescore networks using either:
  - Multiple trajectory files, or  
  - A single combined dataset  
- Support for discrete and continuous data  
- Flexible labeling options (union of labels, user-provided subsets, etc.)  
- Parallel execution across multiple CPU cores  
- Output results as a simple CSV table  

## Setup Guide

Cloning the repository is enough to utilize the code.

Make sure to have all neccessary dependencies and setup the path as below:

### Python requirements

Typical dependencies include:

- Python ≥ 3.8  
- numpy  
- pandas  
- networkx  
- scipy  

### Environment setup

Modify this line inside `run.py` to match the location of your module:

```python
sys.path.append('/path/to/rescore/')
```

## Input Requirements (Change in `run.py`)

### 1. Bayesian Network DOT File

Provide the path to a `.dot` file:

```python
dotfile = '/path/to/dotfile'
```

---

### 2. Data Input Options

You may provide data in one of two ways:

#### Option A — Multiple trajectory files

```python
files = ['Rep1.csv', 'Rep2.csv', 'Rep3.csv']
union = False
filename = None
subsets = None
```

Use when your dataset is split across multiple replicates.  
If `union = True`, labels will be merged across files.

#### Option B — Single combined file

```python
filename = 'filename.csv'
files = None
subsets = None
```

Use when your data is already merged.  
The file must contain a `label` column as the last column denoting the subsets to rescore.  
If not, provide subset indices manually through `subsets`.

---

### 3. Data Type

Specify whether your data is discrete or continuous:

```python
discrete = True
```

---

### 4. Output File

```python
outfile = 'out.csv'
```

---

### 5. Parallelization

Choose the number of CPU cores:

Suggested: use one core per subset.


```python
ncore = 3
```

---

## Run from command line

```bash
python run.py
```

## Example Output (subset of `out.csv`)

Below is a sample of what the output file may look like.  
Source and target columns give the edges in the Bayesian network.  
Remaining columns represent rescored values for each replicate.  
(Values are for show only)

| source | target | Rep1   | Rep2   | Rep3   |
|--------|--------|--------|--------|--------|
| A      | B      | 0.1832 | 0.9344 | 0.5210 |
| C      | D      | 0.7729 | 0.1043 | 0.6637 |
| E      | F      | 0.2285 | 0.5830 | 0.9831 |
| G      | H      | 0.9013 | 0.3711 | 0.2476 |
| I      | J      | 0.4520 | 0.7792 | 0.6184 |

