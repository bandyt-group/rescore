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

sys.path.append('/path/to/rescore/')

## Input Requirements

### 1. Bayesian Network DOT File

Provide the path to a `.dot` file:

dotfile = '/path/to/dotfile'

---

### 2. Data Input Options

You may provide data in one of two ways:

#### Option A — Multiple trajectory files

files = ['Rep1.csv', 'Rep2.csv', 'Rep3.csv']
union = False
filename = None
subsets = None

Use when your dataset is split across multiple replicates.  
If `union = True`, labels will be merged across files.

#### Option B — Single combined file

filename = 'filename.csv'
files = None
subsets = None

Use when your data is already merged.  
The file must contain a `label` column.  
If not, provide subset indices manually through `subsets`.

---

### 3. Data Type

Specify whether your data is discrete or continuous:

discrete = True


