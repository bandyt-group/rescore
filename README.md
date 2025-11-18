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


