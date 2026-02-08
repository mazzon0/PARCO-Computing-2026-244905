#!/bin/bash

module load python-3.10.14

scripts/.venv/bin/python scripts/to_csv.py
scripts/.venv/bin/python scripts/plot.py

module unload python-3.10.14