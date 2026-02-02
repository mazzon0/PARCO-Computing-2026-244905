#!/bin/bash

module load python-3.10.14

python -m venv scripts/.venv
scripts/.venv/bin/pip install -r scripts/requirements.txt

module unload python-3.10.14