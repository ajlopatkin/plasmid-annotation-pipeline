#!/bin/bash
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
plasmid=$1

## blast plasmids
python blast_plasmids.py $plasmid