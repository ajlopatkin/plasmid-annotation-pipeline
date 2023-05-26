#!/bin/bash

out_dir=$1
plasmid=$2

CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

conda activate bakta_env
bakta --db ./../../bakta_database/db --output $out_dir ./../fastas/$plasmid.fasta
conda deactivate