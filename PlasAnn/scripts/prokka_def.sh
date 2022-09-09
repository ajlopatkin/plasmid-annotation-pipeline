#!/bin/bash

out_dir=$1
plasmid=$2

CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

conda activate prokka_env
prokka --kingdom Bacteria --outdir $out_dir --prefix $plasmid --force ./../fastas/$plasmid.fasta
conda deactivate
