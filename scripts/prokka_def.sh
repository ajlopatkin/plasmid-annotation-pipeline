#!/bin/bash

out_dir=$1
plasmid=$2

CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

conda init
conda activate prokka_env
prokka --kingdom Bacteria --outdir $out_dir --prefix $plasmid --force ./../fastas/$plasmid.fasta
conda deactivate

python prokka2go.py -i $out_dir/$plasmid.gbk -o ./../output/plasmids/$plasmid/gb_match/default/prokka2go.txt -p $plasmid
