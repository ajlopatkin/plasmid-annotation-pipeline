#!/bin/bash

out_dir=$1
prefix=$2
gbk_dir=$3
match=$4

CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

conda init
conda activate prokka_env
prokka --kingdom Bacteria --outdir $out_dir --prefix $prefix --force --proteins $gbk_dir.gbk  ./../fastas/$prefix.fasta
conda deactivate

python prokka2go.py -i ./../output/plasmids/$prefix/gb_match/$match/pk_results/$prefix.gbk -o ./../output/plasmids/$prefix/gb_match/$match/prokka2go.txt -p $prefix
