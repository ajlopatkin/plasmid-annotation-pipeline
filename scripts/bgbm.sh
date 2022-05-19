#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh

plasmid=$1

python annotate.py $plasmid

#Use prokka to determine gene accuracy between accession ID's

# Get the Conda source path
#CONDA_BASE=$(conda info --base)

# Activate the conda environment
#source $CONDA_BASE/etc/profile.d/conda.sh

#conda init
#conda activate prokka_env
python prokka_nr.py $plasmid
#conda deactivate

python prokka_process.py $plasmid

#Use abricate to update annotation with antibiotic genes + create gff
#Create final.tsv and final.gff files
#conda activate abricate_env
python abricate_dw.py $plasmid
#conda deactivate
#cat ./plasmids/$plasmid/final.tsv

python abricate_test.py $plasmid

#Create final.gbk file
python gbk_update.py $plasmid

#Create plasmid map
python visual_map.py $plasmid
