#!/bin/bash
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
plasmid=$1

## Find acceptably annotated plasmids 
python annotate.py $plasmid

##Use prokka to determine gene accuracy between accession ID's
python prokka_nr.py $plasmid
conda activate prokka_env
python prokka_process.py $plasmid
conda deactivate

##Use abricate to update annotation with antibiotic genes + create gff
conda activate abricate_env
python abricate_dw.py $plasmid
python abricate_test.py $plasmid

##Create final.gbk file
python gbk_update.py $plasmid

##Create plasmid map
python visual_map.py $plasmid

##deactivate conda
conda deactivate