#!/bin/bash
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
plasmid=$1

#Find acceptably annotated plasmids 
python annotate.py $plasmid

#Use abricate to update annotation with antibiotic genes
conda activate abricate_env
python abricate_dw.py $plasmid
python abricate_test.py $plasmid
conda deactivate

# Use prokka to get initial gbk annotation
python prokka_default.py $plasmid

# Label the plasmid based on genes
python combining_hit_info.py $plasmid


python prokka2go_LB.py -i ./../output/plasmids/$plasmid/the_final.gbk -o ./../output/plasmids/$plasmid/prokka2go.txt -p $plasmid

python update_tsv.py $plasmid

python abr_genes.py $plasmid

python update_gff.py $plasmid


#Create plasmid map
python visual_map.py $plasmid
