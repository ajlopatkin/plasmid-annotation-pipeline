#!/bin/bash
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
plasmid=$1
comb_run=$2

#Find acceptably annotated plasmids 
conda activate plasann_env
python annotate.py $plasmid
conda deactivate

# Use abricate to update annotation with antibiotic genes
conda activate abricate_env
python abricate_dw.py $plasmid
python abricate_test.py $plasmid
conda deactivate

# Use prokka to get initial gbk annotation
conda activate prokka_env
python prokka_default.py $plasmid
conda deactivate

# Label the plasmid based on genes
conda activate plasann_env
python combining_hit_info.py -p $plasmid -c $comb_run

python oric_blast.py $plasmid
python orit_blast.py $plasmid
python TnFinder_loc.py $plasmid
python COG_identify.py $plasmid

python prokka2go_LB.py -i ./../output/plasmids/$plasmid/the_final.tsv -o ./../output/plasmids/$plasmid/prokka2go.txt -p $plasmid

python update_tsv.py $plasmid

python abr_genes.py $plasmid

python update_gff.py $plasmid

python final_edits.py $plasmid

#Create plasmid map
python visual_map.py $plasmid
