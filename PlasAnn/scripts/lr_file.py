#!/usr/bin/perl
"""
LR_FILE_UPDATE.PY

DESCRIPTION: Creates csv file for processing in Linear Regression Script
"""
import argparse
import csv
import os
import sys
from glob import glob
import pandas as pd
import numpy as np
from Bio import SeqIO, GenBank, SeqFeature
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
#-----This function uses the psi.clstr results to 

def pres_abs(df):
	clstr = open("./../lr_files/clstr/seqs_25.clstr")
	clstr_file = clstr.readlines()

	genes_df = pd.DataFrame( columns = ["Cluster", "Plasmid", "Gene", "Locus"])
	plasmids = []
	row = []

	for line in clstr_file: 
		if not line.startswith(">Cluster"): 
			indexes = [i for i, ltr in enumerate(line) if ltr == "|"]
			index = line.index("...")
			carrot_ind = line.index(">")
			if len(indexes) > 0: 
				gene = (line[indexes[1]+1: index])
				locus = (line[indexes[0]+1:indexes[1]])
				plasmid = line[carrot_ind+1:indexes[0]]
				plasmids.append(str(plasmid))
				row.extend([plasmid, gene, locus])
				genes_df.loc[len(genes_df.index)] = row

				row = [row[0]]

		else: 
			row = [line[9:-1]]
    
	total_genes = int(row[0]) + 1
	plasmids = np.array(plasmids)
	plasmids = np.unique(plasmids)

	abs_pres = pd.DataFrame(index = plasmids, columns = np.arange(total_genes))

	for column in abs_pres:
		abs_pres[column].values[:] = 0

	cluster = 0
	new_cols = []

	gene_chosen = 0
	gene_name = {}

	for index, row in genes_df.iterrows():

		curr_clus = int(row["Cluster"])
		curr_plas = row["Plasmid"]

		if curr_clus != cluster: 
			pop_gene = max(gene_name, key=gene_name.get)
			if pop_gene == "nan":
				pop_gene = str(cluster) + "|" + prev_loc + "|" + prev_plas

			new_cols.append(pop_gene)
			cluster += 1
			abs_pres.at[curr_plas, cluster] = 1
			gene_chosen = 0
			gene_name = {}

		gene = str(row["Gene"])
		gene_name[gene] = gene_name.get(gene,0) + 1	

		if index == genes_df.index[-1]:
			pop_gene = max(gene_name, key=gene_name.get)
			if pop_gene == "nan":
				pop_gene = str(cluster) + "|" + str(row["Locus"]) + "|" + curr_plas
			# print(pop_gene)
			new_cols.append(pop_gene)
			cluster += 1
			abs_pres.at[curr_plas, curr_clus] = 1
			gene_chosen = 0
			gene_name = {}
			
		if curr_clus == cluster:
			abs_pres.at[curr_plas, curr_clus] = 1

		prev_plas = curr_plas
		prev_loc = str(row["Locus"]) 

	
	abs_pres.columns = new_cols
	abs_pres = pd.merge(df,abs_pres, left_index=True, right_index=True)
	abs_pres.to_csv("./../lr_files/gene_absence_presence.csv")

#---This function grabs all of our protein seqs from the final genbanks to one file to prep for CD-Hit
#----Grabs the lengths of all plasmids while looping through them to add to covariates of dataframes
def our_proteins(plasmids):

    all_sequence = []
    lengths = {}
    files = ""

    for plasmid in plasmids: 
		
        gbk = "./../output/plasmids/" + plasmid + "/the_final.gbk"
        for i in SeqIO.parse(gbk, "genbank"):
            lengths[plasmid] = len(i.seq)
            # plasmid = i.id
            for f in i.features: 
                if (f.type == "CDS" and ("gene" in f.qualifiers)): 
                    locus_tag = str(f.qualifiers["locus_tag"][0])
                    gene = f.qualifiers["gene"][0]
                    prt_id =  plasmid + "|" + locus_tag + "|" + str(gene) 
                    translation = str(f.qualifiers.get("translation", []))
                    translation = translation.strip("[']")
                    translation = translation.strip('\"')
                    curr_rec = SeqRecord(seq = Seq(translation), id = prt_id, description = str(plasmid))
                    all_sequence.append(curr_rec)

    SeqIO.write(all_sequence, "./../lr_files/all_sequences.faa", "fasta")

    if not os.path.exists('./../lr_files/clstr/'):
        os.makedirs('./../lr_files/clstr/')

    os.system("./../scripts/cd-hit/psi-cd-hit/psi-cd-hit.pl -i ./../lr_files/all_sequences.faa -o ./../lr_files/clstr/seqs_25 -c 0.25 -ce 1e-3")
    return lengths


def biologicalFunctions(list,df):

	ac = df
	go_df_p = pd.DataFrame(index = list)
	go_df_c = pd.DataFrame(index = list)
	go_df_f = pd.DataFrame(index = list)
	cog_df = pd.DataFrame(index = list)
	cog_func_df = pd.DataFrame(index = list)

	go_funcs = {}
	cog_functions = {}
	files = ""

	for plasmid in list: 
		file = "./../output/plasmids/" + plasmid + "/the_final.tsv"
		df = pd.read_csv(file, sep='\t', header=0)
		df1 = df[df['COG'].notna()]
		df_go = df.dropna(subset=['GO_ID_P', 'GO_ID_F', 'GO_ID_C'], how='all')

		#---------Grabbing instances of COG:IDs + functions
		for index, row in df1.iterrows():
			cog_id = str(row['COG'])
			cog_func = str(row['COG_FUNC'])
			cog_col = cog_id 
			cog_functions[cog_id] = cog_func

			if cog_id not in cog_df.columns:
				cog_df[cog_col] = 0
				cog_df.at[plasmid, cog_col] = 1
			else: 
				cog_df.at[plasmid, cog_id] += 1

			if cog_func not in cog_func_df.columns: 
				cog_func_df[cog_func] = 0
				cog_func_df.at[plasmid, cog_func] = 1
			else: 
				cog_func_df.at[plasmid, cog_func] += 1

		#---------Grabbing instances of GO:IDs + functions
		for index, row in df_go.iterrows():
			p_ids = ""
			c_ids = ""
			f_ids = ""
			go_ids = {}

			go_functions = str(row['GO_FUNC'])
			go_functions = go_functions.split('||')
			if '' in go_functions: 
				ind = go_functions.index('')
				del go_functions[ind]
			all_func = []
			for id_ in go_functions: 
				id_ = id_.split('|')
				all_func.extend(id_)

			#---------Grabbing each GO ID category if it exists
			all_gos = []
			if not pd.isna(row['GO_ID_P']): 
				go_p_id = str(row['GO_ID_P'])
				p_ids = go_p_id.split('|')
				p_ids.remove('')

				for i in p_ids: 
					if i not in go_df_p.columns: 
						go_df_p[i] = 0
						go_df_p.at[plasmid, i] = 1
					else: 
						go_df_p.at[plasmid, i] +=1
				all_gos.extend(p_ids)

			if not pd.isna(row['GO_ID_C']): 
				go_c_id = str(row['GO_ID_C'])
				c_ids = go_c_id.split('|')
				c_ids.remove('')

				for i in c_ids: 
					if i not in go_df_c.columns:
						go_df_c[i] = 0
						go_df_c.at[plasmid, i] = 1

					else: 
						go_df_c.at[plasmid, i] +=1

				all_gos.extend(c_ids)
				
			if not pd.isna(row['GO_ID_F']): 
				go_f_id = str(row['GO_ID_F'])
				f_ids = go_f_id.split('|')
				f_ids.remove('')

				for i in f_ids: 
					if i not in go_df_f.columns: 
						go_df_f[i] = 0
						go_df_f.at[plasmid, i] = 1
					else: 
						go_df_f.at[plasmid, i] +=1
				all_gos.extend(f_ids)

			#---------Creating dictionary for all GO functions
			if len(all_func) != len(all_gos): 
				line = ("Function missing at " + str(row['locus_tag']) + " on plasmid " + str(plasmid))
				files = files + line + "\n"
			else: 
				count = 0
				while count < len(all_gos): 
					go_funcs[all_gos[count]] = all_func[count]
					count += 1
					
	#---------File documenting where we are still missing IDs
	with open("./missings_go_spots.txt", "w") as f:
		f.write(files)

	#---------Merging each data frame with acquisition cost + covariate information
	go_df_f = pd.merge(ac,go_df_f, left_index=True, right_index=True)
	go_df_c = pd.merge(ac,go_df_c, left_index=True, right_index=True)
	go_df_p = pd.merge(ac,go_df_p, left_index=True, right_index=True)
	cog_df = pd.merge(ac,cog_df, left_index=True, right_index=True)
	cog_func_df = pd.merge(ac,cog_func_df, left_index=True, right_index=True)

	#---------Writing each DF to csv
	cog_df.to_csv('./../lr_files/lr_cog.csv')
	go_df_f.to_csv('./../lr_files/lr_go_f.csv')
	go_df_c.to_csv('./../lr_files/lr_go_c.csv')
	go_df_p.to_csv('./../lr_files/lr_go_p.csv')
	cog_func_df.to_csv('./../lr_files/lr_cog_functions.csv')
	
	if not os.path.exists('./../lr_files/functions/'):
		os.makedirs('./../lr_files/functions/')

	go_funcs = pd.DataFrame([go_funcs])
	go_funcs.to_csv('./../lr_files/functions/go_functions.csv')

	cog_functions = pd.DataFrame([cog_functions])
	cog_functions.to_csv('./../lr_files/functions/cog_functions.csv')

def main(): 
	#-------------Grabbing all plasmids (temporary method easily can use the fastas csv)
	folders = glob("./../output/plasmids/*/", recursive = True)
	print(folders)
	plasmids = [i.split("/")[-2] for i in folders]
	
	print(plasmids)
	if not os.path.exists('./../lr_files/'):
		os.makedirs('./../lr_files/')


	lengths = our_proteins(plasmids)

	#-------------Grabbing Acquisition cost and covariates like inc group
	df = pd.read_csv("./../Summary_Table_by_Plasmid_drug.csv")
	df = df[["Plasmid", "Antibiotic", "INC", "AC_avg"]].copy()

	# Adding new column to experimental info of plasmid length 
	df['length'] = df['Plasmid'].map(lengths)
	df = df.set_index("Plasmid")

	#-------------Creating gene_absence_presence.csv
	pres_abs(df)

	#-------------Running biologicalFunctions to make COG ID + GO:ID csvs
	biologicalFunctions(plasmids, df)

if __name__=="__main__":
	main()
