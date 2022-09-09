#!/usr/bin/python
'''
DESCRIPTION: After running abricate with 4 different databases for each accession ID, check that the databases annotated the same resistance genes, and
write the names of these genes to the final annotation for the plasmid. Also write the incompatibility group found by plasmidfinder to the final gff for the plasmid
'''
from Bio.Seq import Seq
from Bio import SearchIO
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import os
import csv
import argparse

#list of card coordinates based on gene names

def main(plasmid):
	#Identify start location for each locus tag in the final tsv to verify res_gene match locations

	res_genes=[]
	inc_group=[]


	#For each accession ID match, open their abricate files and sort their found resistance genes or inc groups into corresponding lists
	#Print statements commented out; uncomment them for debugging purposes!
	res_ncbi=[]
	res_card=[]
	res_finder=[]
	pla_finder=[]
	found=[res_card, res_finder, res_ncbi]

	directory="./../output/plasmids/" + plasmid + "/ab_results"
	os.system("mkdir " + directory)
	databases=['ncbi', 'card', 'resfinder', 'plasmidfinder']
	for base in databases:
		with open(directory + "/" + base + ".txt") as csvfile:
				res_read=csv.reader(csvfile, delimiter='	')
				for inc_found in res_read:
					if(inc_found[5] not in "GENE"):
						if base in "ncbi":
							res_ncbi.append(inc_found)
						elif base in "card":
							res_card.append(inc_found)
						elif base in "resfinder":
							res_finder.append(inc_found)
						elif base in "plasmidfinder":
							pla_finder.append(inc_found)


	#Append found inc groups to inc_finder, matching the same in_groups and keeping track of their total matches and accession ID gene count information
	inc_finder=[]
	for group in pla_finder:
		if len(inc_finder)==0:
			inc_finder.append([group[5], 1])
		else:
			for inc in inc_finder:
				if group[5] in inc[0] or inc[0] in group[5]:
					inc[1]=inc[1]+1
				else:
					inc_finder.append([group[5], 1])

	#Sort through inc_finder and append valid inc groups to inc_group
	final_inc=[]
	for inc in inc_finder:
		if len(final_inc)==0:
			final_inc.append(inc)
		else:
			for final in final_inc:
				if inc[1] > final[1]:
					final_inc=[inc]
				elif inc[1] == final[1]:
					if inc[0] not in final[0] or final[0] not in inc[0]:
						final_inc.append(inc)
	for inc in final_inc:
		inc_group.append(inc)

	with open("./../output/plasmids/" + plasmid + "/abricate_results.txt", "w") as new_f: 
		new_f.write("abricate_results")
	new_f.close()


	#sort through chosen inc groups and pick the one that was chosen most frequently
	inc=[]
	# print("INC GROUP")
	for group in inc_group:
		# print("INC")
		# print(inc)
		if len(inc)==0:
			inc.append(group)
		else:
			for cmp in inc:
				if int(cmp[1])<int(group[1]):
					inc=[group]
	inc_concensus=""
	if len(inc)>0:
		inc_concensus=str(inc[0][0])

		if not os.path.exists("./../output/plasmids/"+plasmid+"/abr_genes.txt"): 
			with open("./../output/plasmids/"+plasmid+"/abr_genes.txt", 'wt') as writefile: 
				writefile.write(inc_concensus)
			writefile.close()
		else: 
			with open("./../output/plasmids/"+plasmid+"/abr_genes.txt", 'wt') as writefile: 
				writefile.write(inc_concensus)
			writefile.close()

	

if __name__=="__main__":
        parser = argparse.ArgumentParser()
        parser.add_argument('plasmid', type=str)
        plasname=(parser.parse_args()).plasmid
        main(plasname)
