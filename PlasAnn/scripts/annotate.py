#!/usr/bin/python
"""
DESCRIPTION: Through a BLAST Search and input fasta file, program finds matching accession ID's
to given sequence and tests their reliability. Passing accession ID's are compiled into a CSV file

1) Run BLAST
2) Read in BLAST file and parse through each accession ID, downloading their gbk files and determining if
the given accession ID is fit for use
3) Compile good accession ID's into CSV file

"""
from Bio import SearchIO
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from Bio import Entrez
from io import StringIO
import csv
import argparse
import re
import os

def main(plasmid):

	#QBLAST DOES NOT WORK WITH LARGE XML SIZE: IGNORE COMMENTED OUT CODE HERE
	#TO BE FIXED IN FINAL VERSION
	# read in fasta file
	#record = SeqIO.read("./../fastas/" + plasmid + ".fasta", format="fasta")
	# blast it against the nucleotide database
	#print("Commencing Blast Search")
	#os.system("blastn -db nt -query ./../fastas/" + plasmid + ".fasta -out ./plasmids/" + plasmid + "/" + plasmid + "-BLASTTEST.xml -outfmt 5 -max_target_seqs 1000 -remote")
	#os.system("cat ./plasmids/" + plasmid + "/" + plasmid + "-BLAST.xml")
	#result_handle = NCBIWWW.qblast(program="blastn", database="nt", sequence=record.seq, auto_format="XML", hitlist_size=500, alignment_view=5)
	#print("Finished Blast Search")
	#result_handle=open("./plasmids/" + plasmid + "/" + plasmid + "-BLASTTEST.xml")
	#os.system("mkdir ./plasmids/" + plasmid)
	# write blast results to xml file
	#blast_file = "./plasmids/" + plasmid + "/" + plasmid + "-BLAST.xml"
	#with open(blast_file, "w") as out_handle:
	#	out_handle.write(result_handle.read())

	# read in blast file and parse the records
	blast_file = "./../output/plasmids/" + plasmid + "/" + plasmid + "-BLAST.xml"
	with open(blast_file,'r') as f:
		blast_records = list(NCBIXML.parse(f))[0]
		# collect list of accession ids to download
		accession_list = ["default"]
		genenum=[-1]
		print("\n" + str(len(blast_records.alignments)) + " accession ID's found\n")
		os.system("mkdir ./../output/plasmids/" + plasmid + "/gb_match")

		# determine which accession id meets criteria, and store in list
		align_pos=0;
		for alignment in blast_records.alignments:
			print(str(align_pos) + ". Processing " + str(alignment.accession) + ":")
			pct_id =  alignment.hsps[0].identities/alignment.hsps[0].align_length*100
			e_val = alignment.hsps[0].expect
			description = alignment.hit_def

			#can make pct_id into 60
			#if e_val < 1e-3 and "plasmid" in description and "complete sequence" in description and pct_id > 60:
			if e_val < 1e-3 and "plasmid" in description and pct_id > 60:
				#os.system("mkdir ./../output/plasmids/" + plasmid + "/gb_match")
				Entrez.email = 'jc4919@barnard.edu'
				handle = Entrez.efetch(db="nucleotide", id=alignment.accession, rettype="gb", retmode="text")
				gb_file= "./../output/plasmids/" + plasmid + "/gb_match/" + str(alignment.accession) + ".gbk"

				#download gbk file for each accession ID and check that enough genes are annotated
				try:
					with open(gb_file, "w") as bankfile:
						bankfile.write(handle.read())
						genes = SeqIO.read(gb_file, "genbank")

						if len(genes.features) > 0:
							all_genes=0
							id_genes=0

							gbk_idgene=[]


							print("GENES for " + str(alignment.accession) + ": ")
							for item in genes.features:
								if item.type=="CDS":
									#IGNORE: LEFT IN IN CASE SPECIFICATION IS NEEDED
									#if "gene" in item.qualifiers:
									#	if len(gbk_idgene)>0:
									#		if item.qualifiers.get("gene") not in gbk_idgene:
									#			gbk_idgene.append(item.qualifiers.get("gene"))
									#			id_genes=id_genes+1
									#			all_genes=all_genes+1
									#	else:
									#		gbk_idgene.append(item.qualifiers.get("gene"))
									#		all_genes=all_genes+1
									#		id_genes=id_genes+1
									#
									#else:
									#	all_genes=all_genes+1

									all_genes=all_genes+1
									if "gene" in item.qualifiers:
										id_genes=id_genes+1
							print("ALL vs ID: " + str(all_genes) + " vs " + str(id_genes))
							# if all_genes > 20 and id_genes/all_genes >= 0.3:
							# 	print("VALID\n")
							# 	accession_list.append(alignment.accession)
							# 	genenum.append(id_genes)
							# else:
							# 	print("INVALID: Not enough genes identified\n")
							# 	os.remove(gb_file, dir_fd=None)
							if id_genes > 9:
								print("VALID\n")
								accession_list.append(alignment.accession)
								genenum.append(id_genes)
							else:
								print("INVALID: Not enough genes identified\n")
								os.remove(gb_file, dir_fd=None)
						else:
							print("INVALID: No genes identified\n")
							os.remove(gb_file, dir_fd=None)
				except ValueError:
					print("INVALID: Empty genbank\n")
					os.remove(gb_file, dir_fd=None)
			else:
				print("INVALID: Does not meet initial conditions\n")
			align_pos=align_pos+1


		#create csv file to put all of the accession ID's here
		with open("./../output/plasmids/" + plasmid + "/matches.csv", 'w', newline='') as csvfile:
			input=csv.writer(csvfile, delimiter=' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
			input.writerow(['Matches'] + ['Genes'])
			for i in range(len(accession_list)):
				input.writerow([accession_list[i]] + [str(genenum[i])])

if __name__=="__main__":
        parser = argparse.ArgumentParser()
        parser.add_argument('plasmid', type=str)
        plasname=(parser.parse_args()).plasmid
        main(plasname)
