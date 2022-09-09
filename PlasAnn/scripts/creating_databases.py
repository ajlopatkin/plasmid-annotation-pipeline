# creating_databases
from glob import glob
import csv
import os
import re
from Bio import SeqIO
	
incF_backbone = ['RepB/FIB','ori-1/oriV','Ssb', 'stbB', 'stbA', 'SopA','SopB','sopC','TraI','oriT','Tray','TraD',
	'TraT','TraS','GeneX','TraL','TraC','TraE','TraV','TraK','TraB','TraA','TraX','TraF','TraH','TrbI','TraN','TraU',
	'TraW','TrbC','TrbB','Rep E','ori-2/oriS']

incA_C_backbone = ['RepA','Ssb', 'ParA', 'ParB', 'TraI','TraD',
	'TraA','TraL','TraC','TraE','TraV','TraK','TraB','TrhX','TrhF','TraF','TraG','TraN','TraU',
	'TraW','pNDM-1_Dok01_N0219']


incP_backbone = ["TrfA1", "TrfA2", "Ssb", "oriV", "KorB", "TraC",  "IncC", "IncC1", "IncC2"
				, "nic", "oriT", "TraH", "TraI", "TraJ", "TraK", "TraG", "TrbK", "TrbN", 
				"TrbC", "TrbD", "TrbE", "TrbF", "TrbL", "TrbH", "TrbJ", "TrbG", 
				"TrbI", "TrbB", "TrbF", "TrbL", "TrbM", "TrbN", "TrbP"]

incN_backbone = ["repA", "uvp1", "IS6100", "EcoRII", "EcoRII met", "mrr", "kikA", "korB", 
				"traL", "korA", "traM,A", "traB", "traC" , "eex", "traN,E", "traD", "traO", 
				"traF", "traG", "nuc", "fipA", "traJ", "trak", "oriT", "stdA", "stdB", "stdC"
				, "orfD", "ccgAII", "ccgAI", "ardA", "ccgC", "ccgD", "ccgEII", "ardB", 
				"ardR", "mucA", "mucB","mpr", "ardK"]

incI_backbone = ["tnpA", "orf1158", "orf228", "orf192", "tmrB", "aacC2", "blatem-1", "iroN",
				"insC", "blaCTX-M-14", "tnpA", "erm(B)", "erm(B)L", "orf1", "groEL"]

# incL_M_backbone = ["", "", "", "","", "", "", "","", "", "", "", 
# 					"", "", "", "","", "", "", "","", "", "", "",
# 					"", "", "", "","", "", "", "","", "", "", "",
# 					"", "", "", "","", "", "", "","", "", "", "",
# 					]

# outputfiles = {"incF_backbone" : "backbone_genes/incF_backbone.fasta","incA_C_backbone" : "backbone_genes/incA_C_backbone.fasta",
# 				"incP_backbone" : "backbone_genes/incP_backbone.fasta","incL_M_backbone" : "backbone_genes/incL_M_backbone.fasta"}


incAC = glob("%s/*.gb" % "backbone_genes/IncAC/")

open("backbone_genes/incA_C_backbone.fasta",'w').close()
for record in incAC: 
	for i in SeqIO.parse(record, "gb"): 
		for f in i.features:
			if (f.type == "CDS" and ("gene" in f.qualifiers)):
				gene = f.qualifiers.get("gene", [])
				gene = str(gene[0])

				for x in incA_C_backbone:
					y = x
					if gene.lower() == y.lower(): 
						the_sequence = str(f.qualifiers.get("translation", []))
						the_sequence = (the_sequence.strip("[']"))

						with open("backbone_genes/incA_C_backbone.fasta",'a') as new_file:
							new_file.write("> " + gene + "_backbone"+ "\n" + the_sequence + "\n")
						new_file.close()
						print(gene)
						incA_C_backbone.remove(x)
						break

			

open("backbone_genes/incF_backbone.fasta",'w').close()
incF = glob("%s/*.gb" % "backbone_genes/IncF/")
for record in incF:
	for i in SeqIO.parse(record, "gb"): 
		for f in i.features:
			if (f.type == "CDS" and ("gene" in f.qualifiers)):
				gene = f.qualifiers.get("gene", [])
				gene = str(gene[0])

				for x in incF_backbone:
					y = x
					if gene.lower() == y.lower(): 
						the_sequence = str(f.qualifiers.get("translation", []))
						the_sequence = (the_sequence.strip("[']"))

						with open("backbone_genes/incF_backbone.fasta",'a') as new_file:
							new_file.write("> " + gene + "_backbone"+ "\n" + the_sequence + "\n")
						new_file.close()
						incF_backbone.remove(x)
						break 

open("backbone_genes/incLM_backbone.fasta",'w').close()
incLM = glob("%s/*.gb" % "backbone_genes/IncLM/")
for record in incLM: 
	for i in SeqIO.parse(record, "gb"): 
		for f in i.features:
			if (f.type == "CDS" and ("gene" in f.qualifiers)):
				gene = f.qualifiers.get("gene", [])
				gene = str(gene[0])

				for x in incLM_backbone:
					y = x
					if gene.lower() == y.lower(): 
						the_sequence = str(f.qualifiers.get("translation", []))
						the_sequence = (the_sequence.strip("[']"))

						with open("backbone_genes/incLM_backbone.fasta",'a') as new_file:
							new_file.write("> " + gene + "_backbone"+ "\n" + the_sequence + "\n")
						new_file.close()
						incLM_backbone.remove(x)
						break 

open("backbone_genes/incP_backbone.fasta",'w').close()
incP = glob("%s/*.gb" % "backbone_genes/IncP/")
for record in incP: 
	for i in SeqIO.parse(record, "gb"): 
		for f in i.features:
			if (f.type == "CDS" and ("gene" in f.qualifiers)):
				gene = f.qualifiers.get("gene", [])
				gene = str(gene[0])

				for x in incP_backbone:
					y = x
					if gene.lower() == y.lower(): 
						the_sequence = str(f.qualifiers.get("translation", []))
						the_sequence = (the_sequence.strip("[']"))

						with open("backbone_genes/incP_backbone.fasta",'a') as new_file:
							new_file.write(">_" + gene + "_backbone"+ "\n" + the_sequence + "\n")
						new_file.close()
						incP_backbone.remove(x)
						break 


# open("backbone_genes/incF_backbone.fasta",'w').close()
# for i in SeqIO.parse("backbone_genes/NC_002134.gb", "gb"): 
# 	for f in i.features:
# 		if (f.type == "CDS" and ("gene" in f.qualifiers)):
# 			gene = f.qualifiers.get("gene", [])
# 			gene = str(gene[0])

# 			for x in incF_backbone:
# 				y = x
# 				if gene.lower() == y.lower(): 
# 					the_sequence = str(f.qualifiers.get("translation", []))
# 					the_sequence = (the_sequence.strip("[']"))

# 					with open("backbone_genes/incF_backbone.fasta",'a') as new_file:
# 						new_file.write("> " + gene + "_backbone"+ "\n" + the_sequence + "\n")
# 					new_file.close()
# 					incF_backbone.remove(x)
# 					break 

