#!/usr/bin/python
'''
DESCRIPTION: Sort through final annotation for given plasmid and find concensus genes + GO ID to create a concensus annotation
'''
from BCBio import GFF
from Bio import SeqIO
import os
import csv
import argparse
import time

def main(plasmid):
	start=time.time_ns()
	print("PROKKA PROCESS FOR " + plasmid)
        #delete old match csv and create a new matches.csv file that contains more limited list of matches
	matches=[]
	with open("./../output/plasmids/" + plasmid + "/matches.csv") as csvfile:
		input=csv.DictReader(csvfile, delimiter=' ', quotechar='|')
		for row in input:
			matches.append([row['Matches'], row['Genes']])

	#Create merged tsv by adding the final annotations of each accession ID to it
	with open("./../output/plasmids/" + plasmid + "/merge/mergedgo.tsv", 'w') as mergetsv:
		print('Create Merged TSV')
		mergewrite=csv.writer(mergetsv, delimiter='\t')
		with open("./../output/plasmids/" + plasmid + "/gb_match/default/pk_results/" + plasmid + "-go.tsv") as deftsv:
			defread=csv.reader(deftsv, delimiter='\t')
			index=0
			for line in defread:
				new_row=line
				for match in matches:
					opfile="./../output/plasmids/" + plasmid + "/gb_match/" + match[0] + "/pk_results/" + plasmid + "-go.tsv"
					with open(opfile) as aicsv:
						airead=csv.reader(aicsv, delimiter='\t')
						current=0
						for matchline in airead:
							if current==index:
								tmp=matchline
								for item in tmp[3:10]:
									new_row.append(item)
								break
							else:
								current=current+1
				mergewrite.writerow(new_row)
				print(new_row)
				index=index+1
	end=time.time_ns()
	print("TIME")
	print(end-start)
	#For each line (or position on the plasmid), sort through all the genes found for that position and determine which gene is the concensus gene for that location
	#Tie breakers based on accession ID gene count
	print("Concensus Genes from merged tsv")
	with open("./../output/plasmids/" + plasmid + "/final.tsv", 'w') as finaltsv:
		finalwrite=csv.writer(finaltsv, delimiter='\t')
		with open("./../output/plasmids/" + plasmid + "/merge/mergedgo.tsv") as reftsv:
			refread=csv.reader(reftsv, delimiter='\t')
			for row in refread:
				genes=[]
				#tmp=row
				if row[0] in 'locus_tag':
					finalwrite.writerow(row[0:10])
				else:
					print("Finding concensus gene for locus " + row[0])
					start=3
					for match in matches:
						if row[start] not in '':
							if len(genes)==0:
								genes.append([row[start:start+7], int(match[1]), 1])
								print(genes)
							else:
								equal=0
								for gene in genes:
									if row[start].lower() in gene[0][0].lower() or gene[0][0].lower() in row[start].lower():
										if int(match[1]) > int(gene[1]):

											if ("" in gene[0][4] and "" not in row[start+4]) or ("" in gene[0][4] and "" in row[start+4]):
												gene[0]=row[start:start+7]
											else:
												#Prioritize Reviewed Biological Process
												if 'P' in row[start+5] and 'P' not in gene[0][5]:
													gene[0][4]=row[start+4]
													gene[0][5]=row[start+5]
													gene[0][6]=row[start+6]
												elif 'R' in row[start+5] and 'R' not in gene[0][5]:
													gene[0][4]=row[start+4]
													gene[0][5]=row[start+5]
													gene[0][6]=row[start+6]
											gene[1]=int(match[1])

										#If there is no GO ID currently found for the given gene, append the first GO ID found to it
										if "" in gene[0][4] and "" not in row[start+4]:
											gene[0][4]=row[start+4]
											gene[0][5]=row[start+5]
											gene[0][6]=row[start+6]
										gene[2]=gene[2]+1
										equal=1
										break
								if equal==0:
									genes.append([row[start:start+7], int(match[1]), 1])
						start=start+7
						if start>=len(row):
							break
					mode=[]
					for gene in genes:
						if len(mode)==0:
							mode.append(gene)
						else:
							for most in mode:
								if int(most[2])<int(gene[2]):
									mode=[gene]
								elif int(most[2])==int(gene[2]) and int(gene[1])>int(most[1]):
									mode=[gene]
					if len(mode)==0:
						print("No gene at fiven locus!\n")
						finalwrite.writerow(row[0:10])
						print(row[0:10])
					else:
						extra=mode[0][0]
						print("Gene at given locus: " + extra[0])
						finalwrite.writerow(row[0:3] + extra[0:7])
						print(row[0:3] + extra[0:7])

	print('\n')

if __name__=='__main__':
        parser = argparse.ArgumentParser()
        parser.add_argument('plasmid', type=str)
        plasname=(parser.parse_args()).plasmid
        main(plasname)
