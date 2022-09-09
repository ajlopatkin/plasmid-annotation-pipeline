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
	start_check=[]
	with open("./../output/plasmids/" + plasmid + "/the_final.gbk") as start_chk:
		for line in start_chk.readlines():
			if line.startswith(' ' * 5 + 'CDS'):
				coordinates=line[8:].strip().strip("complement()")
				strt_index=0
				end_index=coordinates.index('..')
				start_check.append(int(coordinates[strt_index:end_index]))

	res_genes=[]

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


		#If a resistance gene is found by 2/3 of the resistance gene databases, it is added to the res_gene list
	compare=4
	for list in found:
			#print("CURR LIST OF RESISTANCE GENES IN PLASMID")
			#print(list)
			for resgene in list:
			#	print("CURRENT GENE:" + resgene[5])
				match=1
				for cmp in found[compare%3]:
					cmpgene=cmp[5]
					if '-' in cmp[5]:
						line_index=cmp[5].index('-')
						# print(line_index)
						# print(cmp[5][:5])
						if line_index>=4 and line_index!=len(cmp[5]):
							cmpgene=cmp[5][:line_index]
			#		print("COMPARE TO " + cmp[5])
					if (resgene[5].lower() in cmpgene.lower() or cmpgene.lower() in resgene[5].lower()) and abs(int(resgene[2])-int(cmp[2]))-300:
						resloc=resgene[2]
						cmploc=cmp[2]
						if int(resloc) not in start_check and cmploc in start_check:
							resgene[2]=cmp[2]

						match=match+1
			#			print("Match Found!")
						#found[compare%3].remove(cmp)

				for cmp in found[(compare+1)%3]:
			#		print("COMPARE TO " + cmp[5])
					cmpgene=cmp[5]
					if '-' in cmp[5]:
						line_index=cmp[5].index('-')
						# print(line_index)
						# print(cmp[5][:5])
						if line_index>=4 and line_index!=len(cmp[5]):
							cmpgene=cmp[5][:line_index]

			#		print("COMPARE TO " + cmp[5])
					if (resgene[5].lower() in cmpgene.lower() or cmpgene.lower() in resgene[5].lower()) and abs(int(resgene[2])-int(cmp[2]))-300:
						resloc=resgene[2]
						cmploc=cmp[2]
						if resloc not in start_check and cmploc in start_check:
							resgene[2]=cmp[2]
						match=match+1
						#print("Match Found!")
						#found[(compare+1)%3].remove(cmp)

				if match > 1:
					if len(res_genes)==0:
						res_genes.append([resgene, match])
					else:
						appnd=1
						for curr in res_genes:
							if (curr[0][5] in resgene[5] or resgene[5] in curr[0][5]):
								curr[1]=int(curr[1]) + match
								appnd=0
								break
						if appnd==1:
							res_genes.append([resgene, match])
				# else:
				# 	print("NO MATCHES FOR " + resgene[5])
				# 	print('\n')
			compare=compare+1


				#Sort through res_genes, match the same resistance genes from different accession ID's, and get concensus resistance gene(s)
	final_res=[]
	# print("RES-GENES")
	# print(res_genes)
	for resgene in res_genes:
		if len(final_res)==0:
			final_res.append(resgene)
		else:
			locmatch=0
			for final in final_res:
				if int(resgene[0][2])==int(final[0][2]):
					locmatch=1
					if int(resgene[1]) > int(final[1]):
						final=resgene
						break
				else:
					# print("COMPARE COORDINATES TO CHECK IF REPLACEMENT IS NEEDED")
					if (resgene[0][5] in final[0][5] or final[0][5] in resgene[0][5]) and (final[0][2] not in start_check and resgene[0][2] in start_check):
						final[0][2]=resgene[0][2]
			if locmatch==0:
				final_res.append(resgene)

	# print('\n')
	# print("FINAL RESISTANCE GENE LIST")
	# print('\n')
	#determine locus locations for writing to final tsv
	locus={}
	# print("DETERMINE POSITIONS")

	for final in final_res:
		change=0
		

		for start in start_check:
			if int(final[0][2])==int(start):
				locus[str(start_check.index(int(final[0][2])))]=final
				break
			else:
				if change==0:
					change=abs(int(final[0][2])-int(start))
				elif change<abs(int(final[0][2])-int(start)):
					locus[str(start_check.index(start))]=final
					break
				else:
					change=abs(int(final[0][2])-int(start))
		with open("./../output/plasmids/" + plasmid + "/abricate_results.txt", "a") as new_f: 
			new_f.write("\n"+ str(final))
		new_f.close()


	#write resistance genes into their corresponding positions in the final annotation

	with open("./../output/plasmids/" + plasmid + "/tmp.tsv", 'w') as tmp_file:
		tmpwrite=csv.writer(tmp_file, delimiter='\t')
		with open("./../output/plasmids/" + plasmid + "/the_final.tsv") as final_file:
			finalread=csv.reader(final_file, delimiter='\t')
			
			currlocus=0
			for finaline in finalread:
				if "locus_tag" not in finaline:
					if str(currlocus) in locus:
						finaline[3]=locus[str(currlocus)][0][5]
						finaline[6]=locus[str(currlocus)][0][13]
						finaline[7]=""
						finaline[9]="Antibiotic Resistance"
						# print(finaline[3])
						# print(finaline[6])
						# print(finaline[7])
						
				tmpwrite.writerow(finaline)	
				currlocus+=1
	os.system("mv ./../output/plasmids/" + plasmid + "/tmp.tsv ./../output/plasmids/" + plasmid + "/the_final.tsv")

if __name__ == "__main__": 
    parser = argparse.ArgumentParser()
    parser.add_argument('plasmid', type=str)
    plasname=(parser.parse_args()).plasmid
    main(plasname)
