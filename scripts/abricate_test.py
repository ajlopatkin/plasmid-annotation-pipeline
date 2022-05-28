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
	with open("./../output/plasmids/" + plasmid + "/gb_match/default/pk_results/" + plasmid + ".gbk") as start_chk:
		for line in start_chk.readlines():
			if line.startswith(' ' * 5 + 'CDS'):
				coordinates=line[8:].strip().strip("complement()")
				strt_index=0
				end_index=coordinates.index('..')
				start_check.append(int(coordinates[strt_index:end_index]))

	inc_group=[]
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
			print("CURR LIST OF RESISTANCE GENES IN PLASMID")
			print(list)
			for resgene in list:
				print("CURRENT GENE:" + resgene[5])
				match=1
				for cmp in found[compare%3]:
					cmpgene=cmp[5]
					if '-' in cmp[5]:
						line_index=cmp[5].index('-')
						print(line_index)
						print(cmp[5][:5])
						if line_index>=4 and line_index!=len(cmp[5]):
							cmpgene=cmp[5][:line_index]
					print("COMPARE TO " + cmp[5])
					if (resgene[5].lower() in cmpgene.lower() or cmpgene.lower() in resgene[5].lower()) and abs(int(resgene[2])-int(cmp[2]))-300:
						resloc=resgene[2]
						cmploc=cmp[2]
						if int(resloc) not in start_check and cmploc in start_check:
							resgene[2]=cmp[2]

						match=match+1
						print("Match Found!")
						#found[compare%3].remove(cmp)

				for cmp in found[(compare+1)%3]:
					print("COMPARE TO " + cmp[5])
					cmpgene=cmp[5]
					if '-' in cmp[5]:
						line_index=cmp[5].index('-')
						print(line_index)
						print(cmp[5][:5])
						if line_index>=4 and line_index!=len(cmp[5]):
							cmpgene=cmp[5][:line_index]

					print("COMPARE TO " + cmp[5])
					if (resgene[5].lower() in cmpgene.lower() or cmpgene.lower() in resgene[5].lower()) and abs(int(resgene[2])-int(cmp[2]))-300:
						resloc=resgene[2]
						cmploc=cmp[2]
						if resloc not in start_check and cmploc in start_check:
							resgene[2]=cmp[2]
						match=match+1
						print("Match Found!")
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
				else:
					print("NO MATCHES FOR " + resgene[5])
					print('\n')
			compare=compare+1

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

			#KEPT IN CASE FURTHER PROCESSING IS NEEDED: IGNORE OTHERWISE!
			#with open("./plasmids/" + plasmid + "/tmp.tsv", 'w') as tmp_file:
			#	tmpwrite=csv.writer(tmp_file, delimiter='\t')
			#	with open("./plasmids/" + plasmid + "/final.tsv") as final_file:
			#		finalread=csv.reader(final_file, delimiter='\t')
			#		for base in found:
			#			for res_gene in base:
			#				location=res_gene[6].find('/')
			#				coverage=res_gene[6][location+1:]
			#
			#				for finaline in finalread:
			#					if finaline[2]==coverage:
			#						if finaline[3] not in res_gene[5] and res_gene[5] not in finaline[3]:
			#							finaline[3]=res_gene[5]
			#							finaline[6]=res_gene[13]
			#						else:
			#							if len(res_gene[5])>len(finaline[3]):
			#								finaline[3]=res_gene[5]
			#								finaline[6]=res_gene[13]
			#					tmpwrite.writerow(finaline)
			#	os.system("cat ./plasmids/" + plasmid + "/tmp.tsv")
			#	os.system("mv ./plasmids/" + plasmid + "/tmp.tsv ./plasmids/" + plasmid + "/final.tsv")

	#Sort through res_genes, match the same resistance genes from different accession ID's, and get concensus resistance gene(s)
	final_res=[]
	print("RES-GENES")
	print(res_genes)
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
					print("COMPARE COORDINATES TO CHECK IF REPLACEMENT IS NEEDED")
					if (resgene[0][5] in final[0][5] or final[0][5] in resgene[0][5]) and (final[0][2] not in start_check and resgene[0][2] in start_check):
						final[0][2]=resgene[0][2]
			if locmatch==0:
				final_res.append(resgene)
	
	print('\n')
	print("FINAL RESISTANCE GENE LIST")
	print(final_res)
	print('\n')
	
	#determine locus locations for writing to final tsv
	locus={}
	print("DETERMINE POSITIONS")
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
	
	#write resistance genes into their corresponding positions in the final annotation
	with open("./../output/plasmids/" + plasmid + "/tmp.tsv", 'w') as tmp_file:
		tmpwrite=csv.writer(tmp_file, delimiter='\t')
		with open("./../output/plasmids/" + plasmid + "/final.tsv") as final_file:
			finalread=csv.reader(final_file, delimiter='\t')
			
			currlocus=0
			for finaline in finalread:
				if "locus_tag" not in finaline:
					if str(currlocus) in locus:
						finaline[3]=locus[str(currlocus)][0][5]
						finaline[6]=locus[str(currlocus)][0][13]
						finaline[7]=""
						finaline[9]="Antibiotic Resistance"
					
				tmpwrite.writerow(finaline)	
				currlocus+=1

			#if len(final_res) > 0:
			#	curr=0
			#	for finaline in finalread:
			#		print("#" + str(curr))
			#		curr=curr+1
			#		for resgene in final_res:
			#			print(resgene)
			#			location=resgene[0][6].find('/')
			#			coverage=resgene[0][6][location+1:]
			#			if finaline[2]==coverage:
			#				if finaline[3] not in resgene[0][5] and resgene[0][5] not in finaline[3]:
			#					finaline[3]=resgene[0][5]
			#					finaline[6]=resgene[0][13]
			#					finaline[9]="Antibiotic Resistance"
			#				else:
			#					if len(resgene[0][5])>len(finaline[3]):
			#						finaline[3]=resgene[0][5]
			#						finaline[6]=resgene[0][13]
			#						finaline[9]="Antibiotic Resistance"
			#		print(finaline)
			#		tmpwrite.writerow(finaline)
			#else:
			#	for finaline in finalread:
			#		tmpwrite.writerow(finaline)

		os.system("mv ./../output/plasmids/" + plasmid + "/tmp.tsv ./../output/plasmids/" + plasmid + "/final.tsv")

	#cat final annotation for testing purposes
	print("CATING")
	os.system("cat ./../output/plasmids/" + plasmid + "/final.tsv")

	#sort through chosen inc groups and pick the one that was chosen most frequently
	inc=[]
	print("INC GROUP")
	print(inc_group)
	for group in inc_group:
		print("INC")
		print(inc)
		if len(inc)==0:
			inc.append(group)
		else:
			for cmp in inc:
				if int(cmp[1])<int(group[1]):
					inc=[group]
	inc_concensus=""
	if len(inc)>0:
		inc_concensus=str(inc[0][0])

	#edit final gff based on final tsv annotation
	#write inc-group to the final gff, adding it before the FASTA information
	with open("./../output/plasmids/" + plasmid + "/final.gff", 'w') as final_gff:
		gffwrite=csv.writer(final_gff, delimiter='\t')
		with open("./../output/plasmids/" + plasmid + "/final.tsv") as final_tsv:
			tsvread=csv.reader(final_tsv, delimiter='\t')
			with open("./../output/plasmids/" + plasmid + "/gb_match/default/pk_results/" + plasmid +".gff") as def_gff:
				gffread=csv.reader(def_gff, delimiter='\t')
				end=0
				for records in gffread:
					if '##gff-version' in records[0] or '##sequence-region' in records[0]:
						gffwrite.writerow(records)
					elif '#FASTA' in records[0] or end==1:
						if end==0:
							if len(inc)==0:
								gffwrite.writerow(["##Inc-GROUP: None found"])
							else:
								gffwrite.writerow(["##Inc-GROUP: " + str(inc[0][0])])
						end=1
						gffwrite.writerow(records)
					else:
						for line in tsvread:
							gfeat=[]
							if line[0] not in 'locus_tag':
								features=records[8].split(';')
								if features[0][3:] in line[0]:
									finfeats=""
									current=[]
									for feature in features:
										if 'Name' in feature and line[3] not in '':
											current.append("Name")
											feature='Name=' + line[3]
										elif 'gene' in feature and line[3] not in '':
											current.append("gene")
											feature='gene=' + line[3]
										elif 'product' in feature and line[6] not in '':
											current.append("product")
											feature='product=' + line[6]
										elif 'eC_number' in feature and line[4] not in '':
											current.append("eC_number")
											feature='eC_number=' + line[4]
										elif 'COG' in feature and line[5] not in '':
											current.append("COG")
											feature='COG:' + line[5]
										if len(finfeats)==0:
											finfeats=finfeats + feature
										else:
											finfeats=finfeats + ';' + feature

									if "Name" not in current and line[3] not in '':
										finfeats=finfeats + ';' + "Name=" + line[3]
									if "gene" not in current and line[3] not in '':
										finfeats=finfeats + ';' + "gene=" + line[3]
									if "product" not in current and line[6] not in '':
										finfeats=finfeats + ';' + "product=" + line[6]
									if "eC_number" not in current and line[4] not in '':
										finfeats=finfeats + ';' + "eC_number=" + line[4]
									if "COG" not in current and line[5] not in '':
										finfeats=finfeats + ';' + "COG:" + line[5]

								gffwrite.writerow(records[0:8] + [finfeats])
								break


if __name__=="__main__":
        parser = argparse.ArgumentParser()
        parser.add_argument('plasmid', type=str)
        plasname=(parser.parse_args()).plasmid
        main(plasname)
