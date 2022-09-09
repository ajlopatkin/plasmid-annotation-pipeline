#!/usr/bin/python3

"""
Description: Call on UNIPROT to find a UNIPROT ID for each gene of the given plasmid and then use the GO ID linked
to that UNIPROT ID to label the gene's function

EDITED FROM PROKKA2KEGG by Heyu Lin (heyu.lin(AT)student.unimelb.edu.au)
 
"""

import os
import re
import csv
import gzip
import curses
import argparse
import json
import http.client
import requests
import urllib.request
from io import StringIO
from email.generator import Generator
import requests
import zlib
from goatools.base import download_go_basic_obo
from goatools.obo_parser import GODag
from goatools.utils import get_b2aset
from goatools.godag.go_tasks import get_go2parents
from goatools.obo_parser import GOTerm


http.client._MAXHEADERS=1000

__author__ = "Heyu Lin"
__contact__ = "heyu.lin(AT)student.unimelb.edu.au"

print("PARSING ARGUMENTS")
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', metavar='input_gbk', dest='i',
                    type=str, required=True)
parser.add_argument('-o', '--output', metavar='output', dest='o',
                    type=str, required=True)
parser.add_argument('-p', '--plasmid', metavar='plasmid', dest='p',
                    type=str, required=True)

args = parser.parse_args()

#Download go-basic and sort given GO ID's to use as references
obs_ids=[]
obs_items=[]
val_ids=[]
val_items=[]

fin_dag=download_go_basic_obo("go-basic.obo")
godag=GODag(fin_dag, optional_attrs={'consider', 'replaced_by', 'relationships'}, load_obsolete=True)
#print(godag)

for go_item in godag.values():
	if go_item.is_obsolete:
		obs_items.append(go_item)
		obs_ids.append(go_item.id)
	else:
		val_items.append(go_item)
		val_ids.append(go_item.id)


#Call PSI-BLAST if gene is left w/no functional description
def psi_blast(gene, plasmid):
    
	sequence=""
	
	with open("./../output/plasmids/" + plasmid + "/the_final.gbk", "r") as proteins:
		cds=0
		locus=0
		found=0
		read=0
		prod=0

		for line in proteins.readlines():
			if line.startswith(' ' * 5 + 'CDS'):
				cds=1
				found=0
			if line.startswith(' ' * 21 + '/locus_tag=') and cds==1:
				locus=1
			if line.startswith(' ' * 21 + '/gene=') and cds == 1:
				nmindex=line.find('=')
				if gene.lower() in line[nmindex+2:-2].lower():
					found=1
					print("GENE FOUND")
				else:
					cds=0
					locus=0
			if found==1 and line.startswith(' ' * 21 + '/protein_id='):
				prtindex=line.find('=')
				protein_id=line[prtindex+2:-2]
			if found==1 and line.startswith(' ' * 21 + "sequence") and len(protein_id)==0:
				protein_id=(line.split(':'))[-1][0:-2]
				print(protein_id)
			if line.startswith(' ' * 21 + '/translation='):
				read=1
				prod=0
			if found==1 and read==1 and prod!=1:
				nmindex=line.find('=')
				trans=line[nmindex+2:-2]
				if (not (line.startswith(' ' * 21 + '/')) or line.startswith(' ' * 21 + '/translation=')) and not line.startswith(' ' * 21 + 'sequence') and not line.startswith('ORIGIN'):
					print("HERE")
					print(sequence)
					sequence = sequence + (trans.replace("\n", " ")).replace(" ", "")
					print(sequence)
				elif line.startswith(' ' * 5 + 'CDS') or line.startswith(' ' * 5 + 'gene'):
					print("RESET")
					cds=0
					if line.startswith(' ' * 5 + 'CDS'):
						cds=1
					read=0
					found=0
					locus=0
				elif line.startswith('ORIGIN'):
					break
			if found==1 and (line.startswith(' ' * 21 + "/product") or prod==1):
				prdcti=line.find('=')
				product=line[prdcti+2:-2].replace("\n", "")
				prod=1
				read=0
				print("ADDED TO PRODUCT")
				print(product)
			if len(sequence)>0:
				break
			if len(product)==0:
				product=gene
			else:
				print("PRODUCT PROCESSING")
				print(product)
				product=product.replace("\n", "")				

	if len(sequence)>0:
		job_id=os.system("python -c " + "'" + "from psiblast import *; print(psiblast.serviceRun(" + "'" + "jc4919@barnard.edu" + "'" + ", " + "'" + gene + "'" + ", [ " + "'" + "sequence:" + "'" "'" + sequence + "'" + "]" + ", " + "'" + "U" "))" + "'")
		print("JOB ID " + str(job_id))		
		#job_id=serviceRun("jc4919@barnard.edu", "", [""+ sequence +""])
		#print(job_id)	
		#print(requests.get("https://www.ebi.ac.uk/Tools/services/rest/psiblast/resulttypes/" + job_id))
		#print(results)

#Process go_id function descriptions based on parent descriptions to make labeling easier when mapping 
def process_go_id(go_id, func_desc, general):
	keywords=['conjugation', 'replication', 'pilus', 'plasmid', 'stress', 'antitoxin', 'toxin', 'damage stimulus', 'metabolism', 'metabolic']
	key=0

	print("\nProcessing GO-ID: " + go_id + " (" + func_desc + ")")

	for word in keywords:	
		if word in func_desc:
			key=1

	if key==1:
		return func_desc
	elif "antibiotic" in func_desc or "antimicrobial resistance" in func_desc:
		return ""
	else:
		found=0
		for word in keywords:
			if word in general:
				return func_desc + "-" + word

		addition=""
		#USE PARENT ID'S TO EDIT FUNC DESC
		#if go_id in val_ids:
		#	loc=val_ids.index(go_id)
		#	parents=val_items[loc].parents
		#	for parent in parents:
		#		for word in keywords:
		#			if word in parent.name:
		#				addition=word
		#				if word in keywords[0:-2]:
		#					break
		return func_desc + "-" + addition
	
#Get GO ID for each uniprot id found by calling on UNIPROT 
#Print statements commented out; uncomment for debugging if needed!
def get_go_id(upro_ids):
	print("GET GO ID's FOR GIVEN PLASMID")
	print(upro_ids)
	interest=['conjugation', 'replication', 'pilus', 'plasmid', 'stress', 'antitoxin', 'toxin', 'damage stimulus', 'metabolism', 'metabolic']


	go_ids=[]
	for uid in upro_ids:
		#print(uid[0])
		if len(uid[2])>0:
			poss_ids=[]
			#found=0
			find="U"
			for id in uid[2]:
				#print(id)			
				uni_http=requests.get("https://www.uniprot.org/uniprot/" + id + ".txt")
				uni_http.encoding="ISO-8859-1"
				uni_http=requests.get("https://www.uniprot.org/uniprot/" + id + ".txt")
				uni_txt=uni_http.text.split("\n")

				keywords=[]
				#Find keywords to use for processing	
				for unitem in uni_txt:
					if "KW" in unitem[0:2]:
						keys=unitem[5:].split(';')
						for item in keys:
							#print(item)
							#print(len(item))
							if len(item)>0:
								if '.' in item[len(item)-1]:
									keywords.append(item[:len(item)-1])
								else:
									keywords.append(item)	
					#print(unitem)
		
				#find="U"
				#poss_ids=[]
				found=0
				uid_gos=[]
				for line in uni_txt:
					if 'DR' in line[0:2] and "GO;" in line:
						info=[]
						found=1	
						spl_line=line.split(" ", 5)
						info.append(uid[0])
					
						curr_id=spl_line[4][0:-1]
						if curr_id in obs_ids:
							loc=obs_ids.index(spl_line[4][0:-1])
							#for item in obs_items[loc].alt_ids:
								#print(item)
							if "" in obs_items[loc].alt_ids:
								continue
							else:
								curr_id=""
								for alt_id in obs_items[loc].alt_ids:
									if alt_id not in obs_ids:
										curr_id=alt_id
										break
								if "" in curr_id:
									continue
								#curr_id=obs_items[loc].alt_ids[0]
					
						info.append(curr_id)
						find=""
						if uid[1]:
							find="R"
						else:
							find="U"
						if 'P' in spl_line[5][0:1]:
							find=find + "P"
						elif 'C' in spl_line[5][0:1] and 'P' not in find:
							find=find + "C"
						elif 'F' in spl_line[5][0:1] and 'P' not in find:
							find=find + "F"
						info.append(find)
						funct=spl_line[5][2:-1].split(';', 1)
						info.append(process_go_id(curr_id, funct[0], keywords))
						uid_gos.append(info)
				#print("KEYWORDS")
				#print(keywords)
				#print("uid_gos")
				#print(uid_gos)
				#print("INFO END")
				if found==0:
					if len(keywords)>0:
						#print("KEYWORDS ANALYSIS")
						keyfind=0
						interest=['conjugation', 'replication', 'pilus', 'plasmid', 'stress', 'toxin', 'damage stimulus', 'metabolism', 'metabolic']
						for key in interest:
							if key in keywords:
								keyfind=1
								rev_stat="U"
								if 'R' in find:
									rev_stat="R"
								poss_ids.append([uid[0], "", rev_stat, key])
								break
						if keyfind==0:
							poss_ids.append([uid[0], "", "", ""])
					else:
						poss_ids.append([uid[0], "", "", ""])
				else:
					match_go=[]
					for go_found in uid_gos:
						if len(match_go)==0:
							match_go=go_found
						else:
							if 'P' not in match_go[2] and 'P' in go_found[2]:
								match_go=go_found
								continue
							elif 'P' in match_go[2] and 'P' in go_found[2]:
								key_find=0
								for key in interest:
									if key in match_go[2]:
										key_find=1
										break
									if key in go_found[1]:
										key_find=2

								if key_find==2:
									match_go=go_found
					poss_ids.append(match_go)

			#for go_id in poss_ids:
			#	if find in go_id[2]:
			#			go_ids.append(go_id)
			#			break

			concensus=[]
			for go_find in poss_ids:
				print("GO FIND")
				print(go_find)
				if len(concensus)==0:
					if len(go_find[1])>0:
						temp=go_find
						temp.append(1)
						concensus.append(temp)
					else:
						continue
				else:
					if len(go_find)>0:
						found=0
						for final_go in concensus:
							print(go_find)
							print(final_go)
							if go_find[1] in final_go[1]:
								found=1
								final_go[-1]+=1
								break
						if found!=1:
							temp=go_find
							temp.append(1)
							concensus.append(temp)

			go_app=[]
			for id in concensus:
				if len(go_app)==0:
					go_app=id
				else:
					if id[-1]>go_app[-1]:
						go_app=id
			print("GO APP")
			print(go_app)
			if len(go_app)>0:
				go_ids.append(go_app[0:len(go_app)-1])
			else:
				print("EMPTY GO APP")
				go_ids.append([uid[0], "", "", ""])
								
		else:
			print("NO GO ID B/C NO UNIPROT ID MATCH FOR " + uid[0])
			go_ids.append([uid[0], "", "", ""])
	
	#print(go_ids)
	print("PROCESSED GO ID's FOR GIVEN PLASMID")			
	return go_ids 

#Parse through GBK parser to determine uniprot id's for each gene, calling on UNIPROT where necessary to find the right uniprot id for it
def get_uniprot(gene_name):
	endpoint = f"https://rest.uniprot.org/uniprotkb/search?query=gene:{gene_name}+AND+taxonomy_id:2"
	resp = requests.get(endpoint)
	
	if resp.status_code != 200:
		return ''
		
	resp_json = json.loads(resp.text)['results']
	
	if len(resp_json) == 0:
		return ''
	else:
		return resp_json[0]['primaryAccession']

def gbk_parser(gbk):
	"""
	gbk: gbk genome file generated by Prokka
	"""
	arr = []  # output array containing locus_tag and UniProtKB
	
	print("\nFINDING UNIPROT ID'S FOR EACH GENE")
	with open(gbk) as input:
		cds = 0
		locus = 0
		pattern_locus = re.compile('"(.*)"')
		pattern_uniprotkb = re.compile('RefSeq:(.*)"')
		gene_up=""
		reviewed=False
		current=0
		find=0
		for line in input.readlines():
			#print(line)
			#print(current)
			current=current+1
			if line.startswith(' ' * 5 + 'CDS'):
				cds = 1  # This is a CDS
				if find==1:
					print("UNIPROT MATCH FOUND " + gene_up)
					arr.append([locus_tag, reviewed, gene_up])		
					locus=0
					find=0
				gene_up=""
			if line.startswith(' ' * 21 + '/locus_tag=') and cds == 1:
				locus_tag = pattern_locus.findall(line)[0]
				locus = 1  # locus_tag was read

			if len(gene_up)>1:
				arr.append([locus_tag, reviewed, gene_up])
				cds=0
				locus=0
				find=0
			gene_up=""
			# if line.startswith(' ' * 21 + '/inference=') and locus == 1:
			# 	inf = 1
			# if line.startswith(' ' * 21 + 'sequence:RefSeq:') and inf ==1 : 
			# 	uniprotkb = pattern_uniprotkb.findall(line)[0]
			# 	print(uniprotkb)
			# 	arr.append([locus_tag, True, [uniprotkb]])
			# 	gene_up=""
			# 	cds = 0
			# 	locus = 0
			if line.startswith(' ' * 21 + '/gene=') and cds == 1 and locus == 1:
				gene_up=""

				nmindex=line.find('=')
				gene=line[nmindex+2:-2]
				uniprotkb = get_uniprot(gene)
				arr.append([locus_tag, True, [uniprotkb]])
				cds = 0
				locus = 0
				#if (gene.find('-'))>-1 and 'tra' in gene:
				#      rmv=gene.find('-')
				#      gene=gene[0:rmv]

				print("CURR GENE: " + gene)
				#Call on uniprot to return list of matching uniprot id's for given gene
				uni_http=requests.get("https://www.uniprot.org/uniprot/?query=gene:" + gene + "+AND+taxonomy:2&taxanomy=2&active=yes&limit=20&format=tab&columns=id,reviewed,genes(PREFERRED),lineage(ALL)")
				uni_http.encoding="ISO-8859-1"
				uni_txt=uni_http.text
				index=0
				
				#for match in uni_txt.split('\n'):
				#	print(match)
				#print('\n')
				up_match=[]
				for match in uni_txt.split('\n'):
					#print(match)
					if len(up_match)==3:
						break
					curr=match.split('\t')
					if len(curr)>2:
						if gene.lower() in curr[2].lower() and "Bacteria" in curr[3]:
							if "unreviewed" not in curr[1]:
								reviewed=True
							else:
								reviewed=False
							up_match.append(curr[0])
                        
		#CODE FOR ANALYZING ID'S BASED ON ORDER
                #fir_ur=0
                #for match in uni_txt.split('\n'):
                        #print(match)
                        #if index==20:
                                #break
                        #curr=match.split('\t')
                        #if len(curr)>2:
				#prioritize first reviewed uniprot id found and then the first unreviewed uniprot id 	
                         #       print("analysing")
                          #      if gene.lower() in curr[2].lower() and "Bacteria" in curr[3]:
                           #                  if "unreviewed" not in curr[1]:
                            #                         reviewed=True
                             #                else:
                              #                       reviewed=False
                               #              if "unreviewed" not in curr[1] or fir_ur==0:
                                #                 gene_up=curr[0]
                                 #                fir_ur=1
                                  #           if "unreviewed" not in curr[1]:
                                   #              break
                        #index=index+1
               
				gene_up=up_match
				if len(gene_up)>0:
					find=1

	
			if line.startswith(' ' * 5 + 'CDS') and locus == 1:
				arr.append([locus_tag, True, ''])
				cds = 0
				locus = 0
	return arr

def output(arr, outfile):
	#print("OUTFILE")
	#print(outfile)
	#print(arr)
	with open(outfile, 'w') as fo:
		fo.write("locus_tag" + "\t" + "GO ID" + "\t" + "Type" + "\t" + "Function" + "\n")
		for cds in arr:
			print(cds)
			if len(cds)>2:
				print("WRITING MATCH")
				fo.write(cds[0] + "\t" + cds[1] + "\t" + cds[2] + "\t" + cds[3] + "\n")
			else:
				print("WRITING BLANK")
				fo.write(cds[0] + "\t" + "" + "\t" + "" + "\t" + "" + "\n")


def main():
	print("PROKKA2GO")
	#psi_blast("tir", "pOXA48K8")
	mapping_array = gbk_parser(args.i)
	output(get_go_id(mapping_array), args.o)

if __name__ == '__main__':
	main()
