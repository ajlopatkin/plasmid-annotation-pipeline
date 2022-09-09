"""
LR_FILE.PY

DESCRIPTION: Creates csv file for processing in Linear Regression Script
"""

import argparse
import csv
import os

#CD-HIT runs through the gbk file of the given plasmid and looks for the given gene.
#If found, the function will return a string containing the fasta sequence of the gene
dummy_number = 1
def find_gene_id(plasmid, gene):
	print("\nPLASMID " + plasmid)
	print("LOOKING AT GENE: " + gene)

	with open("./../output/plasmids/" + plasmid + "/matches.csv") as matches:
		match_read=csv.DictReader(matches, delimiter=' ', quotechar='|')
		for match in match_read:
			gene_found=0
			print("Looking at " + match['Matches'])
			if "default" not in match['Matches']:
				with open("./../output/plasmids/" + plasmid + "/gb_match/" + match['Matches'] + "/pk_results/" + plasmid + ".gbk", "r") as proteins:
					for line in proteins.readlines():
						if "CDS" in line:
							if gene_found==1:
								break
						if "/gene" in line and gene in line:
							gene_found=1
						if gene_found==1 and line.startswith(' ' * 21 + '/protein_id='):
							prtindex=line.find('=')
							protein_id=line[prtindex+2:-2]
							return protein_id
						if gene_found==1 and line.startswith(' ' * 21 + "sequence"):
							protein_id=(line.split(':'))[-1][0:-2]
							return protein_id

		return ""

#Obtain GOID for given gene in given plasmid
def get_GOID(plasmid, gene):
	print("RUNNING GO ID SEARCH")
	with open('./../output/plasmids/'+ plasmid + '/final.tsv') as tsv:
		read=csv.reader(tsv, delimiter='\t')
		for row in read:
			if gene in row[3] and len(gene)==len(row[3]):
				if len(row[7])>0:
					print("GO ID " + str(row[7]) + " for " + gene)
					return str(row[7])
		print("NOID for " + gene)
		return "NOID"

#If gene name is cut off, use the gene id to find the full gene name
def find_faa(id):
	gene_line=""
	faa=open("./../output/prtgene.faa")
	for line in faa.readlines():
		if id in line:
			gene_line=line
			break
	if len(gene_line)==0:
		return ""
	index1=gene_line.find("~~~")+3
	index2=gene_line[index1:].find("~~~")
	index2=index1+index2
	return gene_line[index1:index2]

#Load gene sequences for the genes in a given plasmid
def load_genes(plasmid):
	genes={}
	print("Loading genes from " + plasmid)
	with open("./../output/plasmids/" + plasmid + "/final.gbk") as gbk:
		genefound=""
		sequence=""
		product=""
		trans=False
		for line in gbk.readlines():
			if "CDS" in line:
				if len(genefound)>0:
					genes[genefound]=">"+ find_gene_id(plasmid, genefound[0:genefound.index("|")]) +"~~~" + genefound + "~~~" + "" + product + "\n" + sequence + "\n"
				genefound=""
				sequence=""
				product=""
				trans=False
			if "ORIGIN" in line:
				if len(genefound)>0:
					genes[genefound]=">"+ find_gene_id(plasmid, genefound[0:genefound.index("|")]) +"~~~" + genefound + "~~~" + "" + product + "\n" + sequence + "\n"
				break
			if "/product" in line:
				start=line.index('"')+2
				end=line[start:].index('"')
				product=line[start:start+end]
			if "/gene" in line:
				start=line.index('"')+2
				end=line[start:].index('"')
				genefound=line[start:start+end]
				if len(genefound)>0:
					genefound=genefound + "|" + get_GOID(plasmid, genefound)
			if ("/translation" in line or trans==True) and len(genefound)>0:
				if "/translation" in line:
					line=line[line.index('"')+2:]
					trans=True
				if '"' in line:
					line=line[:line.index('"')]
				sequence = str(sequence) + str(line.replace("\n", " ")).replace(" ", "")
	print(genes)
	print('\n')
	return genes

#Parse given acquisition cost data file for given plasmid
#Will be updated to take in file path information
def parse_sample(plasmid):
	print("PARSE SAMPLE " + plasmid)
	labels=[' ']
	pls_info=[]

	antb=0
	inc=0
	AC=0
	err=0

	#CHANGE FILE NAME HERE IF NEEDED
	with open('./../Summary_Table_by_Plasmid.csv', 'r', encoding='utf-8-sig') as csvfile:
		read=csv.reader(csvfile, delimiter=',')
		for line in read:
			if "Plasmid" in line:
				#labels=line
				for i in range(0, len(line)):
					label=line[i].lower()
					if "antibiotic" in label or "drug" in label:
						antb=i
					elif "inc" in label and len(label)==3:
						inc=i
					elif "ac_avg" in label or "normalized_ac" in label:
						AC=i
					elif "normalized_ac_error" in label or "ac_std" in label:
						err=i
				labels=["Number", "Plasmid", line[antb], line[inc], "normalized_AC", line[err]]
			elif plasmid == line[0]:
				pls_info.append([line[antb], line[inc], line[AC], line[err]])

	pls_info.append(labels)

	return pls_info

#Parse final.tsv file for the given plasmid.
#Return list containing list of genes and list containing number
#of genes with the given function
def parse_ann(plasmid):
	genes=[]
	num_genes=0
	num_ar=0
	num_rep=0
	num_con=0
	num_met=0
	num_str=0
	num_txn=0
	num_atxn=0


	try:
		print("ANN PARSE FOR " + plasmid)
		with open('./../output/plasmids/' + plasmid + '/final.tsv') as anntsv:
			read=csv.reader(anntsv, delimiter='\t')
			for line in read:
				if "locus_tag" not in line:
					num_genes=num_genes+1
					if "antibiotic" in str(line[-1]).lower():
						num_ar=num_ar+1
					elif "replication" in str(line[-1]).lower():
						num_rep +=1
					elif "conjugation" in str(line[-1]).lower():
						num_con +=1
					elif "metabolism" in str(line[-1]).lower():
						num_met +=1
					elif "stress" in str(line[-1]).lower():
						num_str +=1
					elif "antitoxin" in str(line[-1]).lower():
						num_atxn +=1
					elif "toxin" in str(line[-1]).lower():
						num_txn +=1
					if len(line[3])>0:
						if len(line[7])>0:
							genes.append("" + line[3] + "|" + line[7] + "")
						else:
							genes.append("" + line[3] + "|NOID")
					else:
						genes.append("")
	except:
		with open('./../output/plasmids/' + plasmid + '/final_KEGG.tsv') as anntsv:
			read=csv.reader(anntsv, delimiter='\t')
			for line in read:
				if "locus_tag" not in line:
					num_genes=num_genes+1
					if len(line[3])>0:
						genes.append("" + line[3] + "|" + line[7] + "")
					else:
						genes.append("")
	print("PLASMID " + plasmid + " has " + str(num_genes) + " total genes and " + str(num_ar) + " ar genes!")
	genes.append([num_genes, num_ar, num_rep, num_con, num_met, num_str, num_txn, num_atxn])
	return genes

#Get length of the given plasmid
def get_length(plasmid):
	with open('./../output/plasmids/'+ plasmid + '/final.gff') as gff:
		read=csv.reader(gff, delimiter='\t')

		index=0
		for row in read:
			if index==1:
				return row[0].split(" ")[-1]
			index +=1

#Creates file by obtaining information from the annotation and acquisition cost files for each plasmid
def create_file(plasmids):
	with open('./../output/PlasmidRegression.csv', 'w') as reg_csv:
		write_csv=csv.writer(reg_csv)
		plasmid_genes=[]
		plasmid_acqu=[]

		num_genes=[]
		ar_genes=[]
		desc_num=[]

		label_num=0
		labels=[]
		acqu_labels=[]
		label=0
		curr=1
		file=[]
		new_plasmids=[]
		for plasmid in plasmids:
			print("PROCESSING PLASMID: " + plasmid)
			ann_parse=parse_ann(plasmid)
			sample_info=parse_sample(plasmid)
			acqu_labels=sample_info[-1]


			if label==0:
				labels.append(acqu_labels[0])
				labels.append(acqu_labels[1])
				labels.append("Length")
				labels.append("Number-of-Genes")
				labels.append("Antibiotic-Genes")
				labels.append("Replication-Genes")
				labels.append("Conjugation-Genes")
				labels.append("Metabolism-Genes")
				labels.append("Stress-Genes")
				labels.append("Toxin-Genes")
				labels.append("Antitoxin-Genes")
				label=1
			sample_info.remove(acqu_labels)

			if int(ann_parse[-1][0])>label_num:
				label_num=ann_parse[-1][0]
			num_genes.append(ann_parse[-1][0])
			ar_genes.append(ann_parse[-1][1])
			desc_num.append(ann_parse[-1][2:])
			ann_parse.remove(ann_parse[-1])

			for sample in sample_info:
				plasmid_genes.append([plasmid, ann_parse])
				plasmid_acqu.append([plasmid, sample])
				new_plasmids.append(plasmid)

		#Run through all found genes and only analyse genes with gene id's that have not already been identified
		#Append unique genes into .faa file
		plasmids=new_plasmids
		all_genes={}

		for plas in plasmids:
			all_genes[plas]=load_genes(plas)

		unique_ids=[]
		all_ids=[]
		seen=[]

		print("\nFINDING UNIQUE GENES")
		for plas in plasmid_genes:
			print("Plasmid: " + plas[0])
			for gene in plas[1]:
				print("Looking at gene: " + gene)

				if gene not in seen and len(gene)>0:
					print(all_genes[plas[0]])
					try:
						result=all_genes[plas[0]][gene]
						end=result.find('~~~')
						id=result[1:end]

						if id not in unique_ids:
							unique_ids.append(id)
							sequ_match=open("./../output/prtgene.faa", "a")
							sequ_match.write(result)
							sequ_match.close()
						seen.append(gene)
						all_ids.append(id)
					except:
						print("Given gene already seen or duplicate in the final annotation!")
						
		#Run PSI-CD-Hit at 30% homology to make sure all homologous genes have been found
		os.system("./cd-hit/psi-cd-hit/psi-cd-hit.pl -i ./../output/prtgene.faa -o ./../output/psiresult.txt -c 0.3 -ce 1e-6")

		#Run through clstr file to find unique genes to be added to PlasmidRegression.csv
		#Also analyse what genes they are homologous with to
		clstr=open("./../output/psiresult.txt.clstr")
		clstr_replc=[]
		hom_genes=[]
		curr=[]
		last_item=0
		curr_process=0
		clstr_file=clstr.readlines()
		for line in clstr_file:
			if (line.startswith(">Cluster")) or clstr_file.index(line)==len(clstr_file)-1:
				curr_process=curr_process+1
				if clstr_file.index(line)==len(clstr_file)-1:
                                	index1=line.find("~~~")+3
                                	index2=(line[index1:].find("~~~"))
                                	if index2==-1:
                                        	id1=line.find(">")+1
                                        	id2=line.find("~~~")
                                        	curr.append(find_faa(line[id1:id2]))
                                	else:
                                        	index2=index2+index1
                                        	curr.append(line[index1:index2])
				if len(curr)==1:
					hom_genes.append(curr[0])
					clstr_replc.append(curr[0])
					curr=[]
				elif len(curr)>1:
					better=[]
					process=0
					for gene in curr:
						if len(better)==0:
							better=gene
						elif gene in better and ('_' not in gene) and ('_' in better):
							hom_genes.append(better)
							better=gene
							process=process+1
						else:
							hom_genes.append(gene)
							process=process+1

					for i in range(last_item, process+last_item):
						clstr_replc.append(better)
					last_item=process+1
					curr=[]
			else:
				index1=line.find("~~~")+3
				index2=(line[index1:].find("~~~"))
				if index2==-1:
					id1=line.find(">")+1
					id2=line.find("~~~")
					curr.append(find_faa(line[id1:id2]))
				else:
					index2=index2+index1
					curr.append(line[index1:index2])

		clstr.close()

		for i in range(0, len(plasmid_genes)):
			for gene in plasmid_genes[i][1]:
				if len(gene)>0:
					if gene in hom_genes:
						rplc=hom_genes.index(gene)
						plasmid_genes[i][1][plasmid_genes[i][1].index(gene)]=clstr_replc[rplc]
						if clstr_replc[rplc] not in labels:
							labels.append(clstr_replc[rplc])
		for label in acqu_labels[2:]:
			labels.append(label)

		write_csv.writerow(labels)

		num=1
		prev_plasmid = ""
		for plasmid in plasmids:
			unique_plasmids = list(dict.fromkeys(plasmids))
			index = unique_plasmids.index(plasmid)

			curr_line=[num, plasmid, get_length(plasmid), num_genes[index], ar_genes[index]]
			for desc in desc_num[index]:
				curr_line.append(desc)

			for gene in labels[11:]:
				if gene in "Antibiotic":
					break
				if gene in plasmid_genes[num-1][1]:
					curr_line.append(1)
				else:
					curr_line.append(0)

			for info in plasmid_acqu[num-1][1]:
				curr_line.append(info)

			if plasmid == "R1":
				print(curr_line)

			write_csv.writerow(curr_line)
			num +=1

def main():
	with open("./../fasta.csv") as csvfile:
		read=csv.reader(csvfile, delimiter='\t')
		list=[]
		for row in read:
			list.append(row[0])
		create_file(list)

if __name__=="__main__":
	main()
