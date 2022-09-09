# plasmid_annotation
import csv
import os
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from collections import OrderedDict 
import pandas as pd

inc_group = []
pla_finder = []

inc_fold = "output/plasmids/pox38/ab_results/plasmidfinder.txt"

with open(inc_fold) as csvfile:
	inc_read=csv.reader(csvfile, delimiter='	')
	for inc_found in inc_read:
		if(inc_found[5] not in "GENE"):
			pla_finder.append(inc_found)

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
inc = inc_group[0]
inc = inc[0]
inc = inc[0:4]

backbone_fasta = "backbone_genes/" + inc + "_backbone.fasta"


with open("output/plasmids/pox38/gene_locators.tsv", 'wt') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    tsv_writer.writerow(['gene', 'start', 'end'])

blastp_results = []
for gene in SeqIO.parse(backbone_fasta,"fasta"):
	entry = str(">" + gene.description + "\n" + gene.seq)
	f1 = open("test.txt", "w")
	f1.write(entry)
	f1.close()
	f2 = open("test.txt", "r")
	blastp_cline = NcbiblastpCommandline(query = "output/plasmids/pox38/gb_match/default/pk_results/pox38.gbk", subject = 'test.txt', evalue=0.001, outfmt=5, out="output/plasmids/pox38/backbone_"+gene.description+".xml")
	stdout, stderr = blastp_cline()
	f2.close()
	result_handle = open("output/plasmids/pox38/backbone_"+gene.description+".xml")
	blast_record = NCBIXML.read(result_handle)
	## Now we want to loop through the top alignments and find their sequences
	for alignment in blast_record.alignments: 
		# Loop through our HSPS
		for hsp in alignment.hsps:
			start = hsp.query_start
			end = hsp.query_end
			with open("output/plasmids/pox38/gene_locators.tsv", 'a') as out_file:
				csv_writer = csv.writer(out_file, delimiter = '\t')
				csv_writer.writerow([gene.description, start, end])


res_genes=[]

#For each accession ID match, open their abricate files and sort their found resistance genes or inc groups into corresponding lists
#Print statements commented out; uncomment them for debugging purposes!
res_ncbi=[]
res_card=[]
res_finder=[]
pla_finder=[]
found=[res_card, res_finder, res_ncbi]

directory="output/plasmids/pox38/ab_results"
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
		# print("CURR LIST OF RESISTANCE GENES IN PLASMID")
		# print(list)
		for resgene in list:
			# print("CURRENT GENE:" + resgene[5])
			match=1
			for cmp in found[compare%3]:
				cmpgene=cmp[5]
				if '-' in cmp[5]:
					line_index=cmp[5].index('-')
					# print(line_index)
					# print(cmp[5][:5])
					if line_index>=4 and line_index!=len(cmp[5]):
						cmpgene=cmp[5][:line_index]
				# print("COMPARE TO " + cmp[5])
				if (resgene[5].lower() in cmpgene.lower() or cmpgene.lower() in resgene[5].lower()) and abs(int(resgene[2])-int(cmp[2]))-300:
					resloc=resgene[2]
					cmploc=cmp[2]
					# if int(resloc) not in start_check and cmploc in start_check:
					# 	resgene[2]=cmp[2]

					match=match+1
					#found[compare%3].remove(cmp)

			for cmp in found[(compare+1)%3]:
				# print("COMPARE TO " + cmp[5])
				cmpgene=cmp[5]
				if '-' in cmp[5]:
					line_index=cmp[5].index('-')
					# print(line_index)
					# print(cmp[5][:5])
					if line_index>=4 and line_index!=len(cmp[5]):
						cmpgene=cmp[5][:line_index]

				# print("COMPARE TO " + cmp[5])
				if (resgene[5].lower() in cmpgene.lower() or cmpgene.lower() in resgene[5].lower()) and abs(int(resgene[2])-int(cmp[2]))-300:
					resloc=resgene[2]
					cmploc=cmp[2]
					# if resloc not in start_check and cmploc in start_check:
					# 	resgene[2]=cmp[2]
					match=match+1
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


		gene_start = []
		for gene in res_genes: 
			gen = (int(gene[0][2]))
			print(gen)
			print("\n")
			print(gene[0][5])

			

			# if gen not in gene_start: 
			# 	gene_start.append(gen)

			# # username.translate({ ord(c): None for c in "._!" })
			# gen = gen.translate({ord(c): None for c in "._-'()"}).lower()
			# #print(gen[0:].lower())
			# if gen not in gene_name: 
			# 	gene_name.append(gen)
		# print(len(gene_start))
		# for start in gene_start: 
		# 	print(start)
		# 	print("\n")
		# if (gene_start[0]) == (gene_start[3]):
		# 	print("Bruh. they're the same")
		# else: 
		# 	print("we have another problem")
		# the_genes = []
		# for gene in gene_name: 
		# 	for i, genes in enumerate((gene)):
		# 		print(genes[i])
		# 	# if gene not in the_genes: 
			# 	the_genes.append(gene)

		# print(len(the_genes))
		# print(the_genes[0], the_genes[1], the_genes[2])

		# for gene in the_genes:
		# 	print(gene)

		# if (the_genes[0]) == (the_genes[3]): 
		# 	print("Bruh. they're the same")
		# else: 
		# 	print("we have another problem")
		# gene_name = []
		# # start = []
		# # end = []
		# for gene in res_genes:
		# 	gen = (str(gene[0][5]))
		# 	if gen not in gene_name:
		# 		gene_name.append(gen)

		# 	gen = gene[0]
		# 	if gen[5] not in gene_name:
		# 		gene_name.append(gen[5])
		# 	# start.append(gen[2]) 
		# 	# end.append(gen[3])
		# for i in gene_name:
		# 	print(i)

		# gene_names = (OrderedDict.fromkeys(gene_name)
		# end = [i for n, i in enumerate(end) if i not in end[:n]]

		# result = [] 
		# for i in start: 
		#     if i not in result: 
		#         result.append(i) 
		# for i in result: 
		# 	print(i)
		# abr_gene = []

		# for gene in res_genes:
		# 	gen = gene[0]
		# 	g_name = str(gen[5])
		# 	start = (gen[2])
		# 	print(start[0])
		# 	if start not in abr_gene:
		# 		abr_gene.append(start)


			# print(str(gen[2]))
			# print(str(gen[3]))
			# with open("output/plasmids/R100/gene_locators.tsv", 'a') as out_file:
			# 	tsv_writer = csv.writer(out_file, delimiter='\t')
			# 	tsv_writer.writerow([gen[5], gen[2], gen[3]])
		#br_gene = abr_gene.sort()
		# print(len(abr_gene))
		# for gene in abr_gene: 
		# 	print(gene)
		# #data_frame = pd.read_csv()


# 	count = count + 1
# print("There are " + str(count) + " backbone genes so far for " + inc + "!")
