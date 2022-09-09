# update_gff.py
import csv
import os
import argparse


def main(plasmid):
	with open("./../output/plasmids/"+plasmid+"/abr_genes.txt") as f:
	        lines = f.readlines()
	        x = lines
	        inc_group = x[0]
	        inc_group = str(inc_group)
	        inc = inc_group[0:4]

	        f.close()

	with open("./../output/plasmids/"+plasmid+"/the_final.gff", 'w') as final_gff:
			gffwrite=csv.writer(final_gff, delimiter='\t')
			with open("./../output/plasmids/"+plasmid+"/the_final.tsv") as final_tsv:
				tsvread=csv.reader(final_tsv, delimiter='\t')
				with open("./../output/plasmids/"+plasmid+"/gb_match/default/pk_results/"+plasmid+".gff") as def_gff:
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

if __name__ == "__main__": 
	parser = argparse.ArgumentParser()
	parser.add_argument('plasmid', type=str)
	plasname=(parser.parse_args()).plasmid
	main(plasname)
