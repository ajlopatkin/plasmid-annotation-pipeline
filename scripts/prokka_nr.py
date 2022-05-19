#!/usr/bin/python
'''
DESCRIPTION:

1) Creates directories for storing prokka results of each accession ID
2) Runs prokka using each accession ID + runs prokka on default mode
3) Create a new csv file to compile most accurate matches

'''

import os
import csv
import argparse


def main(plasmid):
	#os.system("conda init")
	new_matches=[]
	os.system("mkdir ./../output/plasmids/" + plasmid + "/merge")

        #Open csvfile containing accession + reference ID's and run prokka for each one
	with open("./../output/plasmids/" + plasmid + "/matches.csv") as csvfile:
		reader=csv.DictReader(csvfile, delimiter=' ', quotechar='|')
		for row in reader:
			if str(row['Matches']) in "default":
				os.system("mkdir ./../output/plasmids/" + plasmid + "/gb_match/default")	
				os.system("mkdir ./../output/plasmids/" + plasmid + "/gb_match/default/pk_results")
				outdir="./../output/plasmids/" + plasmid + "/gb_match/default/pk_results"
		
				#Run Prokka
				#OLD CODE
				#os.system("CONDA_BASE=$(conda info --base)")
				#os.system("source $CONDA_BASE/etc/profile.d/conda.sh")
				#os.system("conda init")				

				#os.system("conda init")
				#os.system("conda activate prokka_env")
				#os.system("prokka --kingdom Bacteria --outdir " + outdir + " --prefix " + plasmid + " --force ./../fastas/" + plasmid + ".fasta")
				#os.system("CONDA_BASE=$(conda info --base)")
				#os.system("source $CONDA_BASE/etc/profile.d/conda.sh")
				#os.system("conda deactivate prokka_env")
				#p2k="python prokka2go.py --input ./plasmids/" + plasmid + "/gb_match/default/pk_results/" + plasmid + ".gbk --output ./plasmids/" + plasmid + "/gb_match/" + str(row['Matches']) + "/prokka2go.txt"
				#os.system(p2k)
				
				#NEW CODE FOR RUNNING PROKKA
				new_matches.append([row['Matches'], row['Genes']])
				os.system("/bin/bash prokka_def.sh " + outdir + " " + plasmid)

				#append prokka2go results to initial prokka annotation file
				with open("./../output/plasmids/" + plasmid + "/gb_match/default/pk_results" + "/" + plasmid + "-go.tsv", 'w') as cmptsv:
					cmpwrite=csv.writer(cmptsv, delimiter='\t')
					with open("./../output/plasmids/" + plasmid + "/gb_match/default/pk_results" + "/" + plasmid + ".tsv") as currtsv:
						tsvread=csv.reader(currtsv, delimiter='\t')
						with open("./../output/plasmids/" + plasmid + "/gb_match/default/prokka2go.txt") as kotsv:
							koread=csv.reader(kotsv, delimiter='\t')
							for currline in tsvread:
								for go_id in koread:
									new_row=[]
									if currline[0] in go_id[0]:
										print("COMPARE= " + currline[0] + " AND " + go_id[0])
										if len(go_id)>1:
											new_row.append(currline + [go_id[1]] + [go_id[2]] + [go_id[3]])
										else:
											new_row.append(currline + [""] + [""] + [""])
										cmpwrite.writerow(new_row[0])
										break								

				os.system("cp ./../output/plasmids/" + plasmid + "/gb_match/default/pk_results/" + plasmid + "-go.tsv ./plasmids/" + plasmid + "/merge/mergedgo.tsv")
			else:				
				gbk_dir="./../output/plasmids/" + plasmid + "/gb_match/" + str(row['Matches'])
				os.system("mkdir " + gbk_dir)
				out_dir=gbk_dir + "/pk_results"
				os.system("mkdir " + out_dir)

				#Run prokka
				#os.system("CONDA_BASE=$(conda info --base)")
				#os.system("source $CONDA_BASE/etc/profile.d/conda.sh")
				#os.system("conda init")
				#os.system("conda activate prokka_env")
				#os.system("prokka --kingdom Bacteria --outdir " + out_dir + "  --prefix " + plasmid + " --force --proteins " + gbk_dir + ".gbk  ./../fastas/" + plasmid + ".fasta")
				
				#os.system("CONDA_BASE=$(conda info --base)")
				#os.system("source $CONDA_BASE/etc/profile.d/conda.sh")
				#os.system("conda deactivate prokka_env")
				
				#p2k="python prokka2kegg.py --input ./plasmids/" + plasmid + "/gb_match/" + str(row['Matches']) + "/pk_results/" + plasmid + ".gbk --output ./plasmids/" + plasmid + "/gb_match/" + str(row['Matches']) + "/prokka2go.txt"
				#os.system(p2k)
				try:
					os.system("/bin/bash prokka.sh " + out_dir + " " + plasmid + " " + gbk_dir + " " + str(row['Matches']))
				except:
					print("Prokka not able to process given ID")
					os.system("rm " + gbk_dir + "/pk_results/*")
					os.system("rmdir " + gbk_dir + "/pk_results")
					os.system("rm " + gbk_dir + ".gbk")
					os.system(os.system("rmdir " + gbk_dir))

                        	#Only start comparing genes after reference accession ID has been processed by prokka
				ref_dir="./../output/plasmids/" + plasmid + "/gb_match/default/pk_results"
                                #keep track of total genes processed and conflicting matches between reference and current accession ID
				conflicts=0
				gene_num=0
				
				try:
                                	#open tsv files for reference + current acession IDs and compare genes line by line
					with open(ref_dir + "/" + plasmid + ".tsv") as ref_tsv:
						with open(out_dir + "/" + plasmid + ".tsv") as comp_tsv:
							for ref_line in ref_tsv:
								gene_num=gene_num+1
								ref_out=ref_line.replace('\t', ' ').split(' ')
								for comp_line in comp_tsv:
									comp_out=comp_line.replace('\t', ' ').split(' ')
									if len(ref_out)!=0 and len(comp_out)!=0 and ref_out[0] in comp_out[0]:
										if comp_out[3] not in ref_out[3] and (comp_out[3] not in '' and ref_out[3] not in ''):
											conflicts=conflicts+1
											break
				except IOError:
					print('Prokka failed to produce TSV file')
					os.system("rm " + gbk_dir + "/pk_results/*")
					os.system("rmdir " + gbk_dir + "/pk_results")
					os.system("rm " + gbk_dir + ".gbk")
					os.system("rmdir " + gbk_dir)	
											
                                #if a file has a conflict rate above 5%, delete its gbk file and prokka output
				if gene_num!=0:
					if conflicts/gene_num > 0.05:
						os.system("rm " + gbk_dir + "/pk_results/*")
						os.system("rmdir " + gbk_dir + "/pk_results")
						os.system("rmdir " + gbk_dir)
						print("Conflicts for " + str(row['Matches']) + " exceeded 5%\n")
					else:
						print("Conflicts for " + str(row['Matches']) + " DO NOT exceed 5%")
						print("Add results for " + str(row['Matches']) + " to plasmid-go.tsv\n")
						#os.system("cat ./merge/merged.tsv ./gb_match/" + str(row['Matches']) + "/pk_results/" + plasmid + ".tsv > ./merge/merged.tsv")
						new_matches.append([row['Matches'], row['Genes']])
						
						#Else, merge the current tsv file to the plasmid-go.tsv file 
						with open(out_dir + "/" + plasmid + "-go.tsv", 'w') as cmptsv:
							cmpwrite=csv.writer(cmptsv, delimiter='\t')
							with open(out_dir + "/" + plasmid + ".tsv") as currtsv:
								tsvread=csv.reader(currtsv, delimiter='\t')
								with open(gbk_dir + "/prokka2go.txt") as kotsv:
									koread=csv.reader(kotsv, delimiter='\t')
									for currline in tsvread:
										for go_id in koread:
											if currline[0] in go_id[0]:
												new_row=[]
												if len(go_id)>1:
													new_row.append(currline + [go_id[1]] + [go_id[2]] + [go_id[3]])
												else:
													new_row.append(currline + [""] + [""] + [""] + [""])
												cmpwrite.writerow(new_row[0])
												break

        #delete old match csv and create a new matches.csv file that contains more limited list of matches
	os.system("rm ./../output/plasmids/" + plasmid + "/matches.csv")
	print("Compiling new GBK matches for " + plasmid)
	with open("./../output/plasmids/" + plasmid + "/matches.csv", 'w') as csvfile:
		input=csv.writer(csvfile, delimiter=' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
		input.writerow(['Matches'] + ['Genes'])
		for i in range(len(new_matches)):
			input.writerow([new_matches[i][0]] + [new_matches[i][1]])

if __name__=="__main__":
        parser = argparse.ArgumentParser()
        parser.add_argument('plasmid', type=str)
        plasname=(parser.parse_args()).plasmid
        main(plasname)
