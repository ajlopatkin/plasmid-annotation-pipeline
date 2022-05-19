#!/usr/bin/python
'''
DESCRIPTION: GBK_UPDATE creates a final.gbk file utilizing the default.gbk
created by Prokka and updating the gene information (ie. name, EC_Number, product) based
on the final annotation in final.tsv
'''

from Bio.Seq import Seq
from Bio import SearchIO
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import os
import csv
import argparse


def main(plasmid):
	#Sort through the final annotation and place it in gene list to parse through
	#each locus in the plasmid
    genes=[]
    with open("./../output/plasmids/" + plasmid + "/final.tsv") as geneiter:
        iterate=csv.reader(geneiter, delimiter='\t')
        for line in iterate:
            if "locus_tag" not in line:
                genes.append(line)

	#At every locus, update the gene information with the gene that correspond to that same
	#locus in the final annotation
    with open("./../output/plasmids/" + plasmid + "/gb_match/default/pk_results/" + plasmid + ".gbk") as gbkiter:
        gbk=csv.reader(gbkiter, delimiter='\t')
        locus=0
        gene_seen=0
        with open("./../output/plasmids/" + plasmid + "/final.gbk", 'w') as finaliter:
            final=csv.writer(finaliter, delimiter='\t', quoting=csv.QUOTE_NONE, escapechar='"')
            for line in gbk:
                print(line)
                if "/gene" in line[0]:
                    gene_seen=1
                elif "/locus_tag" in line[0]:
                    if len(genes[locus])>0:
                        final.writerow(['                     /gene=' + '"' + genes[locus][3] + '"'])
                    final.writerow([str(line[0])])
                    if len(genes[locus][4])>0:
                        final.writerow(['                     /EC_number=' + '"' + genes[locus][4] + '"'])
                    if len(genes[locus][6])>0:
                        final.writerow(['                     /product=' + '"' + genes[locus][6] + '"'])
                    print(locus)
                    print(len(genes))
                    locus+=1
                else:
                    if "/product" not in line[0] and "/EC_number" not in line[0]:
                        final.writerow([str(line[0])])


if __name__=="__main__":
        parser = argparse.ArgumentParser()
        parser.add_argument('plasmid', type=str)
        plasname=(parser.parse_args()).plasmid
        main(plasname)
