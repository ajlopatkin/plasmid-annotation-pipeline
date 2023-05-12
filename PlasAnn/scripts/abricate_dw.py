#!/usr/bin/python
'''
abricate_dw.py

DESCRIPTION: Run abricate w/ncbi, card, resfinder, and plasmidfinder databases for each accession ID match
'''
import csv
import os
import argparse

def main(plasmid):
		directory="./../output/plasmids/" + plasmid + "/ab_results"
		os.system("mkdir " + directory)
		databases=['ncbi', 'card', 'resfinder', 'plasmidfinder']
		for base in databases:
			os.system("abricate --db " + base + " ./../fastas/" + plasmid + ".fasta > " + directory + "/" + base + ".txt")
			os.system("abricate --summary " + directory + "/" + base + ".txt > " + directory + "/summary_" + base + ".txt")

if __name__=="__main__":
        parser = argparse.ArgumentParser()
        parser.add_argument('plasmid', type=str)
        plasname=(parser.parse_args()).plasmid
        main(plasname)
