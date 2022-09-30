#!/usr/bin/python
'''
DESCRIPTION:

1) Obtains the default bakta annotation on raw fasta

'''

import os
import csv
import argparse

def main(plasmid):

	# making the directory for our bakta default outputs
	os.system("mkdir ./../output/plasmids/" + plasmid + "/gb_match/default")	
	os.system("mkdir ./../output/plasmids/" + plasmid + "/gb_match/default/bakta_results")
	# setting our output directory	
	outdir="./../output/plasmids/" + plasmid + "/gb_match/default/bakta_results"
	
	# using os.system to run prokka on our plasmid
	os.system("/bin/bash bakta.sh " + outdir + " " + plasmid)

if __name__=="__main__":
        parser = argparse.ArgumentParser()
        parser.add_argument('plasmid', type=str)
        plasname=(parser.parse_args()).plasmid
        main(plasname)