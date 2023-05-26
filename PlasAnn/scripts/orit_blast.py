#blast_orit.py
'''
orit_blast.py

DESCRIPTION: get oriT match
'''
from Bio import SearchIO
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from Bio import Entrez
from io import StringIO
import csv
import argparse
import re
import os

def main(plasmid):
    # create a directory for the output related to each plasmid
    os.system("mkdir ./../output/plasmids/" + plasmid + "/oriT_out")
    
    # blast plasmid fasta against oric.fna
    print("Commencing oriT Blast Search on " + plasmid)
    os.system("blastn -query ./../fastas/" + plasmid + ".fasta -subject ./../bakta_database/db/orit.fna -out ./../output/plasmids/" + plasmid + "/oriT_out/orit.csv -outfmt 10")

if __name__=="__main__":
        parser = argparse.ArgumentParser()
        parser.add_argument('plasmid', type=str)
        plasname=(parser.parse_args()).plasmid
        main(plasname)