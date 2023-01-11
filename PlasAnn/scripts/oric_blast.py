#oric_blast.py

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
    os.system("mkdir ./../output/plasmids/" + plasmid + "/oriC_out")
    
    # blast plasmid fasta against oric.fna
    print("Commencing oriC Blast Search on " + plasmid)
    os.system("blastn -query ./../fastas/" + plasmid + ".fasta -subject ./../bakta_database/db/oric.fna -out ./../output/plasmids/" + plasmid + "/oriC_out/oric.csv -outfmt 10")

if __name__=="__main__":
        parser = argparse.ArgumentParser()
        parser.add_argument('plasmid', type=str)
        plasname=(parser.parse_args()).plasmid
        main(plasname)