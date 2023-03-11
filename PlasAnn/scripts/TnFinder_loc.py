
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
    os.system("mkdir ./../../output/plasmids/" + plasmid + "/TnFinder_out")
    
    # blast plasmid fasta against oric.fna
    print("Commencing TnFinder Search on " + plasmid + "")
    #os.system("cd ./TnComp_finder")
    os.system("python ./TnComp_finder/TnComp_finder.py -f ./../fastas/" + plasmid + ".fasta -o ./../output/plasmids/" + plasmid + "/TnFinder_out")
    #os.system("cd ./..")

if __name__=="__main__":
        parser = argparse.ArgumentParser()
        parser.add_argument('plasmid', type=str)
        plasname=(parser.parse_args()).plasmid
        main(plasname)