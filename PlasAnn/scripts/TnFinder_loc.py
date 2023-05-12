'''
TnFinder_loc.py

DESCRIPTION: 

Get TnFinder transposon information
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
    os.system("mkdir ./../../output/plasmids/" + plasmid + "/TnFinder_out")

    if os.path.exists("./../../output/plasmids/" + plasmid + "/TnFinder_out/blastn"):
        os.remove("./../../output/plasmids/" + plasmid + "/TnFinder_out/blastn/" + plasmid + ".blastn")
        os.rmdir("./../../output/plasmids/" + plasmid + "/TnFinder_out/blastn/")

        if os.path.exists("./../../output/plasmids/" + plasmid + "/TnFinder_out/info.txt"):
            os.remove("./../../output/plasmids/" + plasmid + "/TnFinder_out/info.txt")
    
    # blast plasmid fasta against oric.fna
    print("Commencing TnFinder Search on " + plasmid + "")
    os.system("python ./TnComp_finder/TnComp_finder.py -f ./../fastas/" + plasmid + ".fasta -o ./../output/plasmids/" + plasmid + "/TnFinder_out")

if __name__=="__main__":
        parser = argparse.ArgumentParser()
        parser.add_argument('plasmid', type=str)
        plasname=(parser.parse_args()).plasmid
        main(plasname)