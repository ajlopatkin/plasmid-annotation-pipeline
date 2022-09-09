#blast_plasmids.py

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
	os.system("mkdir ./../output/plasmids/" + plasmid)

	# read in the fasta record
	record = SeqIO.read("./../fastas/" + plasmid + ".fasta", format="fasta")

	# blast the record against the nucleotide database
	print("Commencing Blast Search on " + plasmid)
	os.system("blastn -db nt -query ./../fastas/" + plasmid + ".fasta -out ./../output/plasmids/" + plasmid + "/" + plasmid + "-BLASTTEST.xml -outfmt 5 -max_target_seqs 1000 -remote")

	# copy the blasttest results to the plasmid's output folder
	os.system("cat ./../output/plasmids/" + plasmid + "/" + plasmid + "-BLASTTEST.xml")
	result_handle = NCBIWWW.qblast(program="blastn", database="nt", sequence=record.seq, auto_format="XML", hitlist_size=500, alignment_view=5)
	print("Finished Blast Search")

	# open the xml file
	result_handle=open("./../output/plasmids/" + plasmid + "/" + plasmid + "-BLASTTEST.xml")
	
	# write blast results to xml file
	blast_file = "./../output/plasmids/" + plasmid + "/" + plasmid + "-BLAST.xml"
	with open(blast_file, "w") as out_handle:
		out_handle.write(result_handle.read())
	out_handle.close()

if __name__=="__main__":
        parser = argparse.ArgumentParser()
        parser.add_argument('plasmid', type=str)
        plasname=(parser.parse_args()).plasmid
        main(plasname)
