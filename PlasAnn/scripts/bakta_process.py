#Process tsv file produced by bakta to match locus tags from default tsv.
import os
import csv
import argparse

def main(plasmid):
    with open('./../../../pOXA_bakta/pOXA48K8.tsv') as tsvfile:
        input=csv.DictReader(tsvfile, delimiter='\t', quotechar='|')
        
        #start at negative one to account for header
        reference_num=-1

        with open('./../output/plasmids/pOXA48K8/final.tsv') as reference:
            ref_in=csv.DictReader(reference, delimiter='\t', quotechar='|')
            for row in ref_in:
                reference_num+=1

        #start at negative one to account for header
        totalrows=-1
        for row in input:
            totalrows+=1

        if totalrows!=reference_num:
            print("LOCUS NUM DONT MATCH!")
        else:
            print("LOCUS NUM MATCH")

        

if __name__=="__main__":
        parser = argparse.ArgumentParser()
        parser.add_argument('plasmid', type=str)
        plasname=(parser.parse_args()).plasmid
        main(plasname)