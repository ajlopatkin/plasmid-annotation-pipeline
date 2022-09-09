#update_tsv.py

import csv
import pandas as pd
import argparse

def main(plasmid): 
	update = pd.read_table("./../output/plasmids/"+plasmid+"/prokka2go.txt")
	update.columns = ['locus_tag', 'GO ID','Type','Function']
	final_tsv = pd.read_table("./../output/plasmids/"+plasmid+"/the_final.tsv", index_col=False)
	final_tsv['GO ID'] = update['GO ID']
	final_tsv['Type'] = update['Type']
	final_tsv['Function'] = update['Function']
	print(final_tsv.head())
	final_tsv.to_csv("./../output/plasmids/" + plasmid + "/the_final.tsv", sep="\t", index=False)


	print(update.head())



if __name__=="__main__":
        parser = argparse.ArgumentParser()
        parser.add_argument('plasmid', type=str)
        plasname=(parser.parse_args()).plasmid
        main(plasname)
