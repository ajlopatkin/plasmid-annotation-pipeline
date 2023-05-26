'''
DESCRIPTION: Update final tsv by adding COG and GO functional information
'''
import csv
import pandas as pd
import argparse

def main(plasmid): 
	#Update table to add GO ID's
	update = pd.read_table("./../output/plasmids/"+plasmid+"/prokka2go.txt")
	#'GO_ID_P','GO_ID_C', 'GO_ID_F'
	update.columns = ['locus_tag', 'GO_ID_P', 'GO_ID_C', 'GO_ID_F', 'Function']
	final_tsv = pd.read_table("./../output/plasmids/"+plasmid+"/the_final.tsv", index_col=False)
	final_ver_go=final_tsv.merge(update.set_index('locus_tag'), how='left', on='locus_tag')

	#Update table to add COG descriptions
	update = pd.read_csv("./../output/plasmids/"+plasmid+"/COGClass/classifier_result.tsv", usecols=['QUERY_ID', 'COG_DESCRIPTION'], sep='\t')
	update = update.rename(columns={'QUERY_ID': 'locus_tag', 'COG_DESCRIPTION': 'COG_FUNC'})
	final_ver_cog=final_ver_go.merge(update.set_index('locus_tag'), how='left', on='locus_tag')

	#Update table to add COG ID's
	update_COG = pd.read_csv("./../output/plasmids/"+plasmid+"/COGClass/classifier_result.tsv", usecols=['QUERY_ID', 'COG_ID'], sep='\t')
	update_COG = update_COG.rename(columns={'QUERY_ID': 'locus_tag', 'COG_ID': 'COG'})
	final_ver_cog['COG'] = final_ver_cog.apply(lambda x: x['COG'], axis=1) 
	for index, row in update_COG.iterrows():
		if row['locus_tag'] not in 'locus_tag':
			final_ver_cog.loc[final_ver_cog.locus_tag == row['locus_tag'], 'COG']=row['COG']

	#Save new DF
	final_ver_cog = final_ver_cog.rename(columns={'Function': 'GO_FUNC'})
	final_ver_cog.to_csv("./../output/plasmids/" + plasmid + "/the_final.tsv", sep="\t", index=False)

if __name__=="__main__":
        parser = argparse.ArgumentParser()
        parser.add_argument('plasmid', type=str)
        plasname=(parser.parse_args()).plasmid
        main(plasname)
