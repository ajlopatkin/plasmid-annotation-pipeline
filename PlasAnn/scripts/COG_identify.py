from cogclassifier import cogclassifier
import argparse


def main(plasmid):
    #FROM COG_CLASSIFIER GIT HUB SAMPLE API
    query_fasta_file = "../output/plasmids/" + plasmid + "/gb_match/default/pk_results/" + plasmid + ".faa"
    outdir = "../output/plasmids/" + plasmid + "/COGClass"
    cogclassifier.run(query_fasta_file, outdir)


if __name__=="__main__":
        parser = argparse.ArgumentParser()
        parser.add_argument('plasmid', type=str)
        plasname=(parser.parse_args()).plasmid
        main(plasname)