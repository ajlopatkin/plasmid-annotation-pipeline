# combining_hit_accession_information
from glob import glob
import os
import csv
import sys
import argparse
import shutil
import collections
from collections import OrderedDict
from Bio import SeqIO, GenBank, SeqFeature
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
import pandas as pd

def main(plasmid):

    newpath = "./../output/plasmids/"+plasmid+"/gb_match"

    if not os.path.exists("./../output/plasmids/"+plasmid+"/genbank_genes"):
        os.makedirs("./../output/plasmids/"+plasmid+"/genbank_genes")
    the_path = "./../output/plasmids/"+plasmid+"/genbank_genes/"



    genbanks = glob("%s/*.gbk" % newpath)
    all_genes = []
    for file in genbanks:
        for i in SeqIO.parse(file, "gb"):
            print("We are parsing: " + i.id)
            s = (i.description)
            organism = (" ".join(s.split()[:2]))
            for f in i.features:
                if (f.type == "CDS" and ("gene" in f.qualifiers)):
                    gene = f.qualifiers.get("gene", [])
                    gene = gene[0]
                    product = ""
                    inf = ""

                    if "product" in f.qualifiers: 
                        product = f.qualifiers.get("product", [])
                        product = product[0]

                    if "inference" in f.qualifiers: 
                        inf = f.qualifiers.get("inference", [])
                        inf = inf[0]

                    the_id = str(i.id + "_" + gene +  " " + organism )
                    desc = product + "|" + inf 
                    the_sequence = str(f.qualifiers.get("translation", []))
                    the_sequence = (the_sequence.strip("[']"))

                    curr_rec = SeqRecord(seq = Seq(the_sequence), id = the_id, description = desc)
                    all_genes.append(curr_rec)

    with open("./../output/plasmids/"+plasmid+"/abr_genes.txt") as f:
        lines = f.readlines()
        x = lines
        inc_group = x[0]
        inc_group = str(inc_group)
        inc_group = inc_group[0:4]

        f.close()

    for record in SeqIO.parse("./../backbone_genes/" + inc_group + "_backbone.fasta", "fasta"):
        the_id = (record.id)
        the_sequence = record.seq
        new_rec = SeqRecord(seq = the_sequence, id = the_id)
        all_genes.append(new_rec)

    SeqIO.write(all_genes, "./../output/plasmids/"+plasmid+"/genbank_genes/labelled_genes.fasta", "fasta")

    if not os.path.exists("./../output/plasmids/"+plasmid+"/blast_accession_results"):
        os.makedirs("./../output/plasmids/"+plasmid+"/blast_accession_results")

    with open("./../output/plasmids/"+plasmid+"/blast_genes.tsv", 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(['locus_tag', 'gene', 'product', 'inference'])

    our_file = "./../output/plasmids/"+plasmid+"/gb_match/default/pk_results/"+plasmid+".gbk"
    our_genes = []
    for i in SeqIO.parse(our_file, "gb"):
        for f in i.features:
                if (f.type == "CDS" and ("locus_tag" in f.qualifiers)):
                    locus_tag = str(f.qualifiers["locus_tag"][0])
                    the_sequence = str(f.qualifiers.get("translation", []))
                    the_sequence = (the_sequence.strip("[']"))
                    curr_rec = SeqRecord(seq = Seq(the_sequence), id = locus_tag)
                    our_genes.append(curr_rec)
    SeqIO.write(our_genes, "./../output/plasmids/"+plasmid+"/our_genes.fasta", "fasta")

    for rec in SeqIO.parse("./../output/plasmids/"+plasmid+"/our_genes.fasta", "fasta"):

        locus_tag = (rec.id)
        seq = str(rec.seq)
        entry = "> " + locus_tag + "\n" + seq + "\n"
        f1 = open("hyp_gene.fasta", "w")
        f1.write(entry)
        f1.close()

        if locus_tag == 'NMEEHBII_00034': 
            original = "hyp_gene.fasta"
            target = "maybe_xerD.fasta"
            shutil.copyfile(original, target)

        f2 = open("hyp_gene.fasta", "r")
        blastp_cline = NcbiblastpCommandline(query= "hyp_gene.fasta",subject="./../output/plasmids/"+plasmid+"/genbank_genes/labelled_genes.fasta",
                                         evalue=0.001, outfmt=5, out="./../output/plasmids/"+plasmid+"/blast_accession_results/" + locus_tag + ".xml")
        stdout, stderr = blastp_cline()
        f2.close()
        result_handle = open("./../output/plasmids/"+plasmid+"/blast_accession_results/" + locus_tag + ".xml")
        blast_records = NCBIXML.parse(result_handle)

        all_genes_found = {}
        gene_info = {}
        organism_info = {}
        both = {}
        hits = 0
        max_gene = ""
        for blast_record in blast_records: 

            for alignment in blast_record.alignments: 
                e_val = alignment.hsps[0].expect
                pct_id =  alignment.hsps[0].identities/alignment.hsps[0].align_length*100

                the_desc = alignment.hit_def
                gene_ind = the_desc.find("_")
                prod_ind = the_desc.find(" ")
                inf_ind = the_desc.find("|")
                end_ind = len(the_desc)

                
                if e_val == 0 and pct_id == 100:
                    hits = hits + 1

                    the_gene = (the_desc[gene_ind+1:prod_ind])
                    org = (" ".join(the_desc.split()[1:3]))
                    if "backbone" in the_gene: 
                        the_gene = the_gene.replace("_backbone", "")
                        max_gene = the_gene
                        #print(max_gene)
                        break 
                    if "_" in the_gene:
                        ind = the_gene.index("_")
                        the_gene = the_gene[0:ind]

                    product = (the_desc[prod_ind+1:inf_ind]).replace(org, "")
                    inf_res = (the_desc[inf_ind+1:end_ind])

                    gene_org = (the_gene, org)

                    gene_info[the_gene] = (product, inf_res)
                    all_genes_found[the_gene] = all_genes_found.get(the_gene,0) + 1
                    organism_info[org] = organism_info.get(org,0) + 1
                    both[gene_org] = both.get(gene_org,0) + 1

                elif e_val < 1e-3 and pct_id > 99:

                    hits = hits + 1

                    the_gene = (the_desc[gene_ind+1:prod_ind])

                    org = (" ".join(the_desc.split()[1:3]))
                    if "backbone" in the_gene: 
                        the_gene = the_gene.replace("_backbone", "")
                        max_gene = the_gene
                        print(max_gene)
                        break 
                    if "_" in the_gene:
                        ind = the_gene.index("_")
                        the_gene = the_gene[0:ind]

                    product = (the_desc[prod_ind+1:inf_ind]).replace(org, "")
                    inf_res = (the_desc[inf_ind+1:end_ind])

                    if product is None: 
                        product = ""
                    if inf_res is None: 
                        inf_res = ""

                    gene_info[the_gene] = (product, inf_res)
                    gene_org = (the_gene, org)
                   
                    all_genes_found[the_gene] = all_genes_found.get(the_gene,0) + 1
                    organism_info[org] = organism_info.get(org,0) + 1
                    both[gene_org] = both.get(gene_org,0) + 1
                    #print(both)

        # Now we have 4 dictionaries: 
        # all_genes_found, gene_info, organism_info, both 
        if len(max_gene) > 0 :
            prod = ""
            inf =""
            with open("./../output/plasmids/"+plasmid+"/blast_genes.tsv", 'a') as out_file: 
                tsv_writer = csv.writer(out_file, delimiter = '\t')
                tsv_writer.writerow([locus_tag, max_gene, prod, inf])

        elif len(all_genes_found) == 1:
            #print(all_genes_found)
            #print(organism_info)
            max_gene = max(all_genes_found, key=all_genes_found.get)
            # print(max_gene)
            # print(gene_info)
            information = gene_info.get(max_gene)
            #print(information)
            prod = (information[0])
            inf = (information[1])
            with open("./../output/plasmids/"+plasmid+"/blast_genes.tsv", 'a') as out_file: 
                tsv_writer = csv.writer(out_file, delimiter = '\t')
                tsv_writer.writerow([locus_tag, max_gene, prod, inf])

        elif len(all_genes_found) > 1: 
            genes = []
            for key in all_genes_found: 
                val = all_genes_found[key]
                val = val/hits
                if val > 0.25: 
                    genes.append(key)

            max_org = max(organism_info, key=organism_info.get)

            the_best = {}
            if len(genes) > 1: 
                # for each of our top genes
                for x in genes: 
                    #temp label
                    y = str(x) 
                    #empty dictionary
                    y = {}
                    # both holds (gene, organism ) = number of times this occurs, 
                    for key in both: 
                        # if the gene we are looking at matches the key 
                        if x == key[0]: 
                            val = both[key]
                            y[key[1]] = val

                    total = 0

                    orgs = []
                    for key in y: 
                        total = total + y[key]
                        orgs.append(key)
                    
                    if max_org in orgs: 
                        popular = y[max_org]
                    else: 
                        max_org = max(y, key =y.get)
                        popular = y[max_org]
                        
                    ratio = popular/total
                    the_best[x] = ratio
                
                max_gene = max(the_best, key =the_best.get)

            elif len(genes) == 1: 
                max_gene = genes[0]
                #print(max_gene)
                #print(y)
            #print(the_best)
            #print(all_genes_found)
            #print(organism_info)
            #max_gene = max(the_best, key =the_best.get)


            #max_gene = max_gene[0]
            #new_key = (max_gene, max_org)
            #print(all_genes_found[new_key])
            # if len(max_gene) > 0: 
            # print(max_gene)
             #print(gene_info)
            information = gene_info.get(max_gene)
                # print(information)
             #print(information)
            prod = (information[0])
                # #print(prod)
            inf = (information[1])
            # else: 
            #     max_gene = ""
            #     inf = ""
            #     prod  = ""

            with open("./../output/plasmids/"+plasmid+"/blast_genes.tsv", 'a') as out_file: 
                tsv_writer = csv.writer(out_file, delimiter = '\t')
                tsv_writer.writerow([locus_tag, max_gene, prod, inf])

        else: 
            max_gene = ""
            prod  = ""
            inf = ""

            with open("./../output/plasmids/"+plasmid+"/blast_genes.tsv", 'a') as out_file: 
                tsv_writer = csv.writer(out_file, delimiter = '\t')
                tsv_writer.writerow([locus_tag, max_gene, prod, inf])

    b_data = pd.read_table("./../output/plasmids/"+plasmid+"/blast_genes.tsv")
    locuses = b_data['locus_tag'].values.tolist()
    genes = b_data['gene'].values.tolist()
    product = b_data['product'].values.tolist()
    inference = b_data['inference'].values.tolist()

    original = "./../output/plasmids/"+plasmid+"/gb_match/default/pk_results/"+plasmid+".gbk"
    target = "./../output/plasmids/"+plasmid+"/final.gbk"
    shutil.copyfile(original, target)

    final_data = pd.read_table("./../output/plasmids/"+plasmid+"/gb_match/default/pk_results/"+plasmid+".tsv")
    final_locuses = final_data['locus_tag'].values.tolist()
    # print("The original")
    # print(final_data.to_string())


    final_features = []
    outfile = "./../output/plasmids/"+plasmid+"/the_final.gbk"

    for record in SeqIO.parse("./../output/plasmids/"+plasmid+"/final.gbk", "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                new_qual = {}
                found_locus = 0
                our_locus = feature.qualifiers["locus_tag"][0]

                for locus in locuses:
                    if our_locus == locus:
                        ind = locuses.index(locus)
                        newgene = genes[ind]
                        newprod = product[ind]
                        newinf = inference[ind]

                        fin_ind = final_locuses.index(locus)
                        #print(fin_ind)
                        final_data.loc[fin_ind,'gene']=newgene
                        # final_data.insert(fin_ind,'gene', newgene)
                        # final_data.insert(fin_ind,'EC_number', '')
                        # final_data.insert(fin_ind,'COG', '')


                        final_data.loc[fin_ind, 'product'] = newprod
                        final_data.loc[fin_ind, 'EC_number'] = ''
                        final_data.loc[fin_ind, 'COG'] = ''
                        # feature.qualifiers = SeqFeature.qualifiers()
                        found_locus = 1
                if found_locus == 1:
                    for x in feature.qualifiers:
                        new_qual[x] = feature.qualifiers[x]

                    new_qual['gene'] = [f'{newgene}']
                    new_qual['product'] = [f'{newprod}']
                    new_qual['inference'] = [f'{newinf}']


                    #print(new_qual)
                    feature.qualifiers = new_qual
            final_features.append(feature)
        record.features = final_features
        with open(outfile, "w") as new_gb:
            SeqIO.write(record, new_gb, "genbank")

        final_data.to_csv("./../output/plasmids/" + plasmid + "/the_final.tsv", sep="\t",  index=False)
    
if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('plasmid', type=str)
    plasname=(parser.parse_args()).plasmid
    main(plasname)

