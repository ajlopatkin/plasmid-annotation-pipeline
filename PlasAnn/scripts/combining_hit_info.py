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

#get gene and indicate what choosen gene name it should be replaced with
def appendPlasmidDB(inc, choosen_gene, sequence, gene_info, consider):
    fasta=open("./../ref_dbs/" + inc + ".fasta", 'a+')
    
    #process choosen gene in case it's an orf__ gene to just say orf
    if (choosen_gene.lower()).startswith('orf'):
        choosen_gene='orf'

    #create record for each gene in consider (each gene close in apperance to choosen gene) to identify prefered naming for these genes
    for gene in consider:
        fasta_format=">" + gene + ":" + choosen_gene + "~~~" + gene_info[0] + "||" + gene_info[1] + "\n" + sequence + '\n'
        fasta.write(fasta_format)
    fasta.close()

#return prefered gene name
def findNaming(inc, gene):
    new_name=''
    if os.path.exists("./../ref_dbs/" + inc + ".fasta"):
        fasta=open("./../ref_dbs/" + inc + ".fasta", 'r')
    
        #if gene naming found for given gene, obtain new gene name and information
        for line in fasta:
            if '>' in line:
                if gene in line[line.index('>')+1:line.index(':')]:
                    new_name=line[line.index(':')+1:line.index('~')]
                    prod=line[line.index('~~~')+3:line.index('||')]
                    inf=line[line.index('||')+2:line.index('\n')]
                    return [new_name, [prod, inf]]
    
    #if no naing convention found for given gene, return empty result
    result=[new_name]
    return result

#ask user to break ties; return user gene name choice
def namingAsk(inc, consider, locus):
    choosen_gene=consider[0]
    print("\nSELECT GENE FOR " + locus + " (INC: " + inc + ")")
    for i in range(0, len(consider)):
        print("Option " + str(i+1) + ": " + str(consider[i]))

    while True:
        try:
            option=int(input("\nPlease indicate which gene to use: "))
            if option<1 or option>len(consider)+1:
                print("Enter valid integer input corresponding to the gene option only!")
            else:
                choosen_gene=consider[option-1]
                break
        except:
            print("Enter valid integer input only!")

    return choosen_gene

#NOT CURRENTLY BEING USED
#Get gene sequence from labelled_genes.fasta for given gene
def getGeneSequence(gene, other_genes, plasmid):
    sequence=""
    protein_id=""
    product=""

    sequence=[]
    with open("./../output/plasmids/"+plasmid+"/genbank_genes/labelled_genes.fasta") as fasta:
        write=0
        for line in fasta:
            if '>' in line:
                if write==0 and gene.lower() in line.lower() and 'non-experimental' not in line:
                    write=1
                elif write==1:
                    break
            
            if write==1:
                sequence.append(line)
    return sequence

#def getPercentId(gene_file, plasmid):
#    os.system("tblastn -query ./../fastas/" + plasmid + ".fasta -subject " + gene_file + " -out ./../output/plasmids/" + plasmid + "/temp_gene_file.csv -outfmt 10")

#return inc group for current plasmid
def getINC(plasmid):
    with open("./../output/plasmids/"+plasmid+"/abr_genes.txt") as f:
	    lines = f.readlines()
	    x = lines
	    inc_group = x[0]
	    inc_group = str(inc_group)
	    inc = inc_group[0:4]
	    f.close()
    
    return inc

def main(plasmid, comb_run):
    INC_GROUP=getINC(plasmid)

    newpath = "./../output/plasmids/"+plasmid+"/gb_match"

    if not os.path.exists("./../output/plasmids/"+plasmid+"/genbank_genes"):
        os.makedirs("./../output/plasmids/"+plasmid+"/genbank_genes")
    the_path = "./../output/plasmids/"+plasmid+"/genbank_genes/"

    if not os.path.exists("./../ref_dbs"):
        os.makedirs("./../ref_dbs")

    accepted_match=[]
    with open("./../output/plasmids/" + plasmid + "/matches.csv") as csvfile:
        matchreader = csv.DictReader(csvfile, delimiter=' ', quotechar='|')
        for match in matchreader:
            if os.path.exists("./../output/plasmids/"+plasmid+"/gb_match/" + match['Matches']):
                accepted_match.append(match['Matches'])

    print(accepted_match)
    genbanks = glob("%s/*.gbk" % newpath)
    all_genes = []
    for file in genbanks:
        #fileid=(file.split('/')[-1])[0:(file.split('/')[-1]).index('.')]
        #if fileid not in accepted_match:
        #    continue

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

   #  with open("./../output/plasmids/"+plasmid+"/abr_genes.txt") as f:
#         lines = f.readlines()
#         x = lines
#         inc_group = x[0]
#         inc_group = str(inc_group)
#         inc_group = inc_group[0:4]
# 
#         f.close()
# 
#     for record in SeqIO.parse("./../backbone_genes/" + inc_group + "_backbone.fasta", "fasta"):
#         the_id = (record.id)
#         the_sequence = record.seq
#         new_rec = SeqRecord(seq = the_sequence, id = the_id)
#         all_genes.append(new_rec)

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
            max_gene = max(all_genes_found, key=all_genes_found.get)
            information = gene_info.get(max_gene)

            #if comb_run==1, ask for user input and look at naming convention file
            if int(comb_run)==1:
                result=findNaming(INC_GROUP, max_gene)
                if result[0] not in '':
                    max_gene=result[0]
                    information=result[1]

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
            max_gene = max(all_genes_found, key=all_genes_found.get)

            the_best = {}
            if len(genes) > 1: 
                # for each of our top genes
                for x in genes: 
                    total_rec=0
                    total = 0
                    #temp label
                    y = str(x) 
                    #empty dictionary
                    y = {}
                    # both holds (gene, organism ) = number of times this occurs, 
                    for key in both: 
                        # if the gene we are looking at matches the key 
                        val = both[key]
                        total_rec+=val
                        if x == key[0]: 
                            y[key[1]] = val
                            total+=val

                    #total = 0

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

                    #ratio=total/total_rec
                    #the_best[x]=ratio
            
                #obtain max_gene, or most commonly seen gene at current locus
                max_gene = max(the_best, key =the_best.get)
                information = gene_info.get(max_gene)

                #get list of genes that have a incidence rate that is within 0.2 of max gene
                max_val=the_best[max_gene]
                consider=[max_gene]
                for gene in the_best:
                    if gene!=max_gene:
                        if the_best[gene]==max_val or abs(max_val-the_best[gene])<=0.2:
                            consider.append(gene)

                #process consider genes to determine if different naming should be used
                if len(consider)>1:
                    #if comb_run==1, ask for user input and look at naming convention file
                    if int(comb_run)==1:
                        result=findNaming(INC_GROUP, max_gene)

                        if result[0] in '':
                            max_gene=namingAsk(INC_GROUP, consider, locus_tag)
                            appendPlasmidDB(INC_GROUP, max_gene, seq, gene_info.get(max_gene), consider)
                            information = gene_info.get(max_gene)
                        else:
                            max_gene=result[0]
                            information=result[1]
                else:
                    #if comb_run==1, ask for user input and look at naming convention file
                    if int(comb_run)==1:
                        result=findNaming(INC_GROUP, max_gene)

                        if result[0] not in '':
                            max_gene=result[0]
                            information=result[1]
                            
                prod = (information[0])
                inf = (information[1])
                with open("./../output/plasmids/"+plasmid+"/blast_genes.tsv", 'a') as out_file: 
                    tsv_writer = csv.writer(out_file, delimiter = '\t')
                    tsv_writer.writerow([locus_tag, max_gene, prod, inf])
                        

                    #os.system("mkdir ./../output/plasmids/" + plasmid + "/temp")
                    #for gene in consider:
                    #    fasta_sequ=getGeneSequence(gene, consider, plasmid)

                    #    f=open("./../output/plasmids/" + plasmid + "/temp/" + plasmid + "_" + gene + "_" + locus_tag + ".fasta", 'w')
                    #    for line in fasta_sequ:
                    #        f.write(line)
                    #    f.close()

                #pct_ids=[]
                #for gene in consider:
                #    pct_ids.append(getPercentId("./../output/plasmids/" + plasmid + "/temp/" + plasmid + "_" + gene + "_" + locus_tag + ".fasta", plasmid))
                

                #max_gene=max(genes, key=genes.get)

            elif len(genes) == 1: 
                max_gene = genes[0]
                information = gene_info.get(max_gene)
                
                if int(comb_run)==1:
                    result=findNaming(INC_GROUP, max_gene)

                    if result[0] not in '':
                        max_gene=result[0]
                        information=result[1]


                prod = (information[0])
                inf = (information[1])

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

    #add gene data to tsv and gbk files
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
                    feature.qualifiers = new_qual

            final_features.append(feature)
        record.features = final_features
        with open(outfile, "w") as new_gb:
            SeqIO.write(record, new_gb, "genbank")

        final_data.to_csv("./../output/plasmids/" + plasmid + "/the_final.tsv", sep="\t",  index=False)

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', metavar='plasmid', dest='p',
                        type=str, required=True)
    parser.add_argument('-c', metavar='comb_run', dest='c',
                        type=str, required=True)

    args = parser.parse_args()
    main(args.p, args.c)

