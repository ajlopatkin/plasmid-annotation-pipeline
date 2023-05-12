'''
final_edits.py

DESCRIPTION: 

Apply final edits to TSV to include oriC, oriT, and transposon information
'''
import csv
import argparse
import re
import os
import pandas as pd

#remove if needed
import shutil

#TO BE IMPLEMENTED FURTHER: getDeltaGenes should determine if gene has been cut or shortened; add 'delta' to gene name if so
def getDeltaGenes(plasmid, lengths, cut=0.7):
    delta=[]
    locus=1
    with open("./../output/plasmids/"+plasmid+"/the_final.tsv") as ref_tsv:
        tsvread=csv.reader(ref_tsv, delimiter='\t')
        for gene in tsvread:
            if "locus_tag" not in gene:
                print("LOCUS " + str(locus))
                if int(gene[2])<(int(lengths[locus-1][1])-int(lengths[locus-1][0]))*cut:
                    delta.append(gene[0])
                locus+=1
    print(delta)
    return delta

def getGBKTranslations(plasmid):
    arr = []
    
    with open('../output/plasmids/' + plasmid + '/the_final.gbk') as input:
        trans = 0
        sequence=""
        for line in input.readlines():
            print(sequence)
            if "ORIGIN" in line:
                break
            if line.startswith(' ' * 21 + '/translation='):
                trans=1
                line=line[int(line.index('"')): int(line.index('\n'))]
                sequence=((line.replace(" ", '')).replace('"', ''))
            
            elif trans==1 and '/gene' not in line:
                edited=((line.replace(" ", '')).replace('"', '')).replace('\n', '')
                sequence=sequence + edited
            elif trans==1 and '/gene' in line:
                arr.append(sequence)
                trans=0
                sequence=""
    return arr

#check overlap between gene and oriC or oriT site to determine how to order the gene and oriC/oriT sites in the tsv file
def checkOverlapOriX(gene, oriX, prev_gene):
    if (int(oriX[1])<=int(gene[0]) and int(gene[0])-int(oriX[1])>10) or (int(oriX[0])>=int(gene[1]) and int(oriX[0])-int(gene[1])>10):
        return 0
    return 1

#check overlap between gene and Transposons to determine how to order the gene and Transposons sites in the tsv file
def checkOverlapTransp(gene, transp):
    if (int(transp[1])<=int(gene[0]) and int(gene[0])-int(transp[1])>0) or (int(transp[0])>=int(gene[1]) and int(transp[0])-int(gene[1])>0):
        return 0
    return 1

def checkEnd(gene, transp_end):
    if int(transp_end)<=int(gene[0]) and (int(gene[0])-int(transp_end)>200):
        return 0
    return 1

#Get transposons from TNFinder output
def getTransposons(plasmid):
    try:
        tnfinder_out = pd.read_table("./../output/plasmids/" + plasmid + "/TnFinder_out/blastn/" + plasmid + ".blastn", header=None)
        
        tnfinder_out.columns = ["Accession_ID", "Start_Query", "Stop_Query", "Subject_ID", "Length_Subject", "Start_Sub", "Stop_Sub", "Percent_Identity", "Matching_Bps"]
        tnfinder_out = tnfinder_out.drop(columns = ["Accession_ID"])
        tnfinder_out_filtered = tnfinder_out[tnfinder_out["Percent_Identity"] >= 95]

        tnfinder_out_filtered = tnfinder_out_filtered.reset_index(drop = True)
        counter = 0
        for i in tnfinder_out_filtered["Subject_ID"]:
            i = i.split('|')
            tnfinder_out_filtered.at[counter, "Subject_ID"] = i[1] + "|" + i[2] 
            counter = counter + 1

        tnfinder_out_filtered = tnfinder_out_filtered.sort_values(by = ['Start_Query'])
        tnfinder_out_filtered = tnfinder_out_filtered.reset_index(drop = True)

        # creating dictionary to store all MGE

        match = {}
        matches = []
        seq_match = 0
        match_number = 1
        match_name = "MGE_"

        for index, row in tnfinder_out_filtered.iterrows():
            if seq_match == 0: 
                seq_match = row["Stop_Query"]
                matches.append(row["Start_Query"])
                matches.append(row["Stop_Query"])


            else: 
                compare = row["Start_Query"]
                if (compare - seq_match) <= 2500: 
                    seq_match = row["Stop_Query"]
                    matches.append(seq_match)

                else:
                    key = match_name + str(match_number)
                    match[key] = matches
                    matches = []
                    matches.append(compare)
                    seq_match = row["Stop_Query"]
                    matches.append(seq_match)
                    match_number += 1

        key = match_name + str(match_number)
        match[key] = matches 

        print('***************************************************************')
        print("Here is the number of MGEs we're isolating + their starts + stops")
        print(match)


        # Create new df for found features 
        mge = pd.DataFrame(columns = ['Mobile Genetic Element','Start', 'Stop'])

        for lst in match.copy(): 
            if len(match[lst]) == 2: 
                name = tnfinder_out_filtered.loc[tnfinder_out_filtered['Start_Query'] == match[lst][0], 'Subject_ID'].iloc[0]
                mge.loc[len(mge.index)] = [name, match[lst][0], match[lst][1]] 
                del match[lst]

        ## Now we have our clustered sites
        named_mge = {}
        names = []
        for i in match: 
            last = len(match[i])
        
            for j in match[i]:

                if j == match[i][0]:
                    name_1 = tnfinder_out_filtered.loc[tnfinder_out_filtered['Start_Query'] == j, 'Subject_ID'].iloc[0]
                    names.append(j)
                    named_mge[name_1] = names
                else:
                    named_2 = tnfinder_out_filtered.loc[tnfinder_out_filtered['Stop_Query'] == j, 'Subject_ID'].iloc[0]
                
                    if named_2 in named_mge: 
                        named_mge[named_2].append(j)
                    else:
                        names.append(j)
                        named_mge[named_2] = names

                names = []


        for i in named_mge: 
            stop_seq = (len(named_mge[i])) - 1
            if stop_seq >= 1: 
                mge.loc[len(mge.index)] = [i, named_mge[i][0], named_mge[i][stop_seq]] 

        print('***************************************************************')
        print(mge)
        print('***************************************************************')

        final_transposons={}
        mge_list=mge.values.tolist()
        for row in mge_list:
            final_transposons[row[0]]=[row[1], row[2]]

        return final_transposons
    except:
        print("No transposons found by TnCompFinder")
        return []

            
#IN PROGRESS: return ncRNA from Infernal output
def getncRNA(plasmid):
    RNA=[]
    with open("./../output/plasmids/" + plasmid + "/ncRNA_out/pOXA48K8_ncRNA.out") as oriccsv:
        oricread=csv.reader(oriccsv, delimiter='\t')
        found=0
        for row in oricread:
            print(row)
            print(' ------ inclusion threshold ------' in row)
            if "------ inclusion threshold ------" in row:
                break
            elif found==1 or "rank" in row:
                if found==0:
                    found==1
                else:
                    print(row)
    return RNA

def main(plasmid):

    #get top match for oriC
    oriC=[]
    with open("./../output/plasmids/" + plasmid + "/oriC_out/oric.csv") as oriccsv:
        oricread=csv.reader(oriccsv)
        for row in oricread:
            oriC=row
            break
    
    #get top match for oriT
    oriT=[]
    with open("./../output/plasmids/" + plasmid + "/oriT_out/orit.csv") as oritcsv:
        oritread=csv.reader(oritcsv)
        for row in oritread:
            oriT=row
            break

    #transposons
    transposons=getTransposons(plasmid)

    #get gene lengths specified in gff to add in tsv
    gene_lengths=[]
    with open("./../output/plasmids/"+plasmid+"/the_final.gff") as ref_gff:
        gffread=csv.reader(ref_gff, delimiter='\t')
        
        for records in gffread:
            if '##gff-version' in records[0] or '##sequence-region' in records[0]:
                continue
            elif '#FASTA' in records[0] or '##Inc-GROUP' in records[0]:
                 break
            else:
                gene_lengths.append([records[3], records[4]])

    #keep track of previous gene end for overlap purposes
    prev_end=-1
    with open("./../output/plasmids/"+plasmid+"/the_finalTMP.tsv", 'w') as tmp_tsv:
        tsvwrite=csv.writer(tmp_tsv, delimiter='\t')
        with open("./../output/plasmids/"+plasmid+"/the_final.tsv") as ref_tsv:
            tsvread=csv.reader(ref_tsv, delimiter='\t')

            gene_index=0
            end_transposons={}
            for gene in tsvread:
                if "locus_tag" in gene[0]:
                    tsvwrite.writerow(gene)
                else:
                    #edit orf__ genes to just say orf
                    if ((str(gene[3]).lower()).startswith('orf') and str(gene[3][3:-1]).isnumeric()) or gene[3] in '':
                        gene[3]='orf'
                    
                    #add labels for oriC and oriT
                    if len(oriC)>0:
                        if checkOverlapOriX(gene_lengths[gene_index], oriC[6:8], prev_end):
                            ori_length=str(int(oriC[9])-int(oriC[8]))
                            tsvwrite.writerow(oriC[0:1]+["oriC"]+[ori_length]+['.']+['.'])
                            oriC=[]

                    if len(oriT)>0:
                        if checkOverlapOriX(gene_lengths[gene_index], oriT[6:8], prev_end):
                            ori_length=str(int(oriT[9])-int(oriT[8]))
                            tsvwrite.writerow(oriT[0:1]+["oriT"]+[ori_length]+['.']+['.'])
                            oriT=[]
                    
                    #add labels for transposon start
                    if len(transposons)>0:
                        for transp in transposons:
                            if len(transposons[transp])>0:
                                try:
                                    if checkOverlapTransp(gene_lengths[gene_index], [transposons[transp][0], transposons[transp][1]]):
                                        transp_length=str(transposons[transp][1]-transposons[transp][0])
                                        tsvwrite.writerow([str(transp)]+["Transposon START"]+[transp_length]+['.']+['.'])
                                        end_transposons[transp]=transposons[transp]
                                        transposons[transp]=[]
                                except:
                                    print(gene_lengths)
                                    print(gene_index)
                                    print(transposons[transp])
                  
                    #write gene to tsv
                    tsvwrite.writerow(gene)

                    #add labels for transposon end
                    if len(end_transposons)>0:
                        for transp in end_transposons:
                            if len(end_transposons[transp])>0:
                                if checkOverlapTransp(gene_lengths[gene_index], [end_transposons[transp][0], end_transposons[transp][1]]):
                                    transp_length=str(end_transposons[transp][1]-end_transposons[transp][0])
                                    tsvwrite.writerow([str(transp)]+["Transposon END"]+[transp_length]+['.']+['.'])
                                    end_transposons[transp]=[]
                          
                    
                    prev_end=gene_lengths[gene_index][1]
                    gene_index+=1

    os.rename("./../output/plasmids/"+plasmid+"/the_finalTMP.tsv", "./../output/plasmids/"+plasmid+"/the_final.tsv")
    

if __name__=="__main__":
        parser = argparse.ArgumentParser()
        parser.add_argument('plasmid', type=str)
        plasname=(parser.parse_args()).plasmid
        main(plasname)