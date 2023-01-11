#final edits for tsv file
import csv
import argparse
import re
import os

#check overlap between gene and oriC or oriT site to determine how to order the gene and oriC/oriT sites in the tsv file
def checkOverlap(gene, oriX):
    if (int(oriX[1])<=int(gene[0]) and int(gene[0])-int(oriX[1])>20) or (int(oriX[0])>=int(gene[1]) and int(oriX[0])-int(gene[1])>20):
        return 0
    return 1

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

    #ncRNA will be added later to be included in the final tsv

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

    #keep track of previous gene
    prev_gene="NONE"
    with open("./../output/plasmids/"+plasmid+"/the_finalTMP.tsv", 'w') as tmp_tsv:
        tsvwrite=csv.writer(tmp_tsv, delimiter='\t')
        with open("./../output/plasmids/"+plasmid+"/the_final.tsv") as ref_tsv:
            tsvread=csv.reader(ref_tsv, delimiter='\t')

            gene_index=0
            for gene in tsvread:
                if "locus_tag" in gene[0]:
                    tsvwrite.writerow(gene)
                else:
                    #if gene name repeats, skip locus and continue
                    if gene[3] in prev_gene or prev_gene in gene[3]:
                        gene_index+=1
                        continue
                    else:
                        #edit orf__ genes to just say orf
                        if (str(gene[3]).lower()).startswith('orf') and str(gene[3][3:-1]).isnumeric():
                            gene[3]='orf'
                        tsvwrite.writerow(gene)
                    
                    #add labels for oriC and oriT
                    if len(oriC)>0:
                        if checkOverlap(gene_lengths[gene_index], oriC[6:8]):
                            tsvwrite.writerow(oriC[0:1]+["oriC"]+oriC[8:9]+['.']+['.'])
                            oriC=[]

                    if len(oriT)>0:
                        if checkOverlap(gene_lengths[gene_index], oriT[6:8]):
                            tsvrite.writerow(oriT[0:1]+"oriT"+oriT[8:9]+'.'+'.')
                            oriT=[]
                          
                    prev_gene=str(gene[3])
                    gene_index+=1

    os.rename("./../output/plasmids/"+plasmid+"/the_finalTMP.tsv", "./../output/plasmids/"+plasmid+"/the_final.tsv")
    

if __name__=="__main__":
        parser = argparse.ArgumentParser()
        parser.add_argument('plasmid', type=str)
        plasname=(parser.parse_args()).plasmid
        main(plasname)