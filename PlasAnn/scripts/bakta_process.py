#Process tsv file produced by bakta to match locus tags from default tsv.
import os
import csv
import argparse

#def getCorrectPositions(plasmid):
#    positions=[]
#    with open("./../output/plasmids/" + plasmid + "/gb_match/default/pk_results/" + plasmid + ".gff") as ref_gff:
#        gffread=csv.reader(ref_gff, delimiter=' ')
#        for records in gffread:
#            records=' '.join(records).split()
#            print(records)
#            if "#FASTA" in records[0]:
#                break
#            elif '##gff-version' not in records[0] and '##sequence-region' not in records[0]:
#                print(records)
#                positions.append([records[3], records[4]])

#    return positions

def main(plasmid):
    with open('./../../../pOXA_bakta/pOXA48K8.tsv') as tsvfile:
        in_tsv=list(csv.reader(tsvfile, delimiter='\t', quotechar='|'))
        
        #start at negative one to account for header
        reference_num=-1

        with open('./../output/plasmids/pOXA48K8/final.tsv') as reference:
            ref_in=csv.DictReader(reference, delimiter='\t', quotechar='|')
            for row in ref_in:
                reference_num+=1

            #start at negative two to account for header + null row at the end
            totalrows=-2
            for row in in_tsv:
                totalrows+=1

            if totalrows!=reference_num:
                positions=getCorrectPositions(plasmid)
                with open('./../../../pOXA_bakta/pOXA48K8_NonCDS.tsv', 'w') as newtsv:
                    tsv_writer = csv.writer(newtsv, delimiter = '\t')
                    #in_tsv.seek(0)
                    for row in in_tsv:
                        if 'ncRNA-region' in row:
                            break
                        if '#Annotated' in row[0] or '#Database' in row[0] or '#Sequence' in row[0]:
                            tsv_writer.writerow(row)
                        else:
                            if 'cds' not in row[1]:
                                tsv_writer.writerow(row)

        

if __name__=="__main__":
        parser = argparse.ArgumentParser()
        parser.add_argument('plasmid', type=str)
        plasname=(parser.parse_args()).plasmid
        main(plasname)