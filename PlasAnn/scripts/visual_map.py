"""
VISUAL_MAP.PY

Description: Using gene name and function infomartion for the given plasmid from its final.tsv file, create a plasmid gene map.
Color coded based on the pre-defined functions
"""

from reportlab.lib import colors
from Bio.SeqFeature import SeqFeature, FeatureLocation
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO
import csv
import random
import argparse
import json
import os
import sys
import pandas as pd

# ALT: dna_features_viewer_env
# https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer
# gene plotR: http://genoplotr.r-forge.r-project.org/index.php

#process categories to add categories for oriC, oriT, and transposons
def processCategories(old_cat, plasmid):
    new_cat=[]
    index=0
    with open("./../output/plasmids/" + plasmid + "/the_final.tsv") as plstsv:
        plsread=csv.reader(plstsv, delimiter='\t')
        for row in plsread:
            if "locus_tag" not in row[0]:
                if '.' in row[3] and len(row[3])==1:
                    if 'oriC' in row[1]:
                        new_cat.append('oriC')
                    elif 'oriT' in row[1]:
                        new_cat.append('oriT')
                    else:
                        new_cat.append('Transposons (Mobile Element)')
                else:
                    new_cat.append(old_cat[index]) 
                    index+=1
    return new_cat

#go through TSV and GBK files and get information for each gene
def processTSV(plasmid):
    info=[]
    with open("./../output/plasmids/" + plasmid + "/the_final.tsv") as plstsv:
        plsread=csv.reader(plstsv, delimiter='\t')
        for row in plsread:
            if "locus_tag" not in row[0]:
                info.append(row)

    genes=[]
    transposons=[]
    for site in info:
        if site[3] in '.':
            transposons.append(site)
        else:
            genes.append(site)
    
    length=0
    pls_loc=[]
    with open("./../output/plasmids/" + plasmid + "/the_final.gff") as plsgff:
        gffread=csv.reader(plsgff, delimiter='\t')
        for row in gffread:
            if '##' not in row[0]:
                pls_loc.append(row)
            elif 'sequence-region' in row[0]:
                bps=row[0].split(' ')
                length=int(bps[-1])
            elif '##Inc-GROUP' in row[0]:
                break
    
    return genes, transposons, pls_loc, length

pls_end=0
def main(plasmid):
    print("VISUALIZING " + plasmid)
    #Create diagram w/two tracks: one for genes and one for length
    gd_diagram = GenomeDiagram.Diagram(plasmid, circular=True)
    gd_track_for_features = gd_diagram.new_track(0, name="Annotated Features", tracklines=0)
    gd_track_for_length = gd_diagram.new_track(-1, name="Length", tracklines=0)
    gd_feature_set = gd_track_for_features.new_set()
    gd_length_set = gd_track_for_length.new_set()
    
    genes, transposons, pls_loc, length=processTSV(plasmid)
    
    #Select the category + color for each gene based on the functional descriptions
    category=[]
    mobile=' (Mobile Element)'
    prev=""
    for i in range(0, len(genes)):
        description=""

        if genes[i][3].startswith('tra') and len(genes[i][3])==4:
            description=genes[i][6] + "---" + genes[i][10] + "---" +  genes[i][11]
            color="teal"
            category.append('Conjugation')
        elif genes[i][3].startswith('orf'):
            color=colors.silver
            category.append('N/A')
        elif genes[i][3].startswith('pem'):
            color='purple'
            category.append('Toxin/Antitoxin')

        elif len(genes[i][10])>2 or len(genes[i][11]) or len(genes[i][6]):
            if len(genes[i])>10:
                description=genes[i][6] + "---" + genes[i][10] + "---" +  genes[i][11]
            else:
                description=genes[i][6]

            if "antibiotic" in description.lower() or "antimicrobial resistance" in description.lower():
                color="red"
                category.append('Antibiotic Resistance')

            elif "conjugation" in description.lower() or "conjugal" in description.lower() or "pilus" in description.lower() or "plasmid" in description.lower():
                color="teal"
                category.append('Conjugation')
            elif "metabolism" in description.lower() or "metabolic" in description.lower() or "catabolic" in description.lower() or "anabolic" in description.lower():
                color="lightblue"
                category.append('Metabolism')
            elif "stress" in description.lower() or "damage stimulus" in description.lower() or "repair" in description.lower():
                color="orange"
                category.append('Stress Response')
            elif "antitoxin" in description.lower() or "toxin" in description.lower():
                color="purple"
                category.append('Toxin/Antitoxin')
            elif "replication" in description.lower():
                color="pink"
                category.append('Replication')
            else:
                color="yellow"		
                category.append('Other')	

        else:
            if "conjugal" in genes[i][6]:
                color="teal"
                category.append('Conjugation')	
                
            else:
                color=colors.silver
                category.append('N/A')
        
        if genes[i][3].startswith('tnp') or genes[i][3].startswith('ins') or genes[i][3].startswith('tnp'):
            category[-1]=category[-1]+mobile
            color="lightgreen"

        
        value=None
        start=0
        stop=0
        for loc in pls_loc:
            if genes[i][0] in loc[8]:
                start=int(loc[3])
                stop=int(loc[4])

                if "+" in loc[6]:
                    value=1
                elif "-" in loc[6]:
                    value=-1
    
        feature = SeqFeature(FeatureLocation(start, stop, strand=value))
        gd_feature_set.add_feature(
            feature, 
            color=color, 
            label=True, 
            label_position="middle",
            label_strand=1,
            label_size=12, 
            label_angle=90, 
            name=genes[i][3],  
            sigil="ARROW", 
            arrowhead_length=0.25)
            
        feature = SeqFeature(FeatureLocation(start, stop))
        gd_length_set.add_feature(
            feature,
            color=color,
            label=True,
            label_strand=1,
            label_size=6,
            label_angle=90,
            sigil="BOX"
        )

        pls_end=stop

    category=processCategories(category, plasmid)
    if int(length)!=int(pls_end):
        feature = SeqFeature(FeatureLocation(int(pls_end), int(length)))
        gd_feature_set.add_feature(
            feature,
            color="transparent"
        )

        gd_length_set.add_feature(
            feature,
            color="transparent"
        )

    #Make sure that length of full plasmid is rep. in length diagram by filling in gap between last gene and full plasmid length
    divide=0
    for num in range(0, len(str(length))-1):

        if num==0:
            divide=divide+10
        else:
            divide=divide*10

    tics=int(length/divide)
    for i in range(0, tics + 1):
        feature = SeqFeature(FeatureLocation(int(i*divide), int(i*divide)))
        gd_length_set.add_feature(
            feature,
            label=True,
            label_strand=1,
            label_size=6,
            color="black",
            name=str(i*10000)
        )

    gd_diagram.draw(
        format="circular",
        circular=True,
        pagesize=(20 * cm, 20 * cm),
        start=0,
        end=pls_end,
        circle_core=0.3,
    tracklines=0
    )
    
    os.system("mkdir ./../output/plasmid_maps_go") 
    gd_diagram.write("./../output/plasmid_maps_go/TESTER" + plasmid + ".pdf", "PDF")
      
    file = pd.read_csv("./../output/plasmids/" + plasmid + "/the_final.tsv", sep='\t')  
    file['Category']=category
    file.head() 
    file.to_csv("./../output/plasmids/" + plasmid + "/the_final.tsv", sep="\t", index=False)

if __name__=='__main__':
        parser = argparse.ArgumentParser()
        parser.add_argument('plasmid', type=str)
        plasname=(parser.parse_args()).plasmid
        main(plasname)
