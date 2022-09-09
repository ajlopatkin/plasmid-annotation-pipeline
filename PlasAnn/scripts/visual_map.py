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

# ALT: dna_features_viewer_env
# https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer
# gene plotR: http://genoplotr.r-forge.r-project.org/index.php

def main(plasmid_name):
	#Create diagram w/two tracks: one for genes and one for length
	gd_diagram = GenomeDiagram.Diagram(plasmid_name, circular=True)
	gd_track_for_features = gd_diagram.new_track(0, name="Annotated Features", tracklines=0)
	gd_track_for_length = gd_diagram.new_track(-1, name="Length", tracklines=0)
	gd_feature_set = gd_track_for_features.new_set()
	gd_length_set = gd_track_for_length.new_set()

	print("PLASMID NAME:" + plasmid_name)

	pls_rec=[]
	with open("./../output/plasmids/" + plasmid_name + "/the_final.tsv") as plstsv:
		plsread=csv.reader(plstsv, delimiter='\t')
		for row in plsread:
			if "locus_tag" not in row[0]:
				pls_rec.append(row)

	length=0
	pls_loc=[]
	with open("./../output/plasmids/" + plasmid_name + "/the_final.gff") as plsgff:
		gffread=csv.reader(plsgff, delimiter='\t')
		for row in gffread:
			if '##' not in row[0]:
				pls_loc.append(row)
			elif 'sequence-region' in row[0]:
				bps=row[0].split(' ')
				length=int(bps[-1])
			elif '##Inc-GROUP' in row[0]:
				break

	pls_end=0
	color_list=[]
	
	for hue in colors.getAllNamedColors():
		color_list.append(hue)

	brites=[]	
	colors_used=[]

	#OLD JSON CODE: included in cases json file needs to be used as a database of colors for genes
	if os.path.isfile("./visual/colors.json"):
		print("colors.json exists")
		with open("./visual/colors.json", "r+") as jsonfile:
			if os.stat("./visual/colors.json").st_size>0:		
				json_objects=json.load(jsonfile)
				len(json_objects)
				print(json_objects)
				for data in json_objects:
					print(json_objects[data]["color"])	
					brites.append([data, json_objects[data]["color"]])
					colors_used.append(json_objects[data]["color"])

			else:
				json_color={}
				json_color["NO LABEL"]={
					"color": "silver"
				}				
				#jsonfile.write(json.dumps(json_color))
				brites.append(['NO LABEL', 'silver'])
				colors_used.append('silver')
	else:
		print("create color")
		json_color={}
		json_color["Hypothetical"]={
			"color": "silver"
		}
		json_color["Antibiotic Resistance"]={
			"color": "lightblue"
		}
		json_color["Conjugation"]={
			"color": "lightpurple"
		}
		#with open("./visual/colors.json", "w") as jsonfile:
		#	jsonfile.write(json.dumps(json_color))
		brites.append(['Hypothetical', 'silver'])
		brites.append(['Antibiotic Resistance', 'lightblue'])
		brites.append(['Conjugation', 'lightpurple'])
		colors_used.append('silver')
		colors_used.append('lightblue')
		colors_used.append('lightpurple')

	#Label each gene and color it based on determined gene function. Create features per gene to add to the diagram object
	with open("./visual/colors.json", "w") as jsonfile:
		print("Loading")
		for i in range(0, len(pls_loc)):
			print("Starting processing of colors")
			print(i)
			print(pls_rec[i])
			if "CDS" not in pls_rec[i][1] or len(pls_rec[i][3])==0:
				#Exclude this feature
				continue	
			
			description=""
			#print(pls_rec[i][9])
			if len(pls_rec[i][9])>2:
				#color_jsn=open("colors.json")
				if len(pls_rec[i])>10:
					for d in range(9, len(pls_rec[i]), 2):
						if len(pls_rec[i][d]) > len(description):
							description=pls_rec[i][d]
				else:
					description=pls_rec[i][9]
				print("DESCR")
				print(description)
				print(description.lower())

				#antibiotic=["antibiotic", "antimicrobial resistance"]
				#conjugation=["conjugation", "pilus", "plasmid"]

				if "antibiotic" in description.lower() or "antimicrobial resistance" in description.lower():
					print("ANTIBIOTIC")
					color="lightblue"
					#continue
				elif "replication" in description.lower():
					print("REPLICATION")
					color="pink"
				elif "conjugation" in description.lower()  or "conjugal" in description.lower()  or "transfer" in description.lower() or "replication" in description.lower() or "pilus" in description.lower() or "plasmid" in description.lower() or "dna binding" in description.lower():
					print("CONJUGATION")
					color="teal"
                                        #continue
				elif "metabolism" in description.lower() or "metabolic" in description.lower() or "catabolic" in description.lower() or "anabolic" in description.lower():
					print("METABOLISM")
					color="red"
				elif "stress" in description.lower() or "damage stimulus" in description.lower() or "repair" in description.lower():
					print("STRESS RESPONSE")
					color="orange"
				elif "antitoxin" in description.lower():
					print("ANTITOXIN")
					color="purple"
				elif "toxin" in description.lower():
					print("TOXINS")
					color="lightgreen"
				else:
					print("OTHER")
					color="yellow"
                                        #continue				


			else:
				if "conjugal" in pls_rec[i][6]:
					print("CONJUGATION")
					color="teal"	
				else:
					color=colors.silver

		
			value=None
			if "+" in pls_loc[i][6]:
				value=1
			elif "-" in pls_loc[i][6]:
				value=-1
	
			feature = SeqFeature(FeatureLocation(int(pls_loc[i][3]), int(pls_loc[i][4]), strand=value))
			gd_feature_set.add_feature(
				feature, 
				color=color, 
				label=True, 
				label_position="middle",
				label_strand=1,
				label_size=6, 
				label_angle=90, 
				name=pls_rec[i][3],  
				sigil="ARROW", 
				arrowhead_length=0.25)
			
			feature = SeqFeature(FeatureLocation(int(pls_loc[i][3]), int(pls_loc[i][4])))
			gd_length_set.add_feature(
				feature,
				color=color,
				label=True,
				label_strand=1,
				label_size=6,
				label_angle=90,
				sigil="BOX"
			)

			pls_end=int(pls_loc[i][4])
		

		data={}
		for brite in brites:
			data[brite[0]]={
				"color":brite[1]
			}
		jsonfile.write(json.dumps(data))

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
		print(num)
		if num==0:
			divide=divide+10
		else:
			divide=divide*10
	print(divide)		
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
	gd_diagram.write("./../output/plasmid_maps_go/plasmid_circular_" + plasmid_name + "_REF.pdf", "PDF")
	
	print(str(pls_end))
	print(str(length))
	
if __name__=='__main__':
        parser = argparse.ArgumentParser()
        parser.add_argument('plasmid', type=str)
        plasname=(parser.parse_args()).plasmid
        main(plasname)
