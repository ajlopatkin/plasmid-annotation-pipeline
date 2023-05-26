#!/usr/bin/python3

"""
Description: Call on UNIPROT to find a UNIPROT ID for each gene of the given plasmid and then use the GO ID linked
to that UNIPROT ID to label the gene's function

EDITED FROM PROKKA2KEGG by Heyu Lin (heyu.lin(AT)student.unimelb.edu.au)
 
"""

import os
import re
import csv
import gzip
import curses
import argparse
import json
import http.client
import requests
import urllib.request
from io import StringIO
from email.generator import Generator
import requests
import zlib
from goatools.base import download_go_basic_obo
from goatools.obo_parser import GODag
from goatools.utils import get_b2aset
from goatools.godag.go_tasks import get_go2parents
from goatools.obo_parser import GOTerm


http.client._MAXHEADERS=1000

__author__ = "Heyu Lin"
__contact__ = "heyu.lin(AT)student.unimelb.edu.au"

print("PARSING ARGUMENTS")
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', metavar='input_gbk', dest='i',
                    type=str, required=True)
parser.add_argument('-o', '--output', metavar='output', dest='o',
                    type=str, required=True)
parser.add_argument('-p', '--plasmid', metavar='plasmid', dest='p',
                    type=str, required=True)

args = parser.parse_args()

#Download go-basic and sort given GO ID's to use as references
obs_ids=[]
obs_items=[]
val_ids=[]
val_items=[]

fin_dag=download_go_basic_obo("go-basic.obo")
godag=GODag(fin_dag, optional_attrs={'consider', 'replaced_by', 'relationships'}, load_obsolete=True)
#print(godag)

for go_item in godag.values():
    if go_item.is_obsolete:
        obs_items.append(go_item)
        obs_ids.append(go_item.id)
    else:
        val_items.append(go_item)
        val_ids.append(go_item.id)

#Process go_id function descriptions based on parent descriptions to make labeling easier when mapping 
def process_go_id(go_id, func_desc, general, desc):
    keywords=['conjugation', 'conjugal', 'conjugative']

    print("\nProcessing GO-ID: " + go_id + " (" + func_desc + ")")

    addition=""
    for keyword in general:
        if keyword.lower() in keywords:
            addition=addition + " (" + keyword + ")"
    
    for keyword in keywords:
        if keyword in desc:
            addition=addition + " (" + keyword + ")"

    return func_desc + " " + addition
    
#Get GO ID for each uniprot id found by calling on UNIPROT 
def get_go_id(upro_ids):
    interest=['conjugation', 'replication', 'pilus', 'plasmid', 'stress', 'antitoxin', 'toxin', 'damage stimulus', 'metabolism', 'metabolic']


    go_ids=[]
    for uid in upro_ids:
        if len(uid[1])>0:
            poss_ids=[]
            find="U"
            for id in uid[1]:

                uni_http=requests.get("https://rest.uniprot.org/uniprotkb/"+id+"?format=txt")
                uni_txt=uni_http.text.split("\n")
                keywords=[]
                #Find keywords to use for processing	
                for unitem in uni_txt:
                    if "KW" in unitem[0:2]:
                        keys=unitem[5:].split(';')
                        for item in keys:
                            if len(item)>0:
                                if '.' in item[len(item)-1]:
                                    keywords.append(item[:len(item)-1])
                                else:
                                    keywords.append(item)	

                #Parse through UNIPROT returned entry to find info of interest
                found=0
                uid_gos=[]
                reviewed=False
                desc=""
                for line in uni_txt:
                    if 'ID' in line[0:2]:
                        if "Reviewed" in line:
                            reviewed=True
                    elif 'CC' in line[0:2]:
                        desc=desc+line
                    if 'DR' in line[0:2] and "GO;" in line:
                        info=[]
                        found=1	
                        spl_line=line.split(" ", 5)
                        info.append(uid[0])
                    
                        curr_id=spl_line[4][0:-1]
                        if curr_id in obs_ids:
                            loc=obs_ids.index(spl_line[4][0:-1])
                            if "" in obs_items[loc].alt_ids:
                                continue
                            else:
                                curr_id=""
                                for alt_id in obs_items[loc].alt_ids:
                                    if alt_id not in obs_ids:
                                        curr_id=alt_id
                                        break
                                if "" in curr_id:
                                    continue
                    
                        info.append(curr_id)
                        find=""
                        if reviewed:
                            find="R"
                        else:
                            find="U"
                        if 'P' in spl_line[5][0:1]:
                            find=find + "P"
                        elif 'C' in spl_line[5][0:1] and 'P' not in find:
                            find=find + "C"
                        elif 'F' in spl_line[5][0:1] and 'P' not in find:
                            find=find + "F"
                        info.append(find)
                        funct=spl_line[5][2:-1].split(';', 1)
                        info.append(process_go_id(curr_id, funct[0], keywords, desc))
                        uid_gos.append(info)

                if found==0:
                    if len(keywords)>0:
                        keyfind=0
                        interest=['conjugation', 'conjugal', 'replication', 'pilus', 'stress', 'toxin', 'damage stimulus', 'metabolism', 'metabolic']
                        for key in interest:
                            if key in keywords:
                                keyfind=1
                                rev_stat="U"
                                if 'R' in find:
                                    rev_stat="R"
                                go_ids.append([[uid[0], "", rev_stat, key]])
                                break
                        if keyfind==0:
                            go_ids.append([[uid[0], "", "", ""]])
                    else:
                        go_ids.append([[uid[0], "", "", ""]])
                    
                else:
                    compile_match=[]
                    seen=[]
                    for go_found in uid_gos:
                        if go_found[1] not in seen:
                             seen.append(go_found[1])
                             compile_match.append(go_found)

                    go_ids.append(compile_match)

    print("PROCESSED GO ID's FOR GIVEN PLASMID")
    print(go_ids)			
    return go_ids 

#Calling on UNIPROT to find the right uniprot id for each gene
def get_uniprot(gene_name):
    endpoint = f"https://rest.uniprot.org/uniprotkb/search?query=gene:" + gene_name + "+AND+taxonomy_id:2&format=json"
    resp = requests.get(endpoint)
    
    if resp.status_code != 200:
        print("ERROR NOT 200")
        return ''
    
    resp_json = json.loads(resp.text)['results']
    
    if len(resp_json) == 0:
        return ''
    else:
        return resp_json[0]['primaryAccession']

#Parse TSV and get UniprotID for each gene
def tsv_parser(file):
    found_uid=[]
    with open(file) as input:
        tsv_file = csv.reader(input, delimiter="\t")
        for line in tsv_file:
            if 'locus_tag' not in line:
                if not line[3].startswith('orf'):
                    if len(line[3])>0:
                        uid=get_uniprot(line[3])
                        if len(uid)==0:
                            found_uid.append([line[0], []])
                        else:
                            found_uid.append([line[0], [uid]])
    
    return found_uid

#Write GO IDs with correspdoning GO Functions to the TSV file
def output(arr, outfile):
    with open(outfile, 'w') as fo:
        fo.write("locus_tag" + "\t" + "GO_ID_P"  + "\t" + "GO_ID_C" + "\t" + "GO_ID_F" + "\t" + "Function" + "\n")
        for cds in arr:
            locus_tag=""
            go_ids_P=""
            go_ids_C=""
            go_ids_F=""
            function_P=""
            function_C=""
            function_F=""
            for goid in cds:
                if len(goid[1].replace(" ", ""))>0:
                    if len(locus_tag)==0:
                        locus_tag=goid[0]
                    if "P" in goid[2]:
                        go_ids_P=go_ids_P+goid[1]+"|"
                        if len(function_P)==0:
                            function_P=goid[3]
                        else:
                            function_P=function_P + ' | ' + goid[3]
                    elif "C" in goid[2]:
                        go_ids_C=go_ids_C+goid[1]+"|"
                        if len(function_P)==0:
                            function_C=goid[3]
                        else:
                            function_C=function_C + ' | ' + goid[3]
                    elif "F" in goid[2]:
                        go_ids_F=go_ids_F+goid[1]+"|"
                        if len(function_F)==0:
                            function_F=goid[3]
                        else:
                            function_F=function_F + ' | ' + goid[3]
                
            fo.write(locus_tag + "\t" + go_ids_P + "\t" + go_ids_C + "\t" +  go_ids_F + "\t" + function_P + "|" + function_C + "|" + function_F + "\n")

def main():
    mapping_array=tsv_parser(args.i)
    output(get_go_id(mapping_array), args.o)

if __name__ == '__main__':
    main()
