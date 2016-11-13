# -*- coding: utf-8 -*-

import collections, argparse
import HTSeq
from sets import Set

parser = argparse.ArgumentParser(
	description = "Selection of genes under going complex alternative splicing. The selected genes have only two possible transcripts in order to link gene name with event type/s.",
	epilog =  "Authors: Carmen Bravo, Emma Houben, Carlos de Lannoy"    
	)
		
parser.add_argument("suppa_lists", nargs='*', type=str, help='Suppa event lists files')
parser.add_argument("-gff", "--gene_model", type=str, help='Reference annotation in GFF3 format.')
parser.add_argument("-o", "--output_file", type=str, default='complexEvents_genome.txt',  help='Name of the output file containing gene name, event type and number of transcripts.')
args = parser.parse_args()
suppa_lists = args.suppa_lists
name_gff=args.gene_model
output=args.output_file

SElist = list()
RIlist = list()
A3list = list()
A5list = list()
MXlist = list()
ALlist = list()
AFlist = list()

for file in suppa_lists:
	if file[-13:-11] == 'SE':
		SElist = [x.split('\t')[1] for x in open(file).readlines()[1:]]
	elif file[-13:-11] == 'RI':
		RIlist = [x.split('\t')[1] for x in open(file).readlines()[1:]]
	elif file[-13:-11] == 'A3':
		A3list = [x.split('\t')[1] for x in open(file).readlines()[1:]]
	elif file[-13:-11] == 'A5':
		A5list = [x.split('\t')[1] for x in open(file).readlines()[1:]]
	elif file[-13:-11] == 'MX':
		MXlist = [x.split('\t')[1] for x in open(file).readlines()[1:]]
	elif file[-13:-11] == 'AL':
		ALlist = [x.split('\t')[1] for x in open(file).readlines()[1:]]	
	elif file[-13:-11] == 'AF':
		AFlist = [x.split('\t')[1] for x in open(file).readlines()[1:]]				

versusSE = sum([RIlist, A3list, A5list, MXlist, ALlist, AFlist, list([x for x in SElist if SElist.count(x) > 1])],[])
versusRI = sum([SElist, A3list, A5list, MXlist, ALlist, AFlist, list([x for x in RIlist if RIlist.count(x) > 1])],[])
versusA3 = sum([RIlist, SElist, A5list, MXlist, ALlist, AFlist, list([x for x in A3list if A3list.count(x) > 1])],[])
versusA5 = sum([RIlist, A3list, SElist, MXlist, ALlist, AFlist, list([x for x in A5list if A5list.count(x) > 1])],[])

SElist2 = list()
for item in versusSE:
        while item in SElist:
                SElist2.append(item)
                SElist.remove(item)

RIlist2 = list()
for item in versusRI:
        while item in RIlist:
                RIlist2.append(item)
                RIlist.remove(item)

A3list2 = list()
for item in versusA3:
        while item in A3list:
                A3list2.append(item)
                A3list.remove(item)

A5list2 = list()
for item in versusA5:
        while item in A5list:
                A5list2.append(item)
                A5list.remove(item)

GENE = Set(["gene", "pseudogene", "transposable_element_gene"])
EXON = Set(["exon", "pseudogenic_exon"])

num_lines = sum(1 for line in open(name_gff))

file_gff=open(name_gff,'r')
gff_file=HTSeq.GFF_Reader(file_gff)

count = 0
transcript = set()
lines = 0
gene_dict={}

for feature in gff_file:
	lines += 1
	if feature.type in GENE or lines == num_lines:
		if len(transcript) == 2:
			count += 1
			gene_dict[gene_cand.attr["ID"]] = len(transcript)
		gene_cand = feature
		transcript.clear()
	if feature.type in EXON:
		transcript.add(feature.attr["Parent"])

 
complexEvents = open(output, 'w')

for item in SElist2:
	if item in gene_dict.keys():
		print >> complexEvents, item + "\t" + "SE" + "\t" + str(gene_dict.get(item))		

for item in RIlist2:
    if item in gene_dict.keys():
        print >> complexEvents, item + "\t" + "RI" + "\t" + str(gene_dict.get(item))	

for item in A3list2:
    if item in gene_dict.keys():
		print >> complexEvents, item + "\t" + "A3" + "\t" + str(gene_dict.get(item))	

for item in A5list2:
    if item	in gene_dict.keys():
		print >> complexEvents, item + "\t" + "A5" + "\t" + str(gene_dict.get(item))	 

for item in MXlist:
    if item	in gene_dict.keys():
		print >> complexEvents, item + "\t" + "MX" + "\t" + str(gene_dict.get(item))	

for item in ALlist:
    if item	in gene_dict.keys():
		print >> complexEvents, item + "\t" + "AL" + "\t" + str(gene_dict.get(item))

for item in AFlist:
    if item	in gene_dict.keys():
		print >> complexEvents, item + "\t" + "AF" + "\t" + str(gene_dict.get(item))	 	                                                                                                                                                               

complexEvents.close()
