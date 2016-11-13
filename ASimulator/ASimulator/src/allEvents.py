# -*- coding: utf-8 -*-

import collections, argparse
import HTSeq
from sets import Set

parser = argparse.ArgumentParser(
	description = "Selection of genes under going complex and simple alternative splicing. The selected genes have only two possible transcripts in order to link gene name with event type/s",
	epilog =  "Authors: Carmen Bravo, Emma Houben, Carlos de Lannoy"    
	)
		
parser.add_argument("suppa_lists", nargs='*', type=str, help='Suppa event lists files')
parser.add_argument("-gff", "--gene_model", type=str, help='Reference annotation in GFF3 format.')
parser.add_argument("-o", "--output_file", type=str, default='allEvents_genome.txt',  help='Name of the output file containing gene name, event type and number of transcripts.')
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

 
allEvents = open(output, 'w')

for item in SElist:
	if item in gene_dict.keys():
		print >> allEvents, item + "\t" + "SE" + "\t" + str(gene_dict.get(item))		

for item in RIlist:
    if item in gene_dict.keys():
        print >> allEvents, item + "\t" + "RI" + "\t" + str(gene_dict.get(item))	

for item in A3list:
    if item in gene_dict.keys():
		print >> allEvents, item + "\t" + "A3" + "\t" + str(gene_dict.get(item))	

for item in A5list:
    if item	in gene_dict.keys():
		print >> allEvents, item + "\t" + "A5" + "\t" + str(gene_dict.get(item))	 

for item in MXlist:
    if item	in gene_dict.keys():
		print >> allEvents, item + "\t" + "MX" + "\t" + str(gene_dict.get(item))	

for item in ALlist:
    if item	in gene_dict.keys():
		print >> allEvents, item + "\t" + "AL" + "\t" + str(gene_dict.get(item))

for item in AFlist:
    if item	in gene_dict.keys():
		print >> allEvents, item + "\t" + "AF" + "\t" + str(gene_dict.get(item))	 	                                                                                                                                                               

allEvents.close()