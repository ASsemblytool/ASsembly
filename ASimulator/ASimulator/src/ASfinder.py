# -*- coding: utf-8 -*-


try:
	import HTSeq
except ImportError:
	sys.stderr.write( "Could not import HTSeq. Please install the HTSeq Python framework\n" )
	sys.stderr.write( "available from http://www-huber.embl.de/users/anders/HTSeq\n" )
	sys.exit(1)

import argparse, sys
from sets import Set

parser = argparse.ArgumentParser(
	description = "Program to find regions corresponding to potentially AS genes in a genome annotation file (gff). The output is a bed file containing the regions of annotated AS genes. ",
	epilog =  "Authors: Carmen Bravo, Emma Houben, Carlos de Lannoy"    
	)

parser.add_argument("gene_model", type=str, help='Reference annotation in GFF3 format.')
parser.add_argument("-o", "--output", type=str, required=True,  help='Name of the output file.')

args = parser.parse_args()
name_gff = args.gene_model
output = args.output

if len(sys.argv) == 1:
	parser.print_help()
	sys.exit(1)
	
if name_gff[-3:] != "gff":
	print "Provide a GFF file!"
	sys.exit(1)
	
GENE = Set(["gene", "transposable_element_gene","pseudogene"])
EXON = Set(["exon", "pseudogenic_exon"])
num_lines = sum(1 for line in open(name_gff))

file_gff=open(name_gff,'r')
gff_file=HTSeq.GFF_Reader(file_gff)

transcript = set()
gene_list = list()
chrom_list = list()
start_list = list()
end_list = list()
strand_list = list()
lines = 0
for feature in gff_file:
	lines += 1
	if feature.type in GENE or lines == num_lines:
		if len(transcript) > 1:
			gene_list.append(gene_cand.attr["ID"])
			chrom_list.append(gene_cand.iv.chrom)
			start_list.append(gene_cand.iv.start)
			end_list.append(gene_cand.iv.end)
			strand_list.append(gene_cand.iv.strand)
		gene_cand = feature
		transcript.clear()
	if feature.type in EXON:
		transcript.add(feature.attr["Parent"])
print 'Number of potentially AS genes: ', len(gene_list)

ASregs = open(output, 'w')
i=0
while i < len(gene_list):
        print >> ASregs, str(chrom_list[i]) + "\t" + str(start_list[i]) + "\t" +  str(end_list[i]) + "\t" + gene_list[i] + "\t" + "0" + "\t" + str(strand_list[i])
        i += 1

file_gff.close()
ASregs.close()


