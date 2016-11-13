import collections, argparse
import HTSeq
from sets import Set

parser = argparse.ArgumentParser(
	description = "Given a list of gene ids, generates an incomplete annotation, where only the mRNA model of the first transcript is retained.",
	epilog =  "Authors: Carmen Bravo, Emma Houben, Carlos de Lannoy"    
	)
		
parser.add_argument("gene_list", type=str, help='Gene ids whose mRNA annotation will be removed.')
parser.add_argument("-gff", "--gene_model", type=str, help='Reference annotation in GFF3 format.')
parser.add_argument("-o", "--output_file", type=str, default='GFF_incomplete.gff',  help='Name of the incomplete GFF annotation file.')
args = parser.parse_args()
gene_list = args.gene_list
name_gff=args.gene_model
output=args.output_file

GENE = Set(["gene", "pseudogene", "transposable_element_gene", "chromosome"])
PROTEIN = Set(["protein"])
mRNA = Set(["miRNA", "mRNA_TE_gene", "ncRNA", "mRNA", "exon", "rRNA", "snRNA", "tRNA", "three_prime_UTR", "pseudogenic_exon", "pseudogenic_transcript", "CDS", "five_prime_UTR", "three_prime_UTR"])

file_gff=open(name_gff,'r')
outputfile=open(output, 'w')

with open(gene_list, 'r') as genefile:
    genes=genefile.read().split('\n')
    genes.pop()

for line in file_gff:
	if any(gene in line for gene in genes):		
		if line.split()[2] in GENE:
			print >> outputfile, line.rstrip()
		if line.split()[2] in mRNA:
			attr = line.split()[8]
			if attr.split(";")[0][-1] == "1":
				print >> outputfile, line.rstrip()
		if line.split()[2] in PROTEIN:
			attr = line.split()[8]
			if attr.split(";")[1][-1] == "1":
				print >> outputfile, line.rstrip()
		else:
			continue
	else:
		print >> outputfile, line.rstrip()
	
outputfile.close()

