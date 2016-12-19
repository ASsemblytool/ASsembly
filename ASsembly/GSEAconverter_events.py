# -*- coding: utf-8 -*-

import collections, argparse

parser = argparse.ArgumentParser(
	description = "Conversion of the methods output files to a GSEA accepted format.",
	epilog =  "Authors: Carmen Bravo, Emma Houben, Carlos de Lannoy"    
	)
		
parser.add_argument("output_files", nargs='*', type=str, help='Methods output file/s. The file should not have a header and at least it should have in this order two columns: a first with the gene ids and a second with the given p-values.')
args = parser.parse_args()
output_files = args.output_files

for file in output_files:
	handle = open(file, 'r').readlines()
	ids = [line.split("\t")[0].replace('"', '') for line in handle]
	duplicates = set([x for x in ids if ids.count(x) > 1])
	for name in duplicates:
		count = 0
		for i in range(len(ids)):
			if ids[i] == name:
				count = count + 1
				ids[i] = ids[i] + "-" + str(count)	
	pvals = [line.split("\t")[1] for line in handle]
	corpvals = [float(1) - float(element) for element in pvals]
	gseafile = file.replace('.txt', '') + ".rnk"
	gseafileoutput = open(gseafile, 'w')
	for i in range(len(ids)):
		print >> gseafileoutput, ids[i] + "\t" + str(corpvals[i])	
	gseafileoutput.close()
	
		

