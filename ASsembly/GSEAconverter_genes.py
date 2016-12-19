# -*- coding: utf-8 -*-

import collections, argparse

parser = argparse.ArgumentParser(
	description = "Conversion of the methods output files to a GSEA accepted format. Adapted for Cufflinks and DEXseq.",
	epilog =  "Authors: Carmen Bravo, Emma Houben, Carlos de Lannoy"    
	)
		
parser.add_argument("-f", "--output_files", type=str, required='True', help='Methods output file. The file should not have a header and at least it should have in this order two columns: a first with the gene ids and a second with the given p-values.')
parser.add_argument("-el", "--eventList", type=str, required='True', help='List with the annotated events in the genome.')
parser.add_argument("-m", "--method", type=str, required='True', help='Name of the method that generated the results.')
parser.add_argument("-d", "--directory", type=str, required='True', help='Output directory.')
args = parser.parse_args()
output_files = args.output_files
eventList = args.eventList
method = args.method
directory=args.directory
eventhandle = open(eventList, 'r').readlines()

SEset=set()
RIset=set()
A3set=set()
A5set=set()
MXset=set()
AFset=set()
ALset=set()

for line in eventhandle:
	if line.split("\t")[1].strip("\n") == "SE":
		SEset.add(line.split("\t")[0].replace('"', ''))
	if line.split("\t")[1].strip("\n") == "RI":
		RIset.add(line.split("\t")[0].replace('"', ''))
	if line.split("\t")[1].strip("\n") == "A3":
		A3set.add(line.split("\t")[0].replace('"', ''))
	if line.split("\t")[1].strip("\n") == "A5":
		A5set.add(line.split("\t")[0].replace('"', ''))
	if line.split("\t")[1].strip("\n") == "MX":
		MXset.add(line.split("\t")[0].replace('"', ''))
	if line.split("\t")[1].strip("\n") == "AF":
		AFset.add(line.split("\t")[0].replace('"', ''))
	if line.split("\t")[1].strip("\n") == "AL":
		ALset.add(line.split("\t")[0].replace('"', ''))

SElist=list()
RIlist=list()
A3list=list()
A5list=list()
MXlist=list()
AFlist=list()
ALlist=list()

SEpval=list()
RIpval=list()
A3pval=list()
A5pval=list()
MXpval=list()
AFpval=list()
ALpval=list()


filehandle = open(output_files, 'r').readlines()
for line in filehandle:
	if line.split("\t")[0] in SEset:
		SElist.append(line.split("\t")[0].replace('"', ''))
		SEpval.append(line.split("\t")[1])
	elif line.split("\t")[0] in RIset:
		RIlist.append(line.split("\t")[0].replace('"', ''))
		RIpval.append(line.split("\t")[1])
	elif line.split("\t")[0] in A3set:
		A3list.append(line.split("\t")[0].replace('"', ''))
		A3pval.append(line.split("\t")[1])
	elif line.split("\t")[0] in A5set:
		A5list.append(line.split("\t")[0].replace('"', ''))
		A5pval.append(line.split("\t")[1])
	elif line.split("\t")[0] in MXset:
		MXlist.append(line.split("\t")[0].replace('"', ''))
		MXpval.append(line.split("\t")[1])
	elif line.split("\t")[0] in AFset:
		AFlist.append(line.split("\t")[0].replace('"', ''))
		AFpval.append(line.split("\t")[1])
	elif line.split("\t")[0] in ALset:
		ALlist.append(line.split("\t")[0].replace('"', ''))
		ALpval.append(line.split("\t")[1])
	
SEcorpvals = [float(1) - float(element) for element in SEpval]
RIcorpvals = [float(1) - float(element) for element in RIpval]
A3corpvals = [float(1) - float(element) for element in A3pval]
A5corpvals = [float(1) - float(element) for element in A5pval]
MXcorpvals = [float(1) - float(element) for element in MXpval]
AFcorpvals = [float(1) - float(element) for element in AFpval]
ALcorpvals = [float(1) - float(element) for element in ALpval]

SEfile = directory + "SE" + ".rnk"
SEoutput = open(SEfile, 'w')
for i in range(len(SElist)):
	print >> SEoutput, SElist[i] + "\t" + str(SEcorpvals[i])	
SEoutput.close()

RIfile = directory + "RI" + ".rnk"
RIoutput = open(RIfile, 'w')
for i in range(len(RIlist)):
	print >> RIoutput, RIlist[i] + "\t" + str(RIcorpvals[i])	
RIoutput.close()

A3file = directory + "A3" + ".rnk"
A3output = open(A3file, 'w')
for i in range(len(A3list)):
	print >> A3output, A3list[i] + "\t" + str(A3corpvals[i])	
A3output.close()

A5file = directory + "A5" + ".rnk"
A5output = open(A5file, 'w')
for i in range(len(A5list)):
	print >> A5output, A5list[i] + "\t" + str(A5corpvals[i])	
A5output.close()

MXfile = directory + "MX" + ".rnk"
MXoutput = open(MXfile, 'w')
for i in range(len(MXlist)):
	print >> MXoutput, MXlist[i] + "\t" + str(MXcorpvals[i])	
MXoutput.close()

ALfile = directory + "AL" + ".rnk"
ALoutput = open(ALfile, 'w')
for i in range(len(ALlist)):
	print >> ALoutput, ALlist[i] + "\t" + str(ALcorpvals[i])	
ALoutput.close()

AFfile = directory + "AF" + ".rnk"
AFoutput = open(AFfile, 'w')
for i in range(len(AFlist)):
	print >> AFoutput, AFlist[i] + "\t" + str(AFcorpvals[i])	
AFoutput.close()
	
		

