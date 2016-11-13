import sys, argparse, math, random
import numpy as np
import itertools
from sets import Set
from scipy.stats import nbinom
from rpy2.robjects import r
import rpy2.robjects as robjects

try:
	import HTSeq
except ImportError:
	sys.stderr.write( "Could not import HTSeq. Please install the HTSeq Python framework\n" )
	sys.stderr.write( "available from http://www-huber.embl.de/users/anders/HTSeq\n" )
	sys.exit(1)

parser = argparse.ArgumentParser(
	description = "Generated negative binomial relationship gene-level read counts for a synthetic groups based on real experimental data",
	epilog =  "Adaptation from Liu et al., 2014. Authors: Carmen Bravo, Emma Houben, Carlos de Lannoy"    
	)

parser.add_argument("gene_model", type=str, help='Reference annotation in GFF3 format.')
parser.add_argument("-g1", "--group1", nargs='+', type=str, required=True, help='First group bam files separated by space.')
parser.add_argument("-g2", "--group2", nargs='+', type=str, required=True, help='Second group bam files separated by space.')
parser.add_argument("-as", "--ASgenes", type=str, default= 'Random', help = "File of with gene IDs that will be alternatively spliced (if expressed). If not specified, genes are selected randomly.")
parser.add_argument("-n", "--num-reps", type=int, default = 3, dest='nreps', help = "Number of replicates. Default is 3." )
parser.add_argument("-l", "--num-target-gene", type=int, default = 2000, dest='ntarg', help = "Number of genes that will be differentially spliced if random selection is chosen. Default is 2000.")
parser.add_argument("-o", "--output_dir", type=str, required=True,  help='Name of the output directory.')
parser.add_argument("-of", "--output_file", type=str, default="DiffASgenes.txt", help='Name of the output file containing the genes that will be differentially spliced.')

if len(sys.argv) == 1:
	parser.print_help()
	sys.exit(1)

GENE = Set(["gene", "transposable_element_gene","pseudogene"])
EXON = Set(["exon", "pseudogenic_exon"])

## Parse inputs
args = parser.parse_args()
name_gtf = args.gene_model
NREPS=args.nreps
NTARG = args.ntarg
output = args.output_dir
outputfile = args.output_file
num_lines = sum(1 for line in open(name_gtf))

## Define global variables
def countbam(_bam_file, _genes, _dic, _idx):
	'''
	This function process a single bam file and count the reads sitting on each gene
	Results are save in the global variable.
	'''
	num_reads = 0
	for alnmt in _bam_file:
		if alnmt.aligned: 
			intersection_set = None
			#print list(_genes[alnmt.iv].steps())
			for iv1, step_set in _genes[alnmt.iv].steps():
				if intersection_set is None:
					intersection_set = step_set.copy()
				else:
					intersection_set.intersection_update(step_set)
			if len(intersection_set) == 1:
				_dic[ list(intersection_set)[0] ][_idx] += 1
		num_reads += 1
		if num_reads % 1000000 == 0:
			sys.stderr.write("%d reads processed.\n" % num_reads)
	return _dic

def meanVar(_files, _gff_file , _output):

	NFILE=len(_files)
	if NFILE == 1:
		sys.stderr.write("Need at least two samples for each group.\n")
		sys.exit(1)
		
	## Dictionary of gene counts
	_dict_counts = dict() 
	_genes = HTSeq.GenomicArrayOfSets("auto",stranded=False)
	idx=0
	count = 0
	transcript= set()
	cur_line = None
	lines = 0
	for feature in _gff_file:
		lines += 1
		if feature.type in GENE or lines == num_lines:
			if len(transcript) >1:
				_dict_counts[ cur_line.name ] = [0]*NFILE
				_genes[cur_line.iv] += cur_line.name
				count +=1
			cur_line = feature
			transcript.clear()
		if feature.type in EXON:
			transcript.add(feature.attr["Parent"])
	print "Number of genes", count
	_file_raw_count = open(_output+'.rawcounts','w')
	_file_nb_count = open(_output+'.nbcounts','w')
	## This loop read through the input list and call countbam for each input file  
	for f in _files:
		bam_file=HTSeq.BAM_Reader(f)
		_dict_counts=countbam(bam_file, _genes,_dict_counts, idx)
		idx += 1
		sys.stderr.write("Library %d has generated.\n" % idx)
	## Print raw counts in file specified by <out>
	for key, value in sorted(_dict_counts.iteritems()):
		_file_raw_count.write(key+"\t"+"\t".join(map(str,value))+"\n")
	_file_raw_count.close()
	## Calculate group mean and variance
	list_mean = list()
	list_var = list()
	for key, value in sorted(_dict_counts.iteritems()):
		list_mean.append(np.mean(np.array(value)))
		list_var.append(np.var(np.array(value)))
	
	## Computer loess esimates	
	## The following code is using rpy2 module
	a = robjects.FloatVector(list_mean)
	b = robjects.FloatVector(list_var)
	df = robjects.DataFrame({"mean": a, "var": b})
	non0_df=df.rx(df.rx2("mean").ro > 0, True) ## subsetting if mean > 0
	loess_fit = r.loess("var ~ mean", data=non0_df, degree=2)

	var_pred = r.predict(loess_fit, a)
	# This loop overwrite global variable dict_counts for recoding new count data
	count_idx = 0
	for key, value in sorted(_dict_counts.iteritems()):
		n = math.pow(list_mean[count_idx],2)/(var_pred[count_idx]-list_mean[count_idx])
		n = int(n) # n: number of failures
		if n<=0:
			_dict_counts[key] = [0]*NREPS
		else:
			p = n/float(n+list_mean[count_idx]) # p: prob of success
			_dict_counts[key] = nbinom.rvs(n, p, size=NREPS).tolist()
		count_idx += 1
	for key, value in sorted(_dict_counts.iteritems()):
		_file_nb_count.write(key+"\t"+"\t".join(map(str,value))+"\n")
	_file_nb_count.close()
	_file_raw_count.close()
	return _dict_counts
def main():
	## first group
	group1_f = args.group1
	file_gtf=open(name_gtf,'r')
	gff_file=HTSeq.GFF_Reader(file_gtf)

	##### Sanity Check
	bamfile=HTSeq.BAM_Reader(group1_f[0])
	is_chr_bam = None
	set_chr_gff = set()
	for almnt in bamfile:
		is_chr_bam = almnt.iv.chrom
		break
	for feature in gff_file:
		set_chr_gff.add(feature.iv.chrom)
	
	if  is_chr_bam not in set_chr_gff:
		sys.stderr.write("Error: Chromosome id in bam files and GFF file does not agree!!\n")
		sys.exit(1)
	#####
	
	file_gtf=open(name_gtf,'r')
	gff_file=HTSeq.GFF_Reader(file_gtf)	
	outputname= output + "/group1"
	counts1=meanVar(group1_f, gff_file,outputname)	
	file_gtf.close()
	
	## second group
	group2_f = args.group2
	file_gtf=open(name_gtf,'r')
	gff_file=HTSeq.GFF_Reader(file_gtf)
	outputname= output + "/group2"
	counts2=meanVar(group2_f, gff_file,outputname)
	file_gtf.close()
	
	merged_count = {k1: v1+v2 for (k1,v1) in counts1.iteritems() for (k2,v2) in counts2.iteritems() if k1 == k2}
	non0_exp_list = list()
	for k,v in sorted(merged_count.iteritems()):
		value=np.array(v)
		if len(value[value ==0]) == 0:
			non0_exp_list.append(k)

	num_non0 = len(non0_exp_list)
	name_as = args.ASgenes
	out = open(outputfile,'w')
	
	if name_as == 'Random':
		sys.stderr.write("Randomly selected %d out of %d genes as the final target genes that are to undergo differential AS\n" % (NTARG,num_non0))
		l=random.sample(xrange(num_non0), NTARG)
		for i in l:
			out.write(non0_exp_list[i]+"\n")
		out.close()	
	else:	
		file_asg = open(name_as, 'r')
		numAS = 0
		for line in file_asg:
			if line.split()[0] in non0_exp_list:
				out.write(line)
				numAS += 1	
		out.close()
		file_asg.close()
	sys.stderr.write("Selected %d out of %d genes according to event type selection as the final target genes that are to undergo differential AS\n" % (numAS,num_non0))
if __name__ == '__main__': main()
