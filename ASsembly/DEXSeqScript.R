#!/bin/Rscript
#' DEXSeq script for AS detection
#' 
#' DEXSeq script for iterative application in a compound alternative splicing detection.
#' Written in the framework of the course "Integrated Bioinformatics Project" at KU Leuven.
#' 
#' @author: Carmen Bravo, Emma Houben, Carlos de Lannoy

# test if there is at least one argument: if not, return an error
args = commandArgs(TRUE)
if (length(args)!=4) {
  stop("Supply 4 arguments in the following order: path to bam and gff files, number of nodes, regex for treatment files, regex for control files.", call.=FALSE)
}

wd = args[1]
nNodes = as.numeric(args[2])
treatment = args[3]
control = args[4]



# Package management -----

# Install packages if necessary (suppose this doesn't apply once we go run it for real)

list.of.packages = c("rPython")
new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

if(!any(installed.packages()[,"Package"]=="BiocInstaller") | !any(installed.packages()[,"Package"]=="BiocParallel"))
  source("https://bioconductor.org/biocLite.R")

list.of.biocPackages = c("DEXSeq")
new.packages = list.of.biocPackages[!(list.of.biocPackages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocInstaller::biocLite(new.packages)

library(DEXSeq)
library(rPython)
library(BiocParallel)

# Analysis ----

# Define parameter required for parallel computing
BPPARAM = MulticoreParam(workers=nNodes)

# Load count files created by dexseq_count.py.

countFiles = list.files(paste0(wd,"/int/DEXSeq_output/HTSeqCount_files"), pattern="_count.txt$", full.names=TRUE)
# Load the GFF-file (=flattened GTF-file)
flattenedFile = list.files(paste0(wd,"/int/DEXSeq_output"), pattern="gff$", full.names=TRUE)

# Create table containing info on condition and library type per count file
samples = list.files(paste0(wd,"/int/DEXSeq_output/HTSeqCount_files"))
samples = basename(samples)
toRemove = paste(sub(paste0(control,"|",treatment),"",samples),collapse="|")
condition = sub(toRemove,"",samples)
sampleTable = data.frame(
  row.names = samples,
  condition = condition,
  libType = rep("paired-end",length(samples)))

# create dexseq data object
dxd = DEXSeqDataSetFromHTSeq(
  countFiles,
  sampleData=sampleTable,
  design= ~ sample + exon + condition:exon,
  flattenedfile=flattenedFile )

# Normalize for sequencing depth
dxd = estimateSizeFactors( dxd )

# Estimate dispersion
dxd = estimateDispersions( dxd , BPPARAM=BPPARAM)

# Test for Differential Exon Usage (DEU)
dxd = testForDEU( dxd , BPPARAM=BPPARAM)

# # Extra step, may be omitted if not useful: estimate fold change
# dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition")

# create dexseq result object (a dataframe)
dxr1 = DEXSeqResults( dxd )


# Retrieve the gene names and p.adj, order list on p.adj
# (adjustment is Benjamini-Hochberg)
# exonNames = paste(dxr1$groupID,dxr1$featureID,sep="_")
toGSEA = data.frame(geneNames =  dxr1$groupID,padj = dxr1$padj)
toGSEA = toGSEA[order(toGSEA[,2]),]
toGSEA = toGSEA[!is.na(toGSEA[,2]),]

# Write data frame
write.table(toGSEA,file=paste0(wd,"/int/toGSEA/DEXSeq_toGSEA.txt"),sep="\t",quote=F,row.names=F,col.names=F)
