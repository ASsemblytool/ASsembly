#!/bin/Rscript
#' ASpli script for AS detection
#' 
#' ASpli script for iterative application in a compound alternative splicing detection.
#' Written in the framework of the course "Integrated Bioinformatics Project" at KU Leuven.
#' 
#' @author: Carmen Bravo, Emma Houben, Carlos de Lannoy
#' 

# Test if there are 6 arguments, if not return an error.
args = commandArgs(TRUE)
if (length(args)=6) {
  stop("Supply 6 arguments in the follwoing order: path to bam folder, path to gtf file, regex for treatment files, regex for control files, readlength and maximum intron size. \nA sixth argument can be added to find the GTF file in different location.", call.=FALSE)
} else {
  wd = args[1]
  pathToGTF = args[2]
  patcase = args[3]
  patcontr = args[4]
  readlength = args[5]
  maxintronsize = args[6]
}

# Package management -----

# Install packages if necessary (suppose this doesn't apply once we go run it for real)
if(!any(installed.packages()[,"Package"]=="BiocInstaller") | !any(installed.packages()[,"Package"]=="BiocParallel"))
  source("https://bioconductor.org/biocLite.R")

list.of.packages = c("ASpli")
new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocInstaller::biocLite(new.packages)

library(ASpli)

# Analysis ----

# Gene annotation is handled using TxDb objects.
gtfFile = list.files(pathToGTF, pattern=".gtf",full.names=TRUE)
TxDb <- makeTxDbFromGFF(file=gtfFile,format="gtf")

# BinGenome method allows to extract genome coordinations at the gene, exon, intron and junction level and returns an object of class ASpliFeatures.
features <- binGenome(TxDb)

# Create targets object using the information of the available files:
bamFiles = list.files(paste0(wd,"/Bam_files"), pattern="\\.bam$", full.names=TRUE)
# Find number of treatment files and number of control files
lengthControl=length(grep(patcontr, bamFiles))
lengthCase=length(grep(patcase, bamFiles))
targets <- data.frame(bam=bamFiles, condition=c(rep(patcontr,lengthControl),rep(patcase,lengthCase)), row.names=sub("\\.bam$", "", bamFiles))

# Load bam files from targets data.frame:
bam <- loadBAM(targets)

# Overlap features and read alignments
# param l = read length of sequenced library
# param maxISize = maximum intron expected size. junctions longer than this size will be discarded.
# Variables should be given by user!
counts <- readCounts(features, bam, l=readlength, maxISize=maxintronsize)

# Alternative splicing analysis and discovery using junction
pair <- c(patcontr,patcase)
as <- AsDiscover(counts=counts, targets=targets, features=features, bam=bam, l=readlength, pair=pair)

# Differential expression
group <- c(rep(patcontr, lengthControl),rep(patcase, lengthCase))
du <- DUreport(counts, targets, pair, group)

## Writing all information away:
setwd(paste0(wd,"/int/"))
writeAll(counts=counts, du=du, as=as, output.dir="ASpli_output")

# Retrieve gene names and p-value, order based on p-value
summary<-list.files(paste0(wd,"/ASpli_output/"), pattern="summary", full.names=TRUE)
mycols <- rep("NULL",12)
mycols[c(2,3,11)] <- NA
toGSEA = read.table(summary, colClasses=mycols, row.names=NULL)
toGSEA = toGSEA[order(toGSEA[,3]),]
toGSEA = as.data.frame(toGSEA)
toGSEA <- toGSEA[toGSEA$event != '-',]
dir.create(paste0(wd,"/int/GSEA/ASpli/"))
IR_ASpliToGSEA <- toGSEA[toGSEA$event == 'IR' | toGSEA$event == 'IR*',]
ES_ASpliToGSEA <- toGSEA[toGSEA$event == 'ES'| toGSEA$event == 'ES*',]
A3SS_ASpliToGSEA <- toGSEA[toGSEA$event == 'Alt3ss' | toGSEA$event == 'Alt3ss*',]
A5SS_ASpliToGSEA <- toGSEA[toGSEA$event == 'Alt5ss' | toGSEA$event == 'Alt5ss*',]
write.table(IR_ASpliToGSEA,file=paste0(wd,"/int/GSEA/ASpli/RI.txt"),sep="\t",quote=F,row.names=F)
write.table(ES_ASpliToGSEA,file=paste0(wd,"/int/GSEA/ASpli/SE.txt"),sep="\t",quote=F,row.names=F)
write.table(A3SS_ASpliToGSEA,file=paste0(wd,"/int/GSEA/ASpli/A3.txt"),sep="\t",quote=F,row.names=F)
write.table(A5SS_ASpliToGSEA,file=paste0(wd,"/int/GSEA/ASpli/A5.txt"),sep="\t",quote=F,row.names=F)