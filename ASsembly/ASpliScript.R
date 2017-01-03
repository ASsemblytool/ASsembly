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
if (length(args) != 8) {
  stop("Supply 8 arguments in the follwoing order: path to bam folder, path to gtf file, regex for treatment files, regex for control files, readlength and maximum intron size, output path and number of cores.", call.=FALSE)
} else {
  wd = args[1]
  gtfFile = args[2]
  patcase = args[3]
  patcontr = args[4]
  readlength = as.numeric(args[5])
  maxintronsize = as.numeric(args[6])
  outputPath = args[7]
  nCores = as.numeric(args[8])
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
#gtfFile = list.files(pathToGTF, pattern=".gtf",full.names=TRUE)
TxDb <- makeTxDbFromGFF(file=gtfFile,format="gtf")

# BinGenome method allows to extract genome coordinations at the gene, exon, intron and junction level and returns an object of class ASpliFeatures.
features <- binGenome(TxDb)

# Create targets object using the information of the available files:
bamFiles = list.files(paste0(wd,"/bam_files"), pattern="\\.bam$", full.names=TRUE)
# Find number of treatment files and number of control files
lengthControl=length(grep(patcontr, bamFiles))
lengthCase=length(grep(patcase, bamFiles))

bamFilesNoPath=basename(bamFiles)
toRemove = paste(sub(paste0(patcontr,"|",patcase),"",bamFilesNoPath),collapse="|")
condition = sub(toRemove,"",bamFilesNoPath)

targets <- data.frame(bam=bamFiles, 
                      condition=condition, 
                      row.names=bamFilesNoPath)

# Load bam files from targets data.frame:
bam <- loadBAM(targets,cores=nCores)

# Overlap features and read alignments
# param l = read length of sequenced library
# param maxISize = maximum intron expected size. junctions longer than this size will be discarded.
# Variables should be given by user!
print("Start readCounts()")
counts <- readCounts(features, bam, l=readlength, maxISize=maxintronsize, cores=nCores)

# Alternative splicing analysis and discovery using junction
pair <- c(patcontr,patcase)
print("Start AsDiscover()")
as <- AsDiscover(counts=counts, targets=targets, features=features, bam=bam, l=readlength, pair=pair, cores=nCores)

# Differential expression
print("Start DUreport()")
group <- c(rep(patcontr, lengthControl),rep(patcase, lengthCase))
du <- DUreport(counts, targets, pair, group)

## Writing all information away:
oldwd = getwd()
setwd(outputPath) # Issue with writeAll: pastes output.dir argument behind current path. Workaround applied.
writeAll(counts=counts, du=du, as=as, output.dir="raw")
setwd(oldwd)

# Retrieve gene names and p-value, order based on p-value
summary<-list.files(paste0(outputPath,"/raw"), pattern="summary", full.names=TRUE)
mycols <- rep("NULL",12)
mycols[c(2,3,11)] <- NA
toGSEA = read.table(summary, colClasses=mycols, row.names=NULL)
toGSEA = toGSEA[order(toGSEA[,3]),]
toGSEA = as.data.frame(toGSEA)
toGSEA <- toGSEA[toGSEA$event != '-',]
dir.create(paste0(wd,"/int/GSEA/ASpli/"), recursive = TRUE)
IR_ASpliToGSEA <- toGSEA[toGSEA$event == 'IR' | toGSEA$event == 'IR*',]
ES_ASpliToGSEA <- toGSEA[toGSEA$event == 'ES'| toGSEA$event == 'ES*',]
A3SS_ASpliToGSEA <- toGSEA[toGSEA$event == 'Alt3ss' | toGSEA$event == 'Alt3ss*',]
A5SS_ASpliToGSEA <- toGSEA[toGSEA$event == 'Alt5ss' | toGSEA$event == 'Alt5ss*',]
write.table(IR_ASpliToGSEA[,c(2,3,1)],file=paste0(wd,"/int/GSEA/ASpli/RI.txt"),sep="\t",quote=F,row.names=F,col.names=F)
write.table(ES_ASpliToGSEA[,c(2,3,1)],file=paste0(wd,"/int/GSEA/ASpli/SE.txt"),sep="\t",quote=F,row.names=F,col.names=F)
write.table(A3SS_ASpliToGSEA[,c(2,3,1)],file=paste0(wd,"/int/GSEA/ASpli/A3.txt"),sep="\t",quote=F,row.names=F,col.names=F)
write.table(A5SS_ASpliToGSEA[,c(2,3,1)],file=paste0(wd,"/int/GSEA/ASpli/A5.txt"),sep="\t",quote=F,row.names=F,col.names=F)
