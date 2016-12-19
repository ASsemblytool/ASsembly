# ASsembly (v1.0): Assembly of Differential Alternative Splicing Methods 
# Authors: Carmen Bravo, Carlos de Lannoy, Emma Houben
# Email: assemblytoolkuleuven@gmail.com
# Github: https://github.com/ASsemblytool/ASsembly

# Script for running GSEA analysis

# Installing packages
if("devtools" %in% rownames(installed.packages()) == FALSE) {install.packages("devtools", repos = "http://cran.us.r-project.org")}
library(devtools)
if("fgsea" %in% rownames(installed.packages()) == FALSE) {install_github("ctlab/fgsea")}
library(fgsea)
if("ggplot2" %in% rownames(installed.packages()) == FALSE) {install.packages("ggplot2", repos = "http://cran.us.r-project.org")}
library(ggplot2)
if("optparse" %in% rownames(installed.packages()) == FALSE) {install.packages("optparse", repos = "http://cran.us.r-project.org")}
library("optparse")

# Parse arguments
option_list = list(
  make_option(c("-d", "--directory"), type="character", default=NULL, 
              help="Directory where rank files are found. Output will be produced here.", metavar="character"),
  make_option("--rgmt", type="character", default=NULL, 
              help="Reference gmt file containing true events", metavar="character"),
  make_option("--se", type="character", default=NULL, 
              help="Ranking file with skipping exon events", metavar="character"),
  make_option("--ri", type="character", default=NULL, 
              help="Ranking file with retained intron events", metavar="character"),
  make_option("--a3", type="character", default=NULL, 
              help="Ranking file with A3SS events", metavar="character"),
  make_option("--a5", type="character", default=NULL, 
              help="Ranking file with A5SS events", metavar="character"),
  make_option("--mx", type="character", default=NULL, 
              help="Ranking file with mutually exclusive exon events", metavar="character"),
  make_option("--af", type="character", default=NULL, 
              help="Ranking file with alternative first exon events", metavar="character"),
  make_option("--al", type="character", default=NULL, 
              help="Ranking file with alternative last exon events", metavar="character"),
  make_option(c("-m", "--method"), type="character", default=NULL, 
              help="Method which produced the rankings", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# GSEA analyses
gmt.file <- opt$rgmt
geneSet <- gmtPathways(gmt.file)
setwd(opt$directory)
method <- opt$method
df <- data.frame()
rows <- vector()
columns <- c("Total", "Leading Edge", "TE in LE", "Threshold", "ES score", "p-value")

dir.create("plots", recursive = T, showWarnings = F)
dir.create("LeadingEdge_genes", recursive = T, showWarnings = F)

if (!is.null(opt$se)&("SE" %in% names(geneSet))){
  rnk.file <- opt$se
  ranks <- read.table(rnk.file, header=TRUE, colClasses = c("character", "numeric"), col.names = c("ID", "pval"))
  ranks <- setNames(ranks$pval, ranks$ID)
  fgseaRes <- fgsea(geneSet, ranks, minSize=1, maxSize=1000, nperm=1000)
  All <- length(ranks)
  Allgenes <- names(ranks)
  LEgenes <- unlist(fgseaRes$leadingEdge[which(fgseaRes$pathway == "SE")])
  LE <- fgseaRes$size[which(fgseaRes$pathway == "SE")]
  TELE <- sum(LEgenes %in% geneSet$SE) 
  ES <- fgseaRes$ES[which(fgseaRes$pathway == "SE")]
  pval <- 1-ranks[LE]
  df <- rbind(df, c(All, LE, TELE, LE/All, ES, pval))
  rows <- c(rows, "SE")
  
  graphs <- paste("plots/SE", method, "GSEA", ".pdf", sep="")
  plotEnrichment(geneSet[["SE"]], ranks) + labs(title=paste(method, "SE", sep= " "))
  ggsave(graphs)
  
  write(names(ranks[1:LE]), "LeadingEdge_genes/LESEgenes.txt")
}

if (!is.null(opt$ri)&("RI" %in% names(geneSet))){
  rnk.file <- opt$ri
  ranks <- read.table(rnk.file, header=TRUE, colClasses = c("character", "numeric"), col.names = c("ID", "pval"))
  ranks <- setNames(ranks$pval, ranks$ID)

  fgseaRes <- fgsea(geneSet, ranks, minSize=1, maxSize=1000, nperm=1000)
  All <- length(ranks)
  Allgenes <- names(ranks)
  LEgenes <- unlist(fgseaRes$leadingEdge[which(fgseaRes$pathway == "RI")])
  LE <- fgseaRes$size[which(fgseaRes$pathway == "RI")]
  TELE <- sum(LEgenes %in% geneSet$RI) 
  ES <- fgseaRes$ES[which(fgseaRes$pathway == "RI")]
  pval <- 1-ranks[LE]
  df <- rbind(df, c(All, LE, TELE, LE/All, ES, pval))
  rows <- c(rows, "RI")
  
  graphs <- paste("plots/RI", method, "GSEA", ".pdf", sep="")
  plotEnrichment(geneSet[["RI"]], ranks) + labs(title=paste(method, "RI", sep= " "))
  ggsave(graphs)
  
  write(names(ranks[1:LE]), "LeadingEdge_genes/LERIgenes.txt")
}

if (!is.null(opt$a3)&("A3" %in% names(geneSet))){
  rnk.file <- opt$a3
  ranks <- read.table(rnk.file, header=TRUE, colClasses = c("character", "numeric"), col.names = c("ID", "pval"))
  ranks <- setNames(ranks$pval, ranks$ID)

  fgseaRes <- fgsea(geneSet, ranks, minSize=1, maxSize=1000, nperm=1000)
  All <- length(ranks)
  Allgenes <- names(ranks)
  LEgenes <- unlist(fgseaRes$leadingEdge[which(fgseaRes$pathway == "A3")])
  LE <- fgseaRes$size[which(fgseaRes$pathway == "A3")]
  TELE <- sum(LEgenes %in% geneSet$A3) 
  ES <- fgseaRes$ES[which(fgseaRes$pathway == "A3")]
  pval <- 1-ranks[LE]
  df <- rbind(df, c(All, LE, TELE, LE/All, ES, pval))
  rows <- c(rows, "A3")
  
  graphs <- paste("plots/A3", method, "GSEA", ".pdf", sep="")
  plotEnrichment(geneSet[["A3"]], ranks) + labs(title=paste(method, "A3", sep= " "))
  ggsave(graphs)
  
  write(names(ranks[1:LE]), "LeadingEdge_genes/LEA3genes.txt")
}

if (!is.null(opt$a5)&("A5" %in% names(geneSet))){
  rnk.file <- opt$a5
  ranks <- read.table(rnk.file, header=TRUE, colClasses = c("character", "numeric"), col.names = c("ID", "pval"))
  ranks <- setNames(ranks$pval, ranks$ID)

  fgseaRes <- fgsea(geneSet, ranks, minSize=1, maxSize=1000, nperm=1000)
  All <- length(ranks)
  Allgenes <- names(ranks)
  LEgenes <- unlist(fgseaRes$leadingEdge[which(fgseaRes$pathway == "A5")])
  LE <- fgseaRes$size[which(fgseaRes$pathway == "A5")]
  TELE <- sum(LEgenes %in% geneSet$A5) 
  ES <- fgseaRes$ES[which(fgseaRes$pathway == "A5")]
  pval <- 1-ranks[LE]
  df <- rbind(df, c(All, LE, TELE, LE/All, ES, pval))
  rows <- c(rows, "A5")
  
  graphs <- paste("plots/A5", method, "GSEA", ".pdf", sep="")
  plotEnrichment(geneSet[["A5"]], ranks) + labs(title=paste(method, "A5", sep= " "))
  ggsave(graphs)
  
  write(names(ranks[1:LE]), "LeadingEdge_genes/LEA5genes.txt")
}

if (!is.null(opt$mx)&("MX" %in% names(geneSet))){
  rnk.file <- opt$mx
  ranks <- read.table(rnk.file, header=TRUE, colClasses = c("character", "numeric"), col.names = c("ID", "pval"))
  ranks <- setNames(ranks$pval, ranks$ID)
  
  fgseaRes <- fgsea(geneSet, ranks, minSize=1, maxSize=1000, nperm=1000)
  All <- length(ranks)
  Allgenes <- names(ranks)
  LEgenes <- unlist(fgseaRes$leadingEdge[which(fgseaRes$pathway == "MX")])
  LE <- fgseaRes$size[which(fgseaRes$pathway == "MX")]
  TELE <- sum(LEgenes %in% geneSet$MX) 
  ES <- fgseaRes$ES[which(fgseaRes$pathway == "MX")]
  pval <- 1-ranks[LE]
  df <- rbind(df, c(All, LE, TELE, LE/All,ES, pval))
  rows <- c(rows, "MX")
  
  graphs <- paste("plots/MX", method, "GSEA", ".pdf", sep="")
  plotEnrichment(geneSet[["MX"]], ranks) + labs(title=paste(method, "MX", sep= " "))
  ggsave(graphs)
  
  write(names(ranks[1:LE]), "LeadingEdge_genes/LEMXgenes.txt")
}

if (!is.null(opt$al)&("AL" %in% names(geneSet))){
  rnk.file <- opt$al
  ranks <- read.table(rnk.file, header=TRUE, colClasses = c("character", "numeric"), col.names = c("ID", "pval"))
  ranks <- setNames(ranks$pval, ranks$ID)
  
  fgseaRes <- fgsea(geneSet, ranks, minSize=1, maxSize=1000, nperm=1000)
  All <- length(ranks)
  Allgenes <- names(ranks)
  LEgenes <- unlist(fgseaRes$leadingEdge[which(fgseaRes$pathway == "AL")])
  LE <- fgseaRes$size[which(fgseaRes$pathway == "AL")]
  TELE <- sum(LEgenes %in% geneSet$AL) 
  ES <- fgseaRes$ES[which(fgseaRes$pathway == "AL")]
  pval <- 1-ranks[LE]
  df <- rbind(df, c(All, LE, TELE, LE/All, ES, pval))
  rows <- c(rows, "AL")
  
  graphs <- paste("plots/AL", method, "GSEA", ".pdf", sep="")
  plotEnrichment(geneSet[["AL"]], ranks) + labs(title=paste(method, "AL", sep= " "))
  ggsave(graphs)
  
  write(names(ranks[1:LE]), "LeadingEdge_genes/LEALgenes.txt")
}

if (!is.null(opt$af)&("AF" %in% names(geneSet))){
  rnk.file <- opt$af
  ranks <- read.table(rnk.file, header=TRUE, colClasses = c("character", "numeric"), col.names = c("ID", "pval"))
  ranks <- setNames(ranks$pval, ranks$ID)
  
  fgseaRes <- fgsea(geneSet, ranks, minSize=1, maxSize=1000, nperm=1000)
  All <- length(ranks)
  Allgenes <- names(ranks)
  LEgenes <- unlist(fgseaRes$leadingEdge[which(fgseaRes$pathway == "AF")])
  LE <- fgseaRes$size[which(fgseaRes$pathway == "AF")]
  TELE <- sum(LEgenes %in% geneSet$AF) 
  ES <- fgseaRes$ES[which(fgseaRes$pathway == "AF")]
  pval <- 1-ranks[LE]
  df <- rbind(df, c(All, LE, TELE, LE/All, ES, pval))
  rows <- c(rows, "AF")
  
  graphs <- paste("plots/AF", method, "GSEA", ".pdf", sep="")
  plotEnrichment(geneSet[["AF"]], ranks) + labs(title=paste(method, "AF", sep= " "))
  ggsave(graphs)
  
  write(names(ranks[1:LE]), "LeadingEdge_genes/LEAFgenes.txt")
}

rownames(df) <- rows
colnames(df) <- columns
write.table(df, paste(method, "Statistics", ".txt", sep=""), sep="\t", col.names = NA)

