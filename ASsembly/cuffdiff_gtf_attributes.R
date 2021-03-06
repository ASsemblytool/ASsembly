#!/usr/bin/env Rscript
'
USAGE: cuffdiff_gtf_attributes --input=<inputGTF> [--output=outputGTF] | --help

OPTIONS:
   -i --input=<inputGTF>        the inputGTF
   -o --output=[outputGTF]      the outputGTF

EXAMPLE:
   > wget ftp://ftp.ensembl.org/pub/release-82/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.82.gtf.gz
   > gunzip Saccharomyces_cerevisiae.R64-1-1.82.gtf.gz 
   > cuffdiff_gtf_attributes --input  Saccharomyces_cerevisiae.R64-1-1.82.gtf

PURPOSE:

   Provide p_id and tss_id attributes to a GTF formatted file,
   <input.gtf>.  Input file gets overwritten unless <output.gtf> is
   provided.

   These attributes are necessary for cuffdiff to perform all the
   differential splicing/coding/expression contrasts; if the
   attributes are missing, these contrasts are skipped.

   Output is grouped by gene and sorted on chromosome, genestart,
   gene_id, transcript_start, transcript_id.  This is gratuitous, but
   facilitates eyeball spot-check of results

NOTES:
   The documented way of adding these attributes is to create a
   .combined.gtf file using `cuffcompare`, but this method
   unfortunately (unnecessarily!?!) resets the gene_id and
   transcript_id to newly generated unique values (i.e.: gene_id
   "XLOC_000002"; transcript_id "TCONS_00000002").  Thus, this script,
   which preserves the ids.

   Annotating the p_id, required for cuffdiff differential _coding_
   analysis, depends upon the presence of CDS features in the GTF.
   Alas, there are GTF variants which mark the UTR but not the CDS;
   this are not accomodated for here.

VERSION:
   2015-11-13

DEPENDS:
   Developed and tested with R version 3.2.2 with current packages:
   docopt, data.table, rtracklayer and their dependencies.

AUTHOR:
   malcolm_cook@stowers.org (malcolm.cook@gmail.com)

' -> doc

suppressMessages({
    library(docopt)
    })

opts <- docopt(doc)

with(opts,{
    stopifnot(file.exists(input))
})


suppressMessages({
    library(rtracklayer)
    library(data.table)
})

setAs("GRanges", "data.table", function(from) {
    ## PURPOSE: allow to coerce Granges to data.table
    if (length(from) == 0L) {
        return(data.table())
    }
    as.data.table(as.data.frame(from)
                  ##,keep.rownames=TRUE
                  )
})

cuffdiff_gtf_attributes<-function(input.gtf,output.gtf) {
    gtf.dt<-as(import(input.gtf,'gtf'),'data.table') #,asRangedData=FALSE (has been deprecated)
    setkey(gtf.dt,gene_id,transcript_id,type,start) 
    gtf.dt[,tss:=ifelse('-' == strand,max(end),min(start))
           ## transcripts having the same tss will be assigned the
           ## same tss_id, below.
	  ,by=list(gene_id,transcript_id)]
    gtf.dt[,tss_id:=paste0('tss_',.GRP)
	  ,by=list(gene_id,tss)]
    gtf.dt[,cds_path:=.SD[type=='CDS',paste(start,end,sep='-',collapse='^')]
	   ## cds_path identifies the dna pre-image of a protein
	   ## (possibly in reverse order, but that is OK).  Trascripts
	   ## having the same cds_path will be assigned the same p_id
	   ## (below).
	  ,by=list(gene_id,transcript_id)]
    gtf.dt[cds_path!=''
          ,p_id:=paste0('p_',.GRP)
	  ,by=list(gene_id,transcript_id,cds_path)]
    ## Set up export order:
    gtf.dt[,transcript_start:=ifelse('-' == strand,max(end),min(start))
	  ,by=list(transcript_id)]
    gtf.dt[,gene_start:=ifelse('-' == strand,max(end),min(start))
	  ,by=list(gene_id)]
    setkey(gtf.dt,seqnames,gene_start,gene_id,transcript_start,transcript_id,source,start) # establish the order in which they will get exported.   Not required, but pretty.
    gtf.dt[,tss:=NULL]                  # So as not to be exported since it is longer needed.
    gtf.dt[,cds_path:=NULL]             # Likewise.
    gtf.dt[,transcript_start:=NULL]     # Likewise.
    gtf.dt[,gene_start:=NULL]           # Likewise.
    gtf.gr<-as(gtf.dt,'GRanges')
    export(gtf.gr,output.gtf,'gtf')
    output.gtf
}

with(opts,{
    if (is.null(output)) output<-input
    message(cuffdiff_gtf_attributes(input,output))
})
