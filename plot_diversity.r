#!/usr/bin/env Rscript
#
# plot_diversity - bar plots of OTU richness and alpha-diversity metrics
#
# Version: 0.3 (5/19/2016)
#
# Part of rRNA_pipeline - FASTQ filtering, and swarm OTU classification of 16/18S barcodes
#
# Original version: 4/11/2016 John P. McCrow (jmccrow [at] jcvi.org)
# J. Craig Venter Institute (JCVI)
# La Jolla, CA USA
#
args <- commandArgs(TRUE)
intable = paste("", args[1], sep="")
outpdf = paste("", args[2], sep="")

plot_diversity = function(x, fontsize=1, ...) {
    xmat = as.matrix(x[,-(1:3)])
    sample = colnames(xmat)
    xmat.p = sweep(xmat,2,colSums(xmat),`/`)
    
    richness = apply(xmat, 2, function(x) length(x[x>0]))
    shannon = apply(xmat.p, 2, function(x) {y=x[x>0]; z=y*log(y); -sum(z)})
    simpson = apply(xmat.p, 2, function(x) 1-sum(x^2))

    par(mar=c(10,4,5,2))
    barplot(richness, las=2, main="Richness (total OTUs)", cex.names=0.8*fontsize, space=0)
    barplot(shannon, las=2, main="Shannon diversity", cex.names=0.8*fontsize, space=0)
    barplot(simpson, las=2, main="Simpson diversity", cex.names=0.8*fontsize, space=0)
}

if(intable == 'NA' || outpdf == 'NA') {
    "Usage: [input OTU taxonomy counts table] [output pdf]"
} else {
    dat = read.table(intable, sep="\t", header=T, comment="")
    pdf(outpdf)
    plot_diversity(dat)
    dev.off()
}
