#!/usr/bin/env Rscript
#
# plot_taxa_groups - plots of taxonomic group counts
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
intitle = paste("", args[3], sep="")

plot_heatmap = function(x, fontsize=1, min.taxa.frac=0, ...) {
    xmat = as.matrix(x[,-1])
    rownames(xmat) = x[,1]
    xmat = sweep(xmat,2,colSums(xmat),`/`)
    rsum = apply(xmat,1,sum)
    xmat = xmat[rsum >= min.taxa.frac,]
    rsum = apply(xmat,1,sum)
    heatmap(sqrt(xmat), scale="none", Rowv=rsum^2, margins=c(10,10), cex.axis=0.8*fontsize)
}

plot_taxa_groups = function(x, fontsize=1, min.taxa.frac=0, title="", ...) {
    xmat = as.matrix(x[,-1])
    xmat = sweep(xmat,2,colSums(xmat),`/`)
    rmean = apply(xmat, 1, mean)
    rsd = apply(xmat, 1, sd)
    par(mar=c(12,5,5,2))
    h = rmean[rmean >= min.taxa.frac]
    lab = x[rmean >= min.taxa.frac,1]
    se = (rsd/sqrt(nrow(xmat)))[rmean >= min.taxa.frac]
    maxh = max(h+se)
    bp = barplot(h, names.arg=lab, ylim=c(0,maxh), las=2, cex.axis=0.8, cex.names=0.8*fontsize, space=0, ylab="Fraction", main=title)
    arrows(bp, h+se, bp, h-se, angle=90, code=1, length=0)
}

if(intable == 'NA' || outpdf == 'NA') {
    "Usage: [input group counts table] [output pdf]"
} else {
    if(intitle == 'NA') { intitle = "" }
    dat = read.table(intable, sep="\t", header=T, comment="")
    pdf(outpdf)
    plot_taxa_groups(dat, title=intitle)
    plot_heatmap(dat, min.taxa.frac=0.005)
    dev.off()
}
