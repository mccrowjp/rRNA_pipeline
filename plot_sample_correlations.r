#!/usr/bin/env Rscript
#
# plot_sample_correlations - plot hierarchical tree of sample correlations 
#
# Version: 0.2 (5/3/2016)
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

plot_sample_corr = function(x, fontsize=1, title="", ...) {
    xmat = as.matrix(x[,-(1:3)])
    xmat.h = hclust(as.dist(2-cor(sqrt(xmat))))
    plot(xmat.h, cex=0.7*fontsize, xlab="Sample", ylab="Correlation sqrt(reads)", sub="", main=title)
}

if(intable == 'NA' || outpdf == 'NA') {
    "Usage: [input OTU taxonomy counts table] [output pdf] ([title])"
} else {
    if(intitle == 'NA') { intitle = "" }
    dat = read.table(intable, sep="\t", header=T, comment="")
    pdf(outpdf)
    plot_sample_corr(dat, title=intitle)
    dev.off()
}
