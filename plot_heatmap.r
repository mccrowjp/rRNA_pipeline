#!/usr/bin/env Rscript
#
# plot_heatmap - heatmap of normalized OTU sample abundance
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

plot_heatmap = function(x, fontsize=1, ...) {
    xmat = as.matrix(x[,-(1:3)])
    sample = colnames(xmat)
    xmat = xmat + 0.1
    xmat = sweep(xmat,2,colSums(xmat),`/`)
    xmat.n = log(xmat, 10)
    heatmap(t(xmat.n), labCol="", scale="none", Colv=colSums(xmat), margins=c(2,10), col=colorRampPalette(c("white","red","darkred"))(64))
}

if(intable == 'NA' || outpdf == 'NA') {
    "Usage: [input OTU taxonomy counts table] [output pdf]"
} else {
    dat = read.table(intable, sep="\t", header=T, comment="")
    pdf(outpdf)
    plot_heatmap(dat)
    dev.off()
}
