#!/usr/bin/env Rscript
#
# plot_OTU_purity - scatter plot with marginal histograms of OTU purity
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

plotxy_marginhists = function(x, y, margins=c(5,5), border=1, h.axis=TRUE, h1.axis=h.axis, h2.axis=h.axis, h.b=20, h1.b=h.b, h2.b=h.b, h.col=8, h1.col=h.col, h2.col=h.col, log="", title="", ...) {
    lmat = matrix(c(2,4,1,3), ncol=2)
    layout(lmat, widths=c(4/5, 1/5), heights=c(1/5, 4/5))
    
    par(mar=c(0,0,0,0))
    plot(1:1,axes=F,type="n",xlab="",ylab="")
    text(1,1,title, font=2)
    
    df = data.frame(x,y)
    if(grepl("x", log)) { df = df[df[,1]>0,] }
    if(grepl("y", log)) { df = df[df[,2]>0,] }
    x = df[,1]
    y = df[,2]
    
    hx = hist(x, plot=F, breaks=h1.b)
    hy = hist(y, plot=F, breaks=h2.b)
    if(grepl("x", log)) { hx = hist(log(x,10), plot=F, breaks=h1.b) }
    if(grepl("y", log)) { hy = hist(log(y,10), plot=F, breaks=h2.b) }
    
    par(mar=c(0,margins[2],1,1))
    barplot(hx$counts, axes=F, ylim=c(0, max(hx$counts)), space=0, col=h1.col, border=border)
    if(h1.axis) { axis(2) }
    
    par(mar=c(margins[1],0,1,1))
    barplot(hy$counts, axes=F, xlim=c(0, max(hy$counts)), space=0, col=h2.col, border=border, horiz=TRUE)
    if(h2.axis) { axis(1) }
    
    par(mar=c(margins[1],margins[2],1,1))
    plot(x, y, log=log, ...)
}

if(intable == 'NA' || outpdf == 'NA') {
    "Usage: [input purity table] [output pdf] ([title])"
} else {
    if(intitle == 'NA') { intitle = "" }
    dat = read.table(intable, sep="\t", header=T, comment="")
    if(nrow(dat) > 0) {
        pdf(outpdf)
            plotxy_marginhists(dat[,3], dat[,5]*100, log="x", pch=20, col=4, ylim=c(0, 100), xlab="OTU reads", ylab="OTU purity (%)", title=intitle, h2.b=100, h1.b=50)
        dev.off()
    }
}
