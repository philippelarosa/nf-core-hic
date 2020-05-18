#!/usr/bin/env Rscript

library(HiTC)
library(rtracklayer)
library(optparse)

option_list <- list(make_option(c("-i", "--matrix"), type="character", default=NULL, help="Hi-C matrix file in HiC-Pro format."),
                    make_option(c("-b", "--bed"), type="character", default=NULL, help="HiC-C Bed file in HiC-Pro format"),
                    make_option(c("-o", "--odir"), type="character", default='./', help="Output directory", metavar="path"))

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

## Load Hi-C data
d <- importC(con=opt$matrix, xgi.bed=opt$bed, rm.trans=TRUE)
d <- HTClist(lapply(d, forceSymmetric))

## Remove empty matrix if any
d <- d[which(sapply(d, function(x){sum(intdata(x))})>0)]

logs <- lapply(d, function(x){
    xdata <- intdata(x)
    if (round(length(xdata@x)/length(xdata),3) <= 0.1)
        warning(seqlevels(x), ": data are too sparse ! increase the resolution !")
    invisible(NULL)
})

d <- sortSeqlevels(d)

## Call compartments
d.intra <- d[isIntraChrom(d)]
cc <- lapply(d.intra, pca.hic, npc=2, method="mean", asGRangesList=TRUE)
cc <- cc[which(!sapply(cc, is.null))]

## Export BEDgraph file for PC1 and PC2
if (length(cc) > 0){
    pc1export <- unlist(GRangesList(lapply(cc, function(k){k$PC1})))
    score(pc1export)[which(is.na(score(pc1export)))] <- 0
    tline <- new("GraphTrackLine", type="bedGraph", name=gsub(".mat[rix]*","",basename(in.mat)), description="A/B compartments (PC1)")
    export(cp1export, con=file.path(odir, gsub(".mat[rix]*","_pc1_cc.bedgraph",basename(in.mat))), format="bedGraph", trackLine=tline)
    
    pc2export <- unlist(GRangesList(lapply(cc, function(k){k$PC2})))
    score(pc2export)[which(is.na(score(pc2export)))] <- 0
    tline <- new("GraphTrackLine", type="bedGraph", name=gsub(".mat[rix]*","",basename(in.mat)), description="A/B compartments (PC2)")
    export(pc2export, con=file.path(odir, gsub(".mat[rix]*","_pc2_cc.bedgraph",basename(in.mat))), format="bedGraph", trackLine=tline)
}
