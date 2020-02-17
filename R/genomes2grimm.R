#!/usr/bin/env Rscript

## call with
## Rscript --vanilla genomes2grimm.R simout/sim1/sim1_1.RData outpfx
## [.RData containing focalgenome and compgenome]
## or
## Rscript --vanilla genomes2grimm.R genome1 genome2 outpfx
## [one focalgenome and one compgenome, or two focalgenomes,
##  where the latter will be treated as compgenome]

## Example:
## Rscript --vanilla ./R/genomes2grimm.R ./simout/sim1/sim1_1.RData sim1_1

## ------------------------------------------------------------------------
## convert two genome files to GRIMM format
## [GRIMM format == UniMoG format with chromosome ends '$' instead of '|']
## ------------------------------------------------------------------------

args<-commandArgs(trailingOnly=TRUE)

if(length(args)<2) {
    stop("Require >=2 arguments: genomes.RData outpfx or genome1 genome2 outpfx")
}else if(length(args)==2){
    if(!grepl("\\.RData$",args[1])){
        stop("Require RData with one argument")
    }
    load(args[1])
    if(!exists("focalgenome")){
        stop("Couldn't load focalgenome")
    }
    if(!exists("compgenome")){
        stop("Couldn't load compgenome")
    }
    outprefix<-args[2]
}else if(length(args)==3){
    genome<-read.table(file=args[1],header=TRUE,as.is=TRUE)
    if(sum(is.element(c("scaff","marker","strand"),colnames(genome)))==3){
        focalgenome<-genome
    }else if(sum(is.element(c("car","marker","orientation"),
                            colnames(genome)))==3){
        compgenome<-genome
    }else{
        stop("Coudn't identify first genome")
    }
    genome<-read.table(file=args[2],header=TRUE,as.is=TRUE)
    if(sum(is.element(c("scaff","marker","strand"),colnames(genome)))==3){
        if(!exists("focalgenome")){
            focalgenome<-genome
        }else{
            ## transform
            compgenome<-data.frame(marker=genome$marker,
                                   orientation=genome$strand,
                                   car=genome$scaff,
                                   stringsAsFactors=FALSE)
        }
    }
    if(sum(is.element(c("car","marker","orientation"),colnames(genome)))==3){
        if(!exists("compgenome")){
            compgenome<-genome
        }else{
            stop("Require one focalgenome as input")
        }
    }else{
        stop("Coudn't identify second genome")
    }
    outprefix<-args[3]
}



## subset files so that they contain common set of markers
keep1<-focalgenome$marker %in% compgenome$marker
keep2<-compgenome$marker %in% focalgenome$marker
focalgenome<-focalgenome[keep1,]
compgenome<-compgenome[keep2,]

## need to re-assign new marker names to be integers from 1:ngenes
## (class or marker can be character or integer)
ids1<-match(focalgenome$marker,compgenome$marker)
ids2<-1:nrow(compgenome)
focalgenome$marker<-ids1
compgenome$marker<-ids2

## replace strand orientation "+" by ""
focalgenome$strand<-gsub("+","",focalgenome$strand,fixed=TRUE)
compgenome$orientation<-gsub("+","",compgenome$orientation,fixed=TRUE)


outfile<-">compgenome"
for (s in unique(compgenome$car)){
    pos<-which(compgenome$car==s)
    if(length(pos)==0){
        next
    }
    line<-paste0(compgenome$orientation[pos],
                 compgenome$marker[pos],collapse=" ")
    outfile<-c(outfile,paste(line,"$",sep=" "))
}

outfile<-c(outfile,">focalgenome")
for (s in unique(focalgenome$scaff)){
    pos<-which(focalgenome$scaff==s)
    if(length(pos)==0){
        next
    }
    line<-paste0(focalgenome$strand[pos],
                 focalgenome$marker[pos],collapse=" ")
    outfile<-c(outfile,paste(line,"$",sep=" "))
}


write.table(outfile,file=paste0(outprefix,"_grimm.txt"),
            quote=FALSE,col.names=FALSE,row.names=FALSE)


## ------------------------------------------------------------------------
