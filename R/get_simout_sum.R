#!/usr/bin/env Rscript

## call with
## Rscript --vanilla get_simout_sum.R mysim nrep myfrag

args<-commandArgs(trailingOnly=TRUE)

if(length(args)<2) {
    stop("Require at least two arguments: mysim nrep myfrag")
}

## read simulation settings
## ------------------------

mysim<-args[1]
nrep<-as.integer(args[2])
if(length(args)>2){
    myfrag<-args[3]
    frag<-TRUE
}else{
    frag<-FALSE
}

if(frag==TRUE){
    rearrsout<-paste0(mysim,"_",myfrag,"_rearrs.RData")
    precrecout<-paste0(mysim,"_",myfrag,"_PrecRec.RData")
}else{
    rearrsout<-paste0(mysim,"_rearrs.RData")
    precrecout<-paste0(mysim,"_PrecRec.RData")
}


## rearrs
## ------

## get dimensions
if(frag==TRUE){
    myin<-paste(mysim,1,myfrag,"rearrs.txt",sep="_")
}else{
    myin<-paste(mysim,1,"rearrs.txt",sep="_")
}
tmp<-read.table(myin,header=TRUE,as.is=TRUE)

## set up output
rs.sum<-array(dim=c(nrep,ncol(tmp)-1),
              dimnames=list(as.character(1:nrep),colnames(tmp)[-1]))

for(i in 1:nrep){
    if(frag==TRUE){
        myin<-paste(mysim,i,myfrag,"rearrs.txt",sep="_")
    }else{
        myin<-paste(mysim,i,"rearrs.txt",sep="_")
    }
    rs.sum[i,]<-colSums(read.table(myin,header=TRUE)[,-1])
}

save(rs.sum,file=rearrsout)



## PrecRec
## -------

## get dimensions
if(frag==TRUE){
    myin<-paste(mysim,1,myfrag,"PrecRec.txt",sep="_")
}else{
    myin<-paste(mysim,1,"PrecRec.txt",sep="_")
}
tmp<-read.table(myin,header=TRUE,as.is=TRUE)

## set up output
pr.sum<-array(dim=c(nrow(tmp),ncol(tmp)-1,nrep),
              dimnames=list(tmp[,1],colnames(tmp)[-1],as.character(1:nrep)))

for(i in 1:nrep){
    if(frag==TRUE){
        myin<-paste(mysim,i,myfrag,"PrecRec.txt",sep="_")
    }else{
        myin<-paste(mysim,i,"PrecRec.txt",sep="_")
    }
    pr.sum[,,i]<-as.matrix(read.table(myin,header=TRUE)[,-1])
}

save(pr.sum,file=precrecout)

