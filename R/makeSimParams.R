#!/usr/bin/env Rscript

## --------------------------------------------------------------------
## generate RData files with simulations settings
## --------------------------------------------------------------------

myoutdir<-"~/Rearrangements/Simulations/simparams"


## plot size and prob parameters for negative binomial distribution
## ----------------------------------------------------------------


## hist(rnbinom(1000,20,0.92)) ## small rearrs, median ~2, max ~13
## hist(rnbinom(1000,15,0.6)) ## medium rearrs, median ~10, max ~34
## hist(rnbinom(1000,10,0.325)) ## large, median ~20, max ~ 62
## hist(rnbinom(1000,6,0.123)) ## xlarge rearrs, median ~40, max ~163


## --------------------------------------------------------------------
## DEFINE REARRANGEMENT SIZES
## --------------------------------------------------------------------
sm<-c(20,0.92)
md<-c(15,0.6)
lg<-c(10,0.325)
xl<-c(6,0.123)



## --------------------------------------------------------------------
## SIM 1
## --------------------------------------------------------------------
## (medium number small-sized intra-chr and small-sized inter-chr)

## simulation name
mysim<-"sim1"

## number of chromosomes and genes per chromosome
ngenes<-c(1000,1000,1000,1000,1000)

## number of events per type
nrearr<-list(inv=20,synmov=20,nonsynmov=20,
             retra=0,fus=0,fis=0)

## set rearrangement size by size and prob of negative binomial distribution
##   +1 gene to not have zero sizes
## not applicable for retra, fus, fis
rearrsize<-list(inv=sm[1],synmov=sm[1],nonsynmov=sm[1],
                retra=NA,fus=NA,fis=NA)
rearrprob<-list(inv=sm[2],synmov=sm[2],nonsynmov=sm[2],
                retra=NA,fus=NA,fis=NA)

## probability of breakpoints along chromosome
betashape1<-1
betashape2<-1

save(mysim,ngenes,nrearr,rearrsize,rearrprob,betashape1,betashape2,
     file=paste0(myoutdir,"/",mysim,"_params.RData"))



## --------------------------------------------------------------------
## SIM 2
## --------------------------------------------------------------------
## (medium number medium-sized intra-chr and medium-sized inter-chr)

## simulation name
mysim<-"sim2"

## number of chromosomes and genes per chromosome
ngenes<-c(1000,1000,1000,1000,1000)

## number of events per type
nrearr<-list(inv=20,synmov=20,nonsynmov=20,
             retra=0,fus=0,fis=0)

## set rearrangement size by size and prob of negative binomial distribution
##   +1 gene to not have zero sizes
## not applicable for retra, fus, fis
rearrsize<-list(inv=md[1],synmov=md[1],nonsynmov=md[1],
                retra=NA,fus=NA,fis=NA)
rearrprob<-list(inv=md[2],synmov=md[2],nonsynmov=md[2],
                retra=NA,fus=NA,fis=NA)

## probability of breakpoints along chromosome
betashape1<-1
betashape2<-1

save(mysim,ngenes,nrearr,rearrsize,rearrprob,betashape1,betashape2,
     file=paste0(myoutdir,"/",mysim,"_params.RData"))


## --------------------------------------------------------------------
## SIM 3
## --------------------------------------------------------------------
## (medium number large-sized intra-chr and large-sized inter-chr)

## simulation name
mysim<-"sim3"

## number of chromosomes and genes per chromosome
ngenes<-c(1000,1000,1000,1000,1000)

## number of events per type
nrearr<-list(inv=20,synmov=20,nonsynmov=20,
             retra=0,fus=0,fis=0)

## set rearrangement size by size and prob of negative binomial distribution
##   +1 gene to not have zero sizes
## not applicable for retra, fus, fis
rearrsize<-list(inv=lg[1],synmov=lg[1],nonsynmov=lg[1],
                retra=NA,fus=NA,fis=NA)
rearrprob<-list(inv=lg[2],synmov=lg[2],nonsynmov=lg[2],
                retra=NA,fus=NA,fis=NA)

## probability of breakpoints along chromosome
betashape1<-1
betashape2<-1

save(mysim,ngenes,nrearr,rearrsize,rearrprob,betashape1,betashape2,
     file=paste0(myoutdir,"/",mysim,"_params.RData"))



## --------------------------------------------------------------------
## SIM 4
## --------------------------------------------------------------------
## (medium number xlarge-sized intra-chr and xlarge-sized inter-chr)

## simulation name
mysim<-"sim4"

## number of chromosomes and genes per chromosome
ngenes<-c(1000,1000,1000,1000,1000)

## number of events per type
nrearr<-list(inv=20,synmov=20,nonsynmov=20,
             retra=0,fus=0,fis=0)

## set rearrangement size by size and prob of negative binomial distribution
##   +1 gene to not have zero sizes
## not applicable for retra, fus, fis
rearrsize<-list(inv=xl[1],synmov=xl[1],nonsynmov=xl[1],
                retra=NA,fus=NA,fis=NA)
rearrprob<-list(inv=xl[2],synmov=xl[2],nonsynmov=xl[2],
                retra=NA,fus=NA,fis=NA)

## probability of breakpoints along chromosome
betashape1<-1
betashape2<-1

save(mysim,ngenes,nrearr,rearrsize,rearrprob,betashape1,betashape2,
     file=paste0(myoutdir,"/",mysim,"_params.RData"))



## --------------------------------------------------------------------
## SIM 5
## --------------------------------------------------------------------
## (large number small-sized intra-chr and small-sized inter-chr)

## simulation name
mysim<-"sim5"

## number of chromosomes and genes per chromosome
ngenes<-c(1000,1000,1000,1000,1000)

## number of events per type
nrearr<-list(inv=40,synmov=40,nonsynmov=40,
             retra=0,fus=0,fis=0)

## set rearrangement size by size and prob of negative binomial distribution
##   +1 gene to not have zero sizes
## not applicable for retra, fus, fis
rearrsize<-list(inv=sm[1],synmov=sm[1],nonsynmov=sm[1],
                retra=NA,fus=NA,fis=NA)
rearrprob<-list(inv=sm[2],synmov=sm[2],nonsynmov=sm[2],
                retra=NA,fus=NA,fis=NA)

## probability of breakpoints along chromosome
betashape1<-1
betashape2<-1

save(mysim,ngenes,nrearr,rearrsize,rearrprob,betashape1,betashape2,
     file=paste0(myoutdir,"/",mysim,"_params.RData"))



## --------------------------------------------------------------------
## SIM 6
## --------------------------------------------------------------------
## (large number medium-sized intra-chr and medium-sized inter-chr)

## simulation name
mysim<-"sim6"

## number of chromosomes and genes per chromosome
ngenes<-c(1000,1000,1000,1000,1000)

## number of events per type
nrearr<-list(inv=40,synmov=40,nonsynmov=40,
             retra=0,fus=0,fis=0)

## set rearrangement size by size and prob of negative binomial distribution
##   +1 gene to not have zero sizes
## not applicable for retra, fus, fis
rearrsize<-list(inv=md[1],synmov=md[1],nonsynmov=md[1],
                retra=NA,fus=NA,fis=NA)
rearrprob<-list(inv=md[2],synmov=md[2],nonsynmov=md[2],
                retra=NA,fus=NA,fis=NA)

## probability of breakpoints along chromosome
betashape1<-1
betashape2<-1

save(mysim,ngenes,nrearr,rearrsize,rearrprob,betashape1,betashape2,
     file=paste0(myoutdir,"/",mysim,"_params.RData"))


## --------------------------------------------------------------------
## SIM 7
## --------------------------------------------------------------------
## (large number large-sized intra-chr and large-sized inter-chr)

## simulation name
mysim<-"sim7"

## number of chromosomes and genes per chromosome
ngenes<-c(1000,1000,1000,1000,1000)

## number of events per type
nrearr<-list(inv=40,synmov=40,nonsynmov=40,
             retra=0,fus=0,fis=0)

## set rearrangement size by size and prob of negative binomial distribution
##   +1 gene to not have zero sizes
## not applicable for retra, fus, fis
rearrsize<-list(inv=lg[1],synmov=lg[1],nonsynmov=lg[1],
                retra=NA,fus=NA,fis=NA)
rearrprob<-list(inv=lg[2],synmov=lg[2],nonsynmov=lg[2],
                retra=NA,fus=NA,fis=NA)

## probability of breakpoints along chromosome
betashape1<-1
betashape2<-1

save(mysim,ngenes,nrearr,rearrsize,rearrprob,betashape1,betashape2,
     file=paste0(myoutdir,"/",mysim,"_params.RData"))



## --------------------------------------------------------------------
## SIM 8
## --------------------------------------------------------------------
## (large number xlarge-sized intra-chr and xlarge-sized inter-chr)

## simulation name
mysim<-"sim8"

## number of chromosomes and genes per chromosome
ngenes<-c(1000,1000,1000,1000,1000)

## number of events per type
nrearr<-list(inv=40,synmov=40,nonsynmov=40,
             retra=0,fus=0,fis=0)

## set rearrangement size by size and prob of negative binomial distribution
##   +1 gene to not have zero sizes
## not applicable for retra, fus, fis
rearrsize<-list(inv=xl[1],synmov=xl[1],nonsynmov=xl[1],
                retra=NA,fus=NA,fis=NA)
rearrprob<-list(inv=xl[2],synmov=xl[2],nonsynmov=xl[2],
                retra=NA,fus=NA,fis=NA)

## probability of breakpoints along chromosome
betashape1<-1
betashape2<-1

save(mysim,ngenes,nrearr,rearrsize,rearrprob,betashape1,betashape2,
     file=paste0(myoutdir,"/",mysim,"_params.RData"))



## --------------------------------------------------------------------
## SIM 9
## --------------------------------------------------------------------
## (small number small-sized intra-chr and small-sized inter-chr)

## simulation name
mysim<-"sim9"

## number of chromosomes and genes per chromosome
ngenes<-c(1000,1000,1000,1000,1000)

## number of events per type
nrearr<-list(inv=10,synmov=10,nonsynmov=10,
             retra=0,fus=0,fis=0)

## set rearrangement size by size and prob of negative binomial distribution
##   +1 gene to not have zero sizes
## not applicable for retra, fus, fis
rearrsize<-list(inv=sm[1],synmov=sm[1],nonsynmov=sm[1],
                retra=NA,fus=NA,fis=NA)
rearrprob<-list(inv=sm[2],synmov=sm[2],nonsynmov=sm[2],
                retra=NA,fus=NA,fis=NA)

## probability of breakpoints along chromosome
betashape1<-1
betashape2<-1

save(mysim,ngenes,nrearr,rearrsize,rearrprob,betashape1,betashape2,
     file=paste0(myoutdir,"/",mysim,"_params.RData"))



## --------------------------------------------------------------------
## SIM 10
## --------------------------------------------------------------------
## (small number medium-sized intra-chr and medium-sized inter-chr)

## simulation name
mysim<-"sim10"

## number of chromosomes and genes per chromosome
ngenes<-c(1000,1000,1000,1000,1000)

## number of events per type
nrearr<-list(inv=10,synmov=10,nonsynmov=10,
             retra=0,fus=0,fis=0)

## set rearrangement size by size and prob of negative binomial distribution
##   +1 gene to not have zero sizes
## not applicable for retra, fus, fis
rearrsize<-list(inv=md[1],synmov=md[1],nonsynmov=md[1],
                retra=NA,fus=NA,fis=NA)
rearrprob<-list(inv=md[2],synmov=md[2],nonsynmov=md[2],
                retra=NA,fus=NA,fis=NA)

## probability of breakpoints along chromosome
betashape1<-1
betashape2<-1

save(mysim,ngenes,nrearr,rearrsize,rearrprob,betashape1,betashape2,
     file=paste0(myoutdir,"/",mysim,"_params.RData"))


## --------------------------------------------------------------------
## SIM 11
## --------------------------------------------------------------------
## (small number large-sized intra-chr and large-sized inter-chr)

## simulation name
mysim<-"sim11"

## number of chromosomes and genes per chromosome
ngenes<-c(1000,1000,1000,1000,1000)

## number of events per type
nrearr<-list(inv=10,synmov=10,nonsynmov=10,
             retra=0,fus=0,fis=0)

## set rearrangement size by size and prob of negative binomial distribution
##   +1 gene to not have zero sizes
## not applicable for retra, fus, fis
rearrsize<-list(inv=lg[1],synmov=lg[1],nonsynmov=lg[1],
                retra=NA,fus=NA,fis=NA)
rearrprob<-list(inv=lg[2],synmov=lg[2],nonsynmov=lg[2],
                retra=NA,fus=NA,fis=NA)

## probability of breakpoints along chromosome
betashape1<-1
betashape2<-1

save(mysim,ngenes,nrearr,rearrsize,rearrprob,betashape1,betashape2,
     file=paste0(myoutdir,"/",mysim,"_params.RData"))



## --------------------------------------------------------------------
## SIM 12
## --------------------------------------------------------------------
## (small number xlarge-sized intra-chr and xlarge-sized inter-chr)

## simulation name
mysim<-"sim12"

## number of chromosomes and genes per chromosome
ngenes<-c(1000,1000,1000,1000,1000)

## number of events per type
nrearr<-list(inv=10,synmov=10,nonsynmov=10,
             retra=0,fus=0,fis=0)

## set rearrangement size by size and prob of negative binomial distribution
##   +1 gene to not have zero sizes
## not applicable for retra, fus, fis
rearrsize<-list(inv=xl[1],synmov=xl[1],nonsynmov=xl[1],
                retra=NA,fus=NA,fis=NA)
rearrprob<-list(inv=xl[2],synmov=xl[2],nonsynmov=xl[2],
                retra=NA,fus=NA,fis=NA)

## probability of breakpoints along chromosome
betashape1<-1
betashape2<-1

save(mysim,ngenes,nrearr,rearrsize,rearrprob,betashape1,betashape2,
     file=paste0(myoutdir,"/",mysim,"_params.RData"))



## --------------------------------------------------------------------
## SIM 13
## --------------------------------------------------------------------
## (medium number small-sized intra-chr and small-sized inter-chr)
## (plus 1 reciprocal translocation)

## simulation name
mysim<-"sim13"

## number of chromosomes and genes per chromosome
ngenes<-c(1000,1000,1000,1000,1000)

## number of events per type
nrearr<-list(inv=20,synmov=20,nonsynmov=20,
             retra=1,fus=0,fis=0)

## set rearrangement size by size and prob of negative binomial distribution
##   +1 gene to not have zero sizes
## not applicable for retra, fus, fis
rearrsize<-list(inv=sm[1],synmov=sm[1],nonsynmov=sm[1],
                retra=NA,fus=NA,fis=NA)
rearrprob<-list(inv=sm[2],synmov=sm[2],nonsynmov=sm[2],
                retra=NA,fus=NA,fis=NA)

## probability of breakpoints along chromosome
betashape1<-1
betashape2<-1

save(mysim,ngenes,nrearr,rearrsize,rearrprob,betashape1,betashape2,
     file=paste0(myoutdir,"/",mysim,"_params.RData"))



## --------------------------------------------------------------------
## SIM 14
## --------------------------------------------------------------------
## (medium number medium-sized intra-chr and medium-sized inter-chr)
## (plus 1 reciprocal translocation)

## simulation name
mysim<-"sim14"

## number of chromosomes and genes per chromosome
ngenes<-c(1000,1000,1000,1000,1000)

## number of events per type
nrearr<-list(inv=20,synmov=20,nonsynmov=20,
             retra=1,fus=0,fis=0)

## set rearrangement size by size and prob of negative binomial distribution
##   +1 gene to not have zero sizes
## not applicable for retra, fus, fis
rearrsize<-list(inv=md[1],synmov=md[1],nonsynmov=md[1],
                retra=NA,fus=NA,fis=NA)
rearrprob<-list(inv=md[2],synmov=md[2],nonsynmov=md[2],
                retra=NA,fus=NA,fis=NA)

## probability of breakpoints along chromosome
betashape1<-1
betashape2<-1

save(mysim,ngenes,nrearr,rearrsize,rearrprob,betashape1,betashape2,
     file=paste0(myoutdir,"/",mysim,"_params.RData"))


## --------------------------------------------------------------------
## SIM 15
## --------------------------------------------------------------------
## (medium number large-sized intra-chr and large-sized inter-chr)
## (plus 1 reciprocal translocation)

## simulation name
mysim<-"sim15"

## number of chromosomes and genes per chromosome
ngenes<-c(1000,1000,1000,1000,1000)

## number of events per type
nrearr<-list(inv=20,synmov=20,nonsynmov=20,
             retra=1,fus=0,fis=0)

## set rearrangement size by size and prob of negative binomial distribution
##   +1 gene to not have zero sizes
## not applicable for retra, fus, fis
rearrsize<-list(inv=lg[1],synmov=lg[1],nonsynmov=lg[1],
                retra=NA,fus=NA,fis=NA)
rearrprob<-list(inv=lg[2],synmov=lg[2],nonsynmov=lg[2],
                retra=NA,fus=NA,fis=NA)

## probability of breakpoints along chromosome
betashape1<-1
betashape2<-1

save(mysim,ngenes,nrearr,rearrsize,rearrprob,betashape1,betashape2,
     file=paste0(myoutdir,"/",mysim,"_params.RData"))



## --------------------------------------------------------------------
## SIM 16
## --------------------------------------------------------------------
## (medium number xlarge-sized intra-chr and xlarge-sized inter-chr)
## (plus 1 reciprocal translocation)

## simulation name
mysim<-"sim16"

## number of chromosomes and genes per chromosome
ngenes<-c(1000,1000,1000,1000,1000)

## number of events per type
nrearr<-list(inv=20,synmov=20,nonsynmov=20,
             retra=1,fus=0,fis=0)

## set rearrangement size by size and prob of negative binomial distribution
##   +1 gene to not have zero sizes
## not applicable for retra, fus, fis
rearrsize<-list(inv=xl[1],synmov=xl[1],nonsynmov=xl[1],
                retra=NA,fus=NA,fis=NA)
rearrprob<-list(inv=xl[2],synmov=xl[2],nonsynmov=xl[2],
                retra=NA,fus=NA,fis=NA)

## probability of breakpoints along chromosome
betashape1<-1
betashape2<-1

save(mysim,ngenes,nrearr,rearrsize,rearrprob,betashape1,betashape2,
     file=paste0(myoutdir,"/",mysim,"_params.RData"))



## --------------------------------------------------------------------
## SIM 17
## --------------------------------------------------------------------
## (medium number small-sized intra-chr and small-sized inter-chr)
## (plus 1 reciprocal translocation, 1 fission, 1 fusion)

## simulation name
mysim<-"sim17"

## number of chromosomes and genes per chromosome
ngenes<-c(1000,1000,1000,1000,1000)

## number of events per type
nrearr<-list(inv=20,synmov=20,nonsynmov=20,
             retra=1,fus=1,fis=1)

## set rearrangement size by size and prob of negative binomial distribution
##   +1 gene to not have zero sizes
## not applicable for retra, fus, fis
rearrsize<-list(inv=sm[1],synmov=sm[1],nonsynmov=sm[1],
                retra=NA,fus=NA,fis=NA)
rearrprob<-list(inv=sm[2],synmov=sm[2],nonsynmov=sm[2],
                retra=NA,fus=NA,fis=NA)

## probability of breakpoints along chromosome
betashape1<-1
betashape2<-1

save(mysim,ngenes,nrearr,rearrsize,rearrprob,betashape1,betashape2,
     file=paste0(myoutdir,"/",mysim,"_params.RData"))



## --------------------------------------------------------------------
## SIM 18
## --------------------------------------------------------------------
## (medium number medium-sized intra-chr and medium-sized inter-chr)
## (plus 1 reciprocal translocation, 1 fission, 1 fusion)

## simulation name
mysim<-"sim18"

## number of chromosomes and genes per chromosome
ngenes<-c(1000,1000,1000,1000,1000)

## number of events per type
nrearr<-list(inv=20,synmov=20,nonsynmov=20,
             retra=1,fus=1,fis=1)

## set rearrangement size by size and prob of negative binomial distribution
##   +1 gene to not have zero sizes
## not applicable for retra, fus, fis
rearrsize<-list(inv=md[1],synmov=md[1],nonsynmov=md[1],
                retra=NA,fus=NA,fis=NA)
rearrprob<-list(inv=md[2],synmov=md[2],nonsynmov=md[2],
                retra=NA,fus=NA,fis=NA)

## probability of breakpoints along chromosome
betashape1<-1
betashape2<-1

save(mysim,ngenes,nrearr,rearrsize,rearrprob,betashape1,betashape2,
     file=paste0(myoutdir,"/",mysim,"_params.RData"))


## --------------------------------------------------------------------
## SIM 19
## --------------------------------------------------------------------
## (medium number large-sized intra-chr and large-sized inter-chr)
## (plus 1 reciprocal translocation, 1 fission, 1 fusion)

## simulation name
mysim<-"sim19"

## number of chromosomes and genes per chromosome
ngenes<-c(1000,1000,1000,1000,1000)

## number of events per type
nrearr<-list(inv=20,synmov=20,nonsynmov=20,
             retra=1,fus=1,fis=1)

## set rearrangement size by size and prob of negative binomial distribution
##   +1 gene to not have zero sizes
## not applicable for retra, fus, fis
rearrsize<-list(inv=lg[1],synmov=lg[1],nonsynmov=lg[1],
                retra=NA,fus=NA,fis=NA)
rearrprob<-list(inv=lg[2],synmov=lg[2],nonsynmov=lg[2],
                retra=NA,fus=NA,fis=NA)

## probability of breakpoints along chromosome
betashape1<-1
betashape2<-1

save(mysim,ngenes,nrearr,rearrsize,rearrprob,betashape1,betashape2,
     file=paste0(myoutdir,"/",mysim,"_params.RData"))



## --------------------------------------------------------------------
## SIM 20
## --------------------------------------------------------------------
## (medium number xlarge-sized intra-chr and xlarge-sized inter-chr)
## (plus 1 reciprocal translocation, 1 fission, 1 fusion)

## simulation name
mysim<-"sim20"

## number of chromosomes and genes per chromosome
ngenes<-c(1000,1000,1000,1000,1000)

## number of events per type
nrearr<-list(inv=20,synmov=20,nonsynmov=20,
             retra=1,fus=1,fis=1)

## set rearrangement size by size and prob of negative binomial distribution
##   +1 gene to not have zero sizes
## not applicable for retra, fus, fis
rearrsize<-list(inv=xl[1],synmov=xl[1],nonsynmov=xl[1],
                retra=NA,fus=NA,fis=NA)
rearrprob<-list(inv=xl[2],synmov=xl[2],nonsynmov=xl[2],
                retra=NA,fus=NA,fis=NA)

## probability of breakpoints along chromosome
betashape1<-1
betashape2<-1

save(mysim,ngenes,nrearr,rearrsize,rearrprob,betashape1,betashape2,
     file=paste0(myoutdir,"/",mysim,"_params.RData"))



## --------------------------------------------------------------------
## SIM 21
## --------------------------------------------------------------------
## (medium number small-sized intra-chr and small-sized inter-chr)
## (plus 1 fission, 1 fusion)

## simulation name
mysim<-"sim21"

## number of chromosomes and genes per chromosome
ngenes<-c(1000,1000,1000,1000,1000)

## number of events per type
nrearr<-list(inv=20,synmov=20,nonsynmov=20,
             retra=0,fus=1,fis=1)

## set rearrangement size by size and prob of negative binomial distribution
##   +1 gene to not have zero sizes
## not applicable for retra, fus, fis
rearrsize<-list(inv=sm[1],synmov=sm[1],nonsynmov=sm[1],
                retra=NA,fus=NA,fis=NA)
rearrprob<-list(inv=sm[2],synmov=sm[2],nonsynmov=sm[2],
                retra=NA,fus=NA,fis=NA)

## probability of breakpoints along chromosome
betashape1<-1
betashape2<-1

save(mysim,ngenes,nrearr,rearrsize,rearrprob,betashape1,betashape2,
     file=paste0(myoutdir,"/",mysim,"_params.RData"))



## --------------------------------------------------------------------
## SIM 22
## --------------------------------------------------------------------
## (medium number medium-sized intra-chr and medium-sized inter-chr)
## (plus 1 fission, 1 fusion)

## simulation name
mysim<-"sim22"

## number of chromosomes and genes per chromosome
ngenes<-c(1000,1000,1000,1000,1000)

## number of events per type
nrearr<-list(inv=20,synmov=20,nonsynmov=20,
             retra=0,fus=1,fis=1)

## set rearrangement size by size and prob of negative binomial distribution
##   +1 gene to not have zero sizes
## not applicable for retra, fus, fis
rearrsize<-list(inv=md[1],synmov=md[1],nonsynmov=md[1],
                retra=NA,fus=NA,fis=NA)
rearrprob<-list(inv=md[2],synmov=md[2],nonsynmov=md[2],
                retra=NA,fus=NA,fis=NA)

## probability of breakpoints along chromosome
betashape1<-1
betashape2<-1

save(mysim,ngenes,nrearr,rearrsize,rearrprob,betashape1,betashape2,
     file=paste0(myoutdir,"/",mysim,"_params.RData"))



## --------------------------------------------------------------------
## SIM 23
## --------------------------------------------------------------------
## (medium number large-sized intra-chr and large-sized inter-chr)
## (plus 1 fission, 1 fusion)

## simulation name
mysim<-"sim23"

## number of chromosomes and genes per chromosome
ngenes<-c(1000,1000,1000,1000,1000)

## number of events per type
nrearr<-list(inv=20,synmov=20,nonsynmov=20,
             retra=0,fus=1,fis=1)

## set rearrangement size by size and prob of negative binomial distribution
##   +1 gene to not have zero sizes
## not applicable for retra, fus, fis
rearrsize<-list(inv=lg[1],synmov=lg[1],nonsynmov=lg[1],
                retra=NA,fus=NA,fis=NA)
rearrprob<-list(inv=lg[2],synmov=lg[2],nonsynmov=lg[2],
                retra=NA,fus=NA,fis=NA)

## probability of breakpoints along chromosome
betashape1<-1
betashape2<-1

save(mysim,ngenes,nrearr,rearrsize,rearrprob,betashape1,betashape2,
     file=paste0(myoutdir,"/",mysim,"_params.RData"))



## --------------------------------------------------------------------
## SIM 24
## --------------------------------------------------------------------
## (medium number xlarge-sized intra-chr and xlarge-sized inter-chr)
## (plus 1 fission, 1 fusion)

## simulation name
mysim<-"sim24"

## number of chromosomes and genes per chromosome
ngenes<-c(1000,1000,1000,1000,1000)

## number of events per type
nrearr<-list(inv=20,synmov=20,nonsynmov=20,
             retra=0,fus=1,fis=1)

## set rearrangement size by size and prob of negative binomial distribution
##   +1 gene to not have zero sizes
## not applicable for retra, fus, fis
rearrsize<-list(inv=xl[1],synmov=xl[1],nonsynmov=xl[1],
                retra=NA,fus=NA,fis=NA)
rearrprob<-list(inv=xl[2],synmov=xl[2],nonsynmov=xl[2],
                retra=NA,fus=NA,fis=NA)

## probability of breakpoints along chromosome
betashape1<-1
betashape2<-1

save(mysim,ngenes,nrearr,rearrsize,rearrprob,betashape1,betashape2,
     file=paste0(myoutdir,"/",mysim,"_params.RData"))

