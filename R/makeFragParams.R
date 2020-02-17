#!/usr/bin/env Rscript

## --------------------------------------------------------------------
## generate RData files with fragmentation settings (exploratory)
## --------------------------------------------------------------------

myoutdir<-"~/Rearrangements/Simulations/simparams"


## betashapes=c(1.5,8) give a bit more left-shifted distribution,
## i.e., there are still 5 larger chromosomes, while
## betashapes=c(1,1) result in a rather gamma-shaped distribution.
## both are realistic, the first for better quality assemblies
## and the second for lower quality assemblies
## (formerly tested betashapes=c(1.5,4), which is intermediate)

## plot(seq(0,1,0.01),dbeta(seq(0,1,0.01),betashapes[1],betashapes[2]),
##      type="l",xlab="position",ylab="density")


## --------------------------------------------------------------------
## FRAG 1
## --------------------------------------------------------------------
## (10 breaks with betashapes=c(1.5,8) - focal only)

## fragmentation name (always "frag" and some digits)
myfrag<-"frag1"

## which genome(s) to fragment (c("focal","comp"))
tofrag<-"focal"

## artificially fragment a genome into smaller pieces
##  - betashapes control the skew of breakage of a scaffold
##  - nbreaks defines maximum total number of breakages
##  - size1prop controls when simulations are stopped because
##    a proportion x of scaffolds contains only one marker

## probability of breakpoints along chromosome
betashapes<-c(1.5,8)
## number of breaks for genome
nbreaks<-10
## max proportion of size 1 contigs
size1prop<-0.5


save(myfrag,tofrag,betashapes,nbreaks,size1prop,
     file=paste0(myoutdir,"/",myfrag,"_params.RData"))


## --------------------------------------------------------------------
## FRAG 2
## --------------------------------------------------------------------
## (10 breaks with betashapes=c(1.5,8) - comp only)

## fragmentation name (always "frag" and some digits)
myfrag<-"frag2"

## which genome(s) to fragment (c("focal","comp"))
tofrag<-"comp"

## artificially fragment a genome into smaller pieces
##  - betashapes control the skew of breakage of a scaffold
##  - nbreaks defines maximum total number of breakages
##  - size1prop controls when simulations are stopped because
##    a proportion x of scaffolds contains only one marker

## probability of breakpoints along chromosome
betashapes<-c(1.5,8)
## number of breaks for genome
nbreaks<-10
## max proportion of size 1 contigs
size1prop<-0.5


save(myfrag,tofrag,betashapes,nbreaks,size1prop,
     file=paste0(myoutdir,"/",myfrag,"_params.RData"))


## --------------------------------------------------------------------
## FRAG 3
## --------------------------------------------------------------------
## (10 breaks with betashapes=c(1.5,8) - both genomes)

## fragmentation name (always "frag" and some digits)
myfrag<-"frag3"

## which genome(s) to fragment (c("focal","comp"))
tofrag<-c("focal","comp")

## artificially fragment a genome into smaller pieces
##  - betashapes control the skew of breakage of a scaffold
##  - nbreaks defines maximum total number of breakages
##  - size1prop controls when simulations are stopped because
##    a proportion x of scaffolds contains only one marker

## probability of breakpoints along chromosome
betashapes<-c(1.5,8)
## number of breaks for genome
nbreaks<-10
## max proportion of size 1 contigs
size1prop<-0.5


save(myfrag,tofrag,betashapes,nbreaks,size1prop,
     file=paste0(myoutdir,"/",myfrag,"_params.RData"))


## --------------------------------------------------------------------
## FRAG 4
## --------------------------------------------------------------------
## (20 breaks with betashapes=c(1.5,8) - focal only)

## fragmentation name (always "frag" and some digits)
myfrag<-"frag4"

## which genome(s) to fragment (c("focal","comp"))
tofrag<-"focal"

## artificially fragment a genome into smaller pieces
##  - betashapes control the skew of breakage of a scaffold
##  - nbreaks defines maximum total number of breakages
##  - size1prop controls when simulations are stopped because
##    a proportion x of scaffolds contains only one marker

## probability of breakpoints along chromosome
betashapes<-c(1.5,8)
## number of breaks for genome
nbreaks<-20
## max proportion of size 1 contigs
size1prop<-0.5


save(myfrag,tofrag,betashapes,nbreaks,size1prop,
     file=paste0(myoutdir,"/",myfrag,"_params.RData"))


## --------------------------------------------------------------------
## FRAG 5
## --------------------------------------------------------------------
## (20 breaks with betashapes=c(1.5,8) - comp only)

## fragmentation name (always "frag" and some digits)
myfrag<-"frag5"

## which genome(s) to fragment (c("focal","comp"))
tofrag<-"comp"

## artificially fragment a genome into smaller pieces
##  - betashapes control the skew of breakage of a scaffold
##  - nbreaks defines maximum total number of breakages
##  - size1prop controls when simulations are stopped because
##    a proportion x of scaffolds contains only one marker

## probability of breakpoints along chromosome
betashapes<-c(1.5,8)
## number of breaks for genome
nbreaks<-20
## max proportion of size 1 contigs
size1prop<-0.5


save(myfrag,tofrag,betashapes,nbreaks,size1prop,
     file=paste0(myoutdir,"/",myfrag,"_params.RData"))


## --------------------------------------------------------------------
## FRAG 6
## --------------------------------------------------------------------
## (20 breaks with betashapes=c(1.5,8) - both genomes)

## fragmentation name (always "frag" and some digits)
myfrag<-"frag6"

## which genome(s) to fragment (c("focal","comp"))
tofrag<-c("focal","comp")

## artificially fragment a genome into smaller pieces
##  - betashapes control the skew of breakage of a scaffold
##  - nbreaks defines maximum total number of breakages
##  - size1prop controls when simulations are stopped because
##    a proportion x of scaffolds contains only one marker

## probability of breakpoints along chromosome
betashapes<-c(1.5,8)
## number of breaks for genome
nbreaks<-20
## max proportion of size 1 contigs
size1prop<-0.5


save(myfrag,tofrag,betashapes,nbreaks,size1prop,
     file=paste0(myoutdir,"/",myfrag,"_params.RData"))



