#!/usr/bin/env Rscript

## call with (when adding genome fragmentation)
## Rscript --vanilla runRearrvisr.R simout/sim1/sim1.RData 2 simparams/frag1_params.RData
## call with (without adding genome fragmentation)
## Rscript --vanilla runRearrvisr.R simout/sim1/sim1.RData 2
## [where sim1.RData should be the basename of the data set w/o replicate ID]

args<-commandArgs(trailingOnly=TRUE)

if(length(args)<2) {
    stop("Require at least two arguments: out.RData nrep frag_params")
}else if(length(args)==2){
    myfrag<-NULL
}else{
    ##load("simparams/frag5_params.RData")
    load(args[3])
}

filebase<-args[1]
nrep<-as.integer(args[2])


## --------------------------------------------------------------------
## SPECIFY GENOME FRAGMENTATION SETTINGS: makeFragParams.R
## --------------------------------------------------------------------

myscriptsdir<-"~/Rearrangements/Simulations/R"

## load functions
source(paste(myscriptsdir,"rearrsimsOutputUtility.R",sep="/"))

##library(rearrvisr)
library("rearrvisr", lib.loc="~/R/site-library/")


for(x in 1:nrep){

    ## load correct simulation output
    ##load("simout/sim3/sim3_1.RData")
    myin<-sub(".RData",paste0("_",x,".RData"),filebase)
    load(myin)

    ## make output file names
    if(!is.null(myfrag)){
        myout1<-sub(".RData",paste0("_",x,"_",myfrag,"_rearrs.txt"),filebase)
        myout2<-sub(".RData",paste0("_",x,"_",myfrag,"_PrecRec.txt"),filebase)
        myout3<-sub(".RData",paste0("_",x,"_",myfrag,"_FragMeas.txt"),filebase)
        myout4<-sub(".RData",paste0("_",x,"_",myfrag,".RData"),filebase)
    }else{
        myout1<-sub(".RData",paste0("_",x,"_rearrs.txt"),filebase)
        myout2<-sub(".RData",paste0("_",x,"_PrecRec.txt"),filebase)
        myout3<-sub(".RData",paste0("_",x,"_FragMeas.txt"),filebase)
        myout4<-sub(".RData",paste0("_",x,".RData"),filebase)
    }


    ## --------------------------------------------------------------------
    ## OPTIONALLY FRAGMENT GENOME
    ## --------------------------------------------------------------------
    if(!is.null(myfrag)){

        if(file.exists(myout4)){
            ## fragmentation exists already
            ## ----------------------------
            print(paste("using existing genome fragmentation in",myout4))
            load(myout4)

        }else{

            ## compute original chromosome start positions for plots
            ## -----------------------------------------------------
            ## ancestral genome
            mychr<-unique(compgenome$car)
            org.ancchrstart<-numeric()
            cntr<-0
            for(i in 1:length(mychr)){
                org.ancchrstart<-c(org.ancchrstart,cntr+1)
                cntr<-cntr+sum(compgenome$car==mychr[i])
            }
            ## evolved genome
            mychr<-unique(focalgenome$scaff)
            org.evochrstart<-numeric()
            cntr<-0
            for(i in 1:length(mychr)){
                org.evochrstart<-c(org.evochrstart,cntr+1)
                cntr<-cntr+sum(focalgenome$scaff==mychr[i])
            }


            ## fragment genomes
            ## ----------------
            if(is.element("focal",tofrag)){
                focalgenome<-fragmentGenome(focalgenome,"focal",
                                            betashapes,nbreaks,size1prop)
            }
            if(is.element("comp",tofrag)){
                compgenome<-fragmentGenome(compgenome,"comp",
                                           betashapes,nbreaks,size1prop)
            }

        }
    }

    ## --------------------------------------------------------------------
    ## RUN REARRVISR
    ## --------------------------------------------------------------------

    SYNT<-computeRearrs(focalgenome,compgenome,doubled=TRUE)

    ## --------------------------------------------------------------------
    ## STATISTICS
    ## --------------------------------------------------------------------


    ## only test whether marker has been detected, not how many times
    ##  (as some tags can occur multiple times for the same event)

    ## ancestral position translated to evolved position
    newpos<-match(focalgenome$marker,compgenome$marker)
    ## .low and .high estimates including or excluding tags==predThld
    predThld<-0.5



    ## rearrangements between genome segments subset 1
    ## -----------------------------------------------
    ## mat1: nonsynmov, retra, fus, fis
    ## mat2: NM1
    BS1<-precisionRecall(mat1=cbind(rearrs[[3]],rearrs[[4]],
                             rearrs[[5]],rearrs[[6]])[newpos,],
                         mat2=SYNT[[1]],
                         predThld=predThld)

    ## rearrangements between genome segments subset 2
    ## -----------------------------------------------
    ## mat1: nonsynmov, retra, fus, fis
    ## mat2: NM2
    BS2<-precisionRecall(mat1=cbind(rearrs[[3]],rearrs[[4]],
                             rearrs[[5]],rearrs[[6]])[newpos,],
                         mat2=SYNT[[2]],
                         predThld=predThld)


    ## nonsynonymous moves subset 1
    ## ----------------------------
    ## mat1: nonsynmov
    ## mat2: NM1
    NM1<-precisionRecall(mat1=rearrs[[3]][newpos,],
                         mat2=SYNT[[1]],
                         predThld=predThld)

    ## nonsynonymous moves subset 2
    ## ----------------------------
    ## mat1: nonsynmov
    ## mat2: NM2
    NM2<-precisionRecall(mat1=rearrs[[3]][newpos,],
                         mat2=SYNT[[2]],
                         predThld=predThld)

    ## synonymous moves
    ## ----------------
    ## mat1: synmov
    ## mat2: SM
    SM<-precisionRecall(mat1=rearrs[[2]][newpos,],
                        mat2=SYNT[[3]],
                        predThld=predThld)

    ## inversions
    ## ----------
    ## mat1: inv
    ## mat2: IV
    IV<-precisionRecall(mat1=rearrs[[1]][newpos,],
                        mat2=SYNT[[4]],
                        predThld=predThld)



    ## summarize number and type of rearrangement events and breakpoints
    ## -----------------------------------------------------------------
    ordfocal<-unique(focalgenome$scaff)
    rearrsum<-summarizeRearrs(SYNT,focalgenome,compgenome,ordfocal)
    rearrsum<-cbind(segment=ordfocal,rearrsum)


    ## compute measurements of genome fragmentation
    ## --------------------------------------------
    focalFM<-fragMeasure(focalgenome$scaff)
    compFM<-fragMeasure(compgenome$car)


    ## --------------------------------------------------------------------
    ## WRITE OUTPUT
    ## --------------------------------------------------------------------

    precrec<-rbind(unlist(BS1),unlist(BS2),unlist(NM1),unlist(NM2),
                   unlist(SM),unlist(IV))
    precrec<-cbind(rearr.class=c("between.sub1","between.sub2",
                       "nonsyn.sub1","nonsyn.sub2",
                       "syn.moves","inversions"),precrec)

    fragmeas<-rbind(unlist(focalFM),unlist(compFM))
    fragmeas<-cbind(genome=c("focal","comp"),fragmeas)



    write.table(rearrsum,file=myout1,
                row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
    write.table(precrec,file=myout2,
                row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
    write.table(fragmeas,file=myout3,
                row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')

    ## --------------------------------------------------------------------
    ## (RE-)SAVE OUTPUT
    ## --------------------------------------------------------------------

    if(!is.null(myfrag)){
        save(compgenome,focalgenome,SYNT,org.ancchrstart,org.evochrstart,
             file=myout4)
    }else{
        save(compgenome,focalgenome,events,rearrs,brpts,SYNT,file=myout4)
    }

    ## --------------------------------------------------------------------
    ## CLEAN UP
    ## --------------------------------------------------------------------

    if(!is.null(myfrag)){
        rm(org.ancchrstart,org.evochrstart,mychr,cntr,i)
    }
    rm(focalgenome,compgenome,events,rearrs,brpts,SYNT,newpos,predThld,
       BS1,BS2,NM1,NM2,SM,IV,rearrsum,focalFM,compFM,
       precrec,fragmeas,myout1,myout2,myout3,myout4)

    print(paste("Finished",myin),quote=FALSE)
}
