#!/usr/bin/env Rscript

## call with
## Rscript --vanilla runSims.R simparams/sim1_params.RData simout/sim1/sim1.RData 2
## [where sim1.RData should be the basename of the data set w/o replicate ID]

args<-commandArgs(trailingOnly=TRUE)

if(length(args)!=3) {
    stop("Require three arguments: params.RData out.RData nrep")
}

## read simulation settings
## ------------------------

load(args[1])
myout<-args[2]
nrep<-as.integer(args[3])


if(file.exists(sub(".RData",paste0("_",nrep,".RData"),myout))){
    print(paste("Simulations",myout,"already exist"))
    quit(save="no")
}

## --------------------------------------------------------------------
## SPECIFY SIMULATION SETTINGS: makeSimParams.R
## --------------------------------------------------------------------

myscriptsdir<-"~/Rearrangements/Simulations/R"

## load functions
source(paste(myscriptsdir,"rearrsimsUtility.R",sep="/"))
source(paste(myscriptsdir,"rearrsimsOperations.R",sep="/"))


## --------------------------------------------------------------------
## PREPARE SIMULATIONS
## --------------------------------------------------------------------

## genome
nchr<-length(ngenes)
genome0<-vector("list",nchr)
cntr<-0
for(i in 1:nchr){
    genome0[[i]]<-(cntr+1):(cntr+ngenes[i])
    cntr<-cntr+ngenes[i]
}
nallgenes<-cntr


## set up output for rearranged genes
tmp<-matrix(ncol=0,nrow=nallgenes)
rownames(tmp)<-as.character(1:nallgenes)
rearrs0<-list(inv=tmp,synmov=tmp,nonsynmov=tmp,
              retra=tmp,fus=tmp,fis=tmp)

## set up output for breakpoints
tmp<-matrix(0,ncol=2,nrow=nallgenes)
tmp<-cbind(tmp,tmp,tmp)
rownames(tmp)<-as.character(1:nallgenes)
colnames(tmp)<-paste0(rep(c("bl.","ex.","in."),each=2),
                      rep(c("head","tail"),times=3))
brpts0<-list(inv=tmp,synmov=tmp,nonsynmov=tmp,
             retra=tmp,fus=tmp,fis=tmp)
rm(tmp)




## --------------------------------------------------------------------
## SIMULATE
## --------------------------------------------------------------------

for(x in 1:nrep){

    if(file.exists(sub(".RData",paste0("_",x,".RData"),myout))){
        print(paste("Simulation",sub(".RData",paste0("_",x,".RData"),myout),"already exists"))
        next
    }

    genome<-genome0
    rearrs<-rearrs0
    brpts<-brpts0

    ## draw order of events
    events<-integer()
    for(i in 1:length(names(nrearr))){
        events<-c(events,rep(i,nrearr[[i]]))
    }
    events<-sample(events)

    if(length(events)==0){
        stop("No rearrangement events were specified")
    }else{
        for(i in 1:length(events)){
            myevent<-events[i]

            ## re-calculate number of chromosomes and genes
            mynchr<-length(genome)
            myngenes<-integer(mynchr)
            for(n in 1:mynchr){
                myngenes[n]<-length(genome[[n]])
            }


            ## the below requires:
            ## genome,myevent,mynchr,myngenes,betashape1,betashape2,
            ##   rearrsize,rearrprob

            ## need to return:
            ## genome,myblock1,myblock2

            if(myevent==1){ ## === inversion <inv> ===
                G<-do.inv(genome,myevent,mynchr,myngenes,
                          betashape1,betashape2,rearrsize,rearrprob)

            }else if(myevent==2){ ## === syntenic move <synmov> ===
                G<-do.synmov(genome,myevent,mynchr,myngenes,
                             betashape1,betashape2,rearrsize,rearrprob)

            }else if(myevent==3){ ## === nonsyntenic move <nonsynmov> ===
                G<-do.nonsynmov(genome,myevent,mynchr,myngenes,
                                betashape1,betashape2,rearrsize,rearrprob)

            }else if(myevent==4){ ## === reciprocal translocation <retra> ===
                G<-do.retra(genome,mynchr,myngenes,betashape1,betashape2)

            }else if(myevent==5){ ## === fusion <fus> ===
                G<-do.fus(genome,mynchr,myngenes)

            }else if(myevent==6){ ## === fission <fis> ===
                G<-do.fis(genome,mynchr,myngenes,betashape1,betashape2)

            }else{
                stop("Unknown rearrangement type")
            }

            ## double-check that no chromosome has length 0
            my0chr<-numeric()
            for(n in 1:length(G$genome)){
                if(length(G$genome[[n]])==0){
                    my0chr<-c(my0chr,n)
                }
            }
            while(length(my0chr)>0){
                G$genome[[my0chr[1]]]<-NULL
                ## remove index and adjust remaining indices
                my0chr<-(my0chr[-1])-1
            }

            ## store new genome
            genome<-G$genome

            ## store tags and breakpoints
            if(is.element(myevent,c(1))){ ## inv
                ## add tags
                rearrs<-addtags(rearrs,myevent,G$myblock)
                ## add breakpoints of moved block
                brpts<-addbrpts(brpts,myevent,G$myblock,"block")
                ## add breakpoints of excision (==insertion)
                brpts<-addbrpts(brpts,myevent,G$excip,"excision")
            }else if(is.element(myevent,c(2,3))){ ## synmov, nonsynmov
                ## add tags
                ## (>>>note: no tags for leaped block for '2')
                rearrs<-addtags(rearrs,myevent,G$myblock1)
                ## add breakpoints of moved block
                brpts<-addbrpts(brpts,myevent,G$myblock1,"block")
                ## add breakpoints of excision
                brpts<-addbrpts(brpts,myevent,G$excip,"excision")
                ## add breakpoints of insertion
                brpts<-addbrpts(brpts,myevent,G$insep,"insertion")
            }else if(is.element(myevent,c(4,5,6))){
                ## add tags (two columns)
                ## (>>>note: no tags for unmoved blocks for '4')
                rearrs<-addtags(rearrs,myevent,G$myblock1)
                rearrs<-addtags(rearrs,myevent,G$myblock2)
                ## add breakpoints of moved blocks
                ## (>>>note: including telomeres, which are
                ##           actually not true breakpoints)
                brpts<-addbrpts(brpts,myevent,G$myblock1,"block")
                brpts<-addbrpts(brpts,myevent,G$myblock2,"block")
                if(is.element(myevent,c(4))){
                    ## add breakpoints of excisions (==insertions)
                    brpts<-addbrpts(brpts,myevent,G$excip1,"excision")
                    brpts<-addbrpts(brpts,myevent,G$excip2,"excision")
                }
            }
        }
    }

    ## TODO: work with telomeres by adding 'NA' at respective
    ##       end of myblock, but remove it in 'addtags()'
    ##       (this is because for retra(4), fus(5), and fis(6),
    ##        the block side that corresponds to a former chromosome
    ##        end is actually not a breakpoint in 'addbrpts()')



    ## prepare original and evolved genome for further usage
    ## -----------------------------------------------------

    ## prepare ancestral genome
    ancgenome<-1:nallgenes
    ancchrstart<-numeric()
    cntr<-0
    for(i in 1:nchr){
        ancchrstart<-c(ancchrstart,cntr+1)
        cntr<-cntr+ngenes[i]
    }

    ## prepare evolved genome
    evogenome<-numeric()
    evochrstart<-numeric()
    cntr<-0
    for(i in 1:length(genome)){
        evogenome<-c(evogenome,genome[[i]])
        evochrstart<-c(evochrstart,cntr+1)
        cntr<-cntr+length(genome[[i]])
    }
    strand<-rep(1,length(evogenome))
    strand[evogenome<0]<- -1
    evogenome<-abs(evogenome)


    ## --------------------------------------------------------------------
    ## MAKE REARRVISR GENOME FILES
    ## --------------------------------------------------------------------

    ## focalgenome
    evochrsize<-diff(c(evochrstart,length(evogenome)+1))
    markerstart<-numeric()
    for(i in 1:length(evochrsize)){
        ##markerstart<-c(markerstart,seq(10^6,by=10^6,length.out=evochrsize[i]))
        ## ##(10^6 is out of integer range for very large fused chromosomes)
        markerstart<-c(markerstart,seq(10^4,by=10^4,length.out=evochrsize[i]))
    }
    markerstrand<-rep("+",length(evogenome))
    markerstrand[strand<0]<-"-"

    focalgenome<-data.frame(
        marker=as.integer(evogenome),
        scaff=as.character(rep(1:length(evochrsize),evochrsize)),
        start=as.integer(markerstart),
        end=as.integer(markerstart+2),
        strand=markerstrand,
        stringsAsFactors = FALSE)


    ## compgenome
    ancchrsize<-diff(c(ancchrstart,length(ancgenome)+1))
    markerstart<-numeric()
    for(i in 1:length(ancchrsize)){
        markerstart<-c(markerstart,1:(ancchrsize[i]))
    }
    markerstrand<-rep("+",length(ancgenome))

    compgenome<-data.frame(
        marker=as.integer(ancgenome),
        orientation=markerstrand,
        car=as.integer(rep(1:length(ancchrsize),ancchrsize)),
        type1=rep("Q",length(ancgenome)),
        elem1=as.integer(markerstart),
        stringsAsFactors = FALSE)


    ## --------------------------------------------------------------------
    ## COMPUTE BASIC STATS
    ## --------------------------------------------------------------------

    ## count for each marker the number of rearrangement events
    allrearrs<-rowSums(cbind(rearrs[[1]],rearrs[[2]],rearrs[[3]],
                             rearrs[[4]],rearrs[[5]],rearrs[[6]]))
    rearrstats<-character(11)
    for(k in 0:9){
        rearrstats[k+1]<-paste(k,sum(allrearrs==k),sep=":")
    }
    rearrstats[11]<-paste(">=10",sum(allrearrs>=10),sep=":")
    ## tmp<-matrix(unlist(strsplit(rearrstats,":")),ncol=2,byrow=TRUE)
    ## stats<-data.frame(nrearr=tmp[,1],count=as.numeric(tmp[,2]),
    ##                   stringsAsFactors=FALSE)

    ## --------------------------------------------------------------------
    ## SAVE SIMULATION OUTPUT
    ## --------------------------------------------------------------------
    save(compgenome,focalgenome,events,rearrs,brpts,
         file=sub(".RData",paste0("_",x,".RData"),myout))

    ## print counts
    print(paste("Counts",sub(".RData",paste0("_",x),basename(myout)),
                paste(rearrstats,collapse=" ")),quote=FALSE)
    print(paste("Finished",sub(".RData",paste0("_",x,".RData"),myout)),
          quote=FALSE)
}
