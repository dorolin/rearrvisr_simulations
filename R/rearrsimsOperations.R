## --------------------------------------------------------------------
## REARRANGEMENT OPERATIONS for 'rearrsims.R'
## --------------------------------------------------------------------

do.inv<-function(genome,myevent,mynchr,myngenes,
                 betashape1,betashape2,rearrsize,rearrprob){

    ## draw block
    redo<-1
    while(redo==1){
        B<-drawblock(myevent,mynchr,myngenes,
                     betashape1,betashape2,rearrsize,rearrprob)
        ## require that block does not cover full chromosome
        if(length(B$blockend-B$blockstart+1)>=myngenes[B$chr]){
            next
        }
        redo<-0
    }
    ## store excision boundaries
    excip<-c(NA,NA)
    if(B$blockstart>1){
        excip[1]<-genome[[B$chr]][B$blockstart-1]
    }
    if(B$blockend<myngenes[B$chr]){
        excip[2]<-genome[[B$chr]][B$blockend+1]
    }
    ## invert block
    myblock<-genome[[B$chr]][B$blockstart:B$blockend]
    genome[[B$chr]][B$blockstart:B$blockend]<-rev(myblock)*(-1)

    return(list(genome=genome,myblock=myblock,excip=excip))
}

do.synmov<-function(genome,myevent,mynchr,myngenes,
                    betashape1,betashape2,rearrsize,rearrprob){

    ## draw block and additional breakpoint
    redo<-1
    while(redo==1){
        B<-drawblock(myevent,mynchr,myngenes,
                     betashape1,betashape2,rearrsize,rearrprob)
        ## draw additional breakpoint position
        ##   requires that block does not cover full chromosome
        if(length(B$blockend-B$blockstart+1)>=myngenes[B$chr]){
            next
        }
        ## draw breakpoint position *(ngenes+1) to have positions
        ##  bound between 0 and ngenes+1
        mybrpt<-rbeta(1,betashape1,betashape2)*(myngenes[B$chr]+1)
        ##   require at least one gene as gap to block
        if(mybrpt<=(B$blockend+1) & mybrpt>=(B$blockstart-1)){
            next
        }
        redo<-0
    }
    ## move block
    myblock<-genome[[B$chr]][B$blockstart:B$blockend]
    if(mybrpt>B$blockend+1){
        leftp<-floor(mybrpt)
        rightp<-leftp+1
        mynewchr<-numeric()
        if(B$blockstart>1){
            mynewchr<-genome[[B$chr]][1:(B$blockstart-1)]
        }
        mynewchr<-c(mynewchr,
                    genome[[B$chr]][(B$blockend+1):leftp],
                    myblock)
        if(rightp<=myngenes[B$chr]){
            mynewchr<-c(mynewchr,
                        genome[[B$chr]][rightp:myngenes[B$chr]])
        }
    }else if(mybrpt<B$blockstart-1){
        rightp<-ceiling(mybrpt)
        leftp<-rightp-1
        mynewchr<-numeric()
        if(leftp>=1){
            mynewchr<-genome[[B$chr]][1:leftp]
        }
        mynewchr<-c(mynewchr,
                    myblock,
                    genome[[B$chr]][rightp:(B$blockstart-1)])
        if(B$blockend<myngenes[B$chr]){
            mynewchr<-c(mynewchr,
                        genome[[B$chr]][(B$blockend+1):myngenes[B$chr]])
        }
    }else{
        stop("Something went wrong when defining breakpoint position")
    }
    if(length(unique(mynewchr))!=myngenes[B$chr]){
        stop("Something went wrong when moving block")
    }

    ## store excision boundaries
    excip<-c(NA,NA)
    if(B$blockstart>1){
        excip[1]<-genome[[B$chr]][B$blockstart-1]
    }
    if(B$blockend<myngenes[B$chr]){
        excip[2]<-genome[[B$chr]][B$blockend+1]
    }
    ## store insertion boundaries
    insep<-c(NA,NA)
    if(leftp>=1){
        insep[1]<-genome[[B$chr]][leftp]
    }
    if(rightp<=myngenes[B$chr]){
        insep[2]<-genome[[B$chr]][rightp]
    }

    genome[[B$chr]]<-mynewchr

    return(list(genome=genome,myblock1=myblock,myblock2=NULL,
                excip=excip,insep=insep))
}


do.nonsynmov<-function(genome,myevent,mynchr,myngenes,
                       betashape1,betashape2,rearrsize,rearrprob){

    if(mynchr<2){
        warning("nonsyntenic move not possible, skipping event",
                immediate.=TRUE)
        next
    }
    ## draw block
    redo<-1
    while(redo==1){
        B<-drawblock(myevent,mynchr,myngenes,
                     betashape1,betashape2,rearrsize,rearrprob)
        ## require that block does not cover full chromosome
        if(length(B$blockend-B$blockstart+1)>=myngenes[B$chr]){
            next
        }
        redo<-0
    }
    ## draw another chromosome and additional breakpoint
    sinkchr<-sample((1:mynchr)[-B$chr],1,prob=myngenes[-B$chr])
    ## draw breakpoint position *(ngenes+1) to have positions
    ##  bound between 0 and ngenes+1
    mybrpt<-rbeta(1,betashape1,betashape2)*(myngenes[sinkchr]+1)
    ## excise block
    myblock<-genome[[B$chr]][B$blockstart:B$blockend]
    ## insert block
    if(floor(mybrpt)==ceiling(mybrpt)){
        if(floor(mybrpt)==0){
            leftp<-0
        }else if(ceiling(mybrpt)==myngenes[sinkchr]+1){
            leftp<-myngenes[sinkchr]
        }else{
            leftp<-floor(mybrpt)+rbinom(1,1,0.5)
        }
    }else{
        leftp<-floor(mybrpt)
    }
    rightp<-leftp+1
    mynewchr<-numeric()
    if(leftp>=1){
        mynewchr<-genome[[sinkchr]][1:leftp]
    }
    mynewchr<-c(mynewchr,myblock)
    if(rightp<=myngenes[sinkchr]){
        mynewchr<-c(mynewchr,
                    genome[[sinkchr]][rightp:myngenes[sinkchr]])
    }

    ## store excision boundaries
    excip<-c(NA,NA)
    if(B$blockstart>1){
        excip[1]<-genome[[B$chr]][B$blockstart-1]
    }
    if(B$blockend<myngenes[B$chr]){
        excip[2]<-genome[[B$chr]][B$blockend+1]
    }
    ## store insertion boundaries
    insep<-c(NA,NA)
    if(leftp>=1){
        insep[1]<-genome[[sinkchr]][leftp]
    }
    if(rightp<=myngenes[sinkchr]){
        insep[2]<-genome[[sinkchr]][rightp]
    }

    genome[[B$chr]]<-genome[[B$chr]][-c(B$blockstart:B$blockend)]
    genome[[sinkchr]]<-mynewchr

    if(length(unique(mynewchr))!=myngenes[sinkchr]+length(myblock) |
       length(unique(genome[[B$chr]]))!=myngenes[B$chr]-length(myblock)){
        stop("Something went wrong when moving block")
    }

    return(list(genome=genome,myblock1=myblock,myblock2=NULL,
                excip=excip,insep=insep))
}


do.retra<-function(genome,mynchr,myngenes,betashape1,betashape2){

    mycandchr<-(1:mynchr)[myngenes>=2]
    if(length(mycandchr)<2){
        warning("reciprocal translocation not possible, skipping event",
                immediate.=TRUE)
        next
    }
    ## draw chromosomes
    chrs<-sample(mycandchr,2,prob=myngenes[mycandchr])
    ## draw first breakpoint and get block
    ## draw breakpoint position *(ngenes-1)+1 to have positions
    ##  bound between 1 and ngenes
    mybrpt1<-rbeta(1,betashape1,betashape2)*(myngenes[chrs[1]]-1) + 1
    if(mybrpt1<(myngenes[chrs[1]]+1)/2){ ## move left part
        tomove1<-0
    }else if(mybrpt1>(myngenes[chrs[1]]+1)/2){ ## move right part
        tomove1<-1
    }else{
        tomove1<-rbinom(1,1,0.5)
    }
    if(tomove1==0){
        leftp1<-floor(mybrpt1)
        rightp1<-leftp1+1
        myblock1<-genome[[chrs[1]]][1:leftp1]
    }else{
        rightp1<-ceiling(mybrpt1)
        leftp1<-rightp1-1
        myblock1<-genome[[chrs[1]]][rightp1:myngenes[chrs[1]]]
    }
    ## draw second breakpoint and get block
    ## draw breakpoint position *(ngenes-1)+1 to have positions
    ##  bound between 1 and ngenes
    mybrpt2<-rbeta(1,betashape1,betashape2)*(myngenes[chrs[2]]-1) + 1
    if(mybrpt2<(myngenes[chrs[2]]+1)/2){ ## move left part
        tomove2<-0
    }else if(mybrpt2>(myngenes[chrs[2]]+1)/2){ ## move right part
        tomove2<-1
    }else{
        tomove2<-rbinom(1,1,0.5)
    }
    if(tomove2==0){
        leftp2<-floor(mybrpt2)
        rightp2<-leftp2+1
        myblock2<-genome[[chrs[2]]][1:leftp2]
    }else{
        rightp2<-ceiling(mybrpt2)
        leftp2<-rightp2-1
        myblock2<-genome[[chrs[2]]][rightp2:myngenes[chrs[2]]]
    }

    ## store excision boundaries
    excip1<-c(NA,NA)
    if(tomove1==0){
        if(rightp1<=myngenes[chrs[1]]){
            excip1[2]<-genome[[chrs[1]]][rightp1]
        }
    }else{
        if(leftp1>=1){
            excip1[1]<-genome[[chrs[1]]][leftp1]
        }
    }
    excip2<-c(NA,NA)
    if(tomove2==0){
        if(rightp2<=myngenes[chrs[2]]){
            excip2[2]<-genome[[chrs[2]]][rightp2]
        }
    }else{
        if(leftp2>=2){
            excip2[1]<-genome[[chrs[2]]][leftp2]
        }
    }

    ## combine blocks, invert order if tomove1 != tomove2
    if(tomove1==tomove2){
        if(tomove1==0){
            genome[[chrs[1]]]<-c(myblock2,
                                 genome[[chrs[1]]][rightp1:myngenes[chrs[1]]])
            genome[[chrs[2]]]<-c(myblock1,
                                 genome[[chrs[2]]][rightp2:myngenes[chrs[2]]])
        }else{
            genome[[chrs[1]]]<-c(genome[[chrs[1]]][1:leftp1],myblock2)
            genome[[chrs[2]]]<-c(genome[[chrs[2]]][1:leftp2],myblock1)
        }
    }else{
        if(tomove1==0){
            genome[[chrs[1]]]<-c(rev(myblock2)*(-1),
                                 genome[[chrs[1]]][rightp1:myngenes[chrs[1]]])
            genome[[chrs[2]]]<-c(genome[[chrs[2]]][1:leftp2],
                                 rev(myblock1)*(-1))
        }else{
            genome[[chrs[1]]]<-c(genome[[chrs[1]]][1:leftp1],
                                 rev(myblock2)*(-1))
            genome[[chrs[2]]]<-c(rev(myblock1)*(-1),
                                 genome[[chrs[2]]][rightp2:myngenes[chrs[2]]])
        }
    }
    if(length(unique(abs(c(genome[[chrs[1]]],genome[[chrs[2]]]))))!=sum(myngenes[chrs])){
        stop("Something went wrong when moving blocks")
    }

    return(list(genome=genome,myblock1=myblock1,myblock2=myblock2,
                excip1=excip1,excip2=excip2))
}

do.fus<-function(genome,mynchr,myngenes){

    if(mynchr<2){
        warning("fusion not possible, skipping event",
                immediate.=TRUE)
        next
    }
    ## draw chromosomes (prefer short ones)
    chrs<-sample(1:mynchr,2,prob=1/myngenes)
    ## define fusion points (0 head; 1 tail)
    fusp1<-rbinom(1,1,0.5)
    fusp2<-rbinom(1,1,0.5)
    ## save unfused chromosomes for making tags below
    myblock1<-genome[[chrs[1]]]
    myblock2<-genome[[chrs[2]]]
    ## fuse
    if(fusp1==1){
        if(fusp2==0){
            mynewchr<-c(myblock1,myblock2)
        }else{
            mynewchr<-c(myblock1,rev(myblock2)*(-1))
        }
    }else{
        if(fusp2==0){
            mynewchr<-c(rev(myblock1)*(-1),myblock2)
        }else{
            mynewchr<-c(rev(myblock1)*(-1),rev(myblock2)*(-1))
        }
    }
    ## replace first unfused chromosome by fused chromosome
    genome[[chrs[1]]]<-mynewchr
    ## remove second unfused chromosome
    genome[[chrs[2]]]<-NULL

    return(list(genome=genome,myblock1=myblock1,myblock2=myblock2))
}

do.fis<-function(genome,mynchr,myngenes,betashape1,betashape2){

    mycandchr<-(1:mynchr)[myngenes>=2]
    if(length(mycandchr)<1){
        warning("fission not possible, skipping event",
                immediate.=TRUE)
        next
    }
    ## draw chromosome
    chr<-sample(mycandchr,1,prob=myngenes[mycandchr])
    ## draw breakpoint and make blocks
    ## draw breakpoint position *(ngenes-1)+1 to have positions
    ##  bound between 1 and ngenes
    mybrpt<-rbeta(1,betashape1,betashape2)*(myngenes[chr]-1) + 1
    if(mybrpt<(myngenes[chr]+1)/2){
        leftp<-floor(mybrpt)
        rightp<-leftp+1
    }else if(mybrpt>(myngenes[chr]+1)/2){
        rightp<-ceiling(mybrpt)
        leftp<-rightp-1
    }else if(floor(mybrpt)==ceiling(mybrpt)){
        leftp<-floor(mybrpt)+rbinom(1,1,0.5)
        rightp<-leftp+1
    }else{
        leftp<-floor(mybrpt)
        rightp<-leftp+1
    }
    myblock1<-genome[[chr]][1:leftp]
    myblock2<-genome[[chr]][rightp:myngenes[chr]]
    if(length(unique(c(myblock1,myblock2)))!=myngenes[chr]){
        stop("Something went wrong when splitting chromosome")
    }
    ## replace unfissioned chromosome by first fissioned chromosome
    genome[[chr]]<-myblock1
    ## add second fissioned chromosome
    genome[[length(genome)+1]]<-myblock2

    return(list(genome=genome,myblock1=myblock1,myblock2=myblock2))
}
