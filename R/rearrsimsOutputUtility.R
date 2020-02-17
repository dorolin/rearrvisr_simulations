## --------------------------------------------------------------------
## OUTPUT UTILITY FUNCTIONS for 'rearrsims.R'
## --------------------------------------------------------------------

## compute measurements of genome fragmentation
fragMeasure<-function(contigs){
    ##contigs<-focalgenome$scaff
    ##contigs<-compgenome$car

    assembly<-length(contigs)
    sizes<-as.vector(sort(table(contigs),decreasing=TRUE))
    largest<-sizes[1]
    ncontigs<-length(sizes)
    N50<-sizes[which(cumsum(sizes)>=assembly/2)][1]

    return(list(ncontigs=ncontigs,N50=N50,largest=largest))
}

## artificially fragment a genome into smaller pieces
##  - size1prop controls when simulations are stopped because
##    a proportion x of scaffolds contains only one marker
##  - betashapes control the skew of breakage of a scaffold
fragmentGenome<-function(mygenome,genomeclass,betashapes=c(1.5,8),
                         nbreaks=10,size1prop=0.5){
    ##mygenome<-focalgenome
    ##genomeclass<-"focal"
    ##mygenome<-compgenome
    ##genomeclass<-"comp"

    if(size1prop<=0 | size1prop>1){
        stop("size1prop must be in interval (0,1]")
    }

    if(genomeclass=="focal"){
        oldcontigs<-mygenome$scaff
    }else if(genomeclass=="comp"){
        oldcontigs<-mygenome$car
    }

    contigs<-integer(length(oldcontigs))
    mycontigs<-unique(oldcontigs)
    for(i in 1:length(mycontigs)){
        contigs[oldcontigs==mycontigs[i]]<-i
    }
    cntr<-max(contigs)

    while(nbreaks>0){

        if(sum(table(contigs)==1)>length(unique(contigs))*size1prop){
            print(paste0("Stopping contig breaking after ",nbreaks,
                         " breaks as >",size1prop*100,"% are of size 1"))
            break
        }

        ## contig to break
        mycont<-sample(1:cntr,size=1)
        contsize<-sum(contigs==mycont)
        if(contsize<2){
            next
        }

        ## draw breakpoint position *(ngenes-1)+1 to have positions
        ##  bound between 1 and ngenes-1 when using floor
        ## (skew is randomly assigned to head or tail of contig)
        s<-sample(1:2)
        mybreak<-floor(rbeta(1,betashapes[s[1]],
                             betashapes[s[2]])*(contsize-1)+1)

        ## replace contig id
        cntr<-cntr+1
        contigs[contigs==mycont][1:mybreak]<-cntr

        nbreaks<-nbreaks-1
    }
    ## plot(as.integer(sort(table(contigs),decreasing=TRUE)),ylab="Size")

    ## adjust genome representation
    if(genomeclass=="focal"){ ## === FOCAL GENOME ===
        newcontigs<-character(length(oldcontigs))
        for(i in 1:length(mycontigs)){
            tmp<-which(oldcontigs==mycontigs[i])
            subids<-unique(contigs[tmp])
            for(j in 1:length(subids)){
                newcontigs[contigs==subids[j]]<-paste(mycontigs[i],j,sep=".")
            }
        }
        ## compute new start and end positions
        mynewcontigs<-unique(newcontigs)
        chrsize<-integer(length(mynewcontigs))
        for(i in 1:length(mynewcontigs)){
            chrsize[i]<-sum(newcontigs==mynewcontigs[i])
        }
        markerstart<-numeric()
        for(i in 1:length(chrsize)){
            markerstart<-c(markerstart,seq(10^6,by=10^6,length.out=chrsize[i]))
        }

        newgenome<-data.frame(
            marker=as.integer(mygenome$marker),
            scaff=as.character(newcontigs),
            start=as.integer(markerstart),
            end=as.integer(markerstart+2),
            strand=mygenome$strand,
            stringsAsFactors = FALSE)
    }else if(genomeclass=="comp"){ ## === COMPARED GENOME ===
        newcontigs<-integer(length(oldcontigs))
        mynewcontigs<-unique(contigs)
        for(i in 1:length(mynewcontigs)){
            newcontigs[contigs==mynewcontigs[i]]<-i
        }
        ## compute new Q-node element ids
        mynewcontigs<-unique(newcontigs)
        chrsize<-integer(length(mynewcontigs))
        for(i in 1:length(mynewcontigs)){
            chrsize[i]<-sum(newcontigs==mynewcontigs[i])
        }
        QnodePos<-integer()
        for(i in 1:length(chrsize)){
            QnodePos<-c(QnodePos,1:(chrsize[i]))
        }

        newgenome<-data.frame(
            marker=mygenome$marker,
            orientation=mygenome$orientation,
            car=as.integer(newcontigs),
            type1=mygenome$type1,
            elem1=as.integer(QnodePos),
            stringsAsFactors = FALSE)
    }

    return(newgenome)
}

precisionRecall<-function(mat1,mat2,predThld=0.5){
    ##mat1<-rearrs[[1]][newpos,]
    ##mat2<-SYNT[[4]]

    if(!is.matrix(mat1) | !is.matrix(mat2)){
        stop("Require matrix as input")
    }
    if(nrow(mat1) != nrow(mat2)){
        stop("Matrices have different numbers of rows")
    }

    realPositive<-as.vector(apply(mat1,1,function(x) sum(x>0)))

    ## including predThld
    predPositive1<-as.vector(apply(mat2,1,function(x) sum(x>=predThld)))
    TP1<-sum(predPositive1>0 & realPositive>0)
    FP1<-sum(predPositive1>0 & realPositive==0)
    FN1<-sum(predPositive1==0 & realPositive>0)

    ## excluding predThld
    predPositive2<-as.vector(apply(mat2,1,function(x) sum(x>predThld)))
    TP2<-sum(predPositive2>0 & realPositive>0)
    FP2<-sum(predPositive2>0 & realPositive==0)
    FN2<-sum(predPositive2==0 & realPositive>0)

    recall.low<-TP2/(TP2+FN2) ## threshold tags are FN but not TP
    recall.high<-TP1/(TP1+FN1) ## threshold tags are TP but not FN

    precision.low<-TP2/(TP2+FP1) ## threshold tags are FP but not TP
    precision.high<-TP1/(TP1+FP2) ## threshold tags are TP but not FP

    return(list(recall.low=recall.low,recall.high=recall.high,
                precision.low=precision.low,precision.high=precision.high))
}

