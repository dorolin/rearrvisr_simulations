## --------------------------------------------------------------------
## UTILITY FUNCTIONS for 'rearrsims.R'
## --------------------------------------------------------------------

drawblock<-function(myevent,mynchr,myngenes,
                    betashape1,betashape2,rearrsize,rearrprob){

    if(sum(myngenes<1)>0){
        stop("Require at least one gene on chromosome to draw block")
    }
    ## draw chromosome (prob by size)
    mychr<-sample(1:mynchr,1,prob=myngenes)
    ## draw direction
    updown<-rbinom(1,1,0.5)
    ## draw breakpoint position *(ngenes-1)+1 to have positions
    ##  bound between 1 and ngenes
    mypos<-rbeta(1,betashape1,betashape2)*(myngenes[mychr]-1) + 1
    if(updown==0){
        mystart<-ceiling(mypos)
    }else{
        myend<-floor(mypos)
    }
    ## draw size and end position (results in +1 gene)
    ## block is truncated at end of chromosome if too large
    mysize<-rnbinom(1,rearrsize[[myevent]],rearrprob[[myevent]])
    if(updown==0){
        myend<-mystart+mysize
        myend<-min(c(myend,myngenes[mychr]))
    }else{
        mystart<-myend-mysize
        mystart<-max(c(mystart,1))
    }
    if(mystart<1 | myend>myngenes[mychr]){
        stop("Something went wrong when drawing block")
    }
    if(length(myend-mystart+1)<1){
        stop("Something went wrong when drawing block")
    }
    return(list(chr=mychr,blockstart=mystart,blockend=myend))
}


addtags<-function(rearrs,myevent,myblock){

    pos<-match(as.character(abs(myblock)),rownames(rearrs[[myevent]]))
    ## add tags
    tags<-numeric(nrow(rearrs[[myevent]]))
    tags[pos]<-1
    rearrs[[myevent]]<-cbind(rearrs[[myevent]],tags)

    return(rearrs)
}

addbrpts<-function(brpts,myevent,myblock,brtype){

    if(brtype=="block"){
        col1<-1
        col2<-2
        mar1<-myblock[1]
        mar2<-tail(myblock,n=1L)
    }else if(brtype=="excision"){
        col1<-3
        col2<-4
        mar1<-tail(myblock,n=1L)
        mar2<-myblock[1]
    }else if(brtype=="insertion"){
        col1<-5
        col2<-6
        mar1<-tail(myblock,n=1L)
        mar2<-myblock[1]
    }else{
        stop(paste0("unknown brtype \"",brtype,"\""))
    }
    pos1<-match(as.character(abs(mar1)),rownames(rearrs[[myevent]]))
    pos2<-match(as.character(abs(mar2)),rownames(rearrs[[myevent]]))
    ## add breakpoints of moved block
    if(!is.na(pos1)){
        if(mar1>0){ ## breakpoint at marker head
            brpts[[myevent]][pos1,col1]<-brpts[[myevent]][pos1,col1]+1
        }else{ ## breakpoint at marker tail
            brpts[[myevent]][pos1,col2]<-brpts[[myevent]][pos1,col2]+1
        }
    }
    if(!is.na(pos2)){
        if(mar2>0){ ## breakpoint at marker tail
            brpts[[myevent]][pos2,col2]<-brpts[[myevent]][pos2,col2]+1
        }else{ ## breakpoint at marker head
            brpts[[myevent]][pos2,col1]<-brpts[[myevent]][pos2,col1]+1
        }
    }

    return(brpts)
}

