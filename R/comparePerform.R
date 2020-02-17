#!/usr/bin/env Rscript

## compare performance of rearrvisr, GRIMM, and UniMoG in
##  finding correct number of rearrangement events and breakpoints


## run with 'Rscript --vanilla R/comparePerform.R'



paramsdir<-"~/Rearrangements/Simulations/simparams"
simoutdir<-"~/Rearrangements/Simulations/simout"
grimmoutdir<-"~/Rearrangements/Simulations/GRIMM"
unimogoutdir<-"~/Rearrangements/Simulations/UniMoG"

dcj<-"_st" ## standard DCJ in UniMoG: "_st"; uniform sampling DCJ: ""


## simparams.txt 1-12
## -------------------
allsims<-paste0("sim",c(9:12,1:8))
allfrags<-rep("",12)

outname<-"sims-1-12"

## number of rows and columns per page
prow<-3
pcol<-4


## ## simparams.txt 13-24
## ## -------------------
## allsims<-paste0("sim",c(1:4,1:4,1:4))
## allfrags<-rep(c("frag4","frag5","frag6"),each=4)

## outname<-"sims-1-4-frags-4-6"

## ## number of rows and columns per page
## prow<-3
## pcol<-4


## ## simparams.txt 25-32
## ## -------------------
## allsims<-paste0("sim",c(13:16,21:24))
## allfrags<-rep("",8)

## outname<-"sims-13-24"

## ## number of rows and columns per page
## prow<-2
## pcol<-4



## page width (should be 86 or 178 mm)
pwidth<-86
## pointsize (5 or 9 for 86 or 178 mm page works well)
psize<-5


## ---------
## SUMMARIZE
## ---------

nsim<-length(allsims)

MList<-vector("list",nsim)
BList<-vector("list",nsim)
ParamList<-vector("list",nsim)
FragsList<-vector("list",nsim)

for(s in 1:nsim){
    mysim<-allsims[s]
    myfrag<-allfrags[s]

    ## load simulation parameters
    load(paste0(paramsdir,"/",mysim,"_params.RData"))
    ## among others: nrearr, rearrprob, rearrsize

    if(grepl("frag",myfrag)){
        filebase<-paste0(mysim,"_",myfrag)
        load(paste0(paramsdir,"/",myfrag,"_params.RData"))
        ## loads among others: nbreaks, tofrag
        if(is.element("focal",tofrag) & is.element("comp",tofrag)){
            fraggenome<-"f,c"
        }else if(is.element("focal",tofrag)){
            fraggenome<-"f"
        }else if(is.element("comp",tofrag)){
            fraggenome<-"c"
        }else{
            fraggenome<-""
        }
    }else{
        filebase<-mysim
        nbreaks<-0
        fraggenome<-""
    }

    ## load rearrvisr summary
    load(paste0(simoutdir,"/",mysim,"/",filebase,"_rearrs.RData"))
    ## rs.sum
    rs.sum<-as.data.frame(rs.sum)

    ## read GRIMM summary
    grimm<-read.table(paste0(grimmoutdir,"/",mysim,"/",filebase,
                             "_summary.txt"),
                      header=TRUE,as.is=TRUE)
    ## check
    ## remove trivial flips
    trivial<-which(colnames(grimm)=="Trivialchromosomeflip")
    if(length(trivial)>0){
        grimm<-grimm[,-trivial]
    }
    if(sum(grimm$Distance==rowSums(grimm[,-c(1:4),drop=FALSE]))!=nrow(grimm)){
        print(paste("GRIMM output for", mysim, myfrag, "might be missing some events"))
    }


    ## read UniMoG summary
    unimog<-read.table(paste0(unimogoutdir,"/",mysim,"/",filebase,
                              dcj,"_summary.txt"),
                       header=TRUE,as.is=TRUE)
    ## (note that some Distance entries can be NA for uniform sampling DCJ)



    ## bring simulation parameters in a nicer format
    rsize<-numeric(length(rearrprob))
    for(i in 1:length(rearrprob)){
        if(is.element(names(rearrprob)[[i]],c("retra","fus","fis"))){
            next
        }
        rsize[i]<-median(rnbinom(100000,rearrsize[[i]],rearrprob[[i]]))
    }
    simparams<-rbind(nrearr=unlist(nrearr),rsize=rsize)
    simparams[,simparams[1,]==0]<-0
    simparams<-as.data.frame(simparams)



    ## make matrix with counts of rearrangement classes
    M<-matrix(0,nrow=8,ncol=4)
    rownames(M)<-c("other","fusion","fission","nonsyn",
                   "syn","transl","inv","any")
    colnames(M)<-c("true","rearrvisr","grimm","unimog")
    M<-as.data.frame(M)

    ## true
    M$true[4]<-simparams$nonsynmov[1]
    M$true[5]<-simparams$synmov[1]
    M$true[7]<-simparams$inv[1]
    ## rearrvisr
    M$rearrvisr[4]<-median(rs.sum$nonsyntenic)
    M$rearrvisr[5]<-median(rs.sum$syntenic)
    M$rearrvisr[7]<-median(rs.sum$inversions)
    ## grimm
    grimmevents<-colnames(grimm)
    otherevents<-numeric(nrow(grimm))
    for(o in 5:ncol(grimm)){
        if(grimmevents[o]=="Reversal"){
            M$grimm[7]<-median(grimm$Reversal)
        }else if(grimmevents[o]=="Translocation"){
            M$grimm[6]<-median(grimm$Translocation)
        }else if(grimmevents[o]=="Fission"){
            M$grimm[3]<-median(grimm$Fission)
        }else if(grimmevents[o]=="Fusion"){
            M$grimm[2]<-median(grimm$Fusion)
        }else{
            otherevents<-otherevents+grimm[,o]
        }
    }
    M$grimm[1]<-median(otherevents)
    ## unimog
    M$unimog[8]<-median(unimog$Distance,na.rm=TRUE)



    ## make matrix with counts of breakpoints
    B<-matrix(0,nrow=3,ncol=3)
    rownames(B)<-c("any","int","ext")
    colnames(B)<-c("expected","rearrvisr","grimm")
    B<-as.data.frame(B)

    ## prediced (0 to 3 breakpoints per event, according to class)
    B$expected[1]<-sum(simparams$inv[1]*2,simparams$synmov[1]*3,
                       simparams$invsynmov[1]*3,simparams$nonsynmov[1]*3,
                       simparams$invnonsynmov[1]*3,simparams$retra[1]*2,
                       simparams$fus[1]*0,simparams$fis[1]*1)
    ## rearrvisr
    B$rearrvisr[1]<-median(rs.sum$breakpoints)
    ## grimm
    B$grimm[1]<-median(grimm$Int_bpts + grimm$Ext_bpts)
    B$grimm[2]<-median(grimm$Int_bpts)
    B$grimm[3]<-median(grimm$Ext_bpts)


    ## attach to lists
    ## ---------------
    MList[[s]]<-M
    BList[[s]]<-B
    ParamList[[s]]<-simparams
    FragsList[[s]]<-list(nbreaks=nbreaks,fraggenome=fraggenome)


    ## clean up
    ## --------
    if(grepl("frag",myfrag)){
        rm(tofrag,betashapes,size1prop)
    }
    rm(mysim,ngenes,nrearr,rearrsize,rearrprob,betashape1,betashape2,
       myfrag,nbreaks,fraggenome,filebase,rs.sum,grimm,unimog,rsize,
       simparams,M,B)
}


## ------
## PLOTS
## ------



## find maxima per sets of pcol simulations
maxM<-numeric(nsim)
maxB<-numeric(nsim)
tmpM<-0
tmpB<-0
for(s in 1:nsim){
    tmpM<-max(tmpM,colSums(MList[[s]]))
    tmpB<-max(tmpB,as.matrix(BList[[s]][1,]))
    if(s%%pcol == 0){
        maxM[(s-(pcol-1)):s]<-tmpM
        maxB[(s-(pcol-1)):s]<-tmpB
        tmpM<-0
        tmpB<-0
    }
    if(s==nsim & s%%pcol != 0){
        maxM[(s-(s%%pcol)+1):s]<-tmpM
        maxB[(s-(s%%pcol)+1):s]<-tmpB
    }
}


## make simulation names based on simulation settings
## (sum(inv,synmov,nonsynmov))
## (assuming here that all events had same rearrangement size)
filebase<-character(nsim)
for(s in 1:nsim){
    tmp<-paste0("e",sum(ParamList[[s]]$inv[1],
                        ParamList[[s]]$synmov[1],
                        ParamList[[s]]$nonsynmov[1]),
                ",z",ParamList[[s]]$inv[2]+1)
    ## add macro-rearrangements
    if(sum(ParamList[[s]]$retra[1],
           ParamList[[s]]$fus[1],
           ParamList[[s]]$fis[1])>0){
        tmp<-paste0(tmp,",m(",paste0(c(ParamList[[s]]$retra[1],
                                        ParamList[[s]]$fus[1],
                                        ParamList[[s]]$fis[1]),
                                      collapse=","),")")
    }
    ## add genome fragmentation
    if(grepl("frag",allfrags[s])){
        tmp<-paste0(tmp,",s",FragsList[[s]]$nbreaks,
                    "(",FragsList[[s]]$fraggenome,")")
    }
    filebase[s]<-tmp
}




## counts of rearrangement classes
## -------------------------------

rearrcolors<-c("white","#DECD87","#DEAA87","#00AAD4",
##             "other","fusion","fission", "nonsyn" ,
               "#6600FF","#0066FF","#FF0000","gray90")
##             "syn"    ,"transl" ,"inv"    ,"any")
## rearrnames<-c("Other","Fusion","Fission","Nonsynonymous move",
##               "Synonymous move","Translocation","Inversion","Any")
rearrnames<-c("Other","Fus","Fis","N-syn",
              "Syn","Trans","Inv","Any")



panel<-rep(LETTERS[1:(prow*pcol)],times=ceiling(nsim/(prow*pcol)))[1:nsim]



setEPS()
postscript(paste0(simoutdir,"/",outname,"_rearrcnts.eps"),
           width=pwidth/25.4,height=(((pwidth/25.4)/pcol)*prow)*(3/2),
           pointsize=psize,colormodel="srgb",paper="special")
## 1 Inch = 25.4 millimeters
## ("cmyk" better for printing, but looks dark on screen)

nf<-layout(matrix(1:(prow*pcol),prow,pcol,byrow=TRUE),
           widths=c(rep(2,pcol)),
           heights=c(rep(3,prow)+c(rep(0,(prow-1)),0.7)),
           respect=FALSE)
##layout.show(nf)

op<-par(lwd=0.7) ## reduce overall line width

for(s in 1:nsim){

    if(is.element(s,(1:nsim)[rep(c(rep(FALSE,prow-1),TRUE),each=pcol)])){
        par(mar=c(6,3.5,2,0.5))
    }else{
        par(mar=c(2,3.5,2,0.5))
    }

    nr<-as.numeric(sub("e","",unlist(strsplit(filebase[s],","))[1]))

    opb<-par(lwd=0.6) ## set line width of bar borders
    barplot(height=as.matrix(MList[[s]]),beside=FALSE,space=c(0.2,1.0,0.4,0.4),
            names.arg=c("","","",""),col=NA,ylim=c(0,1.02*maxM[s]),
            axes=FALSE,border=NA)
    par(xpd=FALSE)
    abline(h=seq(nr/3,maxM[s],nr/3),col="gray",lty=3,lwd=1.2)
    par(xpd=TRUE)
    barplot(height=as.matrix(MList[[s]]),beside=FALSE,space=c(0.2,1.0,0.4,0.4),
            names.arg=c("","","",""),col=rearrcolors,ylim=c(0,maxM[s]),
            axes=FALSE,add=TRUE)
    par(opb)

    if(is.element(s,(1:nsim)[rep(c(rep(FALSE,prow-1),TRUE),each=pcol)])){
        text(x=c(1.2,3.2,4.6,6.0)-0.5,y=0,
             labels=c("Simulated  ","rearrvisr  ","GRIMM  ","UniMoG  "),
             srt=60,adj=c(1,1),cex=1.4)
    }
    if(s%%pcol == 1){
        axis(2,at=c(0,seq(nr/3,maxM[s],nr/3)),padj= 0.4,cex.axis=1.2)
        mtext("Number of rearrangements",side=2,line=2.0)
    }
    mtext(filebase[s],side=3,line=0.5)
    ##mtext(panel[s],side=3,line=0.5,adj=0,font=2,cex=1.4)


    if(s==1){
        legend(x= -0.3,y=maxM[s],xjust=0,yjust=0.8,
               legend=rev(rearrnames)[1:2],ncol=1,x.intersp=0.6,
               pt.bg=rev(rearrcolors)[1:2],pch=22,col="black",
               bty="n",pt.cex=1.8,cex=1.4,y.intersp=0.9)
    }else if(s==2){
        legend(x= -0.3,y=maxM[s],xjust=0,yjust=0.8,
               legend=rev(rearrnames)[3:5],ncol=1,x.intersp=0.6,
               pt.bg=rev(rearrcolors)[3:5],pch=22,col="black",
               bty="n",pt.cex=1.8,cex=1.4,y.intersp=0.9)
    }else if(s==3){
        legend(x= -0.3,y=maxM[s],xjust=0,yjust=0.8,
               legend=rev(rearrnames)[6:7],ncol=1,x.intersp=0.6,
               pt.bg=rev(rearrcolors)[6:7],pch=22,col="black",
               bty="n",pt.cex=1.8,cex=1.4,y.intersp=0.9)
    }
}

par(op)

dev.off()





## counts of breakpoints
## ---------------------

brptfill<-c("gray90","white","darkgreen")
brptborder<-c("black","black","black")
brptnames<-c("Any","Int","Ext")


setEPS()
postscript(paste0(simoutdir,"/",outname,"_breakpnts.eps"),
           width=pwidth/25.4,height=(((pwidth/25.4)/pcol)*prow)*(3/1.9),
           pointsize=psize*1.1,colormodel="srgb",paper="special")
## 1 Inch = 25.4 millimeters
## ("cmyk" better for printing, but looks dark on screen)


nf<-layout(matrix(1:(prow*pcol),prow,pcol,byrow=TRUE),
           widths=c(rep(2,pcol)),
           heights=c(rep(3,prow)+c(rep(0,(prow-1)),0.7)),
           respect=FALSE)
##layout.show(nf)

op<-par(lwd=0.7) ## reduce overall line width

for(s in 1:nsim){

    if(is.element(s,(1:nsim)[rep(c(rep(FALSE,prow-1),TRUE),each=pcol)])){
        par(mar=c(6,3.5,2,0.8))
    }else{
        par(mar=c(2,3.5,2,0.8))
    }

    nb<-as.numeric(sub("e","",unlist(strsplit(filebase[s],","))[1]))

    opb<-par(lwd=0.6) ## set line width of bar borders
    barplot(height=as.matrix(BList[[s]][1,]),beside=FALSE,
            space=c(0.2,1.0,0.4),names.arg=c("","",""),
            ylim=c(0,1.02*maxB[s]),axes=FALSE,col=NA,border=NA)
    par(xpd=FALSE)
    abline(h=seq(nb/3,maxB[s],nb/3),col="gray",lty=3,lwd=1.2)
    par(xpd=TRUE)
    barplot(height=as.matrix(BList[[s]][1,]),beside=FALSE,
            space=c(0.2,1.0,0.4),names.arg=c("","",""),
            ylim=c(0,maxB[s]),axes=FALSE,col=brptfill[1],
            border=NA,add=TRUE)
    barplot(height=as.matrix(BList[[s]][2:3,]),beside=FALSE,
            space=c(0.2,1.0,0.4),names.arg=c("","",""),
            ylim=c(0,maxB[s]),axes=FALSE,col=brptfill[2:3],
            border=NA,add=TRUE)
    barplot(height=as.matrix(BList[[s]][1,]),beside=FALSE,
            space=c(0.2,1.0,0.4),names.arg=c("","",""),ylim=c(0,maxB[s]),
            axes=FALSE,col=c(NA,NA,NA),border=brptborder,
            add=TRUE)
    par(opb)
    if(is.element(s,(1:nsim)[rep(c(rep(FALSE,prow-1),TRUE),each=pcol)])){
        text(x=c(1.2,3.2,4.6)-0.5,y=0,
             labels=c("Expected  ","rearrvisr  ","GRIMM  "),
             srt=60,adj=c(1,1),cex=1.4)
    }
    if(s%%pcol == 1){
        axis(2,at=c(0,seq(nb/3,maxB[s],nb/3)),padj= 0.4,cex.axis=1.2)
        mtext("Number of breakpoints",side=2,line=2.0)
    }
    mtext(filebase[s],side=3,line=0.5)
    ##mtext(panel[s],side=3,line=0.5,adj=0,font=2,cex=1.4)

    if(s==1){
        legend(x=1.7,y=0,xjust=0.5,yjust= -0.2,
               legend=brptnames,ncol=1,x.intersp=0.6,
               pt.bg=brptfill,pch=22,col=brptborder,
               bty="o",pt.cex=1.8,cex=1.4,y.intersp=0.9,
               bg="white")
    }

}

par(op)

dev.off()

