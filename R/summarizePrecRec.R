#!/usr/bin/env Rscript


## summarize and plot precision and recall from rearrvisr simulations

## run with 'Rscript --vanilla R/summarizePrecRec.R'


paramsdir<-"~/Rearrangements/Simulations/simparams"
simoutdir<-"~/Rearrangements/Simulations/simout"


## simparams.txt 1-12
## -------------------
allsims<-paste0("sim",c(9:12,1:8))
allfrags<-rep("",12)

outname<-"sims-1-12"
Mac<-FALSE
## number of rows and columns per page
prow<-3
pcol<-4


## ## simparams.txt 13-24
## ## -------------------
## allsims<-paste0("sim",c(1:4,1:4,1:4))
## allfrags<-rep(c("frag4","frag5","frag6"),each=4)

## outname<-"sims-1-4-frags-4-6"
## Mac<-FALSE
## ## number of rows and columns per page
## prow<-3
## pcol<-4


## ## simparams.txt 25-32
## ## -------------------
## allsims<-paste0("sim",c(13:16,21:24))
## allfrags<-rep("",8)

## outname<-"sims-13-24"
## Mac<-TRUE
## ## number of rows and columns per page
## prow<-2
## pcol<-4



##myq<-0.05 ## compute myq and (1-myq) quantiles
myq<-0.25

## page width (should be 86 or 178 mm)
pwidth<-178
## pointsize (5 or 9 for 86 or 178 mm page works well)
psize<-9

## --------------------------------------------------------------------

## get dimensions
if(grepl("frag",allfrags[1])){
    myin<-paste0(simoutdir,"/",allsims[1],"/",
                 allsims[1],"_",allfrags[1],"_PrecRec.RData")
}else{
    myin<-paste0(simoutdir,"/",allsims[1],"/",
                 allsims[1],"_PrecRec.RData")
}
load(myin)
## loads pr.sum

nrep<-dim(pr.sum)[3]
rm(pr.sum)
nsim<-length(allsims)

## median over replicates
precrecM<-vector("list",nsim)
## lower quantile
precrecL<-vector("list",nsim)
## upper quantile
precrecH<-vector("list",nsim)
## parameters
ParamList<-vector("list",nsim)
FragsList<-vector("list",nsim)


for(s in 1:nsim){
    if(grepl("frag",allfrags[s])){
        myin<-paste0(simoutdir,"/",allsims[s],"/",
                     allsims[s],"_",allfrags[s],"_PrecRec.RData")
        load(paste0(paramsdir,"/",allfrags[s],"_params.RData"))
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
        myin<-paste0(simoutdir,"/",allsims[s],"/",
                     allsims[s],"_PrecRec.RData")
        nbreaks<-0
        fraggenome<-""
    }
    FragsList[[s]]<-list(nbreaks=nbreaks,fraggenome=fraggenome)

    load(myin)
    ## loads pr.sum

    ## load simulation parameters
    load(paste0(paramsdir,"/",allsims[s],"_params.RData"))
    ## among others: nrearr, rearrprob, rearrsize


    precrecM[[s]]<-apply(pr.sum,c(1,2),function(x) median(x,na.rm=TRUE))
    precrecL[[s]]<-apply(pr.sum,c(1,2),function(x) quantile(x,myq,na.rm=TRUE))
    precrecH[[s]]<-apply(pr.sum,c(1,2),function(x) quantile(x,1-myq,na.rm=TRUE))
    precrecM[[s]]<-as.data.frame(precrecM[[s]])
    precrecL[[s]]<-as.data.frame(precrecL[[s]])
    precrecH[[s]]<-as.data.frame(precrecH[[s]])


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

    ParamList[[s]]<-simparams


    ## clean up
    ## --------
    if(grepl("frag",allfrags[s])){
        rm(myfrag,tofrag,betashapes,size1prop)
    }
    rm(mysim,ngenes,nrearr,rearrsize,rearrprob,betashape1,betashape2,
       nbreaks,fraggenome,pr.sum,rsize,simparams,myin)

}




rearrcolors<-c("#00AAD4", ## between.sub1 (NM1)
               "#00AAD4", ## between.sub2 (NM2)
               "#00AAD4", ## nonsyn.sub1 (NM1)
               "#00AAD4", ## nonsyn.sub2 (NM2)
               "#6600FF", ## syn.moves (SM)
               "#FF0000") ## inversions (IV)


## distinguish inclusion or exclusion of retra, fus, fis
rearrborders<-c("#DECD87", ##
                "#DECD87", ##
                "black",
                "black",
                "black",
                "black")

## distinguish no-sub, sub1, sub2
rearrsymbols<-c(24,
                22,
                24,
                22,
                21,
                21)


## rearrnames<-c("between.sub1",
##               "between.sub2",
##               "nonsyn.sub1",
##               "nonsyn.sub2",
##               "syn.moves",
##               "inversions")
rearrnames<-c("Mac + N-syn I",
              "Mac + N-syn II",
              "N-syn I",
              "N-syn II",
              "Syn",
              "Inv")

## subset precrec classes
## ----------------------
allPR<-FALSE
ss<-rep(T,length(rearrnames))
if(allPR!=TRUE){
    ss<-c(T,T,T,T,T,T)
}
if(Mac!=TRUE){
    ss<-ss & c(F,F,T,T,T,T)
    rearrnames<-gsub("Mac + ","",rearrnames,fixed=TRUE)
}
##rearrnames[ss]
## other order for legend
## ----------------------
lorder<-c(1,2,3,4,5,6)

lrearrnames<-rearrnames[lorder][ss[lorder]]
lrearrcolors<-rearrcolors[lorder][ss[lorder]]
lrearrsymbols<-rearrsymbols[lorder][ss[lorder]]
lrearrborders<-rearrborders[lorder][ss[lorder]]

##lbreak<-ceiling(length(lrearrnames)/2)
lbreak<-length(lrearrnames)


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


panel<-rep(LETTERS[1:(prow*pcol)],times=ceiling(nsim/(prow*pcol)))[1:nsim]



setEPS()
postscript(paste0(simoutdir,"/",outname,"_precrec.eps"),
           width=pwidth/25.4,height=((pwidth/25.4)/pcol)*prow,
           pointsize=psize,colormodel="srgb",paper="special")
## 1 Inch = 25.4 millimeters
## ("cmyk" better for printing, but looks dark on screen)



par(mfrow=c(prow,pcol),mar=c(3.0,3.2,2.0,1.8))

for(s in 1:nsim){
    plot(1,1,xlim=c(0,1),ylim=c(0,1),type="n",axes=FALSE,
         xaxs="i",yaxs="i",xlab="",ylab="")
    axis(1,padj= -0.5,cex.axis=1.2)
    axis(2,padj=0.5,cex.axis=1.2)
    abline(v=seq(0.2,1.0,0.2),col="lightgray",lty=3)
    abline(h=seq(0.2,1.0,0.2),col="lightgray",lty=3)
    if(is.element(s,(1:nsim)[rep(c(rep(FALSE,prow-1),TRUE),each=pcol)])){
        mtext("Recall",side=1,line=1.9)
    }
    if(s%%pcol == 1){
        mtext("Precision",side=2,line=1.9)
    }
    mtext(filebase[s],side=3,line=0.5)
    mtext(panel[s],side=3,line=0.5,adj=0,font=2,cex=1.2)

    par(xpd=TRUE)
    segments(x0=precrecM[[s]]$recall.low[ss],
             x1=precrecM[[s]]$recall.low[ss],
             y0=precrecL[[s]]$precision.low[ss],
             y1=precrecH[[s]]$precision.low[ss],
             col=rearrcolors[ss])
    segments(x0=precrecL[[s]]$recall.low[ss],
             x1=precrecH[[s]]$recall.low[ss],
             y0=precrecM[[s]]$precision.low[ss],
             y1=precrecM[[s]]$precision.low[ss],
             col=rearrcolors[ss])
    points(precrecM[[s]]$recall.low[ss],precrecM[[s]]$precision.low[ss],
           pch=rearrsymbols[ss],bg=rearrcolors[ss],col=rearrborders[ss],
           cex=1.2)
    par(xpd=FALSE)


    if(s==1){
        legend(x=0,y=0,xjust=0,yjust=0,
               legend=lrearrnames[1:lbreak],
               pt.bg=lrearrcolors[1:lbreak],
               pch=lrearrsymbols[1:lbreak],
               col=lrearrborders[1:lbreak],
               bty="n",y.intersp=0.9,ncol=1,
               x.intersp=0.6,pt.cex=1.5,cex=1.4)
    }else if(s==2 & lbreak>length(lrearrnames)){
        legend(x=0,y=0,xjust=0,yjust=0,
               legend=lrearrnames[(lbreak+1):length(lrearrnames)],
               pt.bg=lrearrcolors[(lbreak+1):length(lrearrnames)],
               pch=lrearrsymbols[(lbreak+1):length(lrearrnames)],
               col=lrearrborders[(lbreak+1):length(lrearrnames)],
               bty="n",y.intersp=0.9,ncol=1,
               x.intersp=0.6,pt.cex=1.5,cex=1.4)
    }

}

dev.off()

