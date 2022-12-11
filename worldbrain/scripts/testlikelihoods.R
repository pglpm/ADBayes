nfsamples <- nrow(mcsamples)
subsample <- round(seq(1,nfsamples, length.out=64))
##
ad0 <- cbind(0)
ad1 <- cbind(1)
colnames(ad0) <- colnames(ad1) <- predictands

graphics.off()
pdff(paste0(dirname,'prelim_likelihoods-R',basename,'--',mcmcseed,'-',stage),'a4')
for(v in setdiff(unlist(variate),predictands)){#cat(avar)
    contvar <- varinfo[v,'type'] %in% variatetypes[c('R','L','T','S')]
    rg <- varinfo[v,c('plotmin','plotmax')]
    if(contvar){
        Xgrid <- cbind(seq(rg[1], rg[2], length.out=256))
    }else{
        Xgrid <- seq(varinfo[v,'min'], varinfo[v,'max'], length.out=varinfo[v,'n'])
        Xgrid <- cbind(Xgrid[Xgrid >= rg[1] & Xgrid <= rg[2]])
    }
    colnames(Xgrid) <- v
    ##
    plotsamples0 <- samplesFDistribution(Y=Xgrid, X=ad0, mcsamples=mcsamples, varinfo=varinfo, subsamples=round(seq(1,nrow(mcsamples),length.out=nfsamples)), jacobian=TRUE)
    plotsamples1 <- samplesFDistribution(Y=Xgrid, X=ad1, mcsamples=mcsamples, varinfo=varinfo, subsamples=round(seq(1,nrow(mcsamples),length.out=nfsamples)), jacobian=TRUE)
    ##
    datum0 <- data0[Subgroup_num_==0][[v]]
    datum0 <- datum0[!is.na(datum0)]
    datum1 <- data0[Subgroup_num_==1][[v]]
    datum1 <- datum1[!is.na(datum1)]
##
    par(mfrow=c(1,1))
    if(varinfo[v,'type'] == variatetypes['S']){
        interior <- which(Xgrid < varinfo[v,'max'])
        interiordata0 <- which(datum0 < varinfo[v,'max'])
        interiordata1 <- which(datum1 < varinfo[v,'max'])
    }else{
        interior <- 1:length(Xgrid)
        interiordata0 <- 1:length(datum0)
        interiordata1 <- 1:length(datum1)
    }
    ymax <- max(tquant(apply(plotsamples0[interior,subsample],2,function(x){tquant(x,31/32)}),31/32, na.rm=T), tquant(apply(plotsamples1[interior,subsample],2,function(x){tquant(x,31/32)}),31/32, na.rm=T))
    ##
    ##
    tplot(x=Xgrid[interior], y=plotsamples0[interior,subsample], type='l', col=5, alpha=7/8, lty=1, lwd=2, xlab=paste0(v), ylab=paste0('frequency'), ylim=c(0, ymax), xlim=range(Xgrid), family=family)
    tplot(x=Xgrid[interior], y=plotsamples1[interior,subsample], type='l', col=2, alpha=7/8, lty=1, lwd=2, xlab=paste0(v), ylab=paste0('frequency'), ylim=c(0, ymax), xlim=range(Xgrid), family=family,add=T)
        ##
    tplot(x=(Xgrid[interior]), y=list(rowMeans(plotsamples0, na.rm=T)[interior],rowMeans(plotsamples1, na.rm=T)[interior]), type='l', col=c(1,6), alpha=0.25, lty=c(1,2), lwd=4, add=T)
    legend
    ##
    ##
    if(FALSE){
    histo <- thist(datum0[interiordata0], n=round(length(interiordata0)/64))
    histomax <- max(rowMeans(plotsamples)[interior])/max(histo$density)
    tplot(x=histo$breaks, y=histo$density*histomax, col=6, alpha=15/16, border=NA, lty=1, lwd=1, family=family, ylim=c(0,NA), add=TRUE)
    ##
    histo <- thist(datum1[interiordata1], n=round(length(interiordata1)/64))
    histomax <- max(rowMeans(plotsamples)[interior])/max(histo$density)
    tplot(x=histo$breaks, y=histo$density*histomax, col=5, alpha=15/16, border=NA, lty=1, lwd=1, family=family, ylim=c(0,NA), add=TRUE)
    ##
    }
}
dev.off()



    ##
    if(length(interior) < length(Xgrid)){
    tplot(x=list(Xgrid[-interior]), y=list(plotsamples0[-interior,subsample],plotsamples1[-interior,subsample]), type='l', col=c(5,6), alpha=0.75, lty=1, lwd=2, xlab=paste0(v), ylab=paste0('frequency'), ylim=c(0, ymax), family=family)
        ##
    tplot(x=list(Xgrid[interior]), y=list(rowMeans(plotsamples0, na.rm=T)[interior],rowMeans(plotsamples1, na.rm=T)[interior]), type='l', col=c(1,2), alpha=0.25, lty=c(1,2), lwd=4, add=T)

        
    }



    
            pborder <- sum(datum <= varinfo[v,'min'])/length(datum)
            if(pborder > 0){
                tplot(x=varinfo[v,'min'], y=pborder*ymax, type='p', pch=0, cex=2, col=7, alpha=0, lty=1, lwd=5, family=family, ylim=c(0,NA), add=TRUE)
            }
            ##
            pborder <- sum(datum >= varinfo[v,'max'])/length(datum)
            if(pborder > 0){
                tplot(x=varinfo[v,'max'], y=pborder*ymax, type='p', pch=0, cex=2, col=7, alpha=0, lty=1, lwd=5, family=family, ylim=c(0,NA), add=TRUE)
            }




    interior <- which(Xgrid > varinfo[v,'min'] & Xgrid < varinfo[v,'max'])
        tplot(x=listXgrid[interior], y=plotsamples[interior,subsample], type='l', col=5, alpha=0.75, lty=1, lwd=2, xlab=paste0(v), ylab=paste0('frequency'), ylim=c(0, ymax), family=family)
        if(length(interior) < length(Xgrid)){
            tplot(x=Xgrid[-interior], y=plotsamples[-interior,subsample,drop=F]*ymax, type='p', pch=2, cex=2, col=5, alpha=0.75, lty=1, lwd=2, xlab=paste0(v), ylab=paste0('frequency'), ylim=c(0, ymax), family=family,add=T)
        }
        if(plotmeans){
            tplot(x=Xgrid[interior], y=rowMeans(plotsamples, na.rm=T)[interior], type='l', col=1, alpha=0.25, lty=1, lwd=3, add=T)
            if(length(interior) < length(Xgrid)){
                tplot(x=Xgrid[-interior], y=rowMeans(plotsamples, na.rm=T)[-interior]*ymax, type='p', pch=2, cex=2, col=1, alpha=0.25, lty=1, lwd=3, add=T)
            }
        }
    }
    ##
    if((showdata=='histogram')||(showdata==TRUE && !contvar)){
        datum <- data0[[v]]
        datum <- datum[!is.na(datum)]
        ## fiven <- varinfo[v,c('min','Q1','Q2','Q3','max')]
        fiven <- fivenum(datum)
        if(varinfo[v,'type'] != variatetypes['S']){
            ## histo <- thist(datum, n=(if(contvar){min(max(10,sqrt(ndata)),100)}else{'i'}))#-exp(mean(log(c(round(sqrt(length(datum))), length(Xgrid))))))
            histo <- thist(datum, n=round(ndata/64))#-exp(mean(log(c(round(sqrt(length(datum))), length(Xgrid))))))
            histomax <- max(rowMeans(plotsamples))/max(histo$density)
            tplot(x=histo$breaks, y=histo$density*histomax, col=grey, alpha=0.75, border=NA, lty=1, lwd=1, family=family, ylim=c(0,NA), add=TRUE)
        }else{ # histogram for censored variate
            interior <- which(datum > varinfo[v,'min'] & datum < varinfo[v,'max'])
            histo <- thist(datum[interior], n=round(length(interior)/64))
            interiorgrid <- which(Xgrid > varinfo[v,'min'] & Xgrid < varinfo[v,'max'])
            histomax <- max(rowMeans(plotsamples)[interiorgrid])/max(histo$density)
            tplot(x=histo$breaks, y=histo$density*histomax, col=7, alpha=0.75, border=NA, lty=1, lwd=1, family=family, ylim=c(0,NA), add=TRUE)
            ##
            pborder <- sum(datum <= varinfo[v,'min'])/length(datum)
            if(pborder > 0){
                tplot(x=varinfo[v,'min'], y=pborder*ymax, type='p', pch=0, cex=2, col=7, alpha=0, lty=1, lwd=5, family=family, ylim=c(0,NA), add=TRUE)
            }
            ##
            pborder <- sum(datum >= varinfo[v,'max'])/length(datum)
            if(pborder > 0){
                tplot(x=varinfo[v,'max'], y=pborder*ymax, type='p', pch=0, cex=2, col=7, alpha=0, lty=1, lwd=5, family=family, ylim=c(0,NA), add=TRUE)
            }
        }
        abline(v=fiven,col=paste0(palette()[c(2,4,5,4,2)], '44'),lwd=4)
    }else if((showdata=='scatter')|(showdata==TRUE & contvar)){
        datum <- data0[[v]]
        datum <- datum[!is.na(datum)]
        diffdatum <- c(apply(cbind(c(0,diff(datum)),c(diff(datum),0)),1,min))/2
        scatteraxis(side=1, n=NA, alpha='88', ext=8, x=rnorm(length(datum),mean=datum,sd=diffdatum),col=yellow)
    }
}
dev.off()
