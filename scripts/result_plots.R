## Author: PGL  Porta Mana
## Created: 2021-11-25T14:52:14+0100
## Last-Updated: 2021-12-10T18:24:15+0100
################
## Prediction of population frequencies for Alzheimer study
################
if(file.exists("/cluster/home/pglpm/R")){
    .libPaths(c("/cluster/home/pglpm/R",.libPaths()))
}
#### Custom setup ####
## Colour-blind friendly palettes, from https://personal.sron.nl/~pault/
## library('khroma')
## palette(colour('bright')())
## scale_colour_discrete <- scale_colour_bright
## palette(colour('muted')())
library('data.table')
## library('ggplot2')
## library('ggthemes')
## theme_set(theme_bw(base_size=18))
#library('cowplot')
library('png')
library('foreach')
library('doFuture')
library('doRNG')
registerDoFuture()
print('availableCores:')
print(availableCores())
print('availableCores-multicore:')
print(availableCores('multicore'))
if(file.exists("/cluster/home/pglpm/R")){
    plan(multicore, workers=availableCores()-1)
}else{
    plan(multisession, workers=6)
}
## library('LaplacesDemon')
## library('extraDistr')
## library('mvtnorm')
## options(bitmapType='cairo')
## pdff <- function(filename){pdf(file=paste0(filename,'.pdf'),paper='a4r',height=11.7,width=16.5)} # to output in pdf format
## pngf <- function(filename,res=300){png(file=paste0(filename,'.png'),height=11.7*1.2,width=16.5,units='in',res=res,pointsize=36)} # to output in pdf format
library('nimble')
#### End custom setup ####

## Bernoulli distribution
dbernoulli <- function(x, prob, log=FALSE){
    if(log){
        out <- x*log(prob) + (1-x)*log(1-prob)
        out[is.na(out)] <- 0
    }else{
        out <- x*prob + (1-x)*(1-prob)
    }
    out
}

maincov <- 'Subgroup_num_'
source('functions_mcmc.R')
dirname <- 'newposteriorRd_-V13-D539-K64-I1024'
outfile <- paste0(dirname,'/','results.txt')
frequenciesfile <- paste0(dirname,'/','_frequencies-RnewposteriorRd_3-V13-D539-K64-I1024.rds')
parmList <- readRDS(frequenciesfile)
nclusters <- ncol(parmList$q)
nFsamples <- nrow(parmList$q)
realCovs <- dimnames(parmList$meanR)[[2]]
integerCovs <- dimnames(parmList$probI)[[2]]
binaryCovs <- dimnames(parmList$probB)[[2]]
covNames <- c(realCovs, integerCovs, binaryCovs)
nrcovs <- length(realCovs)
nicovs <- length(integerCovs)
nbcovs <- length(binaryCovs)
ncovs <- length(covNames)
saveinfofile <- 'variates_info.csv'
variateinfo <- fread(saveinfofile, sep=',')
covTypes <- variateinfo$type
covMins <- variateinfo$min
covMaxs <- variateinfo$max
names(covTypes) <- names(covMins) <- names(covMaxs) <- variateinfo$variate
otherCovs <- setdiff(covNames, maincov)
diseasenames <- c('MCI', 'AD')
## 0 is MCI
## 1 is AD
gendernames <- c('male', 'female')
## 0 is male
## 1 is female
##
datafile <- 'data_transformed_shuffled.csv'
alldata <- fread(datafile, sep=',')
alldata <- alldata[Usage_ == 'train']

## probc0 <- samplesF(Y=matrix(0,nrow=1,dimnames=list(NULL,maincov)), X=data.matrix(alldata[2,..otherCovs]), parmList=parmList, inorder=T)
## ##
## probc1 <- samplesF(Y=matrix(1,nrow=1,dimnames=list(NULL,maincov)), X=data.matrix(alldata[2,..otherCovs]), parmList=parmList, inorder=T)
## ##
## probj0 <- samplesF(X=NULL, Y=cbind(matrix(0,nrow=1,dimnames=list(NULL,maincov)),data.matrix(alldata[2,..otherCovs])), parmList=parmList, inorder=T)
## probj1 <- samplesF(X=NULL, Y=cbind(matrix(1,nrow=1,dimnames=list(NULL,maincov)),data.matrix(alldata[2,..otherCovs])), parmList=parmList, inorder=T)
## ##
## probx <- samplesF(X=NULL, Y=data.matrix(alldata[2,..otherCovs]), parmList=parmList, inorder=T)

grids <- foreach(acov=covNames)%do%{
    rg <- range(alldata[[acov]])
    if(acov %in% realCovs){
        rg <- rg+c(-1,1)*IQR(alldata[[acov]],type=8)/2
        Xgrid <- seq(rg[1], rg[2], length.out=256)
    }else if(acov %in% integerCovs){
            rg <- round(c((covMins[acov]+7*rg[1])/8, (covMaxs[acov]+7*rg[2])/8))
            Xgrid <- rg[1]:rg[2]
    }else{Xgrid <- 0:1}
    matrix(Xgrid, ncol=1, dimnames=list(NULL,acov))
}
names(grids) <- covNames
##
xcond <- matrix(0:1,ncol=2,dimnames=list(NULL,rep(maincov,2)))

## Frequencies of each feature given AD state
distsFA <- foreach(acov=otherCovs)%do%{
    dists <- rbind(samplesF(Y=grids[[acov]], X=rbind(xcond[,1]), parmList=parmList, inorder=T),
                   samplesF(Y=grids[[acov]], X=rbind(xcond[,2]), parmList=parmList, inorder=T)
                   )
    dim(dists) <- c(length(grids[[acov]]), 2, ncol(dists))
    dimnames(dists) <- list(NULL, diseasenames, NULL)
    aperm(dists, c(3,1,2))
}
names(distsFA) <- otherCovs

## quantiles
qdistsFA <- foreach(acov=otherCovs)%do%{
    apply(distsFA[[acov]],c(2,3),function(x)c(mean(x,na.rm=T), quant(x=x, probs=c(1,15)/16, na.rm=T)))
}
names(qdistsFA) <- otherCovs



## distsFAG <- foreach(acov=setdiff(covNames, c(maincov, 'Gender_num_')))%do%{
##     dists <- rbind(
##         samplesF(Y=grids[[acov]], X=matrix(c(0,0),nrow=1,dimnames=list(NULL,c(maincov,'Gender_num_'))), parmList=parmList, inorder=T),
##         samplesF(Y=grids[[acov]], X=matrix(c(1,0),nrow=1,dimnames=list(NULL,c(maincov,'Gender_num_'))), parmList=parmList, inorder=T),
##         samplesF(Y=grids[[acov]], X=matrix(c(0,1),nrow=1,dimnames=list(NULL,c(maincov,'Gender_num_'))), parmList=parmList, inorder=T),
##         samplesF(Y=grids[[acov]], X=matrix(c(1,1),nrow=1,dimnames=list(NULL,c(maincov,'Gender_num_'))), parmList=parmList, inorder=T)
##         )
##     dim(dists) <- c(length(grids[[acov]]), 2, 2, ncol(dists))
##     dimnames(dists) <- list(NULL, diseasenames, gendernames, NULL)
##     aperm(dists, c(4,1,2,3))
## }
## names(distsFAG) <- setdiff(covNames, c(maincov, 'Gender_num_'))

## ## quantiles
## qdistsFAG <- foreach(acov=names(distsFAG))%do%{
##     apply(distsFAG[[acov]],c(2,3,4),function(x)c(mean(x,na.rm=T), quant(x=x, probs=c(1,15)/16, na.rm=T)))
## }
## names(qdistsFAG) <- names(distsFAG)



## probability of AD state given features, via Bayes's theorem
bayesAF <- foreach(acov=otherCovs)%do%{
    dist <- distsFA[[acov]]
    zz <- dist[,,'AD']+dist[,,'MCI']
    dist[,,'AD'] <- dist[,,'AD']/zz
    dist[,,'MCI'] <- dist[,,'MCI']/zz
    dist
}
names(bayesAF) <- otherCovs

## quantiles
qbayesAF <- foreach(acov=otherCovs)%do%{
    apply(bayesAF[[acov]],c(2,3),function(x)c(mean(x,na.rm=T), quant(x=x, probs=c(1,15)/16, na.rm=T)))
}
names(qbayesAF) <- otherCovs

## frequencies of AD state given features
distsAF <- foreach(acov=setdiff(covNames, maincov))%do%{
    dists <- rbind(samplesF(X=grids[[acov]], Y=rbind(xcond[,1]), parmList=parmList, inorder=T),
                   samplesF(X=grids[[acov]], Y=rbind(xcond[,2]), parmList=parmList, inorder=T)
                   )
    dim(dists) <- c(length(grids[[acov]]), 2, ncol(dists))
    dimnames(dists) <- list(NULL, diseasenames, NULL)
    aperm(dists, c(3,1,2))
}
names(distsAF) <- otherCovs


## quantiles
qdistsAF <- foreach(acov=otherCovs)%do%{
    apply(distsAF[[acov]],c(2,3),function(x)c(mean(x,na.rm=T), quant(x=x, probs=c(1,15)/16, na.rm=T)))
}
names(qdistsAF) <- otherCovs

## plot of frequencies of features given AD state f(F|AD)
pdff(paste0(dirname,'/','plots_features_given_AD2'))
for(acov in otherCovs){
    agrid <- grids[[acov]]
    ymax <- quant(apply(qdistsFA[[acov]],2,function(x){quant(x,99/100)}),99/100)
    ylim <- c(0,ymax)#max(qdistsFA[[acov]]))
    xlim <- c(NA,NA)
    tpar <- unlist(variateinfo[variate==acov,c('transfM','transfW')])
    if(!any(is.na(tpar))){
        xlabels <- pretty(exp(tpar['transfW']*agrid + tpar['transfM']),n=10)
        xticks <- (log(xlabels)-tpar['transfM'])/tpar['transfW']
    }else{xticks <- NULL
        xlabels <- TRUE}
    if(acov %in% binaryCovs){
        xticks <- 0:1
        xlim <- c(-0.25,1.25)
    }
    for(i in 1:2){
        tplot(x=agrid, y=qdistsFA[[acov]][1,,i], yticks=NULL, xlim=xlim,
              col=i, lty=i, lwd=4, alpha=0.25, ylim=ylim, xticks=xticks, xlabels=xlabels,
              xlab=acov, ylab='frequency of feature for patients with AD/MCI', add=(i==2))
        polygon(x=c(agrid,rev(agrid)), y=c(qdistsFA[[acov]][2,,i], rev(qdistsFA[[acov]][3,,i])), col=paste0(palette()[i],'40'), border=NA)
    }
    legend(x=agrid[1], y=ylim[2]*1.2, legend=c(paste0('distribution for patients with ',diseasenames), '87.5% uncertainty'),
           col=palette()[c(1,2,7)], lty=c(1,2,1), lwd=c(3,3,15), cex=1.5, bty='n', xpd=T
           )
}
dev.off()

## plot of samples of frequencies of features given AD state f(F|AD)
pdff(paste0(dirname,'/','plotssamples_features_given_AD2'))
for(acov in otherCovs){
    agrid <- grids[[acov]]
    ymax <- quant(apply(distsFA[[acov]],2,function(x){quant(x,99/100)}),99/100)
    ylim <- c(0,ymax)#max(qdistsFA[[acov]]))
    xlim <- c(NA,NA)
    tpar <- unlist(variateinfo[variate==acov,c('transfM','transfW')])
    if(!any(is.na(tpar))){
        xlabels <- pretty(exp(tpar['transfW']*agrid + tpar['transfM']),n=10)
        xticks <- (log(xlabels)-tpar['transfM'])/tpar['transfW']
    }else{xticks <- NULL
        xlabels <- TRUE}
    if(acov %in% binaryCovs){
        xticks <- 0:1
        xlim <- c(-0.25,1.25)
    }
    subsam <- seq(1,dim(distsFA[[acov]])[1], length.out=128)
    tplot(x=agrid, y=matrix(
                       rbind(t(distsFA[[acov]][subsam,,1]),t(distsFA[[acov]][subsam,,2])),
                       nrow=length(agrid)),
          yticks=NULL, xlim=xlim,
              col=c(5,2), lty=1, lwd=1, alpha=0.75, ylim=ylim, xticks=xticks, xlabels=xlabels,
              xlab=acov, ylab='frequency of feature for patients with AD/MCI')
    ## for(i in 1:2){
    ##     tplot(x=agrid, y=qdistsFA[[acov]][1,,i], yticks=NULL, xlim=xlim,
    ##           col=c(1,6)[i], lty=i, lwd=3, alpha=0.25, ylim=ylim, xticks=xticks, xlabels=xlabels,
    ##           xlab=acov, ylab='frequency of feature for patients with AD/MCI', add=T)
    ## }
    legend(x=agrid[1], y=ylim[2]*1.2, legend=c(paste0('distributions for patients with ',diseasenames)),
           col=palette()[c(1,2,7)], lty=c(1,1), lwd=c(3,3,15), cex=1.5, bty='n', xpd=T
           )
}
dev.off()

## ## plot of frequencies of features given AD state and gender f(F|AD&G)
## pdff('plots_features_given_ADG2')
## for(acov in names(distsFAG)){
##     agrid <- grids[[acov]]
##    ymax <- quant(apply(qdistsFAG[[acov]],2,function(x){quant(x,99/100)}),99/100)
## ##    ymax <- max(qdistsFAG[[acov]])
##     ylim <- c(0,ymax)#max(qdistsFA[[acov]]))
##     xlim <- c(NA,NA)
##     tpar <- unlist(variateinfo[variate==acov,c('transfM','transfW')])
##     if(!any(is.na(tpar))){
##         xlabels <- pretty(exp(tpar['transfW']*agrid + tpar['transfM']),n=10)
##         xticks <- (log(xlabels)-tpar['transfM'])/tpar['transfW']
##     }else{xticks <- NULL
##         xlabels <- TRUE}
##     if(acov %in% binaryCovs){
##         xticks <- 0:1
##         xlim <- c(-0.25,1.25)
##     }
##     tcols <- matrix(c(1,6,5,2),nrow=2)
##     for(i in 1:2){
##         for(j in 1:2){
##         tplot(x=agrid, y=qdistsFAG[[acov]][1,,i,j], yticks=NULL, xlim=xlim,
##               col=tcols[i,j], lty=i, lwd=4, alpha=0.25, ylim=ylim, xticks=xticks, xlabels=xlabels,
##               xlab=acov, ylab='frequency of feature given AD/MCI & gender', add=(i+j>2))
##         polygon(x=c(agrid,rev(agrid)), y=c(qdistsFAG[[acov]][2,,i,j], rev(qdistsFAG[[acov]][3,,i,j])), col=paste0(palette()[tcols[i,j]],'40'), border=NA)
##     }
##     }
##     legend(x=agrid[1], y=ylim[2]*1.2, legend=c(paste0('distribution for patients with ',diseasenames), '87.5% uncertainty'),
##            col=palette()[c(1,2,7)], lty=c(1,2,1), lwd=c(3,3,15), cex=1.5, bty='n', xpd=T
##            )
##     legend(x=agrid[length(agrid)*4/5], y=ylim[2]*1.2, legend=c('darker: male','lighter: female'), bty='n', xpd=T, cex=1.25)
## }
## dev.off()


## plot of frequencies of AD state given features, using bayes f(F|AD)
pdff(paste0(dirname,'/','plots_predictAD_bayes2'))
for(acov in otherCovs){
    agrid <- grids[[acov]]
    tpar <- unlist(variateinfo[variate==acov,c('transfM','transfW')])
    if(!any(is.na(tpar))){
        xlabels <- pretty(exp(tpar['transfW']*agrid + tpar['transfM']),n=10)
        xticks <- (log(xlabels)-tpar['transfM'])/tpar['transfW']
    }else{xticks <- NULL
    xlabels <- TRUE}
    ylim <- c(0,1)
    xlim <- c(NA,NA)
    if(acov %in% binaryCovs){
        xticks <- 0:1
        xlim <- c(-0.25,1.25)
    }
    tplot(x=agrid, y=qbayesAF[[acov]][1,,'AD'],
          col=7, lty=1, lwd=4, ylim=ylim, xlim=xlim, xticks=xticks, xlabels=xlabels,
          xlab=acov, ylab='probability of AD')
    tpar <- unlist(variateinfo[variate==acov,c('transfM','transfW')])
    if(!any(is.na(tpar))){
        Ogrid <- pretty(exp(tpar['transfW']*agrid + tpar['transfM']),n=10)
        axis(3,at=(log(Ogrid)-tpar['transfM'])/tpar['transfW'],labels=Ogrid,lwd=0,lwd.ticks=1,col.ticks='#bbbbbb80')
    }
    polygon(x=c(agrid,rev(agrid)), y=c(qbayesAF[[acov]][2,,'AD'], rev(qbayesAF[[acov]][3,,'AD'])), col=paste0(palette()[7],'80'), border=NA)
    abline(h=0.5, lty=2, lwd=1, col=2)
legend('topleft', legend=c('87.5% uncertainty on the probability'),
       col=palette()[c(7)], lty=c(1), lwd=c(15), cex=1.5, bty='n'
                       )
}
dev.off()

## plot of frequencies of AD state given features, f(AD|F)
pdff(paste0(dirname,'/','plots_predictAD_direct2'))
for(acov in otherCovs){
    agrid <- grids[[acov]]
    tpar <- unlist(variateinfo[variate==acov,c('transfM','transfW')])
    if(!any(is.na(tpar))){
        xlabels <- pretty(exp(tpar['transfW']*agrid + tpar['transfM']),n=10)
        xticks <- (log(xlabels)-tpar['transfM'])/tpar['transfW']
    }else{xticks <- NULL
    xlabels <- TRUE}
ylim <- c(0,1)
    xlim <- c(NA,NA)
    if(acov %in% binaryCovs){
        xticks <- 0:1
        xlim <- c(-0.25,1.25)
    }
    tplot(x=agrid, y=qdistsAF[[acov]][1,,'AD'],
          col=2, lty=1, lwd=4, ylim=ylim, xlim=xlim, xticks=xticks, xlabels=xlabels,
          xlab=acov, ylab='probability of AD')
    polygon(x=c(agrid,rev(agrid)), y=c(qdistsAF[[acov]][2,,'AD'], rev(qdistsAF[[acov]][3,,'AD'])), col=paste0(palette()[2],'80'), border=NA)
    abline(h=0.5, lty=2, lwd=1, col=2)
legend(x=agrid[1], y=ylim[2]*1, legend=c('87.5% uncertainty on the probability'),
       col=palette()[c(2)], lty=c(1), lwd=c(15), cex=1.5, bty='n'
                       )
}
dev.off()

## plot of samples of frequencies of AD state given features, f(AD|F)
pdff(paste0(dirname,'/','plotssamples_predictAD_direct2'))
for(acov in otherCovs){
    agrid <- grids[[acov]]
    tpar <- unlist(variateinfo[variate==acov,c('transfM','transfW')])
    if(!any(is.na(tpar))){
        xlabels <- pretty(exp(tpar['transfW']*agrid + tpar['transfM']),n=10)
        xticks <- (log(xlabels)-tpar['transfM'])/tpar['transfW']
    }else{xticks <- NULL
    xlabels <- TRUE}
ylim <- c(0,1)
    xlim <- c(NA,NA)
    if(acov %in% binaryCovs){
        xticks <- 0:1
        xlim <- c(-0.25,1.25)
    }
    subsam <- seq(1,dim(distsFA[[acov]])[1], length.out=128)
    tplot(x=agrid, y=t(distsAF[[acov]][,,'AD']),
          col=2, lty=1, lwd=1, alpha=0.875, ylim=ylim, xlim=xlim, xticks=xticks, xlabels=xlabels,
          xlab=acov, ylab='probability of AD')
    tplot(x=agrid, y=qdistsAF[[acov]][1,,'AD'],
          col=6, lty=1, lwd=4, ylim=ylim, xlim=xlim, xticks=xticks, xlabels=xlabels,
          xlab=acov, ylab='probability of AD', add=T)
    abline(h=0.5, lty=2, lwd=1, col=2)
## legend(x=agrid[1], y=ylim[2]*1, legend=c('87.5% uncertainty on the probability'),
##        col=palette()[c(2)], lty=c(1), lwd=c(15), cex=1.5, bty='n'
##                        )
}
dev.off()


## ## Sample of features of future datapoints
## datasamples <- foreach(asample=1:nrow(parmList$q), .combine=rbind, .packages='nimble', .inorder=F)%dorng%{
##     acluster <- rcat(n=1,prob=parmList$q[asample,])
##     sapply(covNames,function(acov){
##         if(acov %in% realCovs){
##             rnorm(n=1,mean=parmList$meanR[asample,acov,acluster],sd=1/sqrt(parmList$tauR[asample,acov,acluster]))
##         }else if(acov %in% integerCovs){
##             rbinom(n=1,prob=parmList$probI[asample,acov,acluster],size=parmList$sizeI[asample,acov,acluster])
##         }else{
##             nimble::rcat(n=1,prob=c(1-parmList$probB[asample,acov,acluster],parmList$probB[asample,acov,acluster]))-1
##         }
##     })
## }

## pADdatasamples <- samplesF(Y=matrix(1,nrow=1,dimnames=list(NULL,maincov)),
##                            X=datasamples[,otherCovs], parmList=parmList, inorder=F)

## pADdata <- rowMeans(pADdatasamples)
## pADdata <- abs(pADdata-0.5)+0.5
## ## > summary(pADdata)
## ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## ##  0.5004  0.5431  0.5784  0.5863  0.6151  0.7269 

## pADdatasamplestest <- samplesF(Y=matrix(1,nrow=1,dimnames=list(NULL,maincov)),
##                            X=testdata[,..otherCovs], parmList=parmList, inorder=F)

## pADdatatest <- rowMeans(pADdatasamplestest)
## pADdatatest <- abs(pADdatatest-0.5)+0.5
## ## > summary(pADdatatest)
## ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## ##  0.5027  0.5469  0.5836  0.5900  0.6109  0.7221 

## psinglefeatures <- sapply(otherCovs,function(acov){
##     rowMeans(samplesF(Y=matrix(1,nrow=1,dimnames=list(NULL,maincov)),
##              X=datasamples[,acov,drop=F], parmList=parmList, inorder=F))})

## sort(colMeans(abs(psinglefeatures-0.5)+0.5), decreasing=T)
## ## >   RAVLT_immediate  AVDEL30MIN_neuro    AVDELTOT_neuro LRHHC_n_long_log_ 
## ##         0.6580801         0.6536124         0.6254173         0.5931097 
## ##    TRABSCOR_neuro   CATANIMSC_neuro    TRAASCOR_neuro            Apoe4_ 
## ##         0.5889707         0.5838083         0.5706103         0.5447360 
## ##          AGE_log_    ANARTERR_neuro       Gender_num_       GDTOTAL_gds 
## ##         0.5421207         0.5364131         0.5274789         0.5243736 

## entropysort <- sort(colMeans(psinglefeatures*log2(1/psinglefeatures)+(1-psinglefeatures)*log2(1/(1-psinglefeatures))), decreasing=F)
##   ## RAVLT_immediate  AVDEL30MIN_neuro    AVDELTOT_neuro    TRABSCOR_neuro 
##   ##       0.8883023         0.8934159         0.9291837         0.9536632 
##   ##  TRAASCOR_neuro LRHHC_n_long_log_   CATANIMSC_neuro          AGE_log_ 
##   ##       0.9636331         0.9658517         0.9660364         0.9915779 
##   ##          Apoe4_    ANARTERR_neuro       Gender_num_       GDTOTAL_gds 
##   ##       0.9920768         0.9939879         0.9951186         0.9978416 
## message(paste0(names(entropysort),collapse='\n'))



###########################
## Exploration on test data
###########################
## Load test file
datafile <- 'data_transformed_shuffled.csv'
testdata <- fread(datafile, sep=',')
testdata <- testdata[Usage_ == 'test']

xcond <- matrix(0:1,ncol=2,dimnames=list(NULL,rep(maincov,2)))
diseasenames <- c('MCI', 'AD')
## subgroup=1 is AD
## subgroup=0 is MCI

## predictive probabilities with uncertainty
ipredictions0 <- samplesF(X=rbind(xcond[,1]), Y=as.matrix(testdata[,..otherCovs]), parmList=parmList, inorder=T)
ipredictions1 <- samplesF(X=rbind(xcond[,2]), Y=as.matrix(testdata[,..otherCovs]), parmList=parmList, inorder=T)

## predictive probabilities via Bayes's theorem
bpredictions <- foreach(adatum=1:nrow(ipredictions0), .combine=rbind)%do%{
    dist <- ipredictions1[adatum,]/(ipredictions0[adatum,]+ipredictions1[adatum,])
}


## direct predictive probabilities 
predictions <- samplesF(Y=rbind(xcond[,2]), X=as.matrix(testdata[,..otherCovs]), parmList=parmList, inorder=T)
## predictionsj <- samplesF(Y=cbind(rbind(xcond[,2]),as.matrix(testdata[,..otherCovs])), parmList=parmList, inorder=T)
## predictionsx <- samplesF(Y=as.matrix(testdata[,..otherCovs]), parmList=parmList, inorder=T)
## predictions <- predictionsj/predictionsx
##


##
sink(outfile)

cat('Average uncertainty in test-set predictions, direct:\n',
    mean(rowMeans(predictions,na.rm=T)*log2(1/rowMeans(predictions,na.rm=T))),
    'bit\n')
## > [1] 0.428 bit
cat('SD of uncertainty in test-set predictions, direct:\n',
    sd(rowMeans(predictions,na.rm=T)*log2(1/rowMeans(predictions,na.rm=T)),na.rm=T),
    'bit\n')
## > [1] 0.132 bit
## mean(dbernoulli(x=testdata[[maincov]], prob=rowMeans(predictions,na.rm=T), log=TRUE))
##
cat('\n')
cat('Average uncertainty in test-set predictions, via Bayes:\n',
    mean(rowMeans(bpredictions,na.rm=T)*log2(1/rowMeans(bpredictions,na.rm=T)),na.rm=T),
    'bit\n')
## > [1] 0.420 bit
cat('SD of uncertainty in test-set predictions, via Bayes:\n',
    sd(rowMeans(bpredictions,na.rm=T)*log2(1/rowMeans(bpredictions,na.rm=T)),na.rm=T),
    'bit\n')
## > [1] 0.108 bit
##mean(dbernoulli(x=testdata[[maincov]], prob=rowMeans(bpredictions,na.rm=T), log=TRUE))

sink()



confm <- function(preds, thr=0.5){
c( TP=sum(rowMeans(preds,na.rm=T)>=thr & testdata[[maincov]]==1),
    FP=sum(rowMeans(preds,na.rm=T)>=thr & testdata[[maincov]]==0),
    TN=sum(rowMeans(preds,na.rm=T)<thr & testdata[[maincov]]==0),
    FN=sum(rowMeans(preds,na.rm=T)<thr & testdata[[maincov]]==1)
  )
}

sink(outfile,append=T)

cat('\n\nConfusion matrix test-set, threshold 0.5, direct prediction:\n')
confm(predictions)
## TP FP TN FN 
## 41 27 47 24 
cat('\nConfusion matrix test-set, threshold 0.5, via Bayes:\n')
confm(bpredictions)
## TP FP TN FN 
## 44 32 42 21 

sink()

## tplot(x=agrid <- seq(0,1,length.out=128),
##       y=t(sapply(agrid,confm)))
## legend('top',legend=names(confm()),lty=1:4,col=palette(),bty='n')

## tplot(x=agrid <- seq(0,1,length.out=256),
##       y=colSums(sapply(agrid,confm)[c('FP','FN'),]))


## tplot(x=agrid <- seq(0,1,length.out=256),
##       y=colSums(sapply(agrid,confm)[c('TP','TN'),]))


## tplot(x=agrid <- seq(0,1,length.out=128),
##       y=t(sapply(agrid,bconfm)))
## legend('top',legend=names(bconfm()),lty=1:4,col=palette(),bty='n')

## tplot(x=agrid <- seq(0,1,length.out=256),
##       y=colSums(sapply(agrid,bconfm)[c('FP','FN'),]))

## tplot(x=agrid <- seq(0,1,length.out=256),
##       y=colSums(sapply(agrid,bconfm)[c('TP','TN'),]))

## plot of predictive probabilities for test data
pdff(paste0(dirname,'/','predictions_testset2'))
for(adatum in 1:nrow(testdata)){
    truev <- testdata[[maincov]][adatum]+1
    aprob <- predictions[adatum,]
    meanprob <- mean(aprob, na.rm=T)
    tcol <- 2
    if((meanprob>=0.5 && truev==2) || (meanprob<=0.5 && truev==1)){tcol <- 1}
    uncprob <- quant(aprob, c(1,15)/16)
    histo <- thist(aprob)
    tplot(x=histo$breaks, y=histo$density,
          xlim=c(0,1), ylim=c(0,NA), col=3, 
          xlab='probability of AD', ylab='density')
    abline(v=meanprob, lty=1, lwd=3, col=3)
    legend(x=0.25,y=max(histo$density)*1.15,legend=c(
                         paste0('probability of AD between [',signif(uncprob[1],2),', ',
                                signif(uncprob[2],2),']')),
           cex=1.5, bty='n', xpd=T)
    ## legend('topleft', legend=c(
    ##                       paste0('probability of AD between\n[',signif(uncprob[1],2),', ',signif(uncprob[2],2),']')
    ##                       ),
    ##        cex=1.5, bty='n')
    legend(if(truev==2){pleg <- 'right'}else{pleg <- 'left'},
           legend=c( paste0('true outcome:\n',diseasenames[truev])),
           col=truev, cex=1.5, bty='n')
}
dev.off()

## plot of predictive probabilities for test data using Bayes's theorem
## pdff('predictions_bayes_testset2')
## for(adatum in 1:nrow(testdata)){
##     truev <- testdata[[maincov]][adatum]+1
##     aprob <- bpredictions[adatum,]
##     meanprob <- mean(aprob)
##     tcol <- 2
##     if((meanprob>=0.5 && truev==2) || (meanprob<=0.5 && truev==1)){tcol <- 1}
##     uncprob <- quant(aprob, c(1,15)/16)
##     histo <- thist(aprob)
##     tplot(x=histo$breaks, y=histo$density,
##           xlim=c(0,1), ylim=c(0,NA), col=7, 
##           xlab='uncertainty over the probability of AD', ylab='density')
##     abline(v=mean(aprob), lty=1, lwd=3, col=tcol)
##     abline(v=0.5, lty=3, lwd=2, col=3)
##     legend('topleft', legend=c(
##                           paste0('true outcome: ',diseasenames[truev]),
##                           paste0('probability of AD: ',signif(meanprob*100,2),'%'),
##                           paste0('87.5% uncertainty:\n[',signif(uncprob[1],2),', ',signif(uncprob[2],2),']')
##                           ),
##            cex=1.5, bty='n')
## }
## dev.off()

## Comparison of direct and indirect predictive probabilities
pdff(paste0(dirname,'/','predictions_testset_compbayes2'))
for(adatum in 1:nrow(testdata)){
    truev <- testdata[[maincov]][adatum]+1
    aprob1 <- predictions[adatum,]
    meanprob1 <- mean(aprob, na.rm=T)
    uncprob1 <- quant(aprob, c(1,15)/16)
    histo1 <- thist(aprob)
    aprob <- bpredictions[adatum,]
    meanprob2 <- mean(aprob, na.rm=T)
    uncprob2 <- quant(aprob, c(1,15)/16)
    histo2 <- thist(aprob)
    ymax <- max(histo1$density,histo2$density)
    tplot(x=histo1$breaks, y=histo1$density,
          xlim=c(0,1), ylim=c(0,ymax), col=3, 
          xlab='probability of AD', ylab='density')
    abline(v=meanprob1, lty=1, lwd=3, col=3)
    tplot(x=histo2$breaks, y=histo2$density,
          xlim=c(0,1), ylim=c(0,NA), col=4, 
          xlab='uncertainty over the probability of AD', ylab='density', add=T)
    abline(v=meanprob2, lty=1, lwd=3, col=4)
    ## legend('topleft', legend=c(
    ##                       paste0('probability of AD between\n[',signif(uncprob[1],2),', ',signif(uncprob[2],2),']')
    ##                       ),
    ##        cex=1.5, bty='n')
    legend(x=0.25,y=ymax*1.2,legend=c(
                         paste0('direct prediction [',signif(uncprob1[1],2),', ',
                                signif(uncprob1[2],2),']'),
                         paste0('prediction via bayes [',signif(uncprob2[1],2),', ',
                                signif(uncprob2[2],2),']')),
           col=c(3,4), lty=1, lwd=3, cex=1.5, bty='n', xpd=T)
    legend(if(truev==2){pleg <- 'right'}else{pleg <- 'left'},
           legend=c( paste0('true outcome:\n',diseasenames[truev])),
           col=truev, cex=1.5, bty='n')
}
dev.off()


#########################################################
## Mutual info for next prediction
#########################################################

YX <- samplesX(parmList=parmList)
attr(YX, 'rng') <- NULL
probYX <- rowMeans(samplesF(Y=YX, parmList=parmList),na.rm=T)
probY <- rowMeans(samplesF(Y=YX[,maincov,drop=F], parmList=parmList), na.rm=T)
probX <- rowMeans(samplesF(Y=YX[,otherCovs,drop=F], parmList=parmList), na.rm=T)
mutualinfo <- mean(log2(probYX/(probY*probX)), na.rm=T)

probD <- rowMeans(samplesF(Y=YX[,c(integerCovs,binaryCovs),drop=F], parmList=parmList), na.rm=T)

entropy <- mean(log2(1/probD))

sink(outfile,append=T)
cat('\n\nMutual information between', maincov, 'and all other features:\n',
mutualinfo, 'bit')
## 0.222 bit
sink()

dropmis <- sapply(otherCovs, function(acov){
    probj <- rowMeans(samplesF(Y=YX[,setdiff(colnames(YX),acov),drop=F], parmList=parmList),na.rm=T)
    probm <- rowMeans(samplesF(Y=YX[,setdiff(otherCovs,acov),drop=F], parmList=parmList), na.rm=T)
    mean(log2(probj/(probY*probm)), na.rm=T)
})

sink(outfile,append=T)
cat('\n\nMutual information between', maincov, 'and all other features minus one (negative values are due to roundoff):\n')
print(cbind(sort(dropmis,decreasing=F)))
## TRABSCOR_neuro    0.200 bit
## ANARTERR_neuro    0.205
## TRAASCOR_neuro    0.210
## AVDEL30MIN_neuro  0.214
## RAVLT_immediate   0.215
## GDTOTAL_gds       0.219
## AVDELTOT_neuro    0.219
## Gender_num_       0.220
## AGE_log_          0.221
## LRHHC_n_long_log_ 0.221
## CATANIMSC_neuro   0.222
## Apoe4_            0.223

cat('\n\nRelative differences between mutual information using all features and those using all features minus one, in %:\n')
## relative difference in mutual information without the variate
print(cbind(sort(1-dropmis/mutualinfo,decreasing=T))*100)
## TRABSCOR_neuro      10 %
## ANARTERR_neuro       8
## TRAASCOR_neuro       5
## AVDEL30MIN_neuro     4
## RAVLT_immediate      3
## GDTOTAL_gds          2
## AVDELTOT_neuro       2
## Gender_num_          1
## AGE_log_             1
## LRHHC_n_long_log_    1
## CATANIMSC_neuro      0
## Apoe4_               0
sink()

## #########################################################
## ## 2D plots
## #########################################################

## buildgrid <- function(X, lgrid=128){
##         if(X %in% realCovs){
##             rgx <- range(alldata[[X]])
##             rgx <- rgx + c(-1,1) * diff(rgx)/4
##             lgrid <- seq(rgx[1], rgx[2], length.out=lgrid)
##         }else if(X %in% integerCovs){
##             rgx <- range(alldata[[X]])
##             rgx <- round(rgx + c(-1,1) * diff(rgx)/4)
##             rgx[1] <- max(rgx[1], variateinfo[variate==X,min])
##             rgx[2] <- min(rgx[2], variateinfo[variate==X,max])
##             if(diff(rgx)<lgrid){lgrid <- rgx[1]:rgx[2]}
##             else{lgrid <- round(seq(rgx[1], rgx[2], length.out=lgrid))}
##         }else{
##             lgrid <- 0:1
##         }
##         lgrid
## }

## acov2 <- 'AVDEL30MIN_neuro'#integerCovs[1]
## acov1 <- maincov
## ##
## xgrid <- buildgrid(acov1, 128)
## ygrid <- buildgrid(acov2, 128)
## ##
## grid2d <- cbind(rep(xgrid,length(ygrid)), rep(ygrid, each=length(xgrid)))
## colnames(grid2d) <- c(acov1, acov2)
## ##
## nfsamples <- 32
## fsamples2d <- samplesF(Y=grid2d, parmList=parmList, nfsamples=nfsamples, inorder=F)
## ##
## ##dim(fsamples2d) <- c(length(xgrid), length(ygrid), nfsamples)
## ##
## asample <- 1
## ## ax <- min(diff(xgrid)[1], diff(ygrid)[1])/2
## ## ay <- min(diff(xgrid)[1], diff(ygrid)[1])/2
## ax <- diff(xgrid)[1]/2
## if(acov1 %in% realCovs){ xticks <- NULL }else{ xticks <- xgrid }
## ay <- diff(ygrid)[1]/2
## if(acov2 %in% realCovs){ yticks <- NULL }else{ yticks <- ygrid }
## pmax <- max(fsamples2d[,asample])
## ##
## pdff('_test2dplot')
## plot2dF(xygrid=grid2d, fsamples=rowMeans(fsamples2d,na.rm=T))
## #plot2dF(xygrid=grid2d, fsamples=apply(fsamples2d,1,function(x)diff(quant(x,c(1,15)/16,na.rm=T))))
## plot2dF(xygrid=grid2d, fsamples=apply(fsamples2d,1,function(x)IQR(x,type=8,na.rm=T)))
## dev.off()




## pdff('_test2dplot')
## par(mfrow=c(4,8))
## for(asample in 1:nfsamples){
## plot2dF(xygrid=grid2d, fsamples=fsamples2d[,asample], ticks=F, labs=F, mar=c(0,0,0,0))
##     }
## dev.off()



## pdff('_test2dplot')
## plot2dF(xygrid=grid2d, fsamples=fsamples2d[,1])
## dev.off()

## tplot(x=NA, y=NA, xlim=extendrange(xgrid), ylim=extendrange(ygrid), xlab=acov1, ylab=acov2, xticks=xticks, yticks=yticks)
## for(i in 1:nrow(grid2d)){
##     rat <- fsamples2d[i,asample]/pmax
##     polygon(x=grid2d[i,1]+c(-1,1,1,-1)*ax,
##             y=grid2d[i,2]+c(-1,-1,1,1)*ay,
##             border=gray(1-rat), col=gray(1-rat))
## }
## #tplot(x=NA, y=NA, xlim=extendrange(xgrid), ylim=extendrange(ygrid), xlab=acov1, ylab=acov2,add=T)
## dev.off()







## ##
## pdff('_test2dplot')
## tplot(x=NA, y=NA, xlim=extendrange(xgrid), ylim=extendrange(ygrid), xlab=acov1, ylab=acov2)
## for(i in 1:nrow(grid2d)){
##     rat <- sqrt(fsamples2d[i,asample]/pmax)
##     polygon(x=grid2d[i,1]+c(-1,1,1,-1)*ax*rat,
##             y=grid2d[i,2]+c(-1,-1,1,1)*ay*rat,
##             border='white', col='black')
## }
## dev.off()


## plot2DsamplesF <- function(X, Y, parmList, xgrid=128, ygrid=128){
##     if(length(xgrid)==1){
##         if(X %in% realCovs){
##             rgx <- range(alldata[[X]])
##             rgx <- rgx + c(-1,1) * diff(rgx)/4
##             xgrid <- seq(rgx[1], rgx[2], length.out=xgrid)
##         }else if(X %in% integerCovs){
##             rgx <- range(alldata[[X]])
##             rgx <- round(rgx + c(-1,1) * diff(rgx)/4)
##             rgx[1] <- max(rgx[1], thminicovs[X])
##             rgx[2] <- min(rgx[2], thmaxicovs[X])
##             if(diff(rgx)<xgrid){xgrid <- rgx[1]:rgx[2]}
##             else{xgrid <- round(seq(rgx[1], rgx[2], length.out=xgrid))}
##         }else{
##             xgrid <- 0:1
##         }
##         ##
##         if(Y %in% realCovs){
##             rgy <- range(alldata[[Y]])
##             rgy <- rgy + c(-1,1) * diff(rgy)/4
##             ygrid <- seq(rgy[1], rgy[2], length.out=ygrid)
##         }else if(Y %in% integerCovs){
##             rgy <- range(alldata[[Y]])
##             rgy <- round(rgy + c(-1,1) * diff(rgy)/4)
##             rgy[1] <- max(rgy[1], thminicovs[Y])
##             rgy[2] <- min(rgy[2], thmaxicovs[Y])
##             if(diff(rgy)<ygrid){ygrid <- rgy[1]:rgy[2]}
##             else{ygrid <- round(seq(rgy[1], rgy[2], length.out=ygrid))}
##         }else{
##             ygrid <- 0:1
##         }
##         ##
        


        
##         if(length(ygrid)==1){
##             rgx <- range(alldata[[X]])
##             rgx <- rgx + c(-1,1) * diff(rgx)/4
##             xgrid <- seq(rgx[1], rgx[2], length.out=xgrid)
##         }
## }
