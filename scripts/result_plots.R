## Author: PGL  Porta Mana
## Created: 2021-11-25T14:52:14+0100
## Last-Updated: 2021-12-02T17:25:57+0100
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
## library('nimble')
#### End custom setup ####

maincov <- 'Subgroup_num_'
source('new_functions_mcmc.R')
frequenciesfile <- 'newIposterior_-V13-D539-K50-I1024/_frequencies-RnewIposterior_2-V13-D539-K50-I1024.rds'
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
subgroupnames <- c('MCI', 'AD')
## subgroup=1 is AD
## subgroup=0 is MCI


datafile <- 'data_transformed_shuffled.csv'
alldata <- fread(datafile, sep=',')
alldata <- alldata[Usage_ == 'train']

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

xcond <- matrix(0:1,ncol=2,dimnames=list(NULL,rep(maincov,2)))

distsFA <- foreach(acov=otherCovs)%do%{
    dists <- rbind(samplesF(Y=grids[[acov]], X=rbind(xcond[,1]), parmList=parmList, inorder=T),
                   samplesF(Y=grids[[acov]], X=rbind(xcond[,2]), parmList=parmList, inorder=T)
                   )
    dim(dists) <- c(length(grids[[acov]]), 2, ncol(dists))
    dimnames(dists) <- list(NULL, subgroupnames, NULL)
    aperm(dists, c(3,1,2))
}
names(distsFA) <- otherCovs

qdistsFA <- foreach(acov=otherCovs)%do%{
    apply(distsFA[[acov]],c(2,3),function(x)c(mean(x,na.rm=T), quant(x=x, probs=c(1,15)/16, na.rm=T)))
}
names(qdistsFA) <- otherCovs

bayesAF <- foreach(acov=otherCovs)%do%{
    dist <- distsFA[[acov]]
    zz <- dist[,,'AD']+dist[,,'MCI']
    dist[,,'AD'] <- dist[,,'AD']/zz
    dist[,,'MCI'] <- dist[,,'MCI']/zz
    dist
}
names(bayesAF) <- otherCovs

qbayesAF <- foreach(acov=otherCovs)%do%{
    apply(bayesAF[[acov]],c(2,3),function(x)c(mean(x,na.rm=T), quant(x=x, probs=c(1,15)/16, na.rm=T)))
}
names(qbayesAF) <- otherCovs




distsAF <- foreach(acov=setdiff(covNames, maincov))%do%{
    dists <- rbind(samplesF(X=grids[[acov]], Y=rbind(xcond[,1]), parmList=parmList, inorder=T),
                   samplesF(X=grids[[acov]], Y=rbind(xcond[,2]), parmList=parmList, inorder=T)
                   )
    dim(dists) <- c(length(grids[[acov]]), 2, ncol(dists))
    dimnames(dists) <- list(NULL, subgroupnames, NULL)
    aperm(dists, c(3,1,2))
}
names(distsAF) <- otherCovs

qdistsAF <- foreach(acov=otherCovs)%do%{
    apply(distsAF[[acov]],c(2,3),function(x)c(mean(x,na.rm=T), quant(x=x, probs=c(1,15)/16, na.rm=T)))
}
names(qdistsAF) <- otherCovs


pdff('plots_features_given_AD')
for(acov in otherCovs){
    agrid <- grids[[acov]]
    ymax <- quant(apply(qdistsFA[[acov]],2,function(x){quant(x,99/100)}),99/100)
ylim <- c(0,ymax)#max(qdistsFA[[acov]]))
for(i in 1:2){
    tplot(x=agrid, y=qdistsFA[[acov]][1,,i],
          col=i, lty=i, lwd=4, alpha=0.25, ylim=ylim,
          xlab=acov, ylab='frequency of feature for patients with AD/MCI', add=(i==2))
        tpar <- unlist(variateinfo[variate==acov,c('transfM','transfW')])
    if(!any(is.na(tpar))){
        Ogrid <- pretty(exp(tpar['transfW']*agrid + tpar['transfM']),n=10)
        axis(3,at=(log(Ogrid)-tpar['transfM'])/tpar['transfW'],labels=Ogrid,lwd=0,lwd.ticks=1,col.ticks='#bbbbbb80')
    }
    polygon(x=c(agrid,rev(agrid)), y=c(qdistsFA[[acov]][2,,i], rev(qdistsFA[[acov]][3,,i])), col=paste0(palette()[i],'40'), border=NA)
}
legend(x=agrid[1], y=ylim[2]*1.2, legend=c(paste0('distribution for patients with ',subgroupnames), '87.5% uncertainty'),
       col=palette()[c(1,2,7)], lty=c(1,2,1), lwd=c(3,3,15), cex=1.5, bty='n', xpd=T
                       )
}
dev.off()



pdff('plots_predictAD_bayes')
for(acov in otherCovs){
    agrid <- grids[[acov]]
    ylim <- c(0,1)
    tplot(x=agrid, y=qbayesAF[[acov]][1,,'AD'],
          col=7, lty=1, lwd=4, ylim=ylim,
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


pdff('plots_predictAD_direct')
for(acov in otherCovs){
    agrid <- grids[[acov]]
ylim <- c(0,1)
    tplot(x=agrid, y=qdistsAF[[acov]][1,,'AD'],
          col=2, lty=1, lwd=4, ylim=ylim,
          xlab=acov, ylab='probability of AD')
    polygon(x=c(agrid,rev(agrid)), y=c(qdistsAF[[acov]][2,,'AD'], rev(qdistsAF[[acov]][3,,'AD'])), col=paste0(palette()[2],'80'), border=NA)
    abline(h=0.5, lty=2, lwd=1, col=2)
legend(x=agrid[1], y=ylim[2]*1, legend=c('87.5% uncertainty on the probability'),
       col=palette()[c(2)], lty=c(1), lwd=c(15), cex=1.5, bty='n'
                       )
}
dev.off()








datafile <- 'data_transformed_shuffled.csv'
testdata <- fread(datafile, sep=',')
testdata <- testdata[Usage_ == 'test']

xcond <- matrix(0:1,ncol=2,dimnames=list(NULL,rep(maincov,2)))
subgroupnames <- c('MCI', 'AD')
## subgroup=1 is AD
## subgroup=0 is MCI

ipredictions0 <- samplesF(X=rbind(xcond[,1]), Y=as.matrix(testdata[,..otherCovs]), parmList=parmList, inorder=T)
ipredictions1 <- samplesF(X=rbind(xcond[,2]), Y=as.matrix(testdata[,..otherCovs]), parmList=parmList, inorder=T)

bpredictions <- foreach(adatum=1:nrow(ipredictions0), .combine=rbind)%do%{
    dist <- ipredictions1[adatum,]/(ipredictions0[adatum,]+ipredictions1[adatum,])
}


predictions <- samplesF(X=rbind(xcond[,2]), Y=as.matrix(testdata[,..otherCovs]), parmList=parmList, inorder=T)

pdff('predictions_bayes_testset')
for(adatum in 1:nrow(testdata)){
    truev <- testdata[[maincov]][adatum]+1
    aprob <- bpredictions[adatum,]
    meanprob <- mean(aprob)
    tcol <- 2
    if((meanprob>=0.5 && truev==2) || (meanprob<=0.5 && truev==1)){tcol <- 1}
    uncprob <- quant(aprob, c(1,15)/16)
    histo <- thist(aprob)
    tplot(x=histo$breaks, y=histo$density,
          xlim=c(0,1), ylim=c(0,NA), col=7, 
          xlab='uncertainty over the probability of AD', ylab='density')
    abline(v=mean(aprob), lty=1, lwd=3, col=tcol)
    abline(v=0.5, lty=3, lwd=2, col=3)
    legend('topleft', legend=c(
                          paste0('true outcome: ',subgroupnames[truev]),
                          paste0('probability of AD: ',signif(meanprob*100,2),'%'),
                          paste0('87.5% uncertainty:\n[',signif(uncprob[1],2),', ',signif(uncprob[2],2),']')
                          ),
           cex=1.5, bty='n')
}
dev.off()

pdff('predictions_testset')
for(adatum in 1:nrow(testdata)){
    truev <- testdata[[maincov]][adatum]+1
    aprob <- predictions[adatum,]
    meanprob <- mean(aprob)
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
           legend=c( paste0('true outcome:\n',subgroupnames[truev])),
           col=truev, cex=1.5, bty='n')
}
dev.off()


pdff('predictions_testset_compbayes')
for(adatum in 1:nrow(testdata)){
    truev <- testdata[[maincov]][adatum]+1
    aprob1 <- predictions[adatum,]
    meanprob1 <- mean(aprob)
    uncprob1 <- quant(aprob, c(1,15)/16)
    histo1 <- thist(aprob)
    aprob <- bpredictions[adatum,]
    meanprob2 <- mean(aprob)
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
           legend=c( paste0('true outcome:\n',subgroupnames[truev])),
           col=truev, cex=1.5, bty='n')
}
dev.off()
