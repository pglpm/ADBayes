## Author: PGL  Porta Mana
## Created: 2021-11-25T14:52:14+0100
## Last-Updated: 2021-11-28T19:53:27+0100
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
## library('coda')
#### End custom setup ####


seed <- 149
baseversion <- 'posterior_'
nclusters <- as.integer(32) #as.integer(2^6)
niter <- as.integer(256) #as.integer(2^11)
niter0 <- as.integer(256) #as.integer(2^10)
thin <- as.integer(2)
nstages <- as.integer(1)
ncheckpoints <- as.integer(4)
maincov <- 'Subgroup_num_'
posterior <- TRUE
##

saveinfofile <- 'variates_info.csv'
variateinfo <- fread(saveinfofile, sep=',')
covNames <- variateinfo$variate
covTypes <- variateinfo$type
covMins <- variateinfo$min
covMaxs <- variateinfo$max
names(covTypes) <- names(covMins) <- names(covMaxs) <- covNames

datafile <- 'data_transformed_shuffled.csv'
alldata <- fread(datafile, sep=',')
alldata <- alldata[Usage_ == 'train']

#################################
## Setup for Monte Carlo sampling
#################################

discreteCovs <- covNames[sapply(covNames, function(x){is.integer(alldata[[x]])})]
continuousCovs <- covNames[sapply(covNames, function(x){is.double(alldata[[x]])})]
covNames <- c(continuousCovs, discreteCovs)
nccovs <- length(continuousCovs)
ndcovs <- length(discreteCovs)
ndata <- as.integer(nrow(alldata))
##
initial.options <- commandArgs(trailingOnly = FALSE)
thisscriptname <- sub('--file=', "", initial.options[grep('--file=', initial.options)])
file.copy(from=thisscriptname, to=paste0('script-R',baseversion,'-V',length(covNames),'-D',ndata,'-K',nclusters,'.Rscript'))


for(obj in c('constants', 'dat', 'inits', 'bayesnet', 'model', 'Cmodel', 'confmodel', 'mcmcsampler', 'Cmcmcsampler')){if(exists(obj)){do.call(rm,list(obj))}}
gc()
##
##
dat <- list(
    X=alldata[, ..continuousCovs],
    Y=alldata[, ..discreteCovs]
)
##
##
source('functions_mcmc.R')
irq2sd <- 1/(2*qnorm(3/4))
alldataRanges <- dataQuantiles <- list()
for(acov in covNames){
        dataQuantiles[[acov]] <- quant(alldata[[acov]], prob=c(0.005,0.995))
        alldataRanges[[acov]] <- range(alldata[[acov]])
}
medianccovs <- apply(alldata[1:ndata,..continuousCovs],2,median)
widthccovs <- irq2sd * apply(alldata[1:ndata,..continuousCovs],2,IQR)
##
mediandcovs <- apply(alldata[1:ndata,..discreteCovs],2,function(x){max(median(x),1)})
widthdcovs <- ceiling(apply(alldata[1:ndata,..discreteCovs],2,IQR))
maxdcovs <- apply(alldata[1:ndata,..discreteCovs],2,max)
thmaxdcovs <- covMaxs[discreteCovs]
matrixmaxdcovs <- matrix(0, nrow=ndcovs, ncol=max(thmaxdcovs), dimnames=list(discreteCovs))
for(acov in discreteCovs){
    matrixmaxdcovs[acov,1:thmaxdcovs[acov]] <- 1/thmaxdcovs[acov]
}
##
print('Creating and saving checkpoints')
checkpointsFile <- paste0('_checkpoints-',ncheckpoints,'-R',baseversion,'-V',length(covNames),'-D',ndata,'-K',nclusters,'.rds')
checkpoints <- rbind(
    c(medianccovs, round(mediandcovs)),
    c(medianccovs+widthccovs, round(mediandcovs+widthdcovs)),
    c(medianccovs-widthccovs, sapply(round(mediandcovs-widthdcovs), function(x){max(0,x)})),
    as.matrix(alldata[sample(1:ndata, size=ncheckpoints), ..covNames])
)
rownames(checkpoints) <- c('Pdatamedians', 'PdatacornerHi', 'PdatacornerLo', paste0('Pdatum',1:ncheckpoints))
saveRDS(checkpoints,file=checkpointsFile)
##
##
constants <- list(
    nClusters=nclusters,
    nData=ndata,
    nCcovs=nccovs,
    nDcovs=ndcovs,
    maxDcovs=ncol(matrixmaxdcovs)
)

##
initsFunction <- function(){
list(
    alphaK=rep(1/nclusters, nclusters),
    meanCmean=medianccovs,
    meanCshape1=rep(1/2, nccovs),
    meanCrate2=1/(widthccovs/2)^2, # dims = inv. variance
    ##
    tauCshape1=rep(1, nccovs),
    tauCrate2=1/(widthccovs/2)^2, # dims = inv. variance
    ##
    probDa1=rep(2, ndcovs),
    probDb1=rep(1, ndcovs),
    sizeDprob1=matrixmaxdcovs,
    ##sizeDsize1=2,
    ##
    ##
    meanCtau1=1/(widthccovs/2)^2, # dims = inv. variance
    meanCrate1=(widthccovs/2)^2, # dims = variance
    tauCrate1=(widthccovs/2)^2, # dims = variance
    ##
    ##
    q=rep(1/nclusters, nclusters),
    meanC=matrix(rnorm(n=nccovs*(nclusters), mean=medianccovs, sd=1/sqrt(rgamma(n=nccovs, shape=1/2, rate=rgamma(n=nccovs, shape=1/2, rate=1/(widthccovs/2)^2)))), nrow=nccovs, ncol=nclusters),
    tauC=matrix(rgamma(n=nccovs*(nclusters), shape=1/2, rate=rgamma(n=nccovs, shape=1/2, rate=1/(widthccovs/2)^2)), nrow=nccovs, ncol=nclusters),
    probD=matrix(rbeta(n=ndcovs*(nclusters), shape1=1, shape2=1), nrow=ndcovs, ncol=nclusters),
    sizeD=matrix(maxdcovs, nrow=ndcovs, ncol=nclusters),
    ## sizeD=matrix(apply(matrix(rcat(n=ndcovs*(nclusters), prob=rbeta(n=ndcovs, shape1=16, shape2=16), size=maxdcovs), nrow=ndcovs, ncol=nclusters), 2, function(x){maxdcovs*(x<maxdcovs)+x*(x>=maxdcovs)}), nrow=ndcovs, ncol=nclusters),
    ##
    C=rep(1,ndata)
)
}

##
##
bayesnet <- nimbleCode({
    q[1:nClusters] ~ ddirch(alpha=alphaK[1:nClusters])
    for(acluster in 1:nClusters){
        for(acov in 1:nCcovs){
            meanC[acov,acluster] ~ dnorm(mean=meanCmean[acov], tau=meanCtau1[acov])
            tauC[acov,acluster] ~ dgamma(shape=tauCshape1[acov], rate=tauCrate1[acov])
        }
        for(acov in 1:nDcovs){
            probD[acov,acluster] ~ dbeta(shape1=probDa1[acov], shape2=probDb1[acov])
            sizeD[acov,acluster] ~ dcat(prob=sizeDprob1[acov,1:maxDcovs])
        }
    }
    ##
    for(acov in 1:nCcovs){
        meanCtau1[acov] ~ dgamma(shape=meanCshape1[acov], rate=meanCrate1[acov])
        meanCrate1[acov] ~ dgamma(shape=meanCshape1[acov], rate=meanCrate2[acov])
        tauCrate1[acov] ~ dgamma(shape=tauCshape1[acov], rate=tauCrate2[acov])
    }
    ##
    if(posterior){
        for(adatum in 1:nData){
            C[adatum] ~ dcat(prob=q[1:nClusters])
        }            ##
        for(adatum in 1:nData){
            for(acov in 1:nCcovs){
                X[adatum,acov] ~ dnorm(mean=meanC[acov,C[adatum]], tau=tauC[acov,C[adatum]])
            }
            for(acov in 1:nDcovs){
                Y[adatum,acov] ~ dbinom(prob=probD[acov,C[adatum]], size=sizeD[acov,C[adatum]])
            }
        }
    }
})
##

timecount <- Sys.time()

if(posterior){
    model <- nimbleModel(code=bayesnet, name='model1', constants=constants, inits=initsFunction(), data=dat)
}else{
    model <- nimbleModel(code=bayesnet, name='model1', constants=constants, inits=initsFunction(), data=list())
    }
Cmodel <- compileNimble(model, showCompilerOutput=FALSE)
gc()
##

    confmodel <- configureMCMC(Cmodel,
                               monitors=c('q','meanC', 'tauC', 'probD', 'sizeD'),
                               monitors2=c('C', 'meanCtau1', 'meanCrate1', 'tauCrate1'))


if(posterior){
    confmodel <- configureMCMC(Cmodel, nodes=NULL,
                               monitors=c('q','meanC', 'tauC', 'probD', 'sizeD'),
                               monitors2=c('C', 'meanCtau1', 'meanCrate1', 'tauCrate1'))
    for(adatum in 1:ndata){
        confmodel$addSampler(target=paste0('C[', adatum, ']'), type='categorical')
    }
    confmodel$addSampler(target=paste0('q[1:', nclusters, ']'), type='conjugate')
    for(acluster in 1:nclusters){
        for(acov in 1:nccovs){
            confmodel$addSampler(target=paste0('meanC[', acov, ', ', acluster, ']'), type='conjugate')
            confmodel$addSampler(target=paste0('tauC[', acov, ', ', acluster, ']'), type='conjugate')
        }
        for(acov in 1:ndcovs){
            confmodel$addSampler(target=paste0('probD[', acov, ', ', acluster, ']'), type='conjugate')
            confmodel$addSampler(target=paste0('sizeD[', acov, ', ', acluster, ']'), type='categorical')
        }
    }
    for(acov in 1:nccovs){
        confmodel$addSampler(target=paste0('meanCtau1[', acov, ']'), type='conjugate')
        confmodel$addSampler(target=paste0('meanCrate1[', acov, ']'), type='conjugate')
        confmodel$addSampler(target=paste0('tauCrate1[', acov, ']'), type='conjugate')
    }
}else{
    confmodel <- configureMCMC(Cmodel,
                               monitors=c('q','meanC', 'tauC', 'probD', 'sizeD'))
}
## confmodel$printSamplers(executionOrder=TRUE)
print(confmodel)

mcmcsampler <- buildMCMC(confmodel)
Cmcmcsampler <- compileNimble(mcmcsampler, resetFunctions = TRUE)
gc()

print('Setup time:')
print(Sys.time() - timecount)

##################################################
## Monte Carlo sampler and plots of MC diagnostics
##################################################
## if(!posterior){
##     nstages <- 1
##     niter0 <- niter <- 64
## }
for(stage in 0:nstages){
    totalruntime <- Sys.time()

    print(paste0('==== STAGE ', stage, ' ===='))
    version <- paste0(baseversion, stage)
    gc()
    if(stage==0){
        mcsamples <- runMCMC(Cmcmcsampler, nburnin=1, niter=niter0+1, thin=1, thin2=niter0, inits=initsFunction, setSeed=seed)
    }else{
        Cmcmcsampler$run(niter=niter*thin, thin=thin, thin2=niter*thin, reset=FALSE, resetMV=TRUE)
    }
    mcsamples <- as.matrix(Cmcmcsampler$mvSamples)
    print('MCMC time:')
    print(Sys.time() - totalruntime)
    ## 7 vars, 6000 data, 100 cl, 2000 iter, slice: 2.48 h
    ## 7 vars, 6000 data, 100 cl, 5001 iter, slice: 6.84 h
    ## 7 vars, 6000 data, 100 cl: rougly 8.2 min/(100 iterations)
    ##
    if(any(is.na(mcsamples))){print('WARNING: SOME NA OUTPUTS')}
    if(any(!is.finite(mcsamples))){print('WARNING: SOME INFINITE OUTPUTS')}
    saveRDS(mcsamples,file=paste0('_mcsamples-R',version,'-V',length(covNames),'-D',ndata,'-K',nclusters,'-I',nrow(mcsamples),'.rds'))
    ## save final state of MCMC chain
    finalstate <- as.matrix(Cmcmcsampler$mvSamples2)
    finalstate <- c(mcsamples[nrow(mcsamples),], finalstate[nrow(finalstate),])
    occupations <- finalstate[grepl('^C\\[', names(finalstate))]
    usedclusters <- length(unique(occupations))
    print(paste0('OCCUPIED CLUSTERS: ', usedclusters, ' OF ', nclusters))
    saveRDS(finalstate2list(finalstate),file=paste0('_finalstate-R',version,'-V',length(covNames),'-D',ndata,'-K',nclusters,'-I',nrow(mcsamples),'.rds'))
    ##
    parmList <- mcsamples2parmlist(mcsamples)
    ## Traces to follow for diagnostics
    ll <- llSamples(dat, parmList)
    if(!posterior && !any(is.finite(ll))){ll <- rep(0, length(ll))}
    momentstraces <- moments12Samples(parmList)
    probCheckpoints <- t(probValuesSamples(checkpoints, parmList))
    ## miqrtraces <- calcSampleMQ(parmList)
    ## medians <- miqrtraces[,,1]
    ## colnames(medians) <- paste0('MEDIAN_', colnames(miqrtraces))
    ## Q1s <- miqrtraces[,,2]
    ## colnames(Q1s) <- paste0('Q1_', colnames(miqrtraces))
    ## Q3s <- miqrtraces[,,3]
    ## colnames(Q3s) <- paste0('Q3_', colnames(miqrtraces))
    ## iqrs <- Q3s - Q1s
    ## colnames(iqrs) <- paste0('IQR_', colnames(miqrtraces))
    traces <- cbind(LL=ll, probCheckpoints, #medians, iqrs, Q1s, Q3s,
                    do.call(cbind, momentstraces))
    badcols <- foreach(i=1:ncol(traces), .combine=c)%do%{if(all(is.na(traces[,i]))){i}else{NULL}}
    if(!is.null(badcols)){traces <- traces[,-badcols]}
    saveRDS(traces,file=paste0('_traces-R',version,'-V',length(covNames),'-D',ndata,'-K',nclusters,'-I',nrow(mcsamples),'.rds'))
    
    ##
    if(nrow(traces)>=1000){
        funMCSE <- function(x){LaplacesDemon::MCSE(x, method='batch.means')$se}
    }else{
        funMCSE <- function(x){LaplacesDemon::MCSE(x)}
    }
    diagnESS <- LaplacesDemon::ESS(traces * (abs(traces) < Inf))
    diagnBMK <- LaplacesDemon::BMK.Diagnostic(traces, batches=2)[,1]
    diagnMCSE <- 100*apply(traces, 2, function(x){funMCSE(x)/sd(x)})
    diagnStat <- apply(traces, 2, function(x){LaplacesDemon::is.stationary(as.matrix(x,ncol=1))})
    diagnBurn <- apply(traces, 2, function(x){LaplacesDemon::burnin(matrix(x[1:(10*trunc(length(x)/10))], ncol=1))})
    ##
    ##
    tracenames <- colnames(traces)
    tracegroups <- list(
        'maxD'=tracenames[grepl('^(Pdat|Dcov|LL)', tracenames)],
        '1D'=tracenames[grepl('^(MEDIAN|Q1|Q3|IQR|MEAN|VAR)_', tracenames)],
        '2D'=tracenames[grepl('^COV_', tracenames)]
    )
    grouplegends <- foreach(agroup=1:length(tracegroups))%do%{
        c( paste0('-- STATS ', names(tracegroups)[agroup], ' --'),
          paste0('min ESS = ', min(diagnESS[tracegroups[[agroup]]])),
          paste0('max BMK = ', max(diagnBMK[tracegroups[[agroup]]])),
          paste0('max MCSE = ', max(diagnMCSE[tracegroups[[agroup]]])),
          paste0('all stationary: ', all(diagnStat[tracegroups[[agroup]]])),
          paste0('burn: ', max(diagnBurn[tracegroups[[agroup]]]))
          )
    }
    colpalette <- sapply(tracenames, function(atrace){
        c(2, 3, 1) %*%
        sapply(tracegroups, function(agroup){atrace %in% agroup})
        })
    ##
    samplesQuantiles <- calcSampleQuantiles(parmList)
    ##
    xlimits <- list()
    for(acov in covNames){
        xlimits[[acov]] <- range(c(alldataRanges[[acov]], samplesQuantiles[,acov,]))
    }
    ##

    ##
    pdff(paste0('mcsummary-R',version,'-V',length(covNames),'-D',ndata,'-K',nclusters,'-I',nrow(mcsamples)))
    matplot(1:2, type='l', col='white', main=paste0('Stats stage ',stage), axes=FALSE, ann=FALSE)
    legendpositions <- c('topleft','bottomleft','topright')
    for(alegend in 1:length(grouplegends)){
        legend(x=legendpositions[alegend], bty='n', cex=1.5,
               legend=grouplegends[[alegend]] )
    }
    legend(x='bottomright', bty='n', cex=1.5,
           legend=c(
               paste0('Occupied clusters: ', usedclusters, ' of ', nclusters),
               paste0('LL: ', signif(mean(ll),3), ' +- ', signif(sd(ll),3))
           ))
    ##
    par(mfrow=c(1,1))
    for(acov in continuousCovs){
        Xgrid <- seq(extendrange(alldata[[acov]])[1], extendrange(alldata[[acov]])[2], length.out=2^8)
        df <- 1/min(diff(Xgrid))/2
        plotsamples <- samplesfX(acov, parmList, Xgrid, nfsamples=min(64,nrow(mcsamples)))
        tplot(Xgrid, plotsamples, type='l', col=paste0(palette()[7], '44'), lty=1, lwd=2, xlab=acov, ylab='probability density', ylim=c(0, max(plotsamples[plotsamples<df])))
    }
    for(acov in discreteCovs){
        if(thmaxdcovs[acov] > 1){
        Xgrid <- seq(max(0,min(alldata[[acov]])-1), max(alldata[[acov]])+1, by=1)
        plotsamples <- samplesfX(acov, parmList, Xgrid, nfsamples=64)
        tplot(Xgrid, plotsamples, type='l', col=paste0(palette()[7], '44'), lty=1, lwd=2, xlab=acov, ylab='probability density', ylim=c(0, max(plotsamples)))
        }else{
            Xgrid <- 0:1
            plotsamples <- samplesfX(acov, parmList, Xgrid, nfsamples=min(512,nrow(mcsamples)))
            histo <- thist(plotsamples[2,])
            tplot(histo$breaks, histo$density, col=7, xlab=paste0('P(',acov,' = 1)'), ylab='probability density', ylim=c(0, max(histo$density)), xlim=c(0,1))
        }
    }
    ##
    par(mfrow = rep(ceiling(sqrt(ndcovs+nccovs)), 2))
    for(addvar in setdiff(covNames, maincov)){
        tplot(x=c(rep(alldataRanges[[maincov]], each=2),
                    alldataRanges[[maincov]][1]),
                y=c(alldataRanges[[addvar]], rev(alldataRanges[[addvar]]),
                    alldataRanges[[addvar]][1]),
                type='l', lwd=2, col=paste0(palette()[2], '88'),
                xlim=xlimits[[maincov]],
                ylim=xlimits[[addvar]],
                xlab=maincov,
                ylab=addvar
                )
        matlines(x=c(rep(dataQuantiles[[maincov]], each=2),
                     dataQuantiles[[maincov]][1]),
                 y=c(dataQuantiles[[addvar]], rev(dataQuantiles[[addvar]]),
                     dataQuantiles[[addvar]][1]),
                 lwd=2, col=paste0(palette()[4], '88'))
    }
    ##
    par(mfrow=c(1,1))
#    matplot(ll, type='l', col=palette()[3], lty=1, main='LL', ylab='LL', ylim=range(ll[abs(ll)<Inf]))
    for(acov in colnames(traces)){
        if(grepl('^[PDV]', acov)){transf <- function(x){log(abs(x)+1e-12)}
        }else{transf <- identity}
        tplot(y=transf(traces[,acov]), type='l', lty=1, col=colpalette[acov],
                main=paste0(acov,
                            '\nESS = ', signif(diagnESS[acov], 3),
                            ' | BMK = ', signif(diagnBMK[acov], 3),
                            ' | MCSE(6.27) = ', signif(diagnMCSE[acov], 3),
                            ' | stat: ', diagnStat[acov],
                            ' | burn: ', diagnBurn[acov]
                            ),
                ylab=acov
              #, ylim=range(c(transf(traces[,acov][abs(transf(traces[,acov]))<Inf])))
              )
    }
    dev.off()

    print('Total runtime:')
    print(Sys.time() - totalruntime)

}
############################################################
## End MCMC
############################################################

