## Author: PGL  Porta Mana
## Created: 2022-09-08T17:03:24+0200
## Last-Updated: 2022-12-05T07:58:22+0100
#########################################
## Inference of exchangeable variates (nonparametric density regression)
## using effectively-infinite mixture of product kernels
## Monte Carlo sampling
#########################################


#### USER INPUTS AND CHOICES ####
baseversion <- '_testfun2' # *** ## Base name of output directory
datafile <- 'data_ep.csv' #***
predictors <- 'predictors.csv'
varinfofile <- 'metadata_realint_noSW.csv' #***
requiredESS <- 1024*2/20 # required effective sample size
nsamples <- 8*ceiling((requiredESS*1.5)/8) # number of samples AFTER thinning
## ndata <- 5 # set this if you want to use fewer data
## shuffledata <- TRUE # useful if subsetting data
posterior <- TRUE # if set to FALSE it samples and plots prior samples
minstepincrease <- 8L
savetempsamples <- FALSE # save temporary MCMC samples
plottempdistributions <- FALSE # plot temporary sampled distributions
showdata <- TRUE # 'histogram' 'scatter' FALSE TRUE
plotmeans <- FALSE # plot frequency averages
##
niter0 <- 1024L * 1L # 3L # iterations burn-in
nclusters <- 64L
alpha <- 1
compoundgamma <- TRUE # use beta-prime distribution for variance instead of gamma
compoundgammapars <- c(1,1)/2
categoryprior <- 1 # choices: 'Haldane' (1/n) or a number
casualinitvalues <- FALSE
## stagestart <- 3L # set this if continuing existing MC = last saved + 1
family <- 'Palatino'
####


#### Packages and setup ####
## load customized plot functions
if(!exists('tplot')){source('~/work/pglpm_plotfunctions.R')}
##
## Read MCMC seed from command line
mcmcseed = as.integer(commandArgs(trailingOnly=TRUE))[1]
if(is.na(mcmcseed) | (!is.na(mcmcseed) & mcmcseed <=0)){mcmcseed <- 1}
cat(paste0('\nMCMC seed = ',mcmcseed,'\n'))
##
set.seed(701+mcmcseed)
##
## Packages
library('data.table')
library('png')
library('foreach')
library('doFuture')
library('doRNG')
registerDoFuture()
cat('\navailableCores: ')
cat(availableCores())
cat('\navailableCores-multicore: ')
cat(availableCores('multicore'))
if(Sys.info()['nodename']=='luca-HP-Z2-G9'){
    ncores <- 1}else{
    ncores <- 6}
cat(paste0('\nusing ',ncores,' cores\n'))
if(ncores>1){
    if(.Platform$OS.type=='unix'){
        plan(multicore, workers=ncores)
    }else{
        plan(multisession, workers=ncores)
    }
}else{
    plan(sequential)
}
library('nimble')
## NB: also requires libraries 'LaplacesDemon' and 'extraDistr'


#### EXTRACT INFORMATION ABOUT THE VARIATES AND THEIR PRIOR PARAMETERS

if(Sys.info()['nodename']=='luca-HP-Z2-G9'){origdir <- '../'}else{origdir <- ''}
source(paste0(origdir,'functions_mcmc.R')) # load functions for post-MCMC calculations
variateinfo <- fread(paste0(origdir,variateinfofile))
## variateinfo <- do.call(rbind, list(
##     ##	 'variate',			'type',		'min',	'max',	'precision')
##     data.table('REANAME',		'real',		-100,	100,	NA)
## ,   data.table('BINNAME',		'binary',	0,	1,	NA)
## ,   data.table('INTNAME',		'integer',	0,	9,	NA)
## ,   data.table('CATNAME',		'category',	1,	10,	NA)
## ))
## colnames(variateinfo) <- c('variate', 'type', 'min', 'max', 'precision')

## Effects of shape parameter:
## 1/8 (broader):
## > testdata <- log10(rinvgamma(n=10^7, shape=1/8, scale=1^2))/2 ; 10^sort(c(quantile(testdata, c(1,7)/8), summary(testdata)))
##         Min.        12.5%      1st Qu.       Median         Mean      3rd Qu. 
## 2.878554e-01 1.937961e+00 3.903828e+00 2.034059e+01 6.641224e+01 3.260648e+02 
##        87.5%         Max. 
## 5.220388e+03 5.266833e+27 
##
## 1/4:
## > testdata <- log10(rinvgamma(n=10^7, shape=1/4, scale=1^2))/2 ; 10^sort(c(quantile(testdata, c(1,7)/8), summary(testdata)))
##         Min.        12.5%      1st Qu.       Median         Mean      3rd Qu. 
## 2.881895e-01 1.274590e+00 1.960271e+00 4.784803e+00 8.283664e+00 1.946374e+01 
##        87.5%         Max. 
## 7.795913e+01 4.670417e+14
##
## 1/2 (narrower):
## > testdata <- log10(rinvgamma(n=10^7, shape=1/2, scale=1^2))/2 ; 10^sort(c(quantile(testdata, c(1,7)/8), summary(testdata)))
##         Min.        12.5%      1st Qu.       Median         Mean      3rd Qu. 
## 2.571008e-01 9.218967e-01 1.229980e+00 2.098062e+00 2.670157e+00 4.440326e+00 
##        87.5%         Max. 
## 8.995125e+00 6.364370e+06 

##
varNames <- variateinfo$variate
varTypes <- variateinfo$type
varMins <- variateinfo$min
varMaxs <- variateinfo$max
names(varTypes) <- names(varMins) <- names(varMaxs) <- varNames
realVars <- varNames[varTypes=='real']
categoryVars <- varNames[varTypes=='categorical']
binaryVars <- varNames[varTypes=='binary']
##
nrvars <- length(realVars)
ncvars <- length(categoryVars)
nbvars <- length(binaryVars)
nvars <- length(varNames)
##
realNums <- sort(unclass(factor(realVars))); names(realNums) <- realVars
categoryNums <- sort(unclass(factor(categoryVars))); names(categoryNums) <- categoryVars
binaryNums <- sort(unclass(factor(binaryVars))); names(binaryNums) <- binaryVars
varNums <- c(realNums, integerNums, categoryNums, binaryNums)
typeNums <- c('real'=0, 'categorical'=1, 'binary'=2, 'integer'=3)

###########################################
## READ DATA AND SETUP SOME HYPERPARAMETERS
###########################################

alldata <- fread(paste0(origdir,datafile), sep=',')
if(!all(varNames %in% names(alldata))){cat('\nERROR: variates missing from datafile')}
alldata <- alldata[, ..varNames]
## shuffle data
if(exists('shuffledata') && shuffledata){alldata <- alldata[sample(1:nrow(alldata))]}
if(!exists('ndata') || is.null(ndata) || is.na(ndata)){ndata <- nrow(alldata)}
alldata <- alldata[1:ndata]

basename <- paste0(baseversion,'-V',length(varNames),'-D',ndata,'-K',nclusters,'-I',nsamples)
##
if(Sys.info()['nodename']=='luca-HP-Z2-G9'){
    dirname <- ''
}else{
    dirname <- paste0(basename,'/')
    dir.create(dirname)
}

## Copy this script to output directory
## if(Sys.info()['nodename']=='luca-HP-Z2-G9'){
##     initial.options <- commandArgs(trailingOnly = FALSE)
##     thisscriptname <- sub('--file=', "", initial.options[grep('--file=', initial.options)])
##     if(mcmcseed==1){file.copy(from=thisscriptname, to=paste0(dirname,'/script-R',baseversion,'-V',length(varNames),'-D',ndata,'-K',nclusters,'.Rscript'))
##     }
## }

#### Get missing metadata from data
## categories
if(any(is.na(varMins[categoryVars]))){
    cat('\nWARNING: missing min for some category variates; using min from data')
    nanames <- categoryVars[is.na(varMins[categoryVars])]
    for(avar in nanames){
        varMins[avar] <- min(alldata[[avar]], na.rm=T)
    }
}
if(any(is.na(varMaxs[categoryVars]))){
    cat('\nWARNING: missing max for some category variates; using max from data')
    nanames <- categoryVars[is.na(varMaxs[categoryVars])]
    for(avar in nanames){
        varMaxs[avar] <- max(alldata[[avar]], na.rm=T)
    }
}



## Normalization and standardization of real variates and calculation of hyperparams
sd2iqr <- 0.5/qnorm(0.75)
##
if(mcmcseed==1){
    graphics.off()
    pdff(paste0(dirname,'densities_variances'),'a4')
}
varinfo <- NULL
for(avar in varNames){
    pars <- c(NA,NA)
    ##
    datum <- alldata[[avar]]
    datamin <- min(datum, na.rm=TRUE)
    datamax <- max(datum, na.rm=TRUE)
    datamedian <- quantile(datum, 0.5, type=8, na.rm=TRUE)
    dataQ1 <- quantile(datum, 0.25, type=8, na.rm=TRUE)
    dataQ2 <- quantile(datum, 0.75, type=8, na.rm=TRUE)
    ##
    if(avar %in% realVars){
        alocation <- median(datum)
        ascale <- IQR(datum) * sd2iqr
        if(ascale==0){ascale <- 1L}
        ## shift and rescale data
        datum <- (datum-alocation)/ascale
        rg <- diff(range(datum, na.rm=T))
        if(is.na(variateinfo[variate==avar, precision])){
            dmin <- min(diff(sort(unique(datum))))
        }else{
            dmin <- variateinfo[variate==avar, precision]/ascale
        }
        ##
        ## ## > ss <- 2; ss2 <- 1/2 ; test <- log10(rinvgamma(1e6, shape=ss, rate=rinvgamma(1e6, shape=ss2, scale=1)))/2; htest <- thist(test); testg <- seq(min(htest$breaks),max(htest$breaks),length.out=256) ; testo <- extraDistr::dbetapr(exp(log(10)*(-2*testg)),shape1=ss,shape2=ss2)*2*log(10)*exp(log(10)*(-2*testg)) ; tplot(list(htest$mids,testg),list(htest$density,testo))
        ## qts <- c(2^-14, 1-2^-14)
        ## fn <- function(p, target){
        ##     sum((log(extraDistr::qbetapr(qts,shape1=p[1],shape2=p[2]))/2 -log(target))^2)
        ##     }
        ## resu <- optim(par=c(2,1/2), fn=fn, target=c(dmin/2,rg*3))
        ## pars <- signif(resu$par, 3)
        ## vals <- sqrt(extraDistr::qbetapr(qts, shape1=pars[1], shape2=pars[2]))
        ## if(abs(vals[1] - dmin/2)/(vals[1] + dmin/2)*200 > 5 |
        ##    abs(vals[2] - rg*3)/(vals[2] + rg*3)*200 > 5
        ##    ){cat(paste0('WARNING ', avar, ': bad parameters'))}
        ##
        ## Parameters and plot
        sgrid <- seq(log10(dmin/2), log10(rg*3), length.out=256)
        vv <- exp(log(10)*2*sgrid)
        ##
        if(compoundgamma){
            pars <- compoundgammapars
            ##pars <- c(1,1)
            test <- thist(log10(rinvgamma(1e6, shape=pars[2], rate=rinvgamma(1e6, shape=pars[1], scale=1)))/2)
            vdistr <- extraDistr::dbetapr(x=vv, shape1=pars[1], shape2=pars[2])*vv*2*log(10)
        }else{
            pars <- c(signif(2^4.5,3), 1/2)
            test <- thist(log10(rinvgamma(1e6, shape=pars[2], rate=pars[1]))/2)
            vdistr <- dinvgamma(x=vv, rate=pars[1], shape=pars[2])*vv*2*log(10)
        }
        ##
if(mcmcseed==1){
    tplot(x=list(test$mids,sgrid), y=list(test$density,vdistr),
          xlab=expression(lg~sigma), ylim=c(0,NA), ylab=NA,
          main=paste0(avar, ': shape1/rate = ', pars[1],', shape2/shape = ',pars[2]),
          xlim=log10(c(dmin/2/10, rg*3*10)))
    abline(v=log10(c(dmin/2, dmin, IQR(datum)*sd2iqr, rg, rg*3)), col=yellow, lwd=3)
}
        ##
        rm('test','sgrid','vdistr')
    }else if((avar %in% integerVars) | (avar %in% binaryVars)){
        alocation <- varMins[avar] # integer and binary start from 0
        ascale <- 1L
        dmin <- 1L
        ## shift and rescale data
        datum <- (datum-alocation)
    }else if(avar %in% categoryVars){
        alocation <- varMins[avar] - 1L # category start from 1
        ascale <- 1L
        dmin <- 1L
        ## shift and rescale data
        datum <- (datum-alocation)
    }
    ##
    varinfo <- rbind(varinfo,
                               cbind(type=typeNums[varTypes[avar]],
                                     index=varNums[avar],
                                     location=alocation, scale=ascale,
                                     precision=dmin,
                                     min=min(datum, na.rm=T), max=max(datum, na.rm=T),
                                     shape1rate=pars[1], shape2shape=pars[2],
                                     datamin=datamin, datamax=datamax,
                                     datamedian=datamedian,
                                     dataQ1=dataQ1, dataQ2=dataQ2,
                                     thmin=varMins[avar], thmax=varMaxs[avar]
                                     ))
}
rownames(varinfo) <- varNames
##
if(mcmcseed==1){
    dev.off()
    fwrite(cbind(data.table(variate=varNames),varinfo), file=paste0(dirname,'varinfo.csv'))
}
##
locvarMins <- varMins-varinfo[,'location']
locvarMaxs <- varMaxs-varinfo[,'location']


#################################
## Setup for Monte Carlo sampling
#################################

if(!exists('stagestart')){stagestart <- 0L}
if(stagestart>0){
    resume <- paste0('_finalstate-R',baseversion,'-V',length(varNames),'-D',ndata,'-K',nclusters,'-I',nsamples,'--',stagestart-1,'-',mcmcseed,'.rds')
}else{
    resume <- FALSE
}
##


for(obj in c('constants', 'dat', 'inits', 'finitemix', 'finitemixnimble', 'Cfinitemixnimble', 'confnimble', 'mcsampler', 'Cmcsampler')){if(exists(obj)){do.call(rm,list(obj))}}
gc()





## Data (standardized for real variates)
dat <- list()
if(nrvars>0){ dat$Real=transform[varinfo[realVars,'transform']](t((t(data.matrix(alldata[, ..realVars])) - varinfo[realVars,'location'])/varinfo[realVars,'scale'])) }
if(ncvars>0){ dat$Category=t((t(data.matrix(alldata[, ..categoryVars])) - varinfo[categoryVars,'location']))}
## if(nbvars>0){ dat$Binary=data.matrix(alldata[, ..binaryVars])}
if(nbvars>0){ dat$Binary=t((t(data.matrix(alldata[, ..binaryVars])) - varinfo[binaryVars,'location']))}


####  CONSTANTS, PRIOR PARAMETERS, INITIAL VALUES
##
## In previous versions some statistics of the data were computed
## to decide on the hyperparameters.
## Now this is not done, because wrong in principle
## and because it can lead to silly hyperparameters
##
## Find max number of categories in data
nCat <- varinfo[categoryVars,'max'] - varinfo[categoryVars,'min'] + 1
if(ncvars > 0){
    ncategories <- max(nCat)
    calphapad <- array(2^(-40), dim=c(ncvars, ncategories), dimnames=list(categoryVars,NULL))
    for(avar in categoryVars){
        if(categoryprior=='Haldane'){
            calphapad[avar,1:ncat[avar]] <- 1/locvarMaxs[avar]
        }else{
            calphapad[avar,1:ncat[avar]] <- categoryprior
        }
    }
}
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################

#####################################################
## READ DATA AND SETUP CONSTANTS AND  HYPERPARAMETERS
#####################################################
variate <- lapply(variatetypes, function(x)rownames(varinfo)[varinfo['type']==x])
len <- lapply(variate,length)


data0 <- fread(paste0(origdir,datafile), sep=',')
if(!all(unlist(variate) %in% names(data0))){cat('\nERROR: variates missing from datafile')}
data0 <- data0[, unlist(variate), with=F]
## shuffle data
if(exists('shuffledata') && shuffledata){data0 <- data0[sample(1:nrow(data0))]}
if(!exists('ndata') || is.null(ndata) || is.na(ndata)){ndata <- nrow(data0)}
data0 <- data0[1:ndata]

##
## data
datapoints = c(
    ## real
    if(len$R > 0){list( Rdata = transfdir(data0[,variate$R,with=F], varinfo) )},
    ## logarithmic
    if(len$L > 0){list( Ldata = transfdir(data0[,variate$L,with=F], varinfo) )},
    ## integer
    if(len$I > 0){list( Idata = transfdir(data0[,variate$I,with=F], varinfo, Iinit=FALSE) )},
    ## Two-bounded
    if(len$T > 0){list(
                   Tdata = transfdir(data0[,variate$I,with=F], varinfo, Tout='data'),
                   Taux = transfdir(data0[,variate$I,with=F], varinfo, Tout='aux')
                  )},
    ## binary
    if(len$B > 0){list( Bdata = transfdir(data0[,variate$B,with=F], varinfo) )},
    ## categorical
    if(len$C > 0){list( Cdata = transfdir(data0[,variate$C,with=F], varinfo) )}
)

##
## calculation of some constants
Imaxn <- max(varinfo[variate$I, 'n']) - 1
Iintervals0 <- t(sapply(variate$I, function(v){
    createbounds(n=varinfo[v, 'n'], nmax=Imaxn)
}))
Iauxinit <- Idata = transfdir(data0[,variate$I,with=F], varinfo, Iinit=TRUE)
##
Tmaxn <- 2
Tintervals0 <- t(sapply(variate$T, function(v){ createbounds(n=varinfo[v, 'n']) }))
## constants
constants <- c(
    list(nclusters = nclusters),
    if(nalpha > 0){ list(nalpha = nalpha) },
    if(len$R > 0){ list(Rn = len$R) },
    if(len$L > 0){ list(Ln = len$L) },
    if(len$T > 0){ list(Tn = len$T, Tintervals0 = Tintervals0, Tmaxn = Tmaxn) },
    if(len$I > 0){ list(In = len$I, Iintervals0 = Iintervals0, Imaxn = Imaxn) },
    if(len$B > 0){ list(Bn = len$B) },
    if(len$C > 0){ list(Cn = len$C, Cmaxn = Cmaxn) },
    if(posterior){ list(ndata = ndata)}
)
##
## hyperparameters and some initial values
initsFunction <- function(){
    c(
        if(nalpha > 1){list( # distribution over concentration parameter
                         probalpha0 = rep(1/nalpha, nalpha),
                         walpha0 = matrix(c(0.5, 1, 2)/nclusters,
                                          nrow=nalpha, ncol=nclusters)
                     )}else{list(
                                walpha0 = rep(1/nclusters, nclusters)
                            )},
        if(len$R > 0){list( # real variate
                     Rmean0 = varinfo[variate$R, 'mean'],
                     Rvar0 = varinfo[variate$R, 'sd']^2,
                     Rshapeout0 = varinfo[variate$R, 'shapeout'],
                     Rshapein0 = varinfo[variate$R, 'shapein'],
                     Rvarscale0 = varinfo[variate$R, 'varscale']^2
                 )},
        if(len$L > 0){list( # logarithmic variate
                     Lmean0 = varinfo[variate$L, 'mean'],
                     Lvar0 = varinfo[variate$L, 'sd']^2,
                     Lshapeout0 = varinfo[variate$L, 'shapeout'],
                     Lshapein0 = varinfo[variate$L, 'shapein'],
                     Lvarscale0 = varinfo[variate$L, 'varscale']^2
                 )},
        if(len$T > 0){list( # real variate
                     Tmean0 = varinfo[variate$T, 'mean'],
                     Tvar0 = varinfo[variate$T, 'sd']^2,
                     Tshapeout0 = varinfo[variate$T, 'shapeout'],
                     Tshapein0 = varinfo[variate$T, 'shapein'],
                     Tvarscale0 = varinfo[variate$T, 'varscale']^2
                 )},
        if(len$I > 0){list( # real variate
                     Imean0 = varinfo[variate$I, 'mean'],
                     Ivar0 = varinfo[variate$I, 'sd']^2,
                     Ishapeout0 = varinfo[variate$I, 'shapeout'],
                     Ishapein0 = varinfo[variate$I, 'shapein'],
                     Ivarscale0 = varinfo[variate$I, 'varscale']^2,
                     Iaux = t((t(datapoints$Idata) + 0.5)/varinfo[variate$I, 'n'])
                 )},
        if(len$B > 0){list( # real variate
                     Bshapeout0 = varinfo[variate$B, 'shapeout'],
                     Bshapein0 = varinfo[variate$B, 'shapein']
                 )},
        if(len$C > 0){list( # real variate
                     Calpha0 = t(sapply(variate$B, function(v){
                         c( rep(varinfo[v, 'shapeout'], varinfo[v, 'max']),
                           rep(2^(-40), Cmaxn-varinfo[v, 'max']) )
                     }))
                 )},
        if((!casualinitvalues) & posterior){list(
                                                W = rep(1/nclusters, nclusters),
                                                K = rep(1, ndata) # all in one cluster at first
                                            )},
        if(casualinitvalues & posterior){list(
                                             W = rdirch(1, alpha=rep(1,nclusters)),
                                             K = sample(1:nclusters, ndata, replace=TRUE)
                                         )}
    )}


##
#### Mathematical representation of long-run frequency distributions
finitemix <- nimbleCode({
    if(nalpha > 1){# distribution over concentration parameter
        Alpha <- dcat(prob=probalpha0[1:nalpha])
        W[1:nclusters] ~ ddirch(alpha=walpha0[Alpha, 1:nclusters])
    }else{
        W[1:nclusters] ~ ddirch(alpha=walpha0[1:nclusters])
    }
    ##
    if(len$R > 0){# real variates
            for(v in 1:Rn){
                Rrate[v] ~ dinvgamma(shape=Rshapein0[v], scale=Rvarscale0[v])
            }
        }
    if(len$L > 0){# logarithmic variates
            for(v in 1:Ln){
                Lrate[v] ~ dinvgamma(shape=Lshapein0[v], scale=Lvarscale0[v])
            }
        }
    if(len$T > 0){# bounded continuous variates
            for(v in 1:Tn){
                Trate[v] ~ dinvgamma(shape=Tshapein0[v], scale=Tscale0[v])
            }
        }
    if(len$I > 0){# integer variates
            for(v in 1:In){
                Irate[v] ~ dinvgamma(shape=Ishapein0[v], scale=Iscale0[v])
            }
        }
    ##
    for(k in 1:nclusters){
        if(len$R > 0){# real variates
            for(v in 1:Rn){
                Rmean[v, k] ~ dnorm(mean=Rmean0[v], var=Rvar0[v])
                Rvar[v, k] ~ dinvgamma(shape=Rshapeout0[v], rate=Rrate[v])
            }
        }
        if(len$L > 0){# logarithmic variates
            for(v in 1:Ln){
                Lmean[v, k] ~ dnorm(mean=Lmean0[v], var=Lvar0[v])
                Lvar[v, k] ~ dinvgamma(shape=Lshapeout0[v], rate=Lrate[v])
            }
        }
        if(len$T > 0){# bounded continuous variates
            for(v in 1:Tn){
                Tmean[v, k] ~ dnorm(mean=Tmeanmean0[v], var=Tmeanvar0[v])
                Tvar[v, k] ~ dinvgamma(shape=Tshapeout0[v], rate=Trate[v])
            }
        }
        if(len$I > 0){# bounded continuous variates
            for(v in 1:In){
                Imean[v, k] ~ dnorm(mean=Imeanmean0[v], var=Imeanvar0[v])
                Ivar[v, k] ~ dinvgamma(shape=Ishapeout0[v], rate=Irate[v])
            }
        }
        if(len$B > 0){# binary variates
            for(v in 1:Bn){
                Bprob[v, k] ~ dbeta(shape1=Bshapein0[v], shape2=Bshapeout0[v])
            }
        }
        if(len$C > 0){# categorical variates
            for(v in 1:Cn){
                Cprob[v, k, 1:nmaxcategories] ~ ddirch(alpha=Calpha0[v, 1:Cnmaxcategories])
            }
        }
    }
    ##
    if(posterior){# cluster occupations
        for(d in 1:ndata){
            K[d] ~ dcat(prob=W[1:nclusters])
            ##
            if(len$R > 0){# real variates
                for(v in 1:Rn){
                    Rdata[d, v] ~ dnorm(mean=Rmean[v, K[d]], var=Rvar[v, K[d]])
                }
            }
            if(len$L > 0){# logarithmic variates
                for(v in 1:Ln){
                    Ldata[d, v] ~ dnorm(mean=Lmean[v, K[d]], var=Lvar[v, K[d]])
                }
            }
            if(len$I > 0){# integer variates
                for(v in 1:In){
                    Idata[d, v] ~ dinterval(t=Iaux[v, K[d]], c=Iintervals0[v, 1:Imaxn])
                    Iaux[d, v] ~ dnorm(mean=Imean[v, K[d]], var=Ivar[v, K[d]])
                }
            }
            if(len$T > 0){# bounded continuous variates
                for(v in 1:Tn){
                    Taux[d, v] ~ dinterval(t=Tdata[v, K[d]], c=Tintervals0[v, 1:Tmaxn])
                    Tdata[d, v] ~ dnorm(mean=Tmean[v, K[d]], var=Tvar[v, K[d]])
                }
            }
            if(len$B > 0){# binary variates
                for(v in 1:Bn){
                    Bdata[d, v] ~ dbern(prob=Bprob[v, K[d]])
                }
            }
            if(len$C > 0){# categorical variates
                for(v in 1:Cn){
                    Cdata[d, v] ~ dcat(prob=Cprob[v, K[d], 1:Cmaxn])
                }
            }
        }
    }
})

##
timecount <- Sys.time()
##
if(posterior){
    finitemixnimble <- nimbleModel(code=finitemix, name='finitemixnimble1', constants=constants,
                         inits=initsFunction(), data=dat,
                         dimensions=c(
                             list(q=nclusters),
                             (if(nrvars > 0){list(meanR=c(nrvars,nclusters), tauR=c(nrvars,nclusters))}),
                             (if(nivars > 0){list(probI=c(nivars,nclusters))}),
                             (if(ncvars > 0){list(probC=c(ncvars,nclusters,ncategories))}),
                             (if(nbvars > 0){list(probB=c(nbvars,nclusters))}),
                             list(C=ndata),
                             if(compoundgamma){list(varRrate=nrvars)})
                         )
}else{
    finitemixnimble <- nimbleModel(code=finitemix, name='finitemixnimble1', constants=constants,
                         inits=initsFunction(), data=list(),
                         dimensions=c(
                             list(q=nclusters),
                             (if(nrvars > 0){list(meanR=c(nrvars,nclusters), tauR=c(nrvars,nclusters))}),
                             (if(nivars > 0){list(probI=c(nivars,nclusters))}),
                             (if(ncvars > 0){list(probC=c(ncvars,nclusters,ncategories))}),
                             (if(nbvars > 0){list(probB=c(nbvars,nclusters))}),
                             if(compoundgamma){list(varRrate=nrvars)})
                         )
}
Cfinitemixnimble <- compileNimble(finitemixnimble, showCompilerOutput=FALSE)
gc()


##
if(posterior){# Samplers for posterior sampling
    confnimble <- configureMCMC(Cfinitemixnimble, nodes=NULL,
                               monitors=c('q',
                                          if(nrvars > 0){c('meanR', 'varR')},
                                          if(nivars > 0){c('probI', 'sizeI')},
                                          if(ncvars > 0){c('probC')},
                                          if(nbvars > 0){c('probB')}
                                          ),
                               monitors2=c('C',
                                           if(compoundgamma & nrvars > 0){c('varRrate')}
                                           )
                               )
    ##
    for(adatum in 1:ndata){
        confnimble$addSampler(target=paste0('C[', adatum, ']'), type='categorical')
    }
    if(compoundgamma & nrvars>0){
        for(avar in 1:nrvars){
            confnimble$addSampler(target=paste0('varRrate[', avar, ']'), type='conjugate')
        }
    }
    for(acluster in 1:nclusters){
        if(nrvars>0){
            for(avar in 1:nrvars){
                confnimble$addSampler(target=paste0('varR[', avar, ', ', acluster, ']'), type='conjugate')
                confnimble$addSampler(target=paste0('meanR[', avar, ', ', acluster, ']'), type='conjugate')
            }
        }
        if(nivars>0){
            for(avar in 1:nivars){
                confnimble$addSampler(target=paste0('probI[', avar, ', ', acluster, ']'), type='conjugate')
                confnimble$addSampler(target=paste0('sizeI[', avar, ', ', acluster, ']'), type='categorical')
            }
        }
        if(ncvars>0){
            for(avar in 1:ncvars){
                confnimble$addSampler(target=paste0('probC[', avar, ', ', acluster, ', 1:', ncategories, ']'), type='conjugate')
            }
        }
        if(nbvars>0){
            for(avar in 1:nbvars){
                confnimble$addSampler(target=paste0('probB[', avar, ', ', acluster, ']'), type='conjugate')
            }
        }
    }
    confnimble$addSampler(target=paste0('q[1:', nclusters, ']'), type='conjugate')
##
}else{# sampler for prior sampling
    confnimble <- configureMCMC(Cfinitemixnimble, 
                               monitors=c('q',
                                          if(nrvars>0){c('meanR', 'varR')},
                                          if(nivars>0){c('probI', 'sizeI')},
                                          if(ncvars>0){c('probC')},
                                          if(nbvars>0){c('probB')}
                                          ),
                               monitors2=c(if(compoundgamma & nrvars > 0){c('varRrate')}
                                           )
                               )
}
##
print(confnimble)

mcsampler <- buildMCMC(confnimble)
Cmcsampler <- compileNimble(mcsampler, resetFunctions = TRUE)
gc()

cat('\nSetup time: ')
print(Sys.time() - timecount)

##################################################
## Monte Carlo sampler and plots of MC diagnostics
##################################################
traces <- mcsamples <- NULL
burnin <- 0
continue <- TRUE
stage <- -1
niter <- niter0
thin <- 1
totaliter <- 0
time0 <- Sys.time()
while(continue){
    calctime <- Sys.time()
    stage <- stage+1

    cat(paste0('\n\n==== STAGE ', stage, ' ===='))
    cat(paste0('\nIterations: ',niter,', thinning: ',thin,'\n'))
    gc()
    if(stage==0){# burn-in stage
        set.seed(mcmcseed+stage+100)
        Cfinitemixnimble$setInits(initsFunction())
        newmcsamples <- Cmcsampler$run(niter=niter+1, thin=thin, thin2=niter, nburnin=1, time=T)
    }else if(is.character(resume)){# continuing previous # must be fixed
        initsc <- readRDS(paste0(dirname,resume))
        inits0 <- initsFunction()
        for(aname in names(inits0)){inits0[[aname]] <- initsc[[aname]]}
        thin <- initsc[['thin']]
        set.seed(mcmcseed+stage+100)
        Cfinitemixnimble$setInits(initsc)
        newmcsamples <- Cmcsampler$run(niter=niter*thin, thin=thin, thin2=niter*thin, nburnin=0)
    }else{# subsequent sampling stages
        cat('\nForecasted computation time: ')
        print(comptime*thin*niter)
        newmcsamples <- Cmcsampler$run(niter=niter*thin, thin=thin, thin2=niter*thin, reset=FALSE, resetMV=TRUE)
    }
    ##
    totaliter <- totaliter + niter*thin
    newmcsamples <- as.matrix(Cmcsampler$mvSamples)
    cat('\nTime MCMC: ')
    print(Sys.time() - calctime)
    ##
    ## ## Check sample-time partition
    ## times <- Cmcsampler$getTimes()
    ## names(times) <- sapply(confnimble$getSamplers(),function(x)x$target)
    ## ##
    ## cbind(sort(times[c('C[1]','q[1:64]','meanR[1, 1]', 'varR[1, 1]', 'probB[1, 1]', 'probC[1, 1, 1:21]', 'varRrate[1]')]))
    ## ##
    ## test <- sapply(c('C','q','meanR', 'varR', 'probB', 'probC', 'varRrate'),
    ##        function(x){
    ##            sum(times[grep(paste0('^',x),names(times))])
    ##        })
    ## names(test) <- c('C','q','meanR', 'varR', 'probB', 'probC', 'varRrate')
    ## cbind(sort(test))
    ##
    if(any(is.na(newmcsamples))){cat('\nWARNING: SOME NA OUTPUTS')}
    if(any(!is.finite(newmcsamples))){cat('\nWARNING: SOME INFINITE OUTPUTS')}
    ##
    ## save final state of MCMC chain
    finalstate <- as.matrix(Cmcsampler$mvSamples2)
    finalstate <- c(newmcsamples[nrow(newmcsamples),], finalstate[nrow(finalstate),])
    ##
    ## Check how many "clusters" were occupied. Warns if too many
    occupations <- finalstate[grepl('^C\\[', names(finalstate))]
    usedclusters <- length(unique(occupations))
    if(usedclusters > nclusters-5){cat('\nWARNING: TOO MANY CLUSTERS OCCUPIED')}
    cat(paste0('\nOCCUPIED CLUSTERS: ', usedclusters, ' OF ', nclusters))
##    saveRDS(finalstate2list(finalstate, realVars=realVars, integerVars=integerVars, categoryVars=categoryVars, binaryVars=binaryVars, compoundgamma=compoundgamma), file=paste0(dirname,'_finalstate-R',basename,'--',mcmcseed,'-',stage,'.rds'))
    ##
    ## SAVE THE PARAMETERS
    ##    parmList <- mcsamples2parmlist(mcsamples, realVars, integerVars, categoryVars, binaryVars)
    ##  saveRDS(parmList,file=paste0(dirname,'_frequencies-R',baseversion,'-V',length(varNames),'-D',ndata,'-K',nclusters,'-I',nrow(parmList$q),'--',mcmcseed,'-',stage,'.rds'))
    ##
    ## Diagnostics
    ## Log-likelihood
    if(posterior){
        ll <- llSamplesmc(dat, newmcsamples)
        flagll <- FALSE
        if(!posterior && !any(is.finite(ll))){
            flagll <- TRUE
            ll <- rep(0, length(ll))}
        condprobsd <- c(logsumsamplesFmc(Y=do.call(cbind,dat)[, mainvar, drop=F],
                                     X=do.call(cbind,dat)[, setdiff(varNames, mainvar),
                                                          drop=F],
                                     mcsamples=newmcsamples,
                                     varinfo=varinfo, inorder=T))
        condprobsi <- c(logsumsamplesFmc(Y=do.call(cbind,dat)[, setdiff(varNames, mainvar),
                                                          drop=F],
                                     X=do.call(cbind,dat)[, mainvar, drop=F],
                                     mcsamples=newmcsamples,
                                     varinfo=varinfo, inorder=T))
        ##
        traces <- rbind(traces,
                        10/log(10)/ndata *
                        cbind(loglikelihood=ll,
                              'mean of direct logprobabilities'=condprobsd,
                              'mean of inverse logprobabilities'=condprobsi)
                        )
        badcols <- foreach(i=1:ncol(traces), .combine=c)%do%{if(all(is.na(traces[,i]))){i}else{NULL}}
        if(!is.null(badcols)){traces <- traces[,-badcols]}
        saveRDS(traces,file=paste0(dirname,'_mctraces-R',basename,'--',mcmcseed,'-',stage,'.rds'))
        ##
        if(nrow(traces)>=1000){
            funMCSE <- function(x){LaplacesDemon::MCSE(x, method='batch.means')$se}
        }else{
            funMCSE <- function(x){LaplacesDemon::MCSE(x)}
        }
        diagnESS <- LaplacesDemon::ESS(traces * (abs(traces) < Inf))
        diagnIAT <- apply(traces, 2, function(x){LaplacesDemon::IAT(x[is.finite(x)])})
        diagnBMK <- LaplacesDemon::BMK.Diagnostic(traces, batches=4)[,1]
        diagnMCSE <- 100*apply(traces, 2, function(x){funMCSE(x)/sd(x)})
        diagnStat <- apply(traces, 2, function(x){LaplacesDemon::is.stationary(as.matrix(x,ncol=1))})
        diagnBurn <- apply(traces, 2, function(x){LaplacesDemon::burnin(matrix(x[1:(10*trunc(length(x)/10))], ncol=1))})
        diagnBurn2 <- proposeburnin(traces, batches=10)
        diagnThin <- proposethinning(traces)
        ##
        cat(paste0('\nESSs (',requiredESS,'): ',paste0(round(diagnESS), collapse=', ')))
        cat(paste0('\nIATs: ',paste0(round(diagnIAT), collapse=', ')))
        cat(paste0('\nBMKs: ',paste0(round(diagnBMK,3), collapse=', ')))
        cat(paste0('\nMCSEs: ',paste0(round(diagnMCSE,2), collapse=', ')))
        cat(paste0('\nStationary: ',paste0(diagnStat, collapse=', ')))
        cat(paste0('\nBurn-in I: ',paste0(diagnBurn, collapse=', ')))
        cat(paste0('\nBurn-in II: ',diagnBurn2))
        cat(paste0('\nProposed thinning: ',paste0(diagnThin, collapse=', ')))
        ##
        #########################################
        #### CHECK IF WE NEED TO SAMPLE MORE ####
        #########################################
        if(stage==0){
            thin <- round(max(diagnThin)*1.5)
            burnin <- niter0-1
            continue <- TRUE
        }else{
            if(min(diagnESS) >= ceiling(requiredESS) &
               #max(diagnMCSE) < 6.27 &
               sum(diagnStat) == 3 &
               diagnBurn2 == 0
               ){
                continue <- FALSE
                burnin <- 0
            }else{
                continue <- TRUE
                ## if(max(diagnThin) > 1){
                ##     thin <- thin*max(diagnThin)
                ##     burnin <- nsamples-1
                ## }else{
                burnin <- min(max(diagnBurn2, minstepincrease), nsamples-1)
            }
        }
        niter <- nsamples - nrow(traces) + burnin
        #########################################
        #### END CHECK                       ####
        #########################################
        mcsamples <- rbind(mcsamples, newmcsamples)
        rm(newmcsamples)
        if(savetempsamples | !continue){
            saveRDS(mcsamples,file=paste0(dirname,'_mcsamples-R',basename,'--',mcmcseed,'-',stage,'.rds'))
        }

        ##
        tracegroups <- list(loglikelihood=1,
                            'main given rest'=2,
                            'rest given main'=3
                            )
        grouplegends <- foreach(agroup=1:length(tracegroups))%do%{
            c( paste0('-- STATS ', names(tracegroups)[agroup], ' --'),
              paste0('min ESS = ', signif(min(diagnESS[tracegroups[[agroup]]]),6)),
              paste0('max IAT = ', signif(max(diagnIAT[tracegroups[[agroup]]]),6)),
              paste0('max BMK = ', signif(max(diagnBMK[tracegroups[[agroup]]]),6)),
              paste0('max MCSE = ', signif(max(diagnMCSE[tracegroups[[agroup]]]),6)),
              paste0('stationary: ', sum(diagnStat[tracegroups[[agroup]]]),'/',length(diagnStat[tracegroups[[agroup]]])),
              ## paste0('burn: ', signif(max(diagnBurn[tracegroups[[agroup]]]),6))
              paste0('burn: ', signif(diagnBurn2,6))
              )
        }
        colpalette <- c(7,2,1)
        names(colpalette) <- colnames(traces)
    ##
    ## Plot various info and traces
        cat('\nPlotting MCMC traces')
        graphics.off()
        pdff(paste0(dirname,'mcmcplottraces-R',basename,'--',mcmcseed,'-',stage),'a4')
    ## Summary stats
        matplot(1:2, type='l', col='white', main=paste0('Stats stage ',stage), axes=FALSE, ann=FALSE)
        legendpositions <- c('topleft','topright','bottomleft','bottomright')
        for(alegend in 1:length(grouplegends)){
            legend(x=legendpositions[alegend], bty='n', cex=1.5,
                   legend=grouplegends[[alegend]] )
        }
        legend(x='center', bty='n', cex=1,
               legend=c(
                   paste0('STAGE ',stage),
                   paste0('Occupied clusters: ', usedclusters, ' of ', nclusters),
                   paste0('LL:  ( ', signif(mean(traces[,1]),3), ' +- ', signif(sd(traces[,1]),3),' ) dHart'),
                   'WARNINGS:',
                   if(any(is.na(mcsamples))){'some NA MC outputs'},
                   if(any(!is.finite(mcsamples))){'some infinite MC outputs'},
                   if(usedclusters > nclusters-5){'too many clusters occupied'},
                   if(flagll){'infinite values in likelihood'}
               ))
        ##
        ## Traces of likelihood and cond. probabilities
        par(mfrow=c(1,1))
        for(avar in colnames(traces)){
            tplot(y=traces[,avar], type='l', lty=1, col=colpalette[avar],
                  main=paste0(avar,
                              '\nESS = ', signif(diagnESS[avar], 3),
                              ' | IAT = ', signif(diagnIAT[avar], 3),
                              ' | BMK = ', signif(diagnBMK[avar], 3),
                              ' | MCSE = ', signif(diagnMCSE[avar], 3),
                              ' | stat: ', diagnStat[avar],
                              ' | burn I: ', diagnBurn[avar],
                              ' | burn II: ', diagnBurn2
                              ),
                  ylab=paste0(avar,'/dHart'), xlab='sample', family=family
                  )
        }
dev.off()
    }
    ##
    ## Samples of marginal frequency distributions
    if(plottempdistributions | !continue){
        if(plotmeans){nfsamples <- totsamples
        }else if(!posterior){nfsamples <- 256
        }else{nfsamples <- 63}
        nfsamples <- min(nfsamples, nrow(mcsamples))
        subsample <- round(seq(1,nfsamples, length.out=63))
        ##
        cat('\nPlotting samples of frequency distributions')
        graphics.off()
        pdff(paste0(dirname,'mcmcdistributions-R',basename,'--',mcmcseed,'-',stage),'a4')
        for(avar in varNames){#cat(avar)
            if(avar %in% realVars){
                rg <- signif((varinfo[avar,c('thmin','thmax')] + 
                              7*varinfo[avar,c('datamin','datamax')])/8, 2)
                if(!is.finite(rg[1])){rg[1] <- diff(varinfo[avar,c('scale','datamin')])}
                if(!is.finite(rg[2])){rg[2] <- sum(varinfo[avar,c('scale','datamax')])}
                Xgrid <- cbind(seq(rg[1], rg[2], length.out=256))
            }else{
                rg <- (varinfo[avar,c('thmin','thmax')] + 
                       7*varinfo[avar,c('datamin','datamax')])/8
                rg <- c(floor(rg[1]), ceiling(rg[2]))
                Xgrid <- cbind(rg[1]:rg[2])
            }
            colnames(Xgrid) <- avar
            plotsamples <- samplesFmc(Y=Xgrid, mcsamples=mcsamples, fromsamples=round(seq(1,nrow(mcsamples),length.out=nfsamples)), inorder=FALSE, varinfo=varinfo)
            fiven <- varinfo[avar,c('datamin','dataQ1','datamedian','dataQ2','datamax')]
            ##
            if(posterior){
                par(mfrow=c(1,1))
                ymax <- tquant(apply(plotsamples[,subsample],2,function(x){tquant(x,99/100)}),99/100, na.rm=T)
                tplot(x=Xgrid, y=plotsamples[,subsample], type='l', col=paste0(palette()[7], '44'), lty=1, lwd=2, xlab=paste0(avar,' (',variateinfo[variate==avar,type],')'), ylab=paste0('frequency',if(avar %in% realVars){' density'}), ylim=c(0, ymax), family=family)
                if(plotmeans){
                    tplot(x=Xgrid, y=rowMeans(plotsamples), type='l', col=paste0(palette()[1], '88'), lty=1, lwd=3, add=T)
                }
                abline(v=fiven,col=paste0(palette()[c(2,4,5,4,2)], '44'),lwd=4)
                ##
                if((showdata=='histogram')|(showdata==TRUE & !(avar %in% realVars))){
                    datum <- alldata[[avar]]
                    histo <- thist(datum, n=(if(avar %in% realVars){min(max(10,sqrt(ndata)),100)}else{'i'}))#-exp(mean(log(c(round(sqrt(length(datum))), length(Xgrid))))))
                    histomax <- (if(avar %in% realVars){max(rowMeans(plotsamples))/max(histo$density)}else{1L})
                    tplot(x=histo$breaks, y=histo$density*histomax, col=grey, alpha=0.75, border=NA, lty=1, lwd=1, family=family, ylim=c(0,NA), add=TRUE)
                }else if((showdata=='scatter')|(showdata==TRUE & (avar %in% realVars))){
                    datum <- alldata[[avar]]
                    scatteraxis(side=1, n=NA, alpha='88', ext=8, x=datum+rnorm(length(datum),mean=0,sd=prod(varinfo[avar,c('precision','scale')])/(if(avar %in% binaryVars){16}else{16})),col=yellow)
                }
            }else{
                par(mfrow=c(8,8),mar = c(0,0,0,0))
                tplot(x=list(Xgrid,Xgrid), y=list(rowMeans(plotsamples),rep(0,length(Xgrid))), type='l', col=c(paste0(palette()[3], 'FF'), '#bbbbbb80'), lty=1, lwd=c(2,1), xlab=NA, ylab=NA, ylim=c(0, NA), family=family,
                      xticks=NA, yticks=NA,
                      mar=c(1,1,1,1))
                abline(v=fiven,col=paste0(palette()[c(2,4,5,4,2)], '44'))
                text(sum(range(Xgrid))/2, par('usr')[4]*0.9, avar)
                ##
                for(aplot in 1:63){
                    tplot(x=list(Xgrid,Xgrid), y=list(plotsamples[,subsample[aplot]], rep(0,length(Xgrid))), type='l', col=c(paste0(palette()[1], 'FF'), '#bbbbbb80'), lty=1, lwd=c(1,1), xlab=NA, ylab=NA, ylim=c(0, NA), family=family,
                          xticks=NA, yticks=NA,
                          mar=c(1,1,1,1))
                    abline(v=fiven,col=paste0(palette()[c(2,4,5,4,2)], '44'))
                    if(aplot==1){ text(Xgrid[1], par('usr')[4]*0.9, variateinfo[variate==avar,type], pos=4)}
                    if(aplot==2){ text(Xgrid[1], par('usr')[4]*0.9, paste(signif(c(rg,diff(rg)),2),collapse=' -- '), pos=4)}
                }
            }
        }
    dev.off()
    }
    ##
    cat('\nTime MCMC+diagnostics: ')
    print(Sys.time() - calctime)
    comptime <- (Sys.time() - time0)/totaliter
    ##
    ## mcsamples <- mcsamples[(burnin+1):nrow(mcsamples),,drop=FALSE]
    ## traces <- traces[(burnin+1):nrow(traces),,drop=FALSE]
    ##
    ## if(savetempsamples | !continue){
    ##         saveRDS(mcsamples,file=paste0(dirname,'_mcsamples-R',basename,'--',mcmcseed,'-',stage,'.rds'))
    ##     }

    if(!continue){
        ## Change name of rds files with final results
        file.rename(from=paste0(dirname,'_mcsamples-R',basename,'--',mcmcseed,'-',stage,'.rds'), to=paste0(dirname,'_mcsamples-R',basename,'--',mcmcseed,'-','F','.rds') )
        file.rename(from=paste0(dirname,'_mctraces-R',basename,'--',mcmcseed,'-',stage,'.rds'), to=paste0(dirname,'_mctraces-R',basename,'--',mcmcseed,'-','F','.rds') )
        ## Change name of pdf files with final plots
        file.rename(from=paste0(dirname,'mcmcplottraces-R',basename,'--',mcmcseed,'-',stage,'.pdf'), to=paste0(dirname,'mcmcplottraces-R',basename,'--',mcmcseed,'-','F','.pdf') )
        file.rename(from=paste0(dirname,'mcmcdistributions-R',basename,'--',mcmcseed,'-',stage,'.pdf'), to=paste0(dirname,'mcmcdistributions-R',basename,'--',mcmcseed,'-','F','.pdf'))
        ##
        cat('\n==== RESULTS SEEM STATIONARY. END ====\n\n')
    }else{
        mcsamples <- mcsamples[(burnin+1):nrow(traces),,drop=F]
        traces <- traces[(burnin+1):nrow(traces),,drop=F]
    }
    ##
    
}

############################################################
## End MCMC
############################################################
plan(sequential)

stop('NONE. End of script')

alldata <- fread('data_ep.csv')
metadata <- fread('metadata.csv')
allnames <- names(alldata)
metadatanew <- metadata
metadatanew$precision <- as.integer(metadatanew$precision)
##
boundarycat <- 0
boundaryrea <- 0
for(i in allnames){
    if(metadata[variate==i,'type']=='integer'){
        print(i)
        if(is.na(metadata[variate==i,'max'])){
            print(paste0('(max)'))
            metadatanew[variate==i,'max'] <- max(alldata[[i]],na.rm=T)
        }
        if(is.na(metadata[variate==i,'min'])){
            print(paste0(': min'))
            metadatanew[variate==i,'min'] <- min(alldata[[i]],na.rm=T)
        }
        rg <- metadatanew[variate==i,'max']-metadatanew[variate==i,'min']+1
        if(rg < boundarycat){
            metadatanew[variate==i,'type'] <- 'categorical'
            print(paste0(': to cat, ',rg))
        }else if(rg>=boundaryrea){
            metadatanew[variate==i,'type'] <- 'real'
            print(paste0(': to real, ',rg))
            metadatanew[variate==i,'precision'] <- 1L
        }else{print(': no changes')}
    }
}

fwrite(x=metadatanew, file=paste0('metadata_noint',boundarycat,'_',boundaryrea,'.csv'))

nosw <- metadatanew$variate[grepl('_SW$',metadatanew$variate)]
metadatanew2 <- metadatanew[!(variate %in% nosw),]

fwrite(x=metadatanew2, file=paste0('metadata_noint',boundarycat,'_',boundaryrea,'_noSW.csv'))


for(i in allnames){if(metadata[variate==i,type]=='integer' & !(min(diff(sort(unique(alldata[[i]]))))==1)){cat(c(i,min(diff(sort(unique(alldata[[i]]))))))}}



traces <- traces1
t(sapply(1:17,function(thin){c(
                                 thin,
                                 LaplacesDemon::IAT(traces[1:1024,1]),
                                 LaplacesDemon::IAT(traces[1025:2048,1]),
                                 rev(thisseq <- seq(from=2048+1,by=thin,length.out=2048))[1],
                                 thisiat <- LaplacesDemon::IAT(traces[thisseq,1]),
                                 thisess <- LaplacesDemon::ESS(traces[thisseq,1]),
                                 thisiat <- LaplacesDemon::IAT(traces[thisseq,1]),
                                 thisess <- LaplacesDemon::ESS(traces[thisseq,1]),
                                 floor(thisiat*thin),round(thisiat)*thin,
                                 NA,
                                 floor(length(thisseq)/thisess*thin),round(length(thisseq)/thisess)*thin
                             )}))


traces <- traces2
thin <-  round(1024/LaplacesDemon::ESS(traces[1:1024,1]))
thisseq <- seq(from=1024+1,by=thin,length.out=1024)
c(thin, range(thisseq), max(thisseq)/1024)
LaplacesDemon::IAT(traces[thisseq,1])
LaplacesDemon::ESS(traces[thisseq,1])
tplot(y=traces[thisseq,1])

traces <- traces2
leng <- 2048
thin <-  round(LaplacesDemon::IAT(traces[1:1024,1]))
thisseq <- seq(from=1024+1,by=thin,length.out=leng)
thisseq <- round(max(thisseq)/4)+thisseq-1024
x <- traces[thisseq,1]
c(thin, range(thisseq), max(thisseq)/1024)
LaplacesDemon::BMK.Diagnostic(cbind(x), batches=2)[,1]
LaplacesDemon::BMK.Diagnostic(cbind(x), batches=4)[,1]
LaplacesDemon::IAT(x)
LaplacesDemon::ESS(x)
LaplacesDemon::MCSE(x, method='batch.means')$se*100/sd(x) #6.27
LaplacesDemon::MCSE(x, method='batch.means')$se*100/sd(x) < 6.27
LaplacesDemon::is.stationary(cbind(x))
tplot(y=traces[thisseq,1])

##traces <- traces1
thin <-  round(LaplacesDemon::IAT(traces[1:1024,1]))
t(sapply(1:16,function(i){
    thisseq <- 1024+seq(from=128*(i-1)+1,by=thin,length.out=128)
    x <- traces[thisseq,1]
    c(thin, range(thisseq), max(thisseq)/128,NA,
    LaplacesDemon::BMK.Diagnostic(cbind(x), batches=2)[,1],
    LaplacesDemon::BMK.Diagnostic(cbind(x), batches=4)[,1],
    LaplacesDemon::IAT(x),
    LaplacesDemon::ESS(x),
    LaplacesDemon::MCSE(x, method='batch.means')$se*100/sd(x), #6.27
    LaplacesDemon::MCSE(x, method='batch.means')$se*100/sd(x) < 6.27,
    LaplacesDemon::is.stationary(cbind(x))
)}))

thisseq <- round(max(thisseq)/2)+thisseq-1024
tplot(y=traces[thisseq,1])





traces <- traces2
t(sapply(1:36,function(i){
    thisseq <- seq(from=1024*(i-1)+1,by=1,length.out=1024)
    x <- traces[thisseq,1]
    c(i,max(thisseq),
      LaplacesDemon::BMK.Diagnostic(cbind(x), batches=2)[,1],
      LaplacesDemon::BMK.Diagnostic(cbind(x), batches=4)[,1],
      LaplacesDemon::IAT(x),
      LaplacesDemon::ESS(x),
      LaplacesDemon::MCSE(x, method='batch.means')$se*100/sd(x),
      LaplacesDemon::is.stationary(cbind(x))
)}))

## q         0.029227
## varRrate  0.707669
## probB     1.131390
## C         6.549711
## meanR    11.425076
## varR     15.414792
## probC    94.969083




xgrid <- seq(-5,-3,length.out=512)
fn <- function(x,dx){pnorm(x)/(2*dx)-1}
testy <- fn(xgrid,1e-4)
tplot(x=xgrid,y=((testy)))

qnorm(2*(1e-4))
## qnorm of twice precision of variate = boundary


mm <- c(0,2,4)
ss <- c(1,2,0.2)
qq <- c(7,4,0.4)
##
f <- function(x){sapply(x,function(y)sum(qq*dnorm(y,mean=mm,sd=ss))/sum(qq))}
##
xgrid <- seq(-3,6,length.out=256)
graphics.off()
pdff('testcog1')
tplot(xgrid,f(xgrid),ylim=c(0,NA),xlabels=NA,xlab='cog score #4',ylab='probability',ylabels=NA,lwd=5,ly=2,yticks=NA)
plotquantiles(xgrid,rbind(
                        f(xgrid)*
                        (1+3*(dnorm(xgrid-min(xgrid))+0.2*dnorm(xgrid-2.5,sd=1)+dnorm(max(xgrid)-xgrid))),
                        f(xgrid)*
                        (1-3*(dnorm(xgrid-min(xgrid))+0.2*dnorm(xgrid-2.5,sd=1)+dnorm(max(xgrid)-xgrid)))
                    ))
dev.off()
me <- sum(xgrid*f(xgrid))/sum(f(xgrid))
sqrt(sum((xgrid-me)^2*f(xgrid))/sum(f(xgrid)))
##
##
mm <- c(0,2,3)
ss <- c(1,2,0.2)*0.5
qq <- c(7,2,0.1)
##
f <- function(x){sapply(x,function(y)sum(qq*dnorm(y,mean=mm,sd=ss))/sum(qq))}
##
xgrid <- seq(-3,6,length.out=256)
graphics.off()
pdff('testcog2')
tplot(xgrid,f(xgrid),ylim=c(0,NA),xlabels=NA,xlab='cog score #4',ylab='probability',ylabels=NA,lwd=5,ly=2,yticks=NA)
plotquantiles(xgrid,rbind(
                        f(xgrid)*
                        (1+7*(dnorm(xgrid-min(xgrid))+0.1*dnorm(xgrid-2.5,sd=1)+dnorm(max(xgrid)-xgrid))),
                        f(xgrid)*
                        (1-7*(dnorm(xgrid-min(xgrid))+0.1*dnorm(xgrid-2.5,sd=1)+dnorm(max(xgrid)-xgrid)))
                    ))
dev.off()
me <- sum(xgrid*f(xgrid))/sum(f(xgrid))
sqrt(sum((xgrid-me)^2*f(xgrid))/sum(f(xgrid)))




tplot(xgrid,dbeta((xgrid-min(xgrid))/diff(range(xgrid)), 0.2,0.2))





bounds <- c(-1, 1)
tra1 <- function(x){
    x[x<=bounds[1]] <- bounds[1]
    x[x>=bounds[2]] <- bounds[2]
    x
}
## slower:
tra2 <- function(x){
    lx <- length(x)
    ind <- rinterval(lx, x, bounds)
    bounds[1]*(ind==0)+bounds[2]*(ind==2)+x*(ind==1)
}


bounds <- ((-2):(2-1))+0.5
##
itra2 <- function(x){
    rinterval(length(x), x, bounds)
}
## slower:
itra0 <- function(x){
    rowSums(sapply(bounds, function(i){x>=i}))
}
## slower:
itra1 <- function(x){
    boundse <- c(-Inf, bounds, Inf)
    rowSums(sapply(2:length(boundse), function(i){
        (x>boundse[i-1] & x<=boundse[i])*i
    }))-2
}
