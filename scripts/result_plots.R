## Author: PGL  Porta Mana
## Created: 2021-11-25T14:52:14+0100
## Last-Updated: 2021-12-02T14:40:38+0100
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
frequenciesfile <- 'newIposterior_-V13-D539-K50-I1024/_frequencies-RnewIposterior_1-V13-D539-K50-I1024.rds'
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

distsAF <- foreach(acov=setdiff(covNames, maincov))%do%{
    dists <- rbind(samplesF(X=grids[[acov]], Y=rbind(xcond[,1]), parmList=parmList, inorder=T),
                   samplesF(X=grids[[acov]], Y=rbind(xcond[,2]), parmList=parmList, inorder=T)
                   )
    dim(dists) <- c(length(grids[[acov]]), 2, ncol(dists))
    dimnames(dists) <- list(NULL, subgroupnames, NULL)
    aperm(dists, c(3,1,2))
}
names(distsAF) <- otherCovs

quantFA <- foreach(acov=otherCovs)%do%{
    apply(distsFA[[acov]],c(2,3),function(x)c(mean(x,na.rm=T), quant(x=x, probs=c(1,7)/8, na.rm=T)))
}
names(quantFA) <- otherCovs


pdff('plots_features_given_AD')
for(acov in otherCovs){
    agrid <- grids[[acov]]
    ymax <- quant(apply(quantFA[[acov]],2,function(x){quant(x,99/100)}),99/100)
ylim <- c(0,ymax)#max(quantFA[[acov]]))
for(i in 1:2){
    tplot(x=agrid, y=quantFA[[acov]][1,,i],
          col=i, lty=i, lwd=3, alpha=0.25, ylim=ylim,
          xlab=acov, ylab='population frequency', add=(i==2))
    polygon(x=c(agrid,rev(agrid)), y=c(quantFA[[acov]][2,,i], rev(quantFA[[acov]][3,,i])), col=paste0(palette()[i],'40'), border=NA)
}
legend(x=agrid[1], y=ylim[2]*1.2, legend=c(subgroupnames, '75% uncertainty'),
       col=palette()[c(1,2,7)], lty=c(1,2,1), lwd=c(2,2,15), cex=1.5, bty='n'
                       )
}
dev.off()

tplot(x=grids[[acov]], y=t(quantFA[[acov]][,,2]),col=2, add=T)


pdff('plots_subgroups')
for(acov in setdiff(covNames, maincov)){
    if(acov %in% realCovs){



        
        rg <- range(datum)+c(-1,1)*IQR(datum)
            Xgrid <- seq(rg[1], rg[2], length.out=256)
            tpar <- unlist(variateinfo[variate==acov,c('transfM','transfW')])
            if(!any(is.na(tpar))){
                Ogrid <- pretty(exp(tpar['transfW']*Xgrid + tpar['transfM']),n=10)
            }
        }else{
            rg <- range(datum)
            rg <- round(c((covMins[acov]+7*rg[1])/8, (covMaxs[acov]+7*rg[2])/8))
            Xgrid <- rg[1]:rg[2]
            tpar <- NA
        }
        Xgrid <- cbind(Xgrid)
        colnames(Xgrid) <- acov
        plotsamples <- samplesF(Y=Xgrid, parmList=parmList, nfsamples=min(64,nrow(mcsamples)), inorder=FALSE)

    
}

dev.off()





pdff('testhistograms')
for(acov in colnames(alldata)){
    print(acov)
if(is.numeric(alldata[[acov]])){
    histo <- thist(alldata[[acov]],round((nrow(alldata))^(2/3)))
 tplot(x=histo$breaks, y=histo$density, main=acov,col=2)
}}
dev.off()



