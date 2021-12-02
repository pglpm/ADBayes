## Author: PGL  Porta Mana
## Created: 2021-11-25T14:52:14+0100
## Last-Updated: 2021-12-02T06:40:41+0100
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
frequenciesfile <- '_frequencies-RposteriorRIB_0-V13-D539-K45-I1024.rds'
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

datafile <- 'data_transformed_shuffled.csv'
alldata <- fread(datafile, sep=',')
alldata <- alldata[Usage_ == 'train']


pdff('test_intcorr_1')
for(acov in integerCovs[c(1:6,8)]){
nn <- covMaxs[acov]
##exn <- function(l,r=1){sapply(1:nn,function(x){prod(rep(x,r))}) %*% exp(l*(1:nn))/sum(exp(l*1:nn))}
dtgamma <- function(x,a,b){x^(a-1)*exp(-b*x)/sum((1:nn)^(a-1)*exp(-b*(1:nn)))}
exn <- function(r=1,a,b){sum(sapply(1:nn,function(x){prod(rep(x,r))}) *
                             ((1:nn)^(a-1) * exp(-b*(1:nn))))/
                             sum((1:nn)^(a-1) * exp(-b*(1:nn)))}
exln <- function(r=1,a,b){sum(sapply(log(1:nn),function(x){prod(rep(x,r))}) *
                             ((1:nn)^(a-1) * exp(-b*(1:nn))))/
                             sum((1:nn)^(a-1) * exp(-b*(1:nn)))}
exnln <- function(r=1,a,b){sum(sapply((1:nn)*log(1:nn),function(x){prod(rep(x,r))}) *
                             ((1:nn)^(a-1) * exp(-b*(1:nn))))/
                             sum((1:nn)^(a-1) * exp(-b*(1:nn)))}
##exq <- function(r=1,s=0,a=1,b=1){beta(a+r,b+s)/beta(a,b)}
exq <- function(r=1,a=1,b=1){beta(a+r,b)/beta(a,b)}
##
##
mm <- medianicovs[acov]
mv <- (widthicovs[acov])^2
vm <- log((widthicovs[acov]/2)^2)
vv <- log(sqrt(10))
fn <- function(par, corr=F){
    par <- exp(par)
    an <- par[1]
    bn <- par[2]
    aq <- par[3]
    bq <- par[4]
    n1 <- exn(r=1,a=an,b=bn)
    q1 <- exq(r=1,a=aq,b=bq)
    mm0 <- n1*q1
    mv0 <- exn(r=2,a=an,b=bn)*exq(r=2,a=aq,b=bq) - mm0*mm0
    vl <- exln(r=1,a=an,b=bn)
    vm0 <- vl + digamma(aq) + digamma(bq) - 2*digamma(aq+bq)
    vv0 <- exln(r=2,a=an,b=bn)-vl*vl + trigamma(aq) + trigamma(bq) - 4*trigamma(aq+bq)
    cv <- (exnln(r=1,a=an,b=bn) + n1*(digamma(aq+1) + digamma(bq) -2*digamma(aq+1+bq)))*q1 - mm0 * vm0
    if(corr){denom <- mv0*vv0}else{denom <- 1}
    (mm0/mm - 1)^2 + (mv0/mv - 1)^2 + (vm0/(vm) - 1)^2 + (vv0/vv - 1)^2 + (cv^2)/denom
}
##
test <- myoptim(par=log(c(1,1/(1), 1,1)), fn=fn) ###########################
test$par <- exp(test$par)
test
##
par(mfrow=c(3,2))
tplot(xgr <- seq(0,1,length.out=256), dbeta(xgr, test$par[3],test$par[4]), main='Q', ylim=c(0,NA))
legend('top',legend=signif(test$par[3:4],3), bty='n')
tplot(1:nn, dtgamma(1:nn, test$par[1],test$par[2]),main='N', ylim=c(0,NA))
legend('top',legend=signif(test$par[1:2],3), bty='n')
##
npts <- 10000
ssamples <- rcat(npts, prob=((1:nn)^(test$par[1]-1) * exp(-test$par[2]*(1:nn))))
qsamples <- rbeta(npts, shape1=test$par[3], shape2=test$par[4])
meansamples <- qsamples*ssamples
logsdsamples <- log(qsamples*(1-qsamples)*ssamples)/2
#par(mfrow=c(2,2))
histm <- thist(meansamples)
histv <- thist(logsdsamples)
tplot(histm$mid, histm$density, main='mean', ylim=c(0,NA), xlim=c(covMins[acov],nn))
abline(v=mm)
abline(v=mm+sqrt(mv), col=2)
abline(v=mm-sqrt(mv), col=2)
legend('topright',legend=nn, bty='n')
tplot(histv$mid, histv$density, main='log-sd', ylim=c(0,NA), xlim=c(min(logsdsamples, vm/2-sqrt(vv)),max(logsdsamples,vm/2+sqrt(vv))))
abline(v=vm/2)
abline(v=vm/2+sqrt(vv), col=2)
abline(v=vm/2-sqrt(vv), col=2)
##par(mfrow=c(1,1))
tplot(x=ssamples*qsamples, y=log(ssamples*qsamples*(1-qsamples))/2, xlab='mean', ylab='log-sd', type='p', pch='.', ylim=c(min(logsdsamples, vm/2-sqrt(vv)),max(logsdsamples,vm/2+sqrt(vv))))
abline(h=vm/2)
abline(h=vm/2+sqrt(vv), col=2)
abline(h=vm/2-sqrt(vv), col=2)
abline(v=mm)
abline(v=mm+sqrt(mv), col=2)
abline(v=mm-sqrt(mv), col=2)
}
dev.off()



pdff('test_intcorr_flat')
for(acov in integerCovs[c(1:6,8)]){
nn <- covMaxs[acov]
##exn <- function(l,r=1){sapply(1:nn,function(x){prod(rep(x,r))}) %*% exp(l*(1:nn))/sum(exp(l*1:nn))}
dtgamma <- function(x,a,b){x^(a-1)*exp(-b*x)/sum((1:nn)^(a-1)*exp(-b*(1:nn)))}
exn <- function(r=1,a,b){sum(sapply(1:nn,function(x){prod(rep(x,r))}) *
                             ((1:nn)^(a-1) * exp(-b*(1:nn))))/
                             sum((1:nn)^(a-1) * exp(-b*(1:nn)))}
exln <- function(r=1,a,b){sum(sapply(log(1:nn),function(x){prod(rep(x,r))}) *
                             ((1:nn)^(a-1) * exp(-b*(1:nn))))/
                             sum((1:nn)^(a-1) * exp(-b*(1:nn)))}
exnln <- function(r=1,a,b){sum(sapply((1:nn)*log(1:nn),function(x){prod(rep(x,r))}) *
                             ((1:nn)^(a-1) * exp(-b*(1:nn))))/
                             sum((1:nn)^(a-1) * exp(-b*(1:nn)))}
##exq <- function(r=1,s=0,a=1,b=1){beta(a+r,b+s)/beta(a,b)}
exq <- function(r=1,a=1,b=1){beta(a+r,b)/beta(a,b)}
##
##
mm <- nn/2
mv <- nn^2/12
vm <- log(nn^2/12)
vv <- log(sqrt(10))
fn <- function(par, corr=F){
    par <- exp(par)
    an <- par[1]
    bn <- par[2]
    aq <- 3+par[3]
    bq <- 1+par[4]
    n1 <- exn(r=1,a=an,b=bn)
    q1 <- exq(r=1,a=aq,b=bq)
    mm0 <- n1*q1
    mv0 <- exn(r=2,a=an,b=bn)*exq(r=2,a=aq,b=bq) - mm0*mm0
    vl <- exln(r=1,a=an,b=bn)
    vm0 <- vl + digamma(aq) + digamma(bq) - 2*digamma(aq+bq)
    vv0 <- exln(r=2,a=an,b=bn)-vl*vl + trigamma(aq) + trigamma(bq) - 4*trigamma(aq+bq)
    cv <- (exnln(r=1,a=an,b=bn) + n1*(digamma(aq+1) + digamma(bq) -2*digamma(aq+1+bq)))*q1 - mm0 * vm0
    if(corr){denom <- mv0*vv0}else{denom <- 1}
    (mm0/mm - 1)^2 + (mv0/mv - 1)^2 + (vm0/(vm) - 1)^2 + (vv0/vv - 1)^2 + (cv^2)/denom
}
##
test <- myoptim(par=log(c(1,1/(1), 1,1)), fn=fn) ###########################
test$par <- exp(test$par) + c(0,0,3,1)
test
##
par(mfrow=c(3,2))
tplot(xgr <- seq(0,1,length.out=256), dbeta(xgr, test$par[3],test$par[4]), main='Q', ylim=c(0,NA))
legend('top',legend=signif(test$par[3:4],3), bty='n')
tplot(1:nn, dtgamma(1:nn, test$par[1],test$par[2]),main='N', ylim=c(0,NA))
legend('top',legend=signif(test$par[1:2],3), bty='n')
##
npts <- 10000
ssamples <- rcat(npts, prob=((1:nn)^(test$par[1]-1) * exp(-test$par[2]*(1:nn))))
qsamples <- rbeta(npts, shape1=test$par[3], shape2=test$par[4])
meansamples <- qsamples*ssamples
logsdsamples <- log(qsamples*(1-qsamples)*ssamples)/2
#par(mfrow=c(2,2))
histm <- thist(meansamples)
histv <- thist(logsdsamples)
tplot(histm$mid, histm$density, main='mean', ylim=c(0,NA), xlim=c(covMins[acov],nn))
abline(v=mm)
abline(v=mm+sqrt(mv), col=2)
abline(v=mm-sqrt(mv), col=2)
legend('topright',legend=nn, bty='n')
tplot(histv$mid, histv$density, main='log-sd', ylim=c(0,NA), xlim=c(min(logsdsamples, vm/2-sqrt(vv)),max(logsdsamples,vm/2+sqrt(vv))))
abline(v=vm/2)
abline(v=vm/2+sqrt(vv), col=2)
abline(v=vm/2-sqrt(vv), col=2)
##par(mfrow=c(1,1))
tplot(x=ssamples*qsamples, y=log(ssamples*qsamples*(1-qsamples))/2, xlab='mean', ylab='log-sd', type='p', pch=1, alpha=0.25, ylim=c(min(logsdsamples, vm/2-sqrt(vv)),max(logsdsamples,vm/2+sqrt(vv))))
abline(h=vm/2)
abline(h=vm/2+sqrt(vv), col=2)
abline(h=vm/2-sqrt(vv), col=2)
abline(v=mm)
abline(v=mm+sqrt(mv), col=2)
abline(v=mm-sqrt(mv), col=2)
}
dev.off()





pdff('test_intcorr_mm')
for(acov in integerCovs[c(1:6,8)]){
nn <- covMaxs[acov]
##exn <- function(l,r=1){sapply(1:nn,function(x){prod(rep(x,r))}) %*% exp(l*(1:nn))/sum(exp(l*1:nn))}
dtgamma <- function(x,a,b){x^(a-1)*exp(-b*x)/sum((1:nn)^(a-1)*exp(-b*(1:nn)))}
exn <- function(r=1,a,b){sum(sapply(1:nn,function(x){prod(rep(x,r))}) *
                             ((1:nn)^(a-1) * exp(-b*(1:nn))))/
                             sum((1:nn)^(a-1) * exp(-b*(1:nn)))}
exln <- function(r=1,a,b){sum(sapply(log(1:nn),function(x){prod(rep(x,r))}) *
                             ((1:nn)^(a-1) * exp(-b*(1:nn))))/
                             sum((1:nn)^(a-1) * exp(-b*(1:nn)))}
exnln <- function(r=1,a,b){sum(sapply((1:nn)*log(1:nn),function(x){prod(rep(x,r))}) *
                             ((1:nn)^(a-1) * exp(-b*(1:nn))))/
                             sum((1:nn)^(a-1) * exp(-b*(1:nn)))}
##exq <- function(r=1,s=0,a=1,b=1){beta(a+r,b+s)/beta(a,b)}
exq <- function(r=1,a=1,b=1){beta(a+r,b)/beta(a,b)}
##
##
mm <- medianicovs[acov]
mv <- (widthicovs[acov])^2
vm <- log((widthicovs[acov]/2)^2)
vv <- log(sqrt(10))
fn <- function(par, corr=F){
    par <- exp(par)
    an <- par[1]
    bn <- par[2]
    aq <- 2+par[3]
    bq <- 1+par[4]
    n1 <- exn(r=1,a=an,b=bn)
    q1 <- exq(r=1,a=aq,b=bq)
    mm0 <- n1*q1
    mv0 <- exn(r=2,a=an,b=bn)*exq(r=2,a=aq,b=bq) - mm0*mm0
    vl <- exln(r=1,a=an,b=bn)
    vm0 <- vl + digamma(aq) + digamma(bq) - 2*digamma(aq+bq)
    vv0 <- exln(r=2,a=an,b=bn)-vl*vl + trigamma(aq) + trigamma(bq) - 4*trigamma(aq+bq)
    cv <- (exnln(r=1,a=an,b=bn) + n1*(digamma(aq+1) + digamma(bq) -2*digamma(aq+1+bq)))*q1 - mm0 * vm0
    if(corr){denom <- mv0*vv0}else{denom <- 1}
    (mm0/mm - 1)^2 + (mv0/mv - 1)^2 + (vm0/(vm) - 1)^2 + (vv0/vv - 1)^2 + (cv^2)/denom
}
##
test <- myoptim(par=log(c(1,1/(mm), 1, 1)), fn=fn) ###### mm is ok
test$par <- c(exp(test$par)) + c(0,0,2,1)
test
##
par(mfrow=c(3,2))
tplot(xgr <- seq(0,1,length.out=256), dbeta(xgr, test$par[3],test$par[4]), main='Q', ylim=c(0,NA))
legend('top',legend=signif(test$par[3:4],3), bty='n')
tplot(1:nn, dtgamma(1:nn, test$par[1],test$par[2]),main='N', ylim=c(0,NA))
legend('top',legend=signif(test$par[1:2],3), bty='n')
##
npts <- 10000
ssamples <- rcat(npts, prob=((1:nn)^(test$par[1]-1) * exp(-test$par[2]*(1:nn))))
qsamples <- rbeta(npts, shape1=test$par[3], shape2=test$par[4])
meansamples <- qsamples*ssamples
logsdsamples <- log(qsamples*(1-qsamples)*ssamples)/2
#par(mfrow=c(2,2))
histm <- thist(meansamples)
histv <- thist(logsdsamples)
tplot(histm$mid, histm$density, main='mean', ylim=c(0,NA), xlim=c(covMins[acov],nn))
abline(v=mm)
abline(v=mm+sqrt(mv), col=2)
abline(v=mm-sqrt(mv), col=2)
legend('topright',legend=nn, bty='n')
tplot(histv$mid, histv$density, main='log-sd', ylim=c(0,NA), xlim=c(min(logsdsamples, vm/2-sqrt(vv)),max(logsdsamples,vm/2+sqrt(vv))))
abline(v=vm/2)
abline(v=vm/2+sqrt(vv), col=2)
abline(v=vm/2-sqrt(vv), col=2)
##par(mfrow=c(1,1))
tplot(x=ssamples*qsamples, y=log(ssamples*qsamples*(1-qsamples))/2, xlab='mean', ylab='log-sd', type='p', pch='.', ylim=c(min(logsdsamples, vm/2-sqrt(vv)),max(logsdsamples,vm/2+sqrt(vv))))
abline(h=vm/2)
abline(h=vm/2+sqrt(vv), col=2)
abline(h=vm/2-sqrt(vv), col=2)
abline(v=mm)
abline(v=mm+sqrt(mv), col=2)
abline(v=mm-sqrt(mv), col=2)
}
dev.off()






pdff('test_intcorr_unifb')
for(acov in integerCovs[c(1:6,8)]){
nn <- covMaxs[acov]
##
##
mm <- medianicovs[acov]
mv <- (widthicovs[acov]/2)^2
vm <- log((widthicovs[acov]/2)^2)
vv <- log(sqrt(10))
test <- list(par=c(2,0,1,1))
test
##
par(mfrow=c(3,2))
tplot(xgr <- seq(0,1,length.out=256), dbeta(xgr, test$par[3],test$par[4]), main='Q', ylim=c(0,NA))
legend('top',legend=signif(test$par[3:4],3), bty='n')
tplot(1:nn, dtgamma(1:nn, test$par[1],test$par[2]),main='N', ylim=c(0,NA))
legend('top',legend=signif(test$par[1:2],3), bty='n')
##
npts <- 10000
ssamples <- rcat(npts, prob=((1:nn)^(test$par[1]-1) * exp(-test$par[2]*(1:nn))))
qsamples <- rbeta(npts, shape1=test$par[3], shape2=test$par[4])
meansamples <- qsamples*ssamples
logsdsamples <- log(qsamples*(1-qsamples)*ssamples)/2
#par(mfrow=c(2,2))
histm <- thist(meansamples)
histv <- thist(logsdsamples)
tplot(histm$mid, histm$density, main='mean', ylim=c(0,NA), xlim=c(covMins[acov],nn))
abline(v=mm)
abline(v=mm+sqrt(mv), col=2)
abline(v=mm-sqrt(mv), col=2)
legend('topright',legend=nn, bty='n')
tplot(histv$mid, histv$density, main='log-sd', ylim=c(0,NA), xlim=c(min(logsdsamples, vm/2-sqrt(vv)),max(logsdsamples,vm/2+sqrt(vv))))
abline(v=vm/2)
abline(v=vm/2+sqrt(vv), col=2)
abline(v=vm/2-sqrt(vv), col=2)
##par(mfrow=c(1,1))
tplot(x=ssamples*qsamples, y=log(ssamples*qsamples*(1-qsamples))/2, xlab='mean', ylab='log-sd', type='p', pch=1, alpha=0.25, ylim=c(min(logsdsamples, vm/2-sqrt(vv)),max(logsdsamples,vm/2+sqrt(vv))))
abline(h=vm/2)
abline(h=vm/2+sqrt(vv), col=2)
abline(h=vm/2-sqrt(vv), col=2)
abline(v=mm)
abline(v=mm+sqrt(mv), col=2)
abline(v=mm-sqrt(mv), col=2)
}
dev.off()


