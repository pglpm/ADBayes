## Author: PGL  Porta Mana
## Created: 2021-03-20T10:07:17+0100
## Last-Updated: 2021-12-17T18:39:42+0100
################
## Data preparation for Alzheimer study - FAQ version
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

## NB: These functions require normalized frequencies
## Function to calculate mutual info from matrix of frequencies
mutualinfo <- function(freqs, base=2){##in bits by default
    freqs1 <- rowSums(freqs)
    freqs2 <- colSums(freqs)
    sum(freqs * log2(freqs/outer(freqs1, freqs2)), na.rm=TRUE)/log2(base)
}
## Function to calculate conditional entropy of rows given cols
condentropy12 <- function(freqs, base=2){##in bits by default
    freqs1 <- rowSums(freqs)
    freqs2 <- colSums(freqs)
    sum(freqs * log2(outer(rep(1,nrow(freqs)), freqs2)/freqs), na.rm=TRUE)/log2(base)
}
## Function to calculate entropy
entropy <- function(freqs, base=2){##in bits by default
    -sum(freqs * log2(freqs), na.rm=TRUE)/log2(base)
}
## function to normalize absolute frequencies
normalize <- function(freqs){freqs/sum(freqs)}

datafile <- '3.13_pycharet_K50.csv'
origdata <- fread(datafile, sep=',')

## Variates to study and their info
covSelect <- 'Usage_'
variateinfo <- rbind(data.table(variate=character(), type=character(), min=numeric(), max=numeric(), transfM=numeric(), transfW=numeric()),
list('Subgroup_num_', NA, 0, 1),
list('AGE', 'double', 0, +Inf),
list('Gender_num_', NA, 0, 1),
list('RAVLT_immediate', NA, 0, 75),
list('AVDEL30MIN_neuro', NA, 0, 15),
list('AVDELTOT_neuro', NA, 0, 15),
list('TRAASCOR_neuro', NA, 0, 150),
list('TRABSCOR_neuro', NA, 0, 300),
list('CATANIMSC_neuro', NA, 0, 60),
list('LRHHC_n_long', 'double', 0, +Inf),
list('Apoe4_', NA, 0, 1),
list('GDTOTAL_gds', NA, 0, 6),
list('ANARTERR_neuro', NA, 0, 50),
fill=T)
covNames <- variateinfo$variate
variateinfo$origvariate <- covNames
##
selectdata <- origdata

## Convert categorical variables to integer
## Make sure that integers are represented as integers (not double)
## Transform positive continuous variates to log scale and standardize
for(acov in covNames){
    ##print(acov)
    datum <- selectdata[[acov]]
    ##
    if(is.character(datum)){# string
        cats <- unique(datum)
        newacov <- paste0(acov, 'num')
        newvalues <- as.integer(factor(datum)) - 1L
        selectdata[[newacov]] <- datum
        message(paste0('* ', acov, ': converting to integer "', newacov, '"'))
        message(paste0(paste0(sort(cats), ' = ', as.integer(factor(sort(cats)))-1L), collapse='\n'))
        variateinfo[origvariate==acov, variate] <- newacov
        variateinfo[origvariate==acov, min] <- min(newvalues, na.rm=T)
        variateinfo[origvariate==acov, max] <- max(newvalues, na.rm=T)
        if(diff(range(newvalues, na.rm=T))==1){
            variateinfo[origvariate==acov, type] <- 'binary'
        }else{
            variateinfo[origvariate==acov, type] <- 'integer'
        }
    } else if(# integer
               ( is.na(variateinfo[origvariate==acov, type]) &
                 all(selectdata[[acov]]==as.integer(selectdata[[acov]]), na.rm=T) ) |
               variateinfo[origvariate==acov, type]=='integer' |
               variateinfo[origvariate==acov, type]=='binary'
           ){
        datum <- as.integer(datum)
        ##
        thmin <- variateinfo[origvariate==acov, min]
        if(!is.na(thmin) & thmin>min(datum, na.rm=T)){
            message(paste0('WARNING ', acov, ': data smaller than theoretical min'))
            thmin <- min(datum, na.rm=T)
        }else if(is.na(thmin)){
            thmin <- min(datum, na.rm=T)
        }
        ##
        thmax <- variateinfo[origvariate==acov, max]
        if(!is.na(thmax) & thmax<max(datum, na.rm=T)){
            message(paste0('WARNING ', acov, ': data larger than theoretical max'))
            thmax <- max(datum, na.rm=T)
        }else if(is.na(thmax)){
            thmax <- max(datum, na.rm=T)
        }
        ##
        if(thmax-thmin==1){
            message(paste0('* ', acov, ': classifying as binary'))
            variateinfo[origvariate==acov, 'type'] <- 'binary'
            newacov <- acov
            if(thmin!=0){
                newacov <- paste0(acov, '_bin')
                message(paste0(acov,'=',thmin, ' becomes ',newacov,'=0'))
                message(paste0(acov,'=',thmax, ' becomes ',newacov,'=1'))
                thmin <- 0L
                thmax <- 1L
            }
            selectdata[[newacov]] <- as.integer(datum - thmin)
            variateinfo[origvariate==acov, 'variate'] <- newacov
        } else {
            message(paste0('* ', acov, ': classifying as integer'))
            selectdata[[acov]] <- datum
            variateinfo[origvariate==acov, 'type'] <- 'integer'            
            selectdata[[acov]] <- datum
        }
        variateinfo[origvariate==acov, 'min'] <- thmin
        variateinfo[origvariate==acov, 'max'] <- thmax
    }else if(# double
              ( is.na(variateinfo[origvariate==acov, type]) &
                !all(selectdata[[acov]]==as.integer(selectdata[[acov]]), na.rm=T) ) |
                            variateinfo[origvariate==acov, type]=='double'
          ){
        message(paste0('* ', acov, ': classifying as double'))
        datum <- as.double(datum)
        if(is.na(variateinfo[origvariate==acov, type])){
            variateinfo[origvariate==acov, 'type'] <- 'double'
        }
        ##
        thmin <- variateinfo[origvariate==acov, min]
        if(!is.na(thmin) & thmin==0){
            datum <- log(datum)
            m <- signif(median(datum, na.rm=T),2)
            s <- signif(IQR(datum, type=8, na.rm=T)/(2*qnorm(3/4)),2)
            newacov <- paste0(acov, '_log')
            selectdata[[newacov]] <- (datum-m)/s
            message(paste0('* ', acov, ': transforming to log-scale and standardizing as "', newacov, '"'))
            message(paste0(acov, ' = exp(W * ', newacov, ' + M)'))
            message(paste0('with  W = ', s, ',  M = ', m))
            variateinfo[origvariate==acov, 'variate'] <- newacov
            variateinfo[origvariate==acov, 'max'] <- +Inf
            variateinfo[origvariate==acov, 'min'] <- -Inf
            variateinfo[origvariate==acov, 'transfM'] <- m
            variateinfo[origvariate==acov, 'transfW'] <- s
        }
    }else{message(paste0('WARNING ', acov, ': cannot categorize'))}
}
    
## Safe the information about the variates
saveinfofile <- 'variatesIng_info.csv'
message(paste0('\nSaving variate info in "', saveinfofile, '"'))
fwrite(variateinfo, saveinfofile, sep=',')
##
message('\nSaving original order in "orig_id" and shuffling')
selectdata$orig_id <- as.integer(1:nrow(selectdata))
##
traindata <- selectdata[get(covSelect)=='train', .SD, .SDcols=variateinfo$variate]
set.seed(123)
traindata[1:nrow(traindata)] <- traindata[sample(nrow(traindata))]
savefile <- 'traindataIng_transformed_shuffled.csv'
message(paste0('\nSaving training data in "', savefile, '"'))
fwrite(traindata, savefile, sep=',')
##
testdata <- selectdata[get(covSelect)=='test', .SD, .SDcols=variateinfo$variate]
set.seed(123)
testdata[1:nrow(testdata)] <- testdata[sample(nrow(testdata))]
savefile <- 'testdataIng_transformed_shuffled.csv'
message(paste0('\nSaving test data in "', savefile, '"'))
fwrite(testdata, savefile, sep=',')
##
alldata <- selectdata[, .SD, .SDcols=variateinfo$variate]
set.seed(123)
alldata[1:nrow(alldata)] <- alldata[sample(nrow(alldata))]
savefile <- 'alldataIng_transformed_shuffled.csv'
message(paste0('\nSaving all data in "', savefile, '"'))
fwrite(alldata, savefile, sep=',')







##
selectdata <- selectdata[get(covSelect)=='train', .SD, .SDcols=variateinfo$variate]
message('\nSaving original order in "orig_id" and shuffling')
selectdata$orig_id <- as.integer(1:nrow(selectdata))
set.seed(123)
selectdata[1:nrow(selectdata)] <- selectdata[sample(nrow(selectdata))]
##
savefile <- 'dataIng_transformed_shuffled.csv'
message(paste0('\nSaving data in "', savefile, '"'))
fwrite(selectdata, savefile, sep=',')
##
## Safe the information about the variates
saveinfofile <- 'variatesIng_info.csv'
message(paste0('\nSaving variate info in "', savefile, '"'))
fwrite(variateinfo, saveinfofile, sep=',')






## Data table having only variates we'll work with
alldata <- origdata[, ..covNames]
##
## construct bins to calculate mutual info
doplots <- TRUE
nSamples <- nrow(alldata)
nCovs <- length(covNames)
nbinsq <- 16
breaksCov <- list()
for(acov in covNames){
    datum <- alldata[[acov]]
    rg <- range(datum, na.rm=TRUE)
    if(is.integer(datum)){
        #breaksCov[[acov]] <- unique(quant(datum, (1:(nbinsq-1))/nbinsq))
        breaksCov[[acov]] <- unique(round((seq(min(datum), max(datum), length.out=nbinsq+1))))[-c(1,nbinsq+1)]
    } else {
        ##breaksCov[[acov]] <- quant(datum, (1:(nbinsq-1))/nbinsq)
        breaksCov[[acov]] <- seq(min(datum), max(datum), length.out=nbinsq+1)[-c(1,nbinsq+1)]
    }
}
names(breaksCov) <- covNames
##
## calculate mutual infos and conditional entropies
mainCov <- 'Subgroup_num_'
minfos <- matrix(NA, 3, nCovs)
rownames(minfos) <- c('mutual_info', 'cond_entropy', 'entropy')
colnames(minfos) <- covNames
rangemaincov <- range(breaksCov[[mainCov]])
##
for(acov in covNames){
    freq <-  as.data.frame(table(findInterval(alldata[[mainCov]], breaksCov[[mainCov]]), findInterval(alldata[[acov]], breaksCov[[acov]])))
    freq[,1] <- as.numeric(freq[,1])
    freq[,2] <- as.numeric(freq[,2])
##
    freq2D <- matrix(0, nrow=length(breaksCov[[mainCov]])+1, ncol=length(breaksCov[[acov]])+1)
    freq2D[cbind(freq[,1], freq[,2])] <- freq[,3]
    freq2D <- freq2D/sum(freq2D)
    ##
    minfos[,acov] <- c(
        mutualinfo(freq2D),
        condentropy12(freq2D),
        entropy(colSums(freq2D))
    )
}
## sort orders according to mutualinfo and condentropies
sortmi <- order(minfos[1,], decreasing=TRUE)
names(sortmi) <- covNames[sortmi]
sortce <- order(minfos[2,], decreasing=FALSE)
names(sortce) <- covNames[sortce]
##
## Plots
if(doplots==TRUE){
    ndata <- nrow(alldata)
    pdff(paste0('histograms_transformed_data'))
    for(acov in covNames[sortmi]){
        datum <- alldata[[acov]]
##        breaks <- breakFeatures[[acov]]
        ## print(ggplot(alldata[,..acov], aes_(x=as.name(acov))) + geom_histogram(breaks=breaks))
        histo <- thist(datum)
        tplot(x=histo$breaks, y=histo$density, xlab=acov, ylab='relative frequency',
              main=paste0(acov, ': entropy = ', signif(minfos[3,acov],3), ' bit'))
    }
    dev.off()
    ##
    set.seed(123)
    rwidth <- 1/8
    randy <- runif(ndata, min=-rwidth, max=rwidth)
    maindatum <- alldata[[mainCov]]
    ylim <- c(-0.5, 1.5)
    pdff(paste0('jointplots_transformed_data'))
    for(acov in covNames[sortmi]){
        datum <- alldata[[acov]]
        mi <- signif(minfos[1,acov],3)
        en <- signif(minfos[3,acov],3)
        conden <- signif(minfos[2,acov],3)
        if(is.integer(datum)){
            xlim <- range(datum)+c(-1,1)*0.5
            xticks <- min(datum):max(datum)
            randx <- runif(ndata, min=-rwidth, max=rwidth)
        }else{
            xlim <- range(datum)+c(-1,1)*diff(range(datum))*0.01
            xticks <- pretty(xlim, n=10)
            randx <- 0
        }
        ## plot.new()
        ## par(mai=c(0.8,0.8,0.8,0))
        ## plot.window(xlim=xlim, ylim=ylim, xaxs='i', yaxs='i')
        ## plot.xy(
        ##     xy.coords(alldata[[acov]]+randx, y=alldata[[mainCov]]+randy),
        ##     type='p',
        ##     pch=16,
        ##     col=paste0(palette()[1],'80')
        ## )
        ## axis(1, at=xticks, lwd=0, lwd.ticks=0.5, col.ticks='#bbbbbb80')
        ## axis(2, at=0:1, lwd=0, lwd.ticks=0.5, col.ticks='#bbbbbb80')
        ## grid(nx=NULL, ny=NULL, col='#bbbbbb80', lty=1, lwd=0.5)
        tplot(x=alldata[[acov]]+randx, y=alldata[[mainCov]]+randy, type='p', pch=20,
              alpha=0.5,
              ylim=c(-1,2),
                xlab=paste0(acov, ' (entropy = ', en, ' bit)'),
              ylab=mainCov,
              main=paste0(acov,
                     ': mutual info = ',mi,' bit, conditional entropy = ',conden, ' bit'))
    }
    dev.off()
}
##
message(paste0(paste0(covNames[sortmi], ': mutual info = ',signif(minfos[1,sortmi],3),' bit'),collapse='\n'))
## message()
## message(paste0(paste0(covNames[sortce], ': cond entropy = ',signif(minfos[2,sortce],3),' bit'),collapse='\n'))


