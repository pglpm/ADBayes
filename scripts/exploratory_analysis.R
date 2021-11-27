## Author: PGL  Porta Mana
## Created: 2021-03-20T10:07:17+0100
## Last-Updated: 2021-11-27T16:58:22+0100
################
## Exploration for MMIV poster
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


## Variates to study
covSelect <- c(
    'Subgroup_' # Diagnostic label
    ,'AGE'
    ,'Gender_num_'
    ,'RAVLT_immediate' # immediate memory function
    ,'AVDEL30MIN_neuro' # delayed memory function
    ,'AVDELTOT_neuro' # recognition memory
    ,'TRAASCOR_neuro' # executive function
    ,'TRABSCOR_neuro' # executive funtion
    ,'CATANIMSC_neuro' # executive function
    ,'LRHHC_n_long' #  total volume of hippocampus volum (both hemispheres)
    ,'Apoe4_' #  binary variable ( no alleles / at least one ellele on ApoE-gene)
    ,'GDTOTAL_gds' #   degrees of depressive symptoms
    ,'ANARTERR_neuro' #   premorbid intelligence 
)
##
covtypes = rep(NA, length(covSelect)) #c('integer', 'double', 'integer', 'integer', 'integer', 'integer', 'integer', 'integer', 'integer', 'double', 'integer', 'integer', 'integer'),
covmin = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
covmax = c(1, +Inf, 1, 75, 15, 15, 150, 300, 60, +Inf, 1, 5, 50)
covtransf <- rep(NA, length(covSelect))
names(covtypes) <- names(covmin) <- names(covmax) <- names(covtransf) <- covSelect
##
## load data and select variates of interest
datafile <- 'RF_K50_final_df_best_blended_top5_K50.csv'
origdata <- fread(datafile, sep=',')[, ..covSelect]
##
## Convert categorical variables to integer
## Make sure that integers are represented as integers (not double)
## Transform positive continuous variates to log scale and standardize
covNames <- covSelect
names(covNames) <- covSelect
for(acov in covSelect){
    datum <- origdata[[acov]]
    ##
    if(is.character(datum)){
        cats <- unique(datum)
        newacov <- paste0('int_', acov)
        origdata[[newacov]] <- as.integer(factor(datum)) - 1L
        message(paste0('* ', acov, ': converting to integer "', newacov, '"'))
        message(paste0(paste0(cats, ' = ', as.integer(factor(cats))-1), collapse='\n'))
        covNames[acov] <- newacov
        covtypes[acov] <- 'integer'
    }
    ##
    else if(is.double(datum) && (
        (is.na(covtypes[acov]) && all(round(datum)==datum)) ||
        (!is.na(covtypes[acov]) && covtypes[acov] == 'integer')
    )){
        message(paste0('* ', acov, ': re-classifying as integer'))
        origdata[[acov]] <- as.integer(datum)
        covtypes[acov] <- 'integer'
    }
    ##
    else if(is.double(datum) && (
        (is.na(covtypes[acov]) && !all(round(datum)==datum)) ||
        (!is.na(covtypes[acov]) && covtypes[acov] == 'double')
    ) && covmin[acov] == 0
    ){
        datum <- log(datum)
        m <- signif(mean(datum),2)
        s <- signif(sd(datum),2)
        newacov <- paste0('log_', acov)
        origdata[[newacov]] <- (datum-m)/s
        message(paste0('* ', acov, ': transforming to log-scale and standardizing as "', newacov, '"'))
        message(paste0(acov, ' = exp(A * log_', acov, ' + B)'))
        message(paste0('with  A = ', s, ',  B = ', m))
        covNames[acov] <- newacov
        covtransf[acov] <- m + s*1i
        covtypes[acov] <- 'double'
    } else if(is.na(covtypes[acov])){covtypes[acov] <- typeof(datum)}
    if(is.na(covmin[acov])){covmin[acov] <- min(datum, na.rm=TRUE)}
    if(is.na(covmax[acov])){covmax[acov] <- max(datum, na.rm=TRUE)}
}
##
message('\nSaving original order in "orig_id" and shuffling')
origdata[['orig_id']] <- as.integer(1:nrow(origdata))
set.seed(123)
origdata[1:nrow(origdata)] <- origdata[sample(nrow(origdata))]
##
savefile <- 'data_transformed_shuffled.csv'
message(paste0('\nSaving data in "', savefile, '"'))
fwrite(origdata, savefile, sep=',')
##
## Safe the information about the variates
variateinfo <- data.table(
    variate = covNames,
    orig_name = covSelect,
    type = covtypes,
    min = covmin,
    max = covmax,
    transfparams = covtransf
)
saveinfofile <- 'variates_info.csv'
message(paste0('\nSaving variate info in "', savefile, '"'))
fwrite(variateinfo, saveinfofile, sep=',')

## Data table having only variates we'll work with
alldata <- origdata[, ..covNames]
##
## construct bins to calculate mutual info
doplots <- TRUE
nSamples <- nrow(alldata)
nCovs <- length(covNames)
nbinsq <- 8
breaksCov <- list()
for(acov in covNames){
    datum <- alldata[[acov]]
    rg <- range(datum, na.rm=TRUE)
    if(is.integer(datum)){
        breaksCov[[acov]] <- unique(quant(datum, (1:(nbinsq-1))/nbinsq))
        ##  breaksCov[[acov]] <- unique(round((seq(min(datum), max(datum), length.out=nbinsq+1))))[-c(1,nbinsq+1)]
    } else {
        breaksCov[[acov]] <- quant(datum, (1:(nbinsq-1))/nbinsq)
        ##  breaksCov[[acov]] <- seq(min(datum), max(datum), length.out=nbinsq+1)[-c(1,nbinsq+1)]
    }
}
names(breaksCov) <- covNames
##
## calculate mutual infos and conditional entropies
mainCov <- 'int_Subgroup_'
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


