library('data.table')
library('png')
library('foreach')

Tprob <- 2^-6

dt <- fread('~/repositories/ADBayes/worldbrain/scripts/ingrid_data_nogds6.csv')
varinfo <- read.csv('~/repositories/ADBayes/worldbrain/scripts/varinfo.csv', row.names=1)
varnames <- rownames(varinfo)


variate <- lapply(variatetypes, function(x)rownames(varinfo)[varinfo['type']==x])
len <- lapply(variate,length)

## real
varinfo[variate$R,'location'] <- apply(dt[,variate$R,with=F],2,median,na.rm=T)
varinfo[variate$R,'scale'] <- apply(dt[,variate$R,with=F],2,IQR,na.rm=T)*sd2iqr
varinfo[variate$R,'plotmin'] <- apply(dt[,variate$R,with=F],2,function(x){
    min(x, na.rm=T) - diff(tquant(x, c(0.25,0.5)))
    })
varinfo[variate$R,'plotmax'] <- apply(dt[,variate$R,with=F],2,function(x){
    max(x, na.rm=T) + diff(tquant(x, c(0.75,0.5)))
    })
##
varinfo[variate$L,'location'] <- apply(log(dt[,variate$L,with=F]),2,median,na.rm=T)
varinfo[variate$L,'scale'] <- apply(log(dt[,variate$L,with=F]),2,IQR,na.rm=T)*sd2iqr
varinfo[variate$L,'plotmin'] <- apply(dt[,variate$L,with=F],2,function(x){
    max(0, min(x, na.rm=T) - diff(tquant(x, c(0.25,0.5))))
    })
varinfo[variate$L,'plotmax'] <- apply(dt[,variate$L,with=F],2,function(x){
    max(x, na.rm=T) + diff(tquant(x, c(0.75,0.5)))
    })
##
varinfo[variate$B,'location'] <- varinfo[variate$B,'min']
varinfo[variate$B,'scale'] <- (varinfo[variate$B,'max']-varinfo[variate$B,'min'])/(varinfo[variate$B,'n']-1L)
varinfo[variate$B,'plotmin'] <- 0
varinfo[variate$B,'plotmax'] <- 1
##
varinfo[variate$I,'location'] <- varinfo[variate$I,'min']
varinfo[variate$I,'scale'] <- (varinfo[variate$I,'max']-varinfo[variate$I,'min'])/(varinfo[variate$I,'n']-1L)
varinfo[variate$I,'plotmin'] <- sapply(variate$I,function(v){
    dat <- dt[[v]]
    max(varinfo[v,'min'], min(dat, na.rm=T) - diff(tquant(dat, c(0.25,0.5))))
    })
varinfo[variate$I,'plotmax'] <- sapply(variate$I,function(v){
    dat <- dt[[v]]
    min(varinfo[v,'max'], max(dat, na.rm=T) + diff(tquant(dat, c(0.5,0.75))))
    })
##
varinfo[variate$T,'scale'] <- (varinfo[variate$T,'max']-varinfo[variate$T,'min'])/(1-2*Tprob)
varinfo[variate$T,'location'] <- varinfo[variate$T,'min']-Tprob*varinfo[variate$T,'scale']
varinfo[variate$T,'n'] <- Tprob
varinfo[variate$T,'plotmin'] <- sapply(variate$T,function(v){
    dat <- dt[[v]]
    max(varinfo[v,'min'], min(dat, na.rm=T) - diff(tquant(dat, c(0.25,0.5))))
    })
varinfo[variate$T,'plotmax'] <- sapply(variate$T,function(v){
    dat <- dt[[v]]
    min(varinfo[v,'max'], max(dat, na.rm=T) + diff(tquant(dat, c(0.5,0.75))))
    })


##
varinfo[varnames,'Q1'] <- apply(dt[,..varnames], 2, tquant, 0.25)
varinfo[varnames,'Q2'] <- apply(dt[,..varnames], 2, tquant, 0.5)
varinfo[varnames,'Q3'] <- apply(dt[,..varnames], 2, tquant, 0.75)


varinfo[variate$R,'mean'] <- 0L
varinfo[variate$L,'mean'] <- 0L
varinfo[variate$T,'mean'] <- 0L
varinfo[variate$I,'mean'] <- 0L
varinfo[variate$B,'mean'] <- NA
varinfo[variate$C,'mean'] <- NA
##
varinfo[variate$R,'sd'] <- 3 #***
varinfo[variate$L,'sd'] <- 3 #***
varinfo[variate$T,'sd'] <- 1
varinfo[variate$I,'sd'] <- 7/8
varinfo[variate$B,'sd'] <- NA
varinfo[variate$C,'sd'] <- NA
##
varinfo[variate$R,'shapeout'] <- 0.5 #***
varinfo[variate$L,'shapeout'] <- 0.5 #***
varinfo[variate$T,'shapeout'] <- 1
varinfo[variate$I,'shapeout'] <- 1
varinfo[variate$B,'shapeout'] <- 1
varinfo[variate$C,'shapeout'] <- 1
##
varinfo[variate$R,'shapein'] <- 0.5 #***
varinfo[variate$L,'shapein'] <- 0.5 #***
varinfo[variate$T,'shapein'] <- 1
varinfo[variate$I,'shapein'] <- 1
varinfo[variate$B,'shapein'] <- 1
varinfo[variate$C,'shapein'] <- NA
##
varinfo[variate$R,'varscale'] <- 1 #***
varinfo[variate$L,'varscale'] <- 1 #***
varinfo[variate$T,'varscale'] <- 1/4
varinfo[variate$I,'varscale'] <- 1/4
varinfo[variate$B,'varscale'] <- NA
varinfo[variate$C,'varscale'] <- NA


varinfo[variate$L,'mean'] <- apply(data.matrix(dt)[,variate$R,drop=F], 2, median, na.rm=T)

predictors <- setdiff(unlist(variate),'Subgroup_num_')

## write.csv(varinfo, 'varinfo.csv')
## write.table(predictors, 'predictors.csv', row.names=F, col.names=F)
## varinfo <- data.matrix(read.csv('varinfo.csv', row.names=1))

summary(dt)

summary(sapply(colnames(dt),function(v){transfdir(data.matrix(dt[,..v]), varinfo)}))

summary(sapply(colnames(dt),function(v){recjacobian(data.matrix(dt[,..v]), varinfo)}))


sapply(colnames(dt),function(v){
    x <- data.matrix(dt[,..v])
    y <- transfinv(transfdir(x,varinfo),varinfo)
    rdiff <- y/x-1
    rdiff[x==0 & y==0] <- 0
    rdiff[x==0 & y!=0] <- Inf
    max(abs(rdiff),na.rm=T)
})

v <- 'GDTOTAL_gds'
varinfo[v,]

transfdir(0:6,v)
transfinv(0:6,v)
