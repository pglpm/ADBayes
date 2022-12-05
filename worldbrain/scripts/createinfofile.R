library('data.table')
library('png')
library('foreach')

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
ddd <- 2^-6
varinfo[variate$T,'scale'] <- (varinfo[variate$T,'max']-varinfo[variate$T,'min'])/(1-2*ddd)
varinfo[variate$T,'location'] <- varinfo[variate$T,'min']-ddd*varinfo[variate$T,'scale']
varinfo[variate$T,'n'] <- ddd
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



## write.csv(varinfo, 'varinfo.csv')
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
