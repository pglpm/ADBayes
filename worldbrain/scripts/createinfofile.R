library('data.table')
library('png')
library('foreach')

dt <- fread('~/repositories/ADBayes/worldbrain/scripts/ingrid_data_nogds6.csv')
varinfo <- read.csv('~/repositories/ADBayes/worldbrain/scripts/varinfo.csv', row.names=1)
varnames <- rownames(varinfo)

variates <- lapply(variatetypes, function(x)rownames(varinfo)[varinfo['type']==x])
n <- lapply(variates,length)

varinfo[variates$R,'location'] <- apply(dt[,variates$R,with=F],2,median,na.rm=T)
varinfo[variates$R,'scale'] <- apply(dt[,variates$R,with=F],2,IQR,na.rm=T)*sd2iqr
##
varinfo[variates$L,'location'] <- apply(log(dt[,variates$L,with=F]),2,median,na.rm=T)
varinfo[variates$L,'scale'] <- apply(log(dt[,variates$L,with=F]),2,IQR,na.rm=T)*sd2iqr
##
varinfo[variates$B,'location'] <- varinfo[variates$B,'min']
varinfo[variates$B,'scale'] <- (varinfo[variates$B,'max']-varinfo[variates$B,'min'])/(varinfo[variates$B,'n']-1L)
##
varinfo[variates$I,'location'] <- varinfo[variates$I,'min']
varinfo[variates$I,'scale'] <- (varinfo[variates$I,'max']-varinfo[variates$I,'min'])/(varinfo[variates$I,'n']-1L)
##
ddd <- 2^-6
varinfo[variates$T,'scale'] <- (varinfo[variates$T,'max']-varinfo[variates$T,'min'])/(1-2*ddd)
varinfo[variates$T,'location'] <- varinfo[variates$T,'min']-ddd*varinfo[variates$T,'scale']

##
varinfo[,'Q1'] <- apply(dt, 2, tquant, 0.25)
varinfo[,'Q2'] <- apply(dt, 2, tquant, 0.5)
varinfo[,'Q3'] <- apply(dt, 2, tquant, 0.75)



## write.csv(varinfo, 'varinfo.csv')
## varinfo <- data.matrix(read.csv('varinfo.csv', row.names=1))

summary(dt)

summary(sapply(colnames(dt),function(v){transfdir(data.matrix(dt[,..v]), v)}))

summary(sapply(colnames(dt),function(v){recjacobian(data.matrix(dt[,..v]), v)}))


sapply(colnames(dt),function(v){
    x <- data.matrix(dt[,..v])
    y <- transfinv(transfdir(x,v),v)
    rdiff <- y/x-1
    rdiff[x==0 & y==0] <- 0
    rdiff[x==0 & y!=0] <- Inf
    max(abs(rdiff),na.rm=T)
})

v <- 'GDTOTAL_gds'
varinfo[v,]

transfdir(0:6,v)
transfinv(0:6,v)
