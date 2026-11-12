library(BASTION)
set.seed(25)
n<-100L
calendar_days<-seq(as.Date('2020-01-01'),by='day',length.out=n)
dates<-as.POSIXct(format(calendar_days,'%Y-%m-%d'),tz='America/New_York')
y<-20+.02*seq_len(n)+sin(seq_len(n)*2*pi/7)+rnorm(n,sd=.3)
X<-cbind(covariate=cos(seq_len(n)*.17))
for(sv in c('const','SV'))for(reg in c(FALSE,TRUE)) {
  data<-y;data[30]<-NA_real_;use<-if(reg)seq_len(n) else setdiff(seq_len(n),45L)
  use<-c(tail(use,20),head(use,-20));T<-n
  args<-list(y=data[use],dates=dates[use],data_unit=list(unit='day'),
    seasonalities=list(hello=list(unit='week'),month=list(unit='month')),
    X=if(reg)X[use,,drop=FALSE] else NULL,Outlier=TRUE,obsSV=sv,
    nchains=2L,nsave=3L,nburn=3L,nskip=1L,verbose=FALSE,save_samples=TRUE,cl=.8)
  set.seed(51);fit<-do.call(fit_BASTION_calendar,args);seed<-.Random.seed;D<-6L;M<-5L+reg
  stopifnot(inherits(fit,'bastion_fit'),identical(dim(fit$components),c(T,as.integer(M))),
    identical(dim(fit$draws$components),c(D,T,as.integer(M))),
    identical(names(fit$components),names(fit$intervals$lower)),
    identical(names(fit$components),names(fit$intervals$upper)),
    identical(as.numeric(fit$time),as.numeric(dates)),identical(attr(fit$time,'tzone'),'America/New_York'),
    identical(fit$observed,replace(rep(NA_real_,n),use,args$y)),is.null(fit$settings),is.null(fit$samples),is.null(fit$summary))
  for(name in c('seasonal_hello','seasonal_month')) {
    d<-fit$draws$seasonality[[name]];mapping<-fit$metadata$mapping$seasonality[[name]]
    J<-nrow(fit$seasonality[[name]]$recurring)
    stopifnot(identical(dim(d$recurring),c(D,J)),identical(dim(d$deviation),c(D,T)),
      identical(dim(fit$seasonality[[name]]$recurring),c(J,3L)),
      all(fit$metadata$mapping$recurring_positions[[name]]>=0),
      all(fit$metadata$mapping$recurring_positions[[name]]<1),
      length(mapping$phase_index)==T,length(mapping$pattern_id)==T)
    recurring<-d$recurring[,mapping$phase_index,drop=FALSE]
    for(g in seq_along(mapping$patterns)) {
      ii<-which(mapping$pattern_id==g)
      recurring[,ii]<-recurring[,ii,drop=FALSE]-rowMeans(d$recurring[,mapping$patterns[[g]],drop=FALSE])
    }
    stopifnot(max(abs(recurring+d$deviation-fit$draws$components[,,name]))<1e-10)
  }
  completed<-apply(fit$draws$components,c(1,2),sum)
  known<-which(!is.na(fit$observed))
  stopifnot(max(abs(completed[,known]-matrix(fit$observed[known],D,length(known),byrow=TRUE)))<1e-10,
    length(fit$imputation)==if(reg)1L else 2L,
    identical(as.numeric(fit$imputation),as.numeric(fit$time[is.na(fit$observed)])),
    is.null(fit$calendar),is.null(fit$draws$imputation),
    is.null(fit$noise$residual),is.null(fit$draws$noise$residual),
    identical(tail(names(fit$components),1),'residual'),
    identical(dim(fit$draws$noise$sd),c(D,if(sv=='SV')T else 1L)))
  if(reg)stopifnot(max(abs(fit$draws$regression$coefficients%*%t(X)-fit$draws$components[,,'regression']))<1e-10)
  if(sv=='SV')stopifnot(identical(dim(fit$draws$noise$parameters),c(D,3L)))
  args$save_samples<-FALSE;set.seed(51);without<-do.call(fit_BASTION_calendar,args)
  stopifnot(identical(seed,.Random.seed),is.null(without$draws));fit$draws<-NULL
  stopifnot(identical(fit,without))
}
# Required dates must fail before initialization or fitting.
a<-tryCatch(fit_BASTION_calendar(y=1:20,data_unit=list(unit='day'),seasonalities=list(w=list(unit='week'))),error=identity)
b<-tryCatch(fit_BASTION_calendar(y=1:20,dates=NULL,data_unit=list(unit='day'),seasonalities=list(w=list(unit='week'))),error=identity)
stopifnot(inherits(a,'error'),inherits(b,'error'),grepl('dates must be supplied',conditionMessage(a)),grepl('dates must be supplied',conditionMessage(b)))
# Long periodic spline and native leap-year union grid.
n<-1100L;dates<-as.POSIXct(seq(as.Date('2019-01-01'),by='day',length.out=n),tz='UTC')
y<-20+sin(seq_len(n)*2*pi/365)+cos(seq_len(n)*.23)*.2
set.seed(16)
fit<-fit_BASTION_calendar(y,dates,list(unit='day'),list(yearly=list(unit='year')),
  nchains=1L,nsave=2L,nburn=2L,nskip=0L,verbose=FALSE,save_samples=TRUE)
stopifnot(nrow(fit$seasonality$seasonal_yearly$recurring)==731L,
  identical(dim(fit$draws$seasonality$seasonal_yearly$recurring),c(2L,731L)),
  length(fit$imputation)==0L,is.null(fit$draws$imputation))
cat('Calendar fitted-object checks passed.\n')
