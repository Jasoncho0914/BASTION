# Base-R regression checks for the public fixed-period fitted object.
library(BASTION)
set.seed(7)
n<-61L
values<-20+.1*seq_len(n)+sin(2*pi*seq_len(n)/10)+rnorm(n,sd=.3)
y<-ts(values,start=c(2000,1),frequency=12)
X<-cbind(covariate=cos(seq_len(n)*.37))
for(sv in c('const','SV'))for(reg in c(FALSE,TRUE))for(outlier in c(FALSE,TRUE)) {
  args<-list(y=y,Ks=list(10,15),X=if(reg)X else NULL,Outlier=outlier,
    obsSV=sv,nchains=2,nsave=4,nburn=4,nskip=2,verbose=FALSE,save_samples=TRUE,
    D_trend=if(reg)1L else 2L,cl=.8)
  set.seed(173);invisible(capture.output(fit<-do.call(fit_BASTION,args)))
  final_seed<-.Random.seed;D<-8L;M<-4L+reg+outlier
  stopifnot(inherits(fit,'bastion_fit'),is.data.frame(fit$components),
    identical(dim(fit$components),c(n,as.integer(M))),
    identical(dim(fit$draws$components),c(D,n,as.integer(M))),
    identical(names(fit$components),names(fit$intervals$lower)),
    identical(names(fit$components),names(fit$intervals$upper)),
    identical(names(fit$components),dimnames(fit$draws$components)[[3L]]),
    identical(fit$time,as.numeric(time(y))),
    identical(fit$draws$index$chain,rep(1:2,each=4)),
    identical(as.numeric(fit$draws$index$iteration),rep(c(7,10,13,16),2)),
    is.null(fit$settings),is.null(fit$summary),is.null(fit$samples),
    is.null(fit$components$remainder))
  for(j in seq_len(M)) {
    block<-matrix(fit$draws$components[,,j],D,n)
    stopifnot(isTRUE(all.equal(fit$components[[j]],colMeans(block),check.attributes=FALSE)),
      isTRUE(all.equal(fit$intervals$lower[[j]],apply(block,2,quantile,.1),check.attributes=FALSE)),
      isTRUE(all.equal(fit$intervals$upper[[j]],apply(block,2,quantile,.9),check.attributes=FALSE)))
  }
  for(k in c(10,15)) {
    name<-paste0('seasonal_',k);s<-fit$draws$seasonality[[name]]
    stopifnot(identical(dim(s$recurring),c(D,as.integer(k))),
      identical(dim(s$deviation),c(D,n)),
      identical(dim(fit$seasonality[[name]]$recurring),c(as.integer(k),3L)),
      identical(dim(fit$seasonality[[name]]$deviation),c(n,3L)),
      is.null(s$recurring_at_observations),
      identical(fit$metadata$mapping$recurring_positions[[name]],seq_len(k)))
    phase<-(seq_len(n)-1L)%%k+1L
    stopifnot(max(abs(s$recurring[,phase]+s$deviation-fit$draws$components[,,name]))<1e-10)
  }
  mu<-apply(fit$draws$components,c(1,2),sum)
  stopifnot(max(abs(mu-matrix(values,D,n,byrow=TRUE)))<1e-10,
    identical(dim(fit$draws$noise$sd),c(D,if(sv=='SV')n else 1L)),
    identical(dim(fit$noise$sd),c(if(sv=='SV')n else 1L,3L)),
    is.null(fit$noise$residual),is.null(fit$draws$noise$residual),
    identical(tail(names(fit$components),1),"residual"))
  if(sv=='SV')stopifnot(identical(dim(fit$draws$noise$parameters),c(D,3L)),
    identical(fit$noise$parameters$parameter,c('mu','phi','innovation_sd')),
    all(abs(fit$draws$noise$parameters[,'phi'])<1),
    all(fit$draws$noise$parameters[,'innovation_sd']>0))
  else stopifnot(is.null(fit$noise$parameters),is.null(fit$draws$noise$parameters))
  if(reg)stopifnot(max(abs(fit$draws$regression$coefficients%*%t(X)-fit$draws$components[,,'regression']))<1e-10)
  else stopifnot(is.null(fit$regression),is.null(fit$draws$regression))
  args$save_samples<-FALSE
  set.seed(173);invisible(capture.output(without<-do.call(fit_BASTION,args)))
  stopifnot(identical(.Random.seed,final_seed),is.null(without$draws))
  fit$draws<-NULL
  stopifnot(identical(fit,without))
  printed<-capture.output(print(fit));stopifnot(grepl('61 observations',printed[1]),
    all(c('- Seasonal_10:', '- Seasonal_15:') %in% printed),
    any(grepl('Recurring: Horseshoe on Second Differences',printed,fixed=TRUE)),
    any(grepl('... 56 more rows',printed,fixed=TRUE)))
}
# One saved draw preserves matrix/array dimensions; missing responses retain
# the fixed fitter's existing interpolation rather than introducing imputation.
y_missing<-values;y_missing[c(1,23,24,n)]<-NA_real_
set.seed(18)
invisible(capture.output(single<-fit_BASTION(y_missing,list(10),nchains=1,
  nsave=1,nburn=2,nskip=0,save_samples=TRUE,verbose=FALSE)))
completed<-approx(seq_len(n),y_missing,xout=seq_len(n),rule=2)$y
stopifnot(identical(dim(single$draws$components),c(1L,n,3L)),
  identical(dim(single$draws$seasonality$seasonal_10$recurring),c(1L,10L)),
  identical(single$observed,y_missing),
  max(abs(apply(single$draws$components,c(1,2),sum)-matrix(completed,1,n)))<1e-10)
cat('Fitted-object checks passed.\n')
