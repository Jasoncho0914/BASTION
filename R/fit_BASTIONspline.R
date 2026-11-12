#' BASTION with a locally adaptive cubic spline trend
#'
#' Fits a locally adaptive cubic spline trend and shared seasonal templates
#' with smooth two-dimensional deviations. All observations enter the likelihood.
#'
#' Data are centered and divided by their sample SD. Priors are defined in
#' standardized units. Output trajectories,
#' templates and noise SDs are restored to data units; deviation variances are
#' multiplied by the squared SD.
#'
#' Template and spline coefficient second differences have horseshoe priors.
#' Separate spline boundary scales are learned. Each seasonal deviation variance
#' has an inverse-Gamma(delta_shape, delta_rate) prior. No knot-spacing rescaling
#' of coefficient differences is performed.
#'
#' In constant-noise mode the shared variance is sampled using the actual number
#' of spline and raw-template prior dimensions. Deviation energies are excluded
#' because they have their own absolute variance parameters. In SV mode sigma_e
#' is fixed at 1 on the standardized scale.
#'
#' Posterior summaries and pointwise credible intervals are returned whether
#' or not samples are retained. Partial final cycles have latent padding.
#'
#' @param y Finite numeric observation vector with positive sample SD.
#' @param Ks List of distinct integer seasonal periods, each >=3 and <length(y).
#' @param knot_spacing Required positive knot spacing in observation-index units.
#' @param Outlier Include horseshoe+ outliers at all observations; default FALSE.
#' @param obsSV Observation noise model: "const" or "SV"; default "const".
#' @param nchains Number of sequential chains; default 2.
#' @param nsave Number of saved draws per chain; default 1000.
#' @param nburn Burn-in sweeps per chain; default 1000.
#' @param nskip Discarded sweeps between saved draws; default 0 (no thinning).
#' @param cl Pointwise credible interval probability; default 0.95.
#' @param delta_shape Positive inverse-Gamma shape, shared as a hyperparameter
#'   across the independent seasonal deviation variances; default 0.001.
#' @param delta_rate Positive inverse-Gamma rate parameter; default 0.001.
#' @param init_evol_variance Starting evolution variance; default 0.01. It is
#'   passed explicitly to both initializers, not fixed during sampling.
#' @param verbose Print chain progress; default TRUE.
#' @param save_samples Return saved draws separately by chain; default FALSE.
#' @param seed Optional integer RNG seed. NULL leaves the caller's RNG stream.
#' @param gaussian_sampler Gaussian backend: 'auto' tries the precision
#'   sampler and falls back to square-root QR on failure; 'precision' or 'sqrt'
#'   selects that backend explicitly.
#' @return List with summary, model metadata, settings, normalization, diagnostics,
#'   elapsed_seconds and, when requested, samples$chains. Diagnostics use coda
#'   when installed. Samples include component trajectories, observation SDs,
#'   templates, spline coefficients and seasonal deviation variances.
#' @noRd
fit_BASTIONspline <- function(y,Ks,knot_spacing,Outlier=FALSE,obsSV="const",
                              nchains=2L,nsave=1000L,nburn=1000L,nskip=0L,cl=.95,
                              delta_shape=.001,delta_rate=.001,
                              init_evol_variance=.01,verbose=TRUE,
                              save_samples=FALSE,seed=NULL,
                              gaussian_sampler=c("auto","precision","sqrt")) {
  gaussian_sampler<-match.arg(gaussian_sampler)
  integer_arg<-function(x,name,minimum) {
    if(length(x)!=1L || !is.numeric(x) || !is.finite(x) || x!=as.integer(x) || x<minimum)
      stop(name," must be an integer >= ",minimum)
    as.integer(x)
  }
  if(!is.numeric(y) || !is.null(dim(y)) || length(y)<4L || any(!is.finite(y)))
    stop("y must be a finite numeric vector with at least four observations")
  n<-length(y)
  if(!is.list(Ks) || !length(Ks))stop("Ks must be a nonempty list")
  periods<-vapply(Ks,function(k)integer_arg(k,"Each seasonal period",3L),integer(1))
  if(any(periods>=n) || anyDuplicated(periods))stop("Seasonal periods must be distinct and < length(y)")
  knot_spacing<-.splinetrend_positive(knot_spacing,1L,"knot_spacing")
  nchains<-integer_arg(nchains,"nchains",1L);nsave<-integer_arg(nsave,"nsave",2L)
  nburn<-integer_arg(nburn,"nburn",0L);nskip<-integer_arg(nskip,"nskip",0L)
  delta_shape<-.season2_positive(delta_shape,1L,"delta_shape")
  delta_rate<-.season2_positive(delta_rate,1L,"delta_rate")
  init_evol_variance<-.season2_positive(init_evol_variance,1L,"init_evol_variance")
  if(length(cl)!=1L || !is.finite(cl) || cl<=0 || cl>=1)stop("cl must be between 0 and 1")
  if(length(obsSV)!=1L || !(obsSV%in%c("const","SV")))stop("obsSV must be 'const' or 'SV'")
  for(flag in c("Outlier","verbose","save_samples")) {
    value<-get(flag)
    if(!is.logical(value) || length(value)!=1L || is.na(value))stop(flag," must be TRUE or FALSE")
  }
  if(Outlier && n<5L)stop("Outlier model requires at least five observations")
  if(!is.null(seed))set.seed(integer_arg(seed,"seed",0L))
  offset<-mean(y);scale_y<-stats::sd(y)
  if(!is.finite(scale_y) || scale_y<=0)stop("y must have positive finite sample SD")
  original_y<-y;y<-(y-offset)/scale_y
  names_season<-paste0("Seasonal",periods)
  component_names<-c(names_season,"Trend",if(Outlier)"Outlier")
  p<-length(periods);trend_col<-p+1L
  nsteps<-nburn+(nskip+1)*nsave
  settings<-list(gaussian_sampler=gaussian_sampler,Ks=as.list(periods),knot_spacing=knot_spacing,Outlier=Outlier,obsSV=obsSV,
    nchains=nchains,nsave=nsave,nburn=nburn,nskip=nskip,cl=cl,
    delta_shape=delta_shape,delta_rate=delta_rate,init_evol_variance=init_evol_variance,seed=seed)
  diag_names<-c("mean_observation_SD","log_trend_HS_global","trend_middle",
    paste0("log_delta_variance_",periods),paste0("log_template_HS_global_",periods),
    paste0("seasonal_middle_",periods),if(obsSV=="SV")c("SV_mu","SV_phi","SV_sigma"))
  chains<-vector("list",nchains)
  start<-proc.time()[["elapsed"]]
  for(chain in seq_len(nchains)) {
    obs<-init_sigmaE_0(y)
    sp<-lapply(periods,function(k)init_Sbeta2(y,obs,k=k,delta_shape=delta_shape,
      delta_rate=delta_rate,init_evol_variance=init_evol_variance,sampler=gaussian_sampler))
    tp<-init_SplineTrend(y,obs,knot_spacing=knot_spacing,init_evol_variance=init_evol_variance,sampler=gaussian_sampler)
    components<-matrix(NA_real_,n,length(component_names),dimnames=list(NULL,component_names))
    for(j in seq_len(p))components[,j]<-as.numeric(sp[[j]]$s_mu)
    components[,trend_col]<-as.numeric(tp$s_mu)
    zp<-NULL
    if(Outlier){zp<-init_Outlier(y,obs);components[,"Outlier"]<-as.numeric(zp$s_mu)}
    residual<-y-rowSums(components)
    sv<-NULL
    if(obsSV=="SV") {
      sv<-dsp_initSV(residual);obs$sigma_e<-1;obs$sigma_et<-sv$sigma_wt
    }
    draws<-list(beta=array(NA_real_,c(nsave,n,ncol(components)),
      dimnames=list(NULL,NULL,component_names)),
      observation_SD=matrix(NA_real_,nsave,n),
      templates=setNames(lapply(periods,function(k)matrix(NA_real_,nsave,k)),names_season),
      delta_variance=matrix(NA_real_,nsave,p,dimnames=list(NULL,names_season)),
      spline_coefficients=matrix(NA_real_,nsave,tp$ncoef),
      diagnostics=matrix(NA_real_,nsave,length(diag_names),dimnames=list(NULL,diag_names)))
    isave<-0L
    for(sweep in seq_len(nsteps)) {
      tryCatch({
        if(obsSV=="SV") {
          sv<-dsp_sampleSVparams(residual,sv);obs$sigma_et<-sv$sigma_wt
        } else {
          rank<-tp$prior_rank+sum(vapply(sp,function(s)s$prior_rank,numeric(1)))
          energy<-tp$n_squared_sum+sum(vapply(sp,function(s)s$n_squared_sum,numeric(1)))
          if(Outlier){rank<-rank+zp$Td;energy<-energy+zp$n_squared_sum}
          obs$sigma_e<-sqrt(extraDistr::rinvgamma(1,(n+rank)/2,(sum(residual^2)+energy)/2))
          obs$sigma_et<-rep(obs$sigma_e,n)
        }
        for(j in seq_len(p)) {
          partial<-residual+components[,j]
          sp[[j]]<-fit_Sbeta2(partial,sp[[j]],obs)
          components[,j]<-as.numeric(sp[[j]]$s_mu)
          residual<-partial-components[,j]
        }
        partial<-residual+components[,trend_col]
        tp<-fit_SplineTrend(partial,tp,obs)
        components[,trend_col]<-as.numeric(tp$s_mu)
        residual<-partial-components[,trend_col]
        if(Outlier) {
          partial<-residual+components[,"Outlier"]
          zp<-fit_Outlier(partial,zp,obs)
          components[,"Outlier"]<-as.numeric(zp$s_mu)
          residual<-partial-components[,"Outlier"]
        }
      },error=function(e) {
        condition<-structure(list(message=sprintf("Chain %d, sweep %d: %s",chain,sweep,conditionMessage(e)),
          call=NULL,chain=chain,sweep=sweep,state=list(seasonal=sp,trend=tp,outlier=zp,sv=sv,obs=obs),
          partial_samples=draws),class=c("BASTIONspline_error","error","condition"))
        stop(condition)
      })
      if(sweep>nburn && (sweep-nburn)%%(nskip+1L)==0L) {
        isave<-isave+1L
        draws$beta[isave,,]<-components*scale_y
        draws$beta[isave,,"Trend"]<-draws$beta[isave,,"Trend"]+offset
        draws$observation_SD[isave,]<-as.numeric(obs$sigma_et)*scale_y
        for(j in seq_len(p)) {
          draws$templates[[j]][isave,]<-sp[[j]]$template*scale_y
          draws$delta_variance[isave,j]<-sp[[j]]$delta_variance*scale_y^2
        }
        draws$spline_coefficients[isave,]<-tp$beta*scale_y+offset
        middle<-ceiling(n/2)
        draws$diagnostics[isave,]<-c(mean(draws$observation_SD[isave,]),
          log(tp$s_evolParams$tauLambda),draws$beta[isave,middle,"Trend"],
          log(draws$delta_variance[isave,]),
          vapply(sp,function(s)log(s$s_evolParamsShape$tauLambda),numeric(1)),
          draws$beta[isave,middle,names_season],
          if(obsSV=="SV")as.numeric(sv$svParams[,1]))
      }
      if(verbose && (sweep%%max(1L,as.integer(nsteps/10L))==0L || sweep==nsteps))
        message(sprintf("Chain %d/%d, sweep %d/%d",chain,nchains,sweep,nsteps))
    }
    draws$last_state_standardized<-list(seasonal=sp,trend=tp,outlier=zp,sv=sv,obs=obs)
    chains[[chain]]<-draws
  }
  alpha<-(1-cl)/2
  summarize<-function(z)data.frame(Mean=colMeans(z),
    CR_lower=apply(z,2,stats::quantile,probs=alpha),
    CR_upper=apply(z,2,stats::quantile,probs=1-alpha))
  summary<-list();means<-matrix(NA_real_,n,length(component_names),dimnames=list(NULL,component_names))
  for(nm in component_names) {
    z<-do.call(rbind,lapply(chains,function(s)matrix(s$beta[,,nm],nrow=nsave)))
    sm<-summarize(z);summary[[paste0(nm,"_sum")]]<-sm;means[,nm]<-sm$Mean
  }
  summary$p_means<-data.frame(y=original_y,means,Remainder=original_y-rowSums(means),check.names=FALSE)
  summary$templates<-setNames(vector("list",p),names_season)
  for(j in seq_len(p)) {
    z<-do.call(rbind,lapply(chains,function(s)s$templates[[j]]))
    summary$templates[[j]]<-data.frame(phase=seq_len(periods[j]),summarize(z))
  }
  z<-do.call(rbind,lapply(chains,function(s)s$delta_variance))
  summary$delta_variance<-data.frame(component=names_season,summarize(z))
  summary$Volatility<-summarize(do.call(rbind,lapply(chains,function(s)s$observation_SD)))
  diagnostic<-NULL
  if(requireNamespace("coda",quietly=TRUE)) {
    mc<-coda::mcmc.list(lapply(chains,function(s)coda::mcmc(s$diagnostics)))
    rh<-rep(NA_real_,length(diag_names))
    if(nchains>1L)rh<-coda::gelman.diag(mc,autoburnin=FALSE,multivariate=FALSE)$psrf[,1]
    diagnostic<-data.frame(parameter=diag_names,Rhat=rh,effective_samples=coda::effectiveSize(mc))
  }
  result<-list(summary=summary,settings=settings,
    model=list(trend="cubic_spline_HS",seasonality="shared_HS_template_smooth_2D_deviation",
      spline_degree=3L,ncoefficients=tp$ncoef,knots=tp$knots,
      coefficient_prior_spacing_rescaled=FALSE),
    normalization=list(center=offset,scale=scale_y,priors_in="standardized data units"),
    diagnostics=diagnostic,elapsed_seconds=proc.time()[["elapsed"]]-start)
  if(save_samples)result$samples<-list(chains=chains)
  class(result)<-c("BASTIONspline","list")
  result
}
