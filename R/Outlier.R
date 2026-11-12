#' @keywords internal
#' @importFrom stats rnorm
sampleOutlier <- function(data, obs_sigma_t2, evol_sigma_t2, Td) {
  .outlier_positive(obs_sigma_t2,"observation variance")
  .outlier_positive(evol_sigma_t2,"outlier variance")
  # Harmonic variance avoids forming reciprocals of tiny variances.
  small<-pmin(obs_sigma_t2,evol_sigma_t2)
  large<-pmax(obs_sigma_t2,evol_sigma_t2)
  postVar<-small/(1+small/large)
  weight<-ifelse(evol_sigma_t2<=obs_sigma_t2,
    (evol_sigma_t2/obs_sigma_t2)/(1+evol_sigma_t2/obs_sigma_t2),
    1/(1+obs_sigma_t2/evol_sigma_t2))
  .outlier_positive(postVar,"conditional outlier variance")
  draw<-stats::rnorm(Td,data*weight,sqrt(postVar))
  if(any(!is.finite(draw)))stop("Nonfinite Gaussian outlier draw")
  draw
}

.outlier_positive <- function(x,label) {
  if(!length(x) || any(!is.finite(x) | x<=0))
    stop("Horseshoe+ outlier sampler: invalid ",label,call.=FALSE)
  invisible(x)
}

.outlier_ig <- function(n,shape,rate,label) {
  .outlier_positive(rate,paste(label,"inverse-gamma rate"))
  x<-extraDistr::rinvgamma(n,shape,rate)
  .outlier_positive(x,label)
  x
}

.outlier_scale <- function(state) {
  variance<-state$tau_t2*state$lambda_t2
  .outlier_positive(variance,"prior variance multiplier")
  state$sigma_wt<-sqrt(variance)
  state
}

# Horseshoe+: lambda | eta ~ C+(0, eta), eta ~ C+(0, 1), tau ~ C+(0, 1).
# IG(shape, rate) has density proportional to x^(-shape-1) exp(-rate/x).
#' @keywords internal
#' @import extraDistr
t_initEvolZeta_ps <- function(zeta,Td) {
  xi<-.outlier_ig(1,.5,1,"global auxiliary")
  tau_t2<-.outlier_ig(1,.5,1/xi,"global variance")
  eta_aux<-.outlier_ig(Td,.5,1,"extra local auxiliary")
  eta_t2<-.outlier_ig(Td,.5,1/eta_aux,"extra local variance")
  v<-.outlier_ig(Td,.5,1/eta_t2,"local auxiliary")
  lambda_t2<-.outlier_ig(Td,.5,1/v,"local variance")
  .outlier_scale(list(xi=xi,tau_t2=tau_t2,v=v,lambda_t2=lambda_t2,
    eta_aux=eta_aux,eta_t2=eta_t2))
}

#' @keywords internal
t_sampleEvolZeta_ps <- function(zeta,Td,evolParams) {
  p<-evolParams
  p$lambda_t2<-.outlier_ig(Td,1,1/p$v+as.numeric(zeta)^2/(2*p$tau_t2),"local variance")
  p$v<-.outlier_ig(Td,1,1/p$lambda_t2+1/p$eta_t2,"local auxiliary")
  p$eta_t2<-.outlier_ig(Td,1,1/p$eta_aux+1/p$v,"extra local variance")
  p$eta_aux<-.outlier_ig(Td,1,1+1/p$eta_t2,"extra local auxiliary")
  p$tau_t2<-.outlier_ig(1,(Td+1)/2,
    1/p$xi+sum(as.numeric(zeta)^2/p$lambda_t2)/2,"global variance")
  p$xi<-.outlier_ig(1,1,1+1/p$tau_t2,"global auxiliary")
  .outlier_scale(p)
}

#' @keywords internal
init_Outlier <- function(data,obserror) {
  Td<-length(data)
  s_mu<-sampleOutlier(data,obserror$sigma_et^2,rep(.01,Td),Td)
  params<-t_initEvolZeta_ps(s_mu/obserror$sigma_et,Td)
  energy<-sum((s_mu/params$sigma_wt)^2)
  if(!is.finite(energy))stop("Nonfinite horseshoe+ outlier prior energy")
  list(s_mu=s_mu,s_evolParams=params,Td=Td,n_squared_sum=energy,colname="Outlier")
}

#' @keywords internal
fit_Outlier <- function(data,zParam,obserror) {
  s_mu<-sampleOutlier(data,obserror$sigma_et^2,
    obserror$sigma_e^2*zParam$s_evolParams$sigma_wt^2,zParam$Td)
  params<-t_sampleEvolZeta_ps(s_mu/obserror$sigma_e,zParam$Td,zParam$s_evolParams)
  energy<-sum((s_mu/params$sigma_wt)^2)
  if(!is.finite(energy))stop("Nonfinite horseshoe+ outlier prior energy")
  zParam$s_mu<-s_mu
  zParam$s_evolParams<-params
  zParam$n_squared_sum<-energy
  zParam
}
