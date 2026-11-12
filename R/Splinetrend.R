# Cubic B-spline trend with locally adaptive HS coefficient differences.
#
# Model: residual y = B beta + noise, noise_t ~ N(0, obs_sigma_t2[t]).
# Priors: beta[1:2] have separate boundary variances; diff(beta,2) uses
# the dsp_package.R horseshoe hierarchy. Prior SDs are multiplied by sigma_e.
# The same dsp_* initial/boundary/evolution helpers are used as in Trend.R.
# The scale updates include squared-innovation numerical offsets.
# Gaussian updates use precision factorization or square-root QR.
#
# The curve is C2 continuous: this model permits local curvature adaptation,
# but not discontinuous level or slope jumps. It is not the daily-state trend.
# knot_spacing is REQUIRED, in units of time (observation positions by default).
# Coefficient differences are NOT rescaled automatically with knot spacing.
# Priors at different spacings therefore need a separate sensitivity assessment.
# Cubic degree 3 and repeated boundary knots (4 copies) define the spline basis.
# init_evol_variance=0.01 is an exposed starting variance, matching Trend.R's
# numeric initialization; it does not replace learned variance parameters.
# No data normalization, detrending, or calendar alignment occurs in this file.
#
# A full coefficient draw uses Q = B' W B + D' V^-1 D, where D has identity
# boundary rows followed by coefficient second differences. Q has half-bandwidth
# at most 3. Basis and Gram assembly cost O(T+J); factorization and draws O(J),
# with O(T+J) storage. B and knot locations are cached in the state.
#
# The outer sigma_e update
# must count prior_rank=J (number of coefficients), NOT Td (observations).
# n_squared_sum is the prior quadratic energy of beta using the current scales.
# The initial variances are on beta[1:2], not on trend values at times 1 and 2.
# Only HS is implemented here; other families are rejected explicitly.
#
# Usage after sourcing dsp_package.R and this file:
#   state <- init_SplineTrend(residual, obserror, knot_spacing=14)
#   state <- fit_SplineTrend(residual, state, obserror)
#   daily_trend_draw <- state$s_mu

.splinetrend_positive <- function(x,n,name) {
  if(!is.numeric(x) || !(length(x)%in%c(1L,n)) || any(!is.finite(x)) || any(x<=0)) {
    stop(name," must be finite and positive, with length 1 or ",n)
  }
  rep_len(as.numeric(x),n)
}

#' Sparse cubic spline design with equally spaced interior knots
#' @keywords internal
build_spline_trend_basis <- function(time,knot_spacing) {
  if(!is.numeric(time) || length(time)<4L || any(!is.finite(time)) || any(diff(time)<=0)) {
    stop("time must have at least four finite, strictly increasing coordinates")
  }
  spacing<-.splinetrend_positive(knot_spacing,1L,"knot_spacing")
  boundary<-range(time)
  interior<-seq(boundary[1],boundary[2],by=spacing)
  interior<-interior[interior>boundary[1] & interior<boundary[2]]
  knots<-c(rep(boundary[1],4L),interior,rep(boundary[2],4L))
  B<-splines::splineDesign(knots=knots,x=time,ord=4L,sparse=TRUE)
  list(B=B,knots=knots,interior_knots=interior,boundary_knots=boundary,
       knot_spacing=spacing,time=as.numeric(time),degree=3L)
}

# First two identity rows + ordinary second differences, constructed by bands.
.splinetrend_prior_precision <- function(evol_sigma_t2,J) {
  variance<-.splinetrend_positive(evol_sigma_t2,J,"evol_sigma_t2")
  precision<-1/variance
  if(any(!is.finite(precision))) stop("Nonfinite prior precision; no clipping applied")
  w<-precision[-c(1L,2L)]
  diagonal<-c(w,0,0)+c(0,4*w,0)+c(0,0,w)
  diagonal[1:2]<-diagonal[1:2]+precision[1:2]
  Matrix::bandSparse(J,k=c(0,1,2),
    diagonals=list(diagonal,c(-2*w,0)+c(0,-2*w),w),symmetric=TRUE)
}

#' Build the sparse spline coefficient posterior precision
#' @keywords internal
build_Q_spline_trend <- function(B,obs_sigma_t2,evol_sigma_t2,
                                 observed=rep(TRUE,nrow(B))) {
  if(!inherits(B,"sparseMatrix") || ncol(B)<4L) stop("B must be a sparse cubic spline basis")
  n<-nrow(B);J<-ncol(B)
  if(!is.logical(observed) || length(observed)!=n || anyNA(observed)) {
    stop("observed must be a logical vector of length nrow(B)")
  }
  vy<-.splinetrend_positive(obs_sigma_t2,n,"obs_sigma_t2")
  weights<-as.numeric(observed)/vy
  prior<-.splinetrend_prior_precision(evol_sigma_t2,J)
  Matrix::forceSymmetric(Matrix::crossprod(B,Matrix::Diagonal(x=weights)%*%B)+prior)
}

#' Joint Gaussian draw of all spline coefficients
#' @keywords internal
sampleSplineTrend_precision <- function(data,B,obs_sigma_t2,evol_sigma_t2) {
  if(!is.numeric(data) || length(data)!=nrow(B) || any(is.infinite(data))) {
    stop("data must have nrow(B) entries, each finite or NA")
  }
  observed<-!is.na(data)
  vy<-.splinetrend_positive(obs_sigma_t2,nrow(B),"obs_sigma_t2")
  Q<-build_Q_spline_trend(B,vy,evol_sigma_t2,observed)
  linear_data<-numeric(length(data))
  linear_data[observed]<-data[observed]/vy[observed]
  linear<-Matrix::crossprod(B,linear_data)
  U<-Matrix::chol(Q,pivot=FALSE)
  beta<-as.numeric(Matrix::solve(U,
    Matrix::solve(Matrix::t(U),linear)+stats::rnorm(ncol(B))))
  list(beta=beta,s_mu=as.matrix(B%*%beta))
}

.splinetrend_finish <- function(state) {
  innovations<-c(state$beta[1:2],diff(state$beta,differences=2))
  relative_sd<-c(state$s_evolParams0$sigma_w0,state$s_evolParams$sigma_wt)
  state$n_squared_sum<-sum((innovations/relative_sd)^2)
  state
}

#' Initialize the separate spline trend option
#' @keywords internal
init_SplineTrend <- function(data,obserror,knot_spacing,evol_error="HS",
                             time=seq_along(data),init_evol_variance=0.01,
                             sampler=c("auto","precision","sqrt")) {
  sampler<-match.arg(sampler)
  if(!identical(evol_error,"HS")) stop("Splinetrend currently supports evol_error='HS' only")
  if(length(time)!=length(data)) stop("time and data lengths differ")
  design<-build_spline_trend_basis(time,knot_spacing)
  sigma<-.splinetrend_positive(obserror$sigma_e,1L,"obserror$sigma_e")
  initial<-.splinetrend_positive(init_evol_variance,1L,"init_evol_variance")
  draw<-sampleSplineTrend(data,design$B,obserror$sigma_et^2,initial,sampler=sampler)
  # Separate initial coefficient variances with a global
  # auxiliary hierarchy. Initialization is explicit about commonSD=FALSE.
  p0<-dsp_initEvol0(draw$beta[1:2]/sigma,commonSD=FALSE)
  pe<-dsp_initEvolParams(diff(draw$beta,differences=2)/sigma,evol_error="HS")
  values<-c(p0$sigma_w0,pe$sigma_wt,pe$tauLambdaj,pe$xiLambdaj,pe$tauLambda,pe$xiLambda)
  if(any(!is.finite(values)) || any(values<=0)) stop("Horseshoe initialization produced invalid scales")
  state<-c(list(model="cubic_spline_hs_trend",Td=length(data),ncoef=ncol(design$B),
    gaussian_sampler=sampler,prior_rank=ncol(design$B),beta=draw$beta,s_mu=draw$s_mu,
    s_evolParams0=p0,s_evolParams=pe,evol_error=evol_error,colname="Trend"),design)
  .splinetrend_finish(state)
}

#' One coefficient and variance-parameter Gibbs sweep
#' @keywords internal
fit_SplineTrend <- function(data,tParam,obserror,evol_error=tParam$evol_error,
                             sampler=if(is.null(tParam$gaussian_sampler)) "auto" else tParam$gaussian_sampler) {
  sampler<-match.arg(sampler,c("auto","precision","sqrt"))
  tParam$gaussian_sampler<-sampler
  if(!identical(tParam$model,"cubic_spline_hs_trend")) stop("Initialize a spline trend state first")
  if(!identical(evol_error,"HS")) stop("Splinetrend currently supports evol_error='HS' only")
  if(length(data)!=tParam$Td) stop("Data length differs from initialized spline state")
  sigma<-.splinetrend_positive(obserror$sigma_e,1L,"obserror$sigma_e")
  relative_sd<-c(tParam$s_evolParams0$sigma_w0,tParam$s_evolParams$sigma_wt)
  variance<-sigma^2*relative_sd^2
  draw<-sampleSplineTrend(data,tParam$B,obserror$sigma_et^2,variance,sampler=sampler)
  tParam$beta<-draw$beta
  tParam$s_mu<-draw$s_mu
  tParam$s_evolParams0<-dsp_sampleEvol0(draw$beta[1:2]/sigma,
    tParam$s_evolParams0,commonSD=FALSE)
  tParam$s_evolParams<-dsp_sampleEvolParams(diff(draw$beta,differences=2)/sigma,
    tParam$s_evolParams,sigma_e=1,evol_error="HS")
  .splinetrend_finish(tParam)
}

# Same coefficient posterior, using the weighted operators rather than Q.
sampleSplineTrend_sqrt <- function(data,B,obs_sigma_t2,evol_sigma_t2,
                                    column_scale=TRUE,qr_order=0L) {
  if(!is.numeric(data)||length(data)!=nrow(B)||any(is.infinite(data)))stop("Invalid spline data")
  observed<-!is.na(data)
  vy<-.splinetrend_positive(obs_sigma_t2,nrow(B),"obs_sigma_t2")
  beta<-as.numeric(.gaussian_sqrt_draw(B[observed,,drop=FALSE],data[observed],
    .gaussian_difference_operator(ncol(B)),vy[observed],evol_sigma_t2,column_scale,qr_order))
  list(beta=beta,s_mu=as.matrix(B%*%beta))
}
sampleSplineTrend <- function(data,B,obs_sigma_t2,evol_sigma_t2,
                               sampler=c("auto","precision","sqrt"),notify=TRUE,
                               column_scale=TRUE,qr_order=0L) {
  .gaussian_sampler_select(
    function()sampleSplineTrend_precision(data,B,obs_sigma_t2,evol_sigma_t2),
    function()sampleSplineTrend_sqrt(data,B,obs_sigma_t2,evol_sigma_t2,column_scale,qr_order),
    sampler,notify,"Spline coefficients")
}
