#' @keywords internal
#' @importFrom Matrix bandSparse
build_Q_season <- function(obs_sigma_t2 = NULL, evol_sigma_t2, Td, k) {
  prec_prior = 1 / evol_sigma_t2
  # initializing indices to reduce runtime
  i0_1 = 1:(k-3)
  i0_2 = 1:(Td-k)
  i0_3 = 2:(k-2)
  diag0 = prec_prior
  diag0[i0_1] = diag0[i0_1] + prec_prior[i0_1 + 2]
  diag0[i0_2] = diag0[i0_2] + prec_prior[i0_2 + k]
  diag0[i0_3] = diag0[i0_3] + 4*prec_prior[i0_3 + 1]

  diag1 = numeric(Td-1)
  diag1[i0_1] = diag1[i0_1] -2*prec_prior[i0_1+2]
  diag1[i0_3] = diag1[i0_3] -2*prec_prior[i0_3+1]

  diag2 = numeric(Td-2)
  diag2[i0_1] = diag2[i0_1] + prec_prior[i0_1 + 2]

  diagk = -prec_prior[-c(1:k)]

  # Add observation variance to main diagonal if provided
  if (!is.null(obs_sigma_t2)) {
    diag0 = diag0 + 1 / obs_sigma_t2
  }

  # Create sparse matrix using the bandSparse function
  Q = Matrix::bandSparse(
    Td,
    k = c(0, 1, 2, k),
    diagonals = list(diag0, diag1, diag2, diagk),
    symmetric = TRUE
  )
  row_indices <- c(rep(1:(k-1), each = (k-1)), 1:(k-1), rep(k, k-1))
  col_indices <- c(rep(1:(k-1), (k-1)), rep(k, k-1), 1:(k-1))
  values <- rep(prec_prior[k], length(row_indices))
  update <- Matrix::sparseMatrix(i = row_indices, j = col_indices, x = values,
                         dims = c(Td, Td))
  Q <- Q + update
  return(Q)
}
#' @keywords internal
#' @importFrom Matrix t chol  solve
#' @importFrom stats rnorm
sampleBeta_season <- function(data, obs_sigma_t2, evol_sigma_t2, Td, k) {
  QHt_Matrix = build_Q_season(obs_sigma_t2,
                              evol_sigma_t2,
                              Td,
                              k)
  linht = data / obs_sigma_t2
  chQht_Matrix <- robust_cholesky(QHt_Matrix)

  mu <- as.matrix(Matrix::solve(chQht_Matrix,
                                Matrix::solve(Matrix::t(chQht_Matrix),
                                              linht) +
                                  stats::rnorm(Td)))
  return(mu)
}
#' @keywords internal
init_Sbeta <- function(data, obserror, evol_error, k) {
  Td = length(data)
  obs_sigma_2T = obserror$sigma_et
  sigma_e = obserror$sigma_e
  s_mu = sampleBeta_season(
    data,
    obs_sigma_t2 = obs_sigma_2T ^ 2,
    evol_sigma_t2 =  sigma_e ^ 2 * rep(0.01, Td),
    Td = Td,
    k = k
  )

  # S 12
  s_evolParams23 = dsp_initEvol0(s_mu[1:2]/sigma_e,
                                 commonSD = TRUE)
  # S(3:km1)
  i1_km1 = c(1:(k-1))
  omega_2 = diff(s_mu[i1_km1], differences = 2)
  s_evolParams3k = dsp_initEvolParams(omega_2/sigma_e,
                                      evol_error = "HS")
  # S k+1,..,T
  omega_kT = matrix(c((s_mu[k] + sum(s_mu[i1_km1])),
                      diff(s_mu, lag = k))
  )
  s_evolParamskT = dsp_initEvolParams(omega_kT/sigma_e,
                                      evol_error = "HS")
  # s_evolParamskT = initEvolParams_HS_sparse(omega_kT/sigma_e,Td-k+1,
  #                                           1/(Td%%k))
  #
  #normalized squared sum
  n_squared_sum = sum((
    c(s_mu[1:2], omega_2, omega_kT) /
      c(
        s_evolParams23$sigma_w0,
        s_evolParams3k$sigma_wt,
        s_evolParamskT$sigma_wt
      )
  ) ^ 2)
  return(
    list(
      s_mu = s_mu,
      s_evolParams23 = s_evolParams23,
      s_evolParams3k = s_evolParams3k,
      s_evolParamskT = s_evolParamskT,
      Td = Td,
      n_squared_sum = n_squared_sum,
      colname = paste0("Seasonal", k)
    )
  )
}
#' @keywords internal
fit_Sbeta <- function(data, sParam, obserror, evol_error, k) {
  Td = sParam$Td
  obs_sigma_2T = obserror$sigma_et
  obs_sigma_e  = obserror$sigma_e
  s_mu = sampleBeta_season(
    data,
    obs_sigma_t2 = obs_sigma_2T ^ 2,
    evol_sigma_t2 = robust_prod(obs_sigma_e^2,
                                c(sParam$s_evolParams23$sigma_w0,
                                  sParam$s_evolParams3k$sigma_wt,
                                  sParam$s_evolParamskT$sigma_wt
                                )^2),
    Td = Td,
    k = k
  )
  s_evolParams23 = dsp_sampleEvol0(s_mu[1:2]/obs_sigma_e,
                                   sParam$s_evolParams23,
                                   commonSD = TRUE)

  # S(3:km1)
  i1_km1 = c(1:(k-1))
  omega_2 = diff(s_mu[i1_km1], differences = 2)
  s_evolParams3k = dsp_sampleEvolParams(
    omega_2 / obs_sigma_e,
    evolParams = sParam$s_evolParams3k,
    evol_error = "HS"
  )

  # S k+1,..,T
  omega_kT = matrix(c((s_mu[k] + sum(s_mu[i1_km1])),
                      diff(s_mu, lag = k)))
  s_evolParamskT = dsp_sampleEvolParams(
    omega_kT /obs_sigma_e,
    evolParams = sParam$s_evolParamskT,
    evol_error = "HS"
  )

  #normalized squared sum
  n_squared_sum = sum((
    c(s_mu[1:2], omega_2, omega_kT) /
      c(
        s_evolParams23$sigma_w0,
        s_evolParams3k$sigma_wt,
        s_evolParamskT$sigma_wt
      )
  ) ^ 2)

  sParam$s_mu = s_mu
  sParam$s_evolParams23 = s_evolParams23
  sParam$s_evolParams3k = s_evolParams3k
  sParam$s_evolParamskT = s_evolParamskT
  sParam$n_squared_sum = n_squared_sum
  return(sParam)
}


# ---- Shared seasonal helpers and smooth two-dimensional deviations ----
# Shared seasonal template + smooth two-dimensional deviations.
#
# Deviation backend: delta_sampler="spectral" (default) or
# "precision" (sparse Cholesky Gaussian block Gibbs draw). Both use the
# same conditional distribution and centering constraints. Set at init_Sbeta2
# or override in fit_Sbeta2; the default backend is spectral.
# The separate sampler argument still controls the TEMPLATE Gaussian backend.
# Precision cost depends on sparse factor fill; the O(T log T) discussion
# below describes the spectral backend, not the precision backend.
#
# s[m,j] = f[j] + delta[m,j],  m=1,...,M; j=1,...,k.
# sum(f)=0; rowSums(delta)=0; colSums(delta)=0.
# The phase direction is CIRCULAR; the cycle direction has OPEN boundaries.
# E(delta) = sum(diff across cycles ^2) + sum(circular phase differences ^2).
# p(delta | lambda) on the doubly centered subspace is proportional to
# lambda^(-r/2) exp(-E(delta)/(2*lambda)), r=(M-1)*(k-1).
# lambda ~ InvGamma(delta_shape, delta_rate), using the reciprocal-Gamma convention.
# lambda is an ABSOLUTE deviation variance, not multiplied by sigma_e^2.
# Default inverse-gamma hyperparameters: delta_shape=delta_rate=0.001.
# The default starting lambda is sigma_e^2 (initialization only, not a prior).
#
# Template prior: boundary scales and horseshoe innovations.
# raw template h[1:2] share an initial half-Cauchy SD; diff(h,2) has the
# dsp_package.R auxiliary horseshoe hierarchy. f=h-mean(h).
# Keeping the uncentered h as an auxiliary variable preserves conjugate scale
# updates: no overlapping-penalty determinant or truncation normalizer is omitted.
# All template prior variances scale with sigma_e^2.
# Squared-innovation offsets are not added.
# Template second differences are ordinary (not circular). Only delta wraps.
#
# One fit_Sbeta2 sweep:
#   1. Gaussian auxiliary-data draw for missing/padded/heteroskedastic observations;
#   2. raw template Gaussian draw, then f=h-mean(h);
#   3. template initial/local/global HS scale draws;
#   4. full delta Gaussian draw using DCT (cycles) and FFT (phases);
#   5. conjugate inverse-Gamma lambda draw.
# No Metropolis, clipping, jitter, or approximate iterative solves.
#
# Complete constant-variance data require no auxiliary random draws. For unequal
# observation variances, exact Gaussian noise-splitting produces a homoskedastic
# auxiliary grid. This is an additional Gibbs block, not a direct draw of the
# marginalized heterogeneous-observation conditional.
#
# A partial final cycle is padded internally as missing data. s_mu is truncated
# to Td; s_full, delta and template retain all k positions. Only FULL latent
# cycles sum to zero; an observed partial cycle need not sum to zero.
#
# Complexity per sweep: O(k + M*k*(log(M)+log(k))) time, O(k+M*k) storage.
# The template uses pentadiagonal Cholesky plus an exact rank-one Gaussian
# covariance correction. Neither template nor deviation updates construct a dense
# precision matrix. Transform initialization is O(M*k).
#
# IMPORTANT OUTER-FITTER CONTRACT:
# * prior_rank=k counts raw-template prior dimensions scaled by sigma_e^2.
# * n_squared_sum contains ONLY raw-template energy; delta has its own lambda.
# * The outer sigma_e update must use observed likelihood + prior_rank, NOT Td
#   prior dimensions or padded auxiliary observations. Resample auxiliary data
#   inside this function after any outer sigma_e update; never reuse it outside.
# * This model has separate template and deviation scales, not a single sigma_wt.
#
# Requires Matrix, and dsp_initEvol0/dsp_initEvolParams from dsp_package.R.

.season2_positive <- function(x, n, name) {
  if (!is.numeric(x) || !(length(x) %in% c(1L,n)) ||
      any(!is.finite(x)) || any(x<=0)) {
    stop(name, " must be finite and positive, with length 1 or ", n)
  }
  rep_len(as.numeric(x),n)
}

.season2_dimensions <- function(Td,k) {
  if (length(Td)!=1L || !is.finite(Td) || Td<1L || Td!=as.integer(Td) ||
      length(k)!=1L || !is.finite(k) || k<3L || k!=as.integer(k)) {
    stop("Require integer Td >= 1 and integer k >= 3")
  }
  c(Td=as.integer(Td), k=as.integer(k), M=as.integer(ceiling(Td/k)))
}

.season2_data <- function(data,Td) {
  if (!is.numeric(data) || length(data)!=Td || any(is.infinite(data))) {
    stop("data must be numeric, length Td, and contain only finite values or NA")
  }
  as.numeric(data)
}

.season2_geometry <- function(M,k) {
  # Laplacian eigenvalues for a path of M cycles and a ring of k phases.
  list(M=M,k=k,
    eigenvalues=outer(4*sin(pi*(0:(M-1))/(2*M))^2,
                     4*sin(pi*(0:(k-1))/k)^2,"+"),
    dct_phase=exp(-1i*pi*(0:(M-1))/(2*M)),
    dct_scale=c(1/sqrt(M),rep(sqrt(2/M),M-1L)))
}

# Orthonormal DCT-II along matrix columns, implemented with base R FFTs.
.season2_dct <- function(x,geometry,inverse=FALSE) {
  n<-geometry$M
  if (n==1L) return(x)
  if (!inverse) {
    z<-stats::mvfft(rbind(x,x[n:1,,drop=FALSE]))[seq_len(n),,drop=FALSE]
    z<-Re(z*geometry$dct_phase)/2
    z*geometry$dct_scale
  } else {
    z<-2*(x/geometry$dct_scale)/geometry$dct_phase
    extended<-rbind(z,matrix(0+0i,1,ncol(x)),Conj(z[n:2,,drop=FALSE]))
    Re(stats::mvfft(extended,inverse=TRUE)[seq_len(n),,drop=FALSE])/(2*n)
  }
}

.season2_transform <- function(x,geometry,inverse=FALSE) {
  if (!inverse) {
    z<-.season2_dct(x,geometry)
    t(stats::mvfft(t(z)))/sqrt(geometry$k)
  } else {
    z<-Re(t(stats::mvfft(t(x),inverse=TRUE))/sqrt(geometry$k))
    .season2_dct(z,geometry,inverse=TRUE)
  }
}

.season2_energy <- function(delta) {
  M<-nrow(delta); k<-ncol(delta)
  across_cycles<-if(M>1L) sum((delta[2:M,,drop=FALSE]-delta[1:(M-1),,drop=FALSE])^2) else 0
  across_cycles+sum((delta[,c(2:k,1L),drop=FALSE]-delta)^2)
}

# Full constrained Gaussian draw, diagonal in the grid spectral basis.
.season2_draw_delta <- function(residual,obs_variance,lambda,geometry) {
  if (geometry$M==1L) return(matrix(0,1L,geometry$k))
  precision<-1/obs_variance+geometry$eigenvalues/lambda
  if (any(!is.finite(precision)) || any(precision<=0)) {
    stop("Invalid deviation precision; no clipping applied")
  }
  variance<-1/precision
  coefficients<-.season2_transform(residual,geometry)*variance/obs_variance+
    .season2_transform(matrix(stats::rnorm(length(residual)),nrow(residual),ncol(residual)),
                      geometry)*sqrt(variance)
  # Remove phase-constant AND cycle-constant modes exactly in spectral space.
  coefficients[1L,]<-0
  coefficients[,1L]<-0
  .season2_transform(coefficients,geometry,inverse=TRUE)
}

# Exact homoskedastic auxiliary data. For an observed variance v>=c,
# z|y,s ~ N(s+c/v*(y-s), c*(v-c)/v). Missing z|s ~ N(s,c).
.season2_auxiliary <- function(data,obs_variance,signal,M,k) {
  Td<-length(data); N<-M*k
  v<-.season2_positive(obs_variance,Td,"obs_sigma_t2")
  base_variance<-min(v)
  full<-as.numeric(t(signal))
  mean<-full
  variance<-rep(base_variance,N)
  seen<-which(!is.na(data))
  mean[seen]<-full[seen]+base_variance/v[seen]*(data[seen]-full[seen])
  variance[seen]<-base_variance*(v[seen]-base_variance)/v[seen]
  z<-mean
  uncertain<-which(variance>0)
  if(length(uncertain)) {
    z[uncertain]<-stats::rnorm(length(uncertain),mean[uncertain],sqrt(variance[uncertain]))
  }
  list(data=matrix(z,M,k,byrow=TRUE),variance=base_variance)
}

.season2_template_precision <- function(k,initial_precision,evolution_precision) {
  q<-.season2_positive(initial_precision,1L,"initial_precision")
  w<-.season2_positive(evolution_precision,k-2L,"evolution_precision")
  diagonal<-c(w,0,0)+c(0,4*w,0)+c(0,0,w)
  diagonal[1:2]<-diagonal[1:2]+q
  off1<-c(-2*w,0)+c(0,-2*w)
  Matrix::bandSparse(k,k=c(0,1,2),diagonals=list(diagonal,off1,w),symmetric=TRUE)
}

.season2_draw_template_precision <- function(residual,obs_variance,sigma_e,
                                   initial_precision,evolution_precision) {
  k<-ncol(residual); M<-nrow(residual)
  a<-M/obs_variance
  Q0<-.season2_template_precision(k,initial_precision,evolution_precision)
  # Posterior Q = A - alpha*11', where A is pentadiagonal and alpha=a/k.
  A<-Q0/sigma_e^2+Matrix::Diagonal(k,x=rep(a,k))
  linear<-colSums(residual)/obs_variance
  linear<-linear-mean(linear)
  U<-Matrix::chol(A,pivot=FALSE) # A=U' U; no fill-reducing permutation needed.
  # Batch the two RHS while retaining O(k) work and storage.
  solved<-as.matrix(Matrix::solve(U,Matrix::solve(Matrix::t(U),
    cbind(linear,rep(1,k)))))
  m0<-solved[,1L]
  v<-solved[,2L] # A^-1 1
  alpha<-a/k
  # Stable identity for d=1-alpha*sum(v): Q0*1 has only its first two
  # entries nonzero. This avoids subtracting two nearly equal numbers.
  first_sum<-sum(v[1:2])
  d<-(initial_precision/sigma_e^2)*(first_sum/k)
  if(!is.finite(d) || d<=0 || !is.finite(first_sum) || first_sum<=0) {
    stop("Nonpositive rank-one denominator; no jitter or clipping applied")
  }
  # Equivalent to m0 + v*alpha*sum(m0)/d, using sum(linear)=0 and
  # Q0*1=(initial_precision)*(1,1,0,...). Avoid cancellation in sum(m0).
  mean_raw<-m0-v*(sum(m0[1:2])/first_sum)
  # Q^-1 = A^-1 + (alpha/d)*v*v'. Independent Gaussian contributions.
  raw<-mean_raw+as.numeric(Matrix::solve(U,stats::rnorm(k)))+
    v*sqrt(alpha/d)*stats::rnorm(1)
  list(raw=raw,template=raw-mean(raw))
}

# Same initial-value half-Cauchy hierarchy as commonSD=TRUE in dsp_package.R,
# without offsets added to squared innovations.
.season2_update_initial <- function(raw,sigma_e,params) {
  rho<-stats::rgamma(1,shape=1.5,
    rate=params$px_sigma_w0[1]+sum((raw[1:2]/sigma_e)^2)/2)
  px<-stats::rgamma(1,shape=1,rate=rho+1)
  params$sigma_w0<-rep(1/sqrt(rho),2L)
  params$px_sigma_w0<-rep(px,2L)
  params
}

# Horseshoe precision/auxiliary hierarchy with conjugate innovation conditionals.
.season2_update_hs <- function(omega,params) {
  n<-length(omega)
  q<-stats::rgamma(n,shape=1,rate=as.numeric(params$xiLambdaj)+omega^2/2)
  xi<-stats::rgamma(n,shape=1,rate=q+params$tauLambda)
  global<-stats::rgamma(1,shape=(n+1)/2,rate=sum(xi)+params$xiLambda)
  global_aux<-stats::rgamma(1,shape=1,rate=global+1)
  values<-c(q,xi,global,global_aux)
  if(any(!is.finite(values)) || any(values<=0)) {
    stop("Nonfinite/zero HS draw; no clipping or innovation offsets applied")
  }
  params$tauLambdaj<-matrix(q,ncol=1L)
  params$xiLambdaj<-matrix(xi,ncol=1L)
  params$tauLambda<-global
  params$xiLambda<-global_aux
  params$sigma_wt<-matrix(1/sqrt(q),ncol=1L)
  params
}

.season2_finish <- function(state) {
  full<-sweep(state$delta,2L,state$template,"+")
  state$s_full<-full
  state$s_mu<-matrix(as.numeric(t(full))[seq_len(state$Td)],ncol=1L)
  raw<-state$template_raw
  state$n_squared_sum<-sum((raw[1:2]/state$s_evolParams0$sigma_w0)^2)+
    sum(diff(raw,differences=2)^2/as.numeric(state$s_evolParamsShape$sigma_wt)^2)
  state$delta_energy<-.season2_energy(state$delta)
  state
}

#' Initialize common locally adaptive template + smooth cycle deviations
#' @keywords internal
init_Sbeta2 <- function(data,obserror,evol_error="HS",k,
                        delta_shape=0.001,delta_rate=0.001,
                        delta_variance_init=NULL,init_evol_variance=0.01,
                        sampler=c("auto","precision","sqrt"),
                        delta_sampler=c("spectral","precision")) {
  sampler<-match.arg(sampler)
  delta_sampler<-match.arg(delta_sampler)
  if(!identical(evol_error,"HS")) stop("Template currently supports evol_error='HS' only")
  dims<-.season2_dimensions(length(data),k)
  Td<-dims[["Td"]]; M<-dims[["M"]]; k<-dims[["k"]]
  data<-.season2_data(data,Td)
  sigma<-.season2_positive(obserror$sigma_e,1L,"obserror$sigma_e")
  shape<-.season2_positive(delta_shape,1L,"delta_shape")
  rate<-.season2_positive(delta_rate,1L,"delta_rate")
  initial<-.season2_positive(init_evol_variance,1L,"init_evol_variance")
  lambda<-if(is.null(delta_variance_init)) sigma^2 else
    .season2_positive(delta_variance_init,1L,"delta_variance_init")
  geometry<-.season2_geometry(M,k)
  if(delta_sampler=="precision") geometry$laplacian<-.season2_grid_laplacian(M,k)
  delta<-matrix(0,M,k)
  auxiliary<-.season2_auxiliary(data,obserror$sigma_et^2,delta,M,k)
  initial_draw<-.season2_draw_template(auxiliary$data,auxiliary$variance,sigma,
    1/initial,rep(1/initial,k-2L),sampler=sampler)
  p0<-dsp_initEvol0(initial_draw$raw[1:2]/sigma,commonSD=TRUE)
  ps<-dsp_initEvolParams(diff(initial_draw$raw,differences=2)/sigma,evol_error="HS")
  values<-c(p0$sigma_w0,p0$px_sigma_w0,ps$tauLambdaj,ps$xiLambdaj,ps$tauLambda,ps$xiLambda)
  if(any(!is.finite(values)) || any(values<=0)) {
    stop("Horseshoe initialization produced invalid scales")
  }
  state<-list(delta_sampler=delta_sampler,gaussian_sampler=sampler,model="shared_template_2d_deviation",Td=Td,k=k,M=M,
    template_raw=initial_draw$raw,template=initial_draw$template,delta=delta,
    delta_variance=lambda,delta_shape=shape,delta_rate=rate,
    delta_rank=(M-1L)*(k-1L),geometry=geometry,
    s_evolParams0=p0,s_evolParamsShape=ps,
    # Raw k-dimensional template prior scales with sigma_e; delta prior does not.
    prior_rank=k,colname=paste0("Seasonal",k))
  .season2_finish(state)
}

#' One full Gibbs sweep, including template HS and deviation variance updates
#' @keywords internal
fit_Sbeta2 <- function(data,sParam,obserror,evol_error="HS",k=sParam$k,
                       sampler=if(is.null(sParam$gaussian_sampler)) "auto" else sParam$gaussian_sampler,
                       delta_sampler=if(is.null(sParam$delta_sampler)) "spectral" else sParam$delta_sampler) {
  sampler<-match.arg(sampler,c("auto","precision","sqrt"))
  sParam$gaussian_sampler<-sampler
  delta_sampler<-match.arg(delta_sampler,c("spectral","precision"))
  sParam$delta_sampler<-delta_sampler
  if(delta_sampler=="precision" && is.null(sParam$geometry$laplacian))
    sParam$geometry$laplacian<-.season2_grid_laplacian(sParam$M,sParam$k)
  if(!identical(sParam$model,"shared_template_2d_deviation")) {
    stop("Seasonal state is incompatible with the smooth-deviation sampler; initialize a matching state")
  }
  if(!identical(evol_error,"HS")) stop("Template currently supports evol_error='HS' only")
  if(length(data)!=sParam$Td || k!=sParam$k) stop("State/data dimensions or k changed")
  data<-.season2_data(data,sParam$Td)
  sigma<-.season2_positive(obserror$sigma_e,1L,"obserror$sigma_e")
  auxiliary<-.season2_auxiliary(data,obserror$sigma_et^2,sParam$s_full,sParam$M,k)
  template<-.season2_draw_template(auxiliary$data-sParam$delta,auxiliary$variance,sigma,
    1/sParam$s_evolParams0$sigma_w0[1]^2,
    as.numeric(sParam$s_evolParamsShape$tauLambdaj),sampler=sampler)
  sParam$template_raw<-template$raw
  sParam$template<-template$template
  sParam$s_evolParams0<-.season2_update_initial(template$raw,sigma,sParam$s_evolParams0)
  sParam$s_evolParamsShape<-.season2_update_hs(diff(template$raw,differences=2)/sigma,
    sParam$s_evolParamsShape)
  residual<-sweep(auxiliary$data,2L,sParam$template,"-")
  # Both backends draw the same full Gaussian conditional.
  delta_draw<-if(delta_sampler=="spectral") .season2_draw_delta else .season2_draw_delta_precision
  sParam$delta<-delta_draw(residual,auxiliary$variance,
    sParam$delta_variance,sParam$geometry)
  energy<-.season2_energy(sParam$delta)
  precision<-stats::rgamma(1,shape=sParam$delta_shape+sParam$delta_rank/2,
    rate=sParam$delta_rate+energy/2)
  lambda<-1/precision
  if(!is.finite(lambda) || lambda<=0) {
    stop("Nonfinite/zero deviation variance draw; no clipping applied")
  }
  sParam$delta_variance<-lambda
  .season2_finish(sParam)
}

# The likelihood acts on P h, while the proper boundary/difference prior acts
# on raw h. Centering h by conditioning would change this target distribution.
.season2_draw_template_sqrt <- function(residual,obs_variance,sigma_e,
    initial_precision,evolution_precision,column_scale=TRUE,qr_order=0L) {
  k<-ncol(residual);M<-nrow(residual)
  initial_precision<-.gaussian_positive(initial_precision,1L,"initial_precision")
  evolution_precision<-.gaussian_positive(evolution_precision,k-2L,"evolution_precision")
  sigma_e<-.gaussian_positive(sigma_e,1L,"sigma_e")
  obs_variance<-.gaussian_positive(obs_variance,1L,"obs_variance")
  # k-by-k centering operator: no T-by-T dense matrix is constructed.
  P<-diag(k)-matrix(1/k,k,k)
  raw<-as.numeric(.gaussian_sqrt_draw(P,colMeans(residual),
    .gaussian_difference_operator(k),obs_variance/M,
    sigma_e^2*c(rep(1/initial_precision,2L),1/evolution_precision),column_scale,qr_order))
  list(raw=raw,template=raw-mean(raw))
}
.season2_draw_template <- function(residual,obs_variance,sigma_e,
    initial_precision,evolution_precision,sampler=c("auto","precision","sqrt"),
    notify=TRUE,column_scale=TRUE,qr_order=0L) {
  .gaussian_sampler_select(
    function().season2_draw_template_precision(residual,obs_variance,sigma_e,initial_precision,evolution_precision),
    function().season2_draw_template_sqrt(residual,obs_variance,sigma_e,initial_precision,evolution_precision,column_scale,qr_order),
    sampler,notify,"Common seasonal template")
}

# Sparse path-by-ring Laplacian in R's column-major vectorization.
# Its quadratic form equals .season2_energy(delta).
.season2_grid_laplacian <- function(M,k) {
  if(length(M)!=1L || !is.finite(M) || M<1L || M!=as.integer(M) ||
     length(k)!=1L || !is.finite(k) || k<3L || k!=as.integer(k))
    stop("Require integer M>=1 and k>=3")
  cycle<-if(M==1L) Matrix::Matrix(0,1L,1L,sparse=TRUE) else {
    D<-Matrix::sparseMatrix(i=rep(seq_len(M-1L),2L),
      j=c(seq_len(M-1L),2:M),x=c(rep(-1,M-1L),rep(1,M-1L)),dims=c(M-1L,M))
    Matrix::crossprod(D)
  }
  D<-Matrix::sparseMatrix(i=rep(seq_len(k),2L),
    j=c(seq_len(k),c(2:k,1L)),x=c(rep(-1,k),rep(1,k)),dims=c(k,k))
  phase<-Matrix::crossprod(D)
  Matrix::forceSymmetric(kronecker(Matrix::Diagonal(k),cycle)+
    kronecker(phase,Matrix::Diagonal(M)))
}

# Exact precision alternative to .season2_draw_delta(), with no transforms.
# Q=I/v + grid_laplacian/lambda; b=vec(residual)/v.
# Draw x~N(Q^-1 b,Q^-1), then double-center. This is exact here because
# Q commutes with the centering projection: both Laplacians annihilate
# constants. This is NOT valid for arbitrary heterogeneous likelihoods.
# Unequal/missing observations use auxiliary Gibbs updates.
# No jitter, clipping, or diagonal repair is used.
.season2_draw_delta_precision <- function(residual,obs_variance,lambda,geometry) {
  M<-geometry$M;k<-geometry$k
  if(!is.matrix(residual) || !identical(dim(residual),c(as.integer(M),as.integer(k))) ||
     any(!is.finite(residual)))stop("Invalid residual dimensions or values")
  v<-.season2_positive(obs_variance,1L,"obs_variance")
  lambda<-.season2_positive(lambda,1L,"lambda")
  if(M==1L)return(matrix(0,1L,k))
  laplacian<-geometry$laplacian
  if(is.null(laplacian))laplacian<-.season2_grid_laplacian(M,k)
  Q<-Matrix::Diagonal(M*k,x=rep(1/v,M*k))+laplacian/lambda
  U<-Matrix::chol(Q,pivot=FALSE)
  x<-as.numeric(Matrix::solve(U,
    Matrix::solve(Matrix::t(U),as.numeric(residual)/v)+stats::rnorm(M*k)))
  delta<-matrix(x,M,k)
  delta<-sweep(delta,1L,rowMeans(delta),"-")
  delta<-sweep(delta,2L,colMeans(delta),"-")
  if(any(!is.finite(delta)))stop("Nonfinite precision deviation draw")
  delta
}


# ---- Shared templates and group-horseshoe seasonal differences ----
# Shared recurring template and locally adaptive seasonal differences.
# eta[m,] = delta[m+1,] - delta[m,], with one local scale per transition.
# eta[m,] is isotropic on the k-1 dimensional zero-phase-sum subspace:
# eta[m,] | tau2,lambda2[m] ~ N_{1-perp}(0,tau2*lambda2[m]*I).
# tau and each lambda are independent half-Cauchy(0,1), in standardized
# fitter units and independent of the observation scale. There is no
# circular within-cycle penalty and no extra delta initial-value prior.
# Row/column sums of delta are zero. The across-cycle mean constraint fixes
# delta[1,] = -sum_{m=1}^{M-1}(M-m)*eta[m,]/M.
# The constrained innovation map has a scale-independent Jacobian, so each
# transition contributes k-1 Gaussian dimensions to the HS updates.
#
# Under the homoskedastic auxiliary construction, the delta
# precision is c^-1 I + D_cycle' diag(1/(tau2*lambda2)) D_cycle.
# Factor once, solve k right-hand sides. This is O(M*k) time/storage.
# Both centering projections commute with this precision (Q*1=c^-1*1),
# so ordinary row/column centering gives the exact constrained Gaussian.
# Precision warnings/errors/nonfinite draws use the square-root
# QR machinery with the identical likelihood and difference operator.
# No jitter, clipping, MH, circular penalty, or iterative approximation.

.season3_innovations <- function(delta) {
  M<-nrow(delta)
  if(M==1L)return(matrix(numeric(),0L,ncol(delta)))
  delta[2:M,,drop=FALSE]-delta[seq_len(M-1L),,drop=FALSE]
}

.season3_draw_delta <- function(residual,obs_variance,innovation_variance,
                                sampler=c("auto","precision","sqrt")) {
  sampler<-match.arg(sampler)
  M<-nrow(residual);k<-ncol(residual)
  if(M==1L)return(structure(matrix(0,1L,k),sampler_backend="zero_rank"))
  c0<-.season2_positive(obs_variance,1L,"auxiliary observation variance")
  v<-.season2_positive(innovation_variance,M-1L,"seasonal innovation variance")
  draw<-.gaussian_sampler_select(
    function(){
      w<-1/v
      Q<-Matrix::bandSparse(M,k=c(0,1),
        diagonals=list(rep(1/c0,M)+c(w,0)+c(0,w),-w),symmetric=TRUE)
      U<-Matrix::chol(Q,pivot=FALSE)
      as.matrix(Matrix::solve(U,Matrix::solve(Matrix::t(U),residual/c0)+
        matrix(stats::rnorm(M*k),M,k)))
    },
    function(){
      D<-Matrix::sparseMatrix(i=rep(seq_len(M-1L),each=2L),
        j=as.numeric(t(cbind(seq_len(M-1L),2:M))),
        x=rep(c(-1,1),M-1L),dims=c(M-1L,M))
      vapply(seq_len(k),function(j)as.numeric(.gaussian_sqrt_draw(
        Matrix::Diagonal(M),residual[,j],D,c0,v)),numeric(M))
    },sampler=sampler,block="Group-horseshoe seasonal deviations")
  backend<-attr(draw,"sampler_backend")
  draw<-sweep(draw,1L,rowMeans(draw),"-")
  draw<-sweep(draw,2L,colMeans(draw),"-")
  attr(draw,"sampler_backend")<-backend
  draw
}

.season3_update_scales <- function(delta,params) {
  eta<-.season3_innovations(delta)
  n<-nrow(eta);d<-ncol(delta)-1L
  energy<-rowSums(eta^2)
  if(n) {
    params$lambda2<-1/stats::rgamma(n,shape=(d+1)/2,
      rate=1/params$nu+energy/(2*params$tau2))
    params$nu<-1/stats::rgamma(n,shape=1,rate=1+1/params$lambda2)
  }
  params$tau2<-1/stats::rgamma(1,shape=(n*d+1)/2,
    rate=1/params$xi+sum(energy/params$lambda2)/2)
  params$xi<-1/stats::rgamma(1,shape=1,rate=1+1/params$tau2)
  values<-c(params$tau2,params$lambda2,params$nu,params$xi)
  if(any(!is.finite(values)) || any(values<=0))
    stop("Invalid group-horseshoe seasonal-deviation scale draw")
  params
}

.season3_finish <- function(state) {
  # Initialize common template and output fields.
  state<-.season2_finish(state)
  state$delta_energy<-sum(.season3_innovations(state$delta)^2)
  state$delta_innovation_variance<-state$delta_hs$tau2*state$delta_hs$lambda2
  # There is no single shared delta variance in this model.
  state$delta_variance<-NA_real_
  state
}

init_Sbeta3 <- function(data,obserror,evol_error="HS",k,
                        init_evol_variance=0.01,
                        sampler=c("auto","precision","sqrt"),
                        delta_sampler=c("auto","precision","sqrt")) {
  sampler<-match.arg(sampler);delta_sampler<-match.arg(delta_sampler)
  if(!identical(evol_error,"HS"))stop("Template currently supports evol_error='HS' only")
  dims<-.season2_dimensions(length(data),k)
  Td<-dims[["Td"]];M<-dims[["M"]];k<-dims[["k"]]
  data<-.season2_data(data,Td)
  sigma<-.season2_positive(obserror$sigma_e,1L,"obserror$sigma_e")
  initial<-.season2_positive(init_evol_variance,1L,"init_evol_variance")
  delta<-matrix(0,M,k)
  auxiliary<-.season2_auxiliary(data,obserror$sigma_et^2,delta,M,k)
  initial_draw<-.season2_draw_template(auxiliary$data,auxiliary$variance,sigma,
    1/initial,rep(1/initial,k-2L),sampler=sampler)
  p0<-dsp_initEvol0(initial_draw$raw[1:2]/sigma,commonSD=TRUE)
  ps<-dsp_initEvolParams(diff(initial_draw$raw,differences=2)/sigma,evol_error="HS")
  values<-c(p0$sigma_w0,p0$px_sigma_w0,ps$tauLambdaj,ps$xiLambdaj,
    ps$tauLambda,ps$xiLambda)
  if(any(!is.finite(values)) || any(values<=0))
    stop("Horseshoe initialization produced invalid scales")
  state<-list(model="shared_template_group_HS_seasonal_differences",Td=Td,k=k,M=M,
    template_raw=initial_draw$raw,template=initial_draw$template,delta=delta,
    delta_rank=(M-1L)*(k-1L),
    delta_hs=list(tau2=1,lambda2=rep(1,M-1L),nu=rep(1,M-1L),xi=1),
    delta_sampler=delta_sampler,delta_sampler_backend="initial_zero",
    gaussian_sampler=sampler,s_evolParams0=p0,s_evolParamsShape=ps,
    prior_rank=k,colname=paste0("Seasonal",k))
  .season3_finish(state)
}

fit_Sbeta3 <- function(data,sParam,obserror,evol_error="HS",k=sParam$k,
                       sampler=sParam$gaussian_sampler,
                       delta_sampler=sParam$delta_sampler) {
  if(!identical(sParam$model,"shared_template_group_HS_seasonal_differences"))
    stop("Initialize a group-horseshoe seasonal-deviation state first")
  if(!identical(evol_error,"HS"))stop("Template currently supports evol_error='HS' only")
  if(length(data)!=sParam$Td || k!=sParam$k)stop("State/data dimensions or k changed")
  sampler<-match.arg(sampler,c("auto","precision","sqrt"))
  delta_sampler<-match.arg(delta_sampler,c("auto","precision","sqrt"))
  data<-.season2_data(data,sParam$Td)
  sigma<-.season2_positive(obserror$sigma_e,1L,"obserror$sigma_e")
  auxiliary<-.season2_auxiliary(data,obserror$sigma_et^2,
    sParam$s_full,sParam$M,k)
  template<-.season2_draw_template(auxiliary$data-sParam$delta,auxiliary$variance,
    sigma,1/sParam$s_evolParams0$sigma_w0[1]^2,
    as.numeric(sParam$s_evolParamsShape$tauLambdaj),sampler=sampler)
  sParam$template_raw<-template$raw;sParam$template<-template$template
  sParam$s_evolParams0<-.season2_update_initial(template$raw,sigma,sParam$s_evolParams0)
  sParam$s_evolParamsShape<-.season2_update_hs(
    diff(template$raw,differences=2)/sigma,sParam$s_evolParamsShape)
  residual<-sweep(auxiliary$data,2L,sParam$template,"-")
  sParam$delta<-.season3_draw_delta(residual,auxiliary$variance,
    sParam$delta_hs$tau2*sParam$delta_hs$lambda2,sampler=delta_sampler)
  sParam$delta_sampler_backend<-attr(sParam$delta,"sampler_backend")
  sParam$delta_hs<-.season3_update_scales(sParam$delta,sParam$delta_hs)
  sParam$gaussian_sampler<-sampler;sParam$delta_sampler<-delta_sampler
  .season3_finish(sParam)
}
