# Square-root Gaussian sampling for the daily-state trend.
# init_Tbeta()/fit_Tbeta() call sampleTrend_auto(): precision factorization
# on success, square-root QR on factorization/solve failure or factor warnings.
# The automatic path does not apply diagonal repair.
#
# D has identity rows for the first D boundary states, followed by Dth differences.
# A = rbind(diag(1/sqrt(obs_variance)), diag(1/sqrt(evol_variance)) %*% D).
# z = c(data/sqrt(obs_variance), rep(0,T)). No cross-product is formed for QR.
# For A*S*P=U*R, draw theta=R^-1(U' z + N(0,I)), then x=S*P*theta.
# S is optional column equilibration; P is the QR column permutation. Both are
# undone explicitly. This is the same Gaussian conditional, up to floating point.
# No priors, local/global scales, innovation floors or model constraints change.
# No jitter, truncation, pseudo-inverse, diagonal repair or iterative tolerance.
#
# Numerical options are explicit: column_scale=TRUE, qr_order=0 (natural columns),
# Fallback is always silent; notify is retained for call compatibility.
# There is no hardcoded condition-number cutoff.
# Auto fallback detects failures/warnings; successful Cholesky is not proof that
# every ill-conditioned input was solved accurately. Call sampleTrend_sqrt()
# directly when the square-root path is desired even without a Cholesky failure.
# Sparse QR does NOT construct dense Q. For fixed D, the operator has O(T) nonzeros;
# ordering and factor fill determine actual QR cost. Stability is tested, not assumed.

.sqrt_trend_inputs <- function(data,obs_sigma_t2,evol_sigma_t2,D,Td) {
  if(length(Td)!=1L || !is.finite(Td) || Td!=as.integer(Td) ||
     length(D)!=1L || !is.finite(D) || !(D%in%1:3) || Td<=D)
    stop("Require integer D in 1:3 and integer Td > D")
  if(!is.numeric(data) || length(data)!=Td || any(!is.finite(data)))
    stop("data must contain Td finite numeric observations")
  positive<-function(x,name) {
    if(!is.numeric(x) || !(length(x)%in%c(1L,Td)) || any(!is.finite(x)) || any(x<=0))
      stop(name," must be finite and positive, with length 1 or Td")
    rep_len(as.numeric(x),Td)
  }
  list(data=as.numeric(data),obs=positive(obs_sigma_t2,"obs_sigma_t2"),
       evol=positive(evol_sigma_t2,"evol_sigma_t2"),D=as.integer(D),Td=as.integer(Td))
}

.sqrt_trend_operator <- function(Td,D) {
  rows<-seq.int(D+1L,Td)
  coefficients<-choose(D,0:D)*(-1)^(D:0)
  Matrix::sparseMatrix(i=c(seq_len(D),rep(rows,each=D+1L)),
    j=c(seq_len(D),as.numeric(t(outer(rows,-D:0,"+")))),
    x=c(rep(1,D),rep(coefficients,length(rows))),dims=c(Td,Td))
}

.sqrt_trend_column_scale <- function(A) {
  # Stable column norms without squaring unscaled large precision roots.
  scale<-numeric(ncol(A))
  for(j in seq_len(ncol(A))) {
    values<-A@x[seq.int(A@p[j]+1L,A@p[j+1L])]
    maximum<-max(abs(values))
    norm<-maximum*sqrt(sum((values/maximum)^2))
    if(!is.finite(norm) || norm<=0)stop("Invalid column norm in square-root operator")
    scale[j]<-1/norm
  }
  scale
}

#' Square-root QR draw of the daily-state trend
#' @keywords internal
sampleTrend_sqrt <- function(data,obs_sigma_t2,evol_sigma_t2,D=1,Td=length(data),
                              column_scale=TRUE,qr_order=0L) {
  inputs<-.sqrt_trend_inputs(data,obs_sigma_t2,evol_sigma_t2,D,Td)
  if(!is.logical(column_scale) || length(column_scale)!=1L || is.na(column_scale))
    stop("column_scale must be TRUE or FALSE")
  if(length(qr_order)!=1L || !is.finite(qr_order) || !(qr_order%in%0:3))
    stop("qr_order must be an integer in 0:3")
  L<-.sqrt_trend_operator(inputs$Td,inputs$D)
  A<-rbind(Matrix::Diagonal(x=1/sqrt(inputs$obs)),
    Matrix::Diagonal(x=1/sqrt(inputs$evol))%*%L)
  A<-as(A,"dgCMatrix")
  rhs<-c(inputs$data/sqrt(inputs$obs),numeric(inputs$Td))
  if(any(!is.finite(A@x)) || any(!is.finite(rhs)))
    stop("Nonfinite square-root operator or RHS; no clipping applied")
  scale<-if(column_scale).sqrt_trend_column_scale(A) else rep(1,inputs$Td)
  if(column_scale)A<-A%*%Matrix::Diagonal(x=scale)
  factor<-Matrix::qr(A,order=as.integer(qr_order))
  R<-as(Matrix::qr.R(factor,complete=FALSE,backPermute=FALSE),"dtCMatrix")
  diagonal<-Matrix::diag(R)
  if(any(!is.finite(diagonal)) || any(diagonal==0))
    stop("Square-root QR is numerically rank deficient; no model modification applied")
  projected<-as.numeric(Matrix::qr.qty(factor,rhs))[seq_len(inputs$Td)]
  permuted<-as.numeric(Matrix::solve(R,projected+stats::rnorm(inputs$Td)))
  q<-factor@q
  permutation<-if(length(q))as.integer(q)+1L else seq_len(inputs$Td)
  unpermuted<-numeric(inputs$Td)
  unpermuted[permutation]<-permuted
  result<-matrix(scale*unpermuted,ncol=1L)
  if(any(!is.finite(result)))stop("Nonfinite square-root draw; no clipping applied")
  attr(result,"sampler_backend")<-"sqrt_qr"
  result
}

#' Precision fast path with square-root QR fallback, without diagonal repair
#' @keywords internal
sampleTrend_auto <- function(data,obs_sigma_t2,evol_sigma_t2,D=1,Td=length(data),
                              notify=TRUE,column_scale=TRUE,qr_order=0L) {
  inputs<-.sqrt_trend_inputs(data,obs_sigma_t2,evol_sigma_t2,D,Td)
  if(!is.logical(notify) || length(notify)!=1L || is.na(notify))stop("notify must be TRUE or FALSE")
  warnings<-character()
  draw<-tryCatch(withCallingHandlers({
    Q<-build_Q_trend(inputs$obs,inputs$evol,D=inputs$D,Td=inputs$Td)
    U<-Matrix::chol(Q,pivot=FALSE)
    if(length(warnings))stop(paste(warnings,collapse="; "))
    if(any(!is.finite(Matrix::diag(U))) || any(Matrix::diag(U)<=0))
      stop("Nonfinite/nonpositive precision factor diagonal")
    result<-as.matrix(Matrix::solve(U,Matrix::solve(Matrix::t(U),inputs$data/inputs$obs)+
      stats::rnorm(inputs$Td)))
    if(any(!is.finite(result)))stop("Nonfinite precision draw")
    result
  },warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")}),
    error=function(e)e)
  if(!inherits(draw,"error")) {
    attr(draw,"sampler_backend")<-"precision"
    return(draw)
  }
  reason<-conditionMessage(draw)
  result<-sampleTrend_sqrt(inputs$data,inputs$obs,inputs$evol,D=inputs$D,Td=inputs$Td,
    column_scale=column_scale,qr_order=qr_order)
  attr(result,"fallback_reason")<-reason
  result
}
