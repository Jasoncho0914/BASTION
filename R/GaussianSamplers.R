# Shared square-root Gaussian regression machinery.
# Targets: y | beta ~ N(B beta, diag(obs_variance));
#         L beta ~ N(0, diag(prior_variance)), including proper boundary rows.
# Factor the stacked WEIGHTED OPERATORS, never their cross-product, for QR.
# Column equilibration and QR column permutations are explicitly undone.
# There is no jitter, diagonal repair, pseudo-inverse or model modification.
# All callers must supply the actual block operator, including centering where
# the observation mean depends on a centered latent parameter.
#
# sampler='auto' tries the preserved precision routine then QR after errors,
# warnings, or nonfinite draws. 'precision' is strict; 'sqrt' skips that routine.
# No condition-number threshold is hardcoded. A successful precision solve can
# still be inaccurate on some ill-conditioned inputs; sqrt can be requested.
# Backend metadata are stored as attributes.
# Fallback is always silent; notify is retained for call compatibility.

.gaussian_positive <- function(x,n,name) {
  if(!is.numeric(x) || !(length(x)%in%c(1L,n)) || any(!is.finite(x)) || any(x<=0))
    stop(name," must be finite and positive, with length 1 or ",n)
  rep_len(as.numeric(x),n)
}

.gaussian_difference_operator <- function(n,D=2L) {
  if(length(n)!=1L || !is.finite(n) || n!=as.integer(n) || n<=D ||
     length(D)!=1L || !is.finite(D) || !(D%in%1:3))stop("Require integer n>D and D in 1:3")
  rows<-seq.int(D+1L,n)
  Matrix::sparseMatrix(i=c(seq_len(D),rep(rows,each=D+1L)),
    j=c(seq_len(D),as.numeric(t(outer(rows,-D:0,"+")))),
    x=c(rep(1,D),rep(choose(D,0:D)*(-1)^(D:0),length(rows))),dims=c(n,n))
}

.gaussian_sqrt_draw <- function(design,response,prior_operator,obs_variance=1,
                                prior_variance=1,column_scale=TRUE,qr_order=0L) {
  B<-as(as(Matrix::Matrix(design,sparse=TRUE),"generalMatrix"),"CsparseMatrix")
  L<-as(as(Matrix::Matrix(prior_operator,sparse=TRUE),"generalMatrix"),"CsparseMatrix")
  n<-nrow(B);p<-ncol(B)
  if(p<1L || ncol(L)!=p || nrow(L)<1L || !is.numeric(response) || length(response)!=n ||
     any(!is.finite(response)) || any(!is.finite(B@x)) || any(!is.finite(L@x)))
    stop("Invalid design, response or prior operator dimensions/values")
  if(!is.logical(column_scale) || length(column_scale)!=1L || is.na(column_scale))
    stop("column_scale must be TRUE or FALSE")
  if(length(qr_order)!=1L || !is.finite(qr_order) || !(qr_order%in%0:3))
    stop("qr_order must be an integer in 0:3")
  vy<-.gaussian_positive(obs_variance,n,"obs_variance")
  vp<-.gaussian_positive(prior_variance,nrow(L),"prior_variance")
  A<-as(rbind(Matrix::Diagonal(x=1/sqrt(vy))%*%B,
    Matrix::Diagonal(x=1/sqrt(vp))%*%L),"dgCMatrix")
  rhs<-c(as.numeric(response)/sqrt(vy),numeric(nrow(L)))
  if(any(!is.finite(A@x)) || any(!is.finite(rhs)))stop("Nonfinite weighted operator/RHS; no clipping applied")
  scale<-rep(1,p)
  if(column_scale) {
    for(j in seq_len(p)) {
      if(A@p[j+1L]==A@p[j])stop("Zero column in square-root operator")
      values<-A@x[seq.int(A@p[j]+1L,A@p[j+1L])]
      maximum<-max(abs(values))
      norm<-maximum*sqrt(sum((values/maximum)^2))
      if(!is.finite(norm) || norm<=0)stop("Invalid square-root column norm")
      scale[j]<-1/norm
    }
    A<-A%*%Matrix::Diagonal(x=scale)
  }
  factor<-Matrix::qr(A,order=as.integer(qr_order))
  R<-as(Matrix::qr.R(factor,complete=FALSE,backPermute=FALSE),"dtCMatrix")
  diagonal<-Matrix::diag(R)
  if(any(!is.finite(diagonal)) || any(diagonal==0))stop("Numerically rank-deficient square-root QR")
  projected<-as.numeric(Matrix::qr.qty(factor,rhs))[seq_len(p)]
  permuted<-as.numeric(Matrix::solve(R,projected+stats::rnorm(p)))
  q<-factor@q
  permutation<-if(length(q))as.integer(q)+1L else seq_len(p)
  draw<-numeric(p);draw[permutation]<-permuted
  draw<-matrix(scale*draw,ncol=1L)
  if(any(!is.finite(draw)))stop("Nonfinite square-root draw; no clipping applied")
  draw
}

.gaussian_sampler_select <- function(precision_draw,sqrt_draw,
                                     sampler=c("auto","precision","sqrt"),
                                     notify=TRUE,block="Gaussian block") {
  sampler<-match.arg(sampler)
  if(!is.logical(notify) || length(notify)!=1L || is.na(notify))stop("notify must be TRUE or FALSE")
  finite<-function(x)if(is.list(x))all(vapply(x,finite,logical(1))) else all(is.finite(x))
  if(sampler=="sqrt") {
    draw<-sqrt_draw();if(!finite(draw))stop("Nonfinite square-root draw")
    attr(draw,"sampler_backend")<-"sqrt_qr";return(draw)
  }
  warning_messages<-character()
  draw<-tryCatch(withCallingHandlers({
    result<-precision_draw()
    if(length(warning_messages))stop(paste(warning_messages,collapse="; "))
    if(!finite(result))stop("Nonfinite precision draw")
    result
  },warning=function(w){warning_messages<<-c(warning_messages,conditionMessage(w));invokeRestart("muffleWarning")}),
    error=function(e)e)
  if(!inherits(draw,"error")) {
    attr(draw,"sampler_backend")<-"precision";return(draw)
  }
  if(sampler=="precision")stop(draw)
  reason<-conditionMessage(draw)
  draw<-sqrt_draw();if(!finite(draw))stop("Nonfinite square-root draw")
  attr(draw,"sampler_backend")<-"sqrt_qr"
  attr(draw,"fallback_reason")<-reason
  draw
}
