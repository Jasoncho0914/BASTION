# Template distance_scale defaults to TRUE with coordinate-matched priors:
# raw slope differences / sqrt(half-span), interior innovation SDs J^(3/2).
# Boundary heights have proper scale priors. distance_scale=FALSE uses ordinary D2.
# Optional deviation_model="none": delta is exactly zero; no variance updates.
# Optional deviation_model="seasonality2": fixed-period chronological blocks,
# with an open-cycle/circular-phase penalty and double centering.
# The calendar template retains its own, independent phase coordinates.
# Optional deviation_model="seasonality3": fixed blocks with group-HS seasonal
# differences only; uses tridiagonal/QR Gaussian draws and conjugate scale updates.
# Optional deviation_model="iid": independent Gaussian deviations before the
# per-cycle zero-sum conditioning; no time/distance penalty, rank N-M.
# Calendar seasonal Gibbs sampler: recurring template on phase grid, deviation
# on chronological collection slots. Missing/boundary slots remain latent.
# r = A h + delta + epsilon; A evaluates and cycle-centers the raw template.
# Template priors use the specified phase-distance scaling.
# Deviation: sum(diff(delta, differences=2)^2 / step_size^3),
# lambda ~ IG(.001,.001),
# with a zero sum within every complete calendar cycle. No phase/cross-cycle
# grid coupling, circular deviation differences, or union-grid delta states.
# step_size is the number of declared data_unit units per collection interval;
# no elapsed-second/DST correction. Missing slots preserve this spacing.
# The intrinsic RW2 prior has rank N-M for M>=2. With M=1 it has rank
# N-2 and retains a flat linear mode, which must be identified by data.
# Exact Gaussian conditioning enforces cycle sums using Q^-1 C' and C Q^-1 C'.
# This is a covariance-weighted conditioning correction, not mean subtraction.
# The unconstrained observation-plus-RW2 precision is pentadiagonal.
# Template updates use Gaussian auxiliaries, banded precision, grouped
# low-rank centering corrections, and automatic square-root QR fallback.
# No FFT/DCT, jitter, clipping, or change to the template HS hierarchy.

.seasoncal_contrast <- function(n) {
  if(n==1L)return(Matrix::Matrix(0,1L,0L,sparse=TRUE))
  Matrix::sparseMatrix(i=c(seq_len(n-1L),rep(n,n-1L)),
    j=rep(seq_len(n-1L),2L),x=c(rep(1,n-1L),rep(-1,n-1L)),dims=c(n,n-1L))
}

# Construct deterministic operators from explicitly supplied coordinates.
build_seasonality_cal <- function(grid,cycle_id,grid_index,distance_scale=TRUE,step_size=1,
    deviation_model=c("RW2","iid","none","seasonality2","seasonality3"),
    deviation_period=NULL,deviation_slot_index=NULL) {
  deviation_model<-match.arg(deviation_model)
  if(!is.numeric(grid) || length(grid)<3L || any(!is.finite(grid)) ||
     any(grid<0 | grid>=1) || any(diff(grid)<=0))
    stop("grid must contain >=3 strictly increasing finite positions in [0,1); do not duplicate phase 0 as 1")
  if(!is.logical(distance_scale) || length(distance_scale)!=1L || is.na(distance_scale))
    stop("distance_scale must be explicitly TRUE or FALSE")
  if(!is.numeric(cycle_id) || !length(cycle_id) || any(!is.finite(cycle_id)) ||
     any(cycle_id<1 | cycle_id!=as.integer(cycle_id)))stop("cycle_id must contain positive integer cycle numbers")
  M<-as.integer(max(cycle_id));J<-length(grid);N<-length(cycle_id)
  if(N<3L)stop("RW2 deviation requires at least three chronological collection slots")
  if(!identical(sort(unique(as.integer(cycle_id))),seq_len(M)))stop("cycle_id must cover successive cycles 1,...,M")
  if(!is.numeric(grid_index) || length(grid_index)!=N || any(!is.finite(grid_index)) ||
     any(grid_index<1 | grid_index>J | grid_index!=as.integer(grid_index)))
    stop("grid_index must map every supplied slot to an integer grid index")
  cycle_id<-as.integer(cycle_id);grid_index<-as.integer(grid_index)
  if(anyDuplicated(paste(cycle_id,grid_index,sep=":")))stop("Each cycle/grid position must occur once; include missing slots as NA")
  slots<-split(seq_len(N),factor(cycle_id,levels=seq_len(M)))
  if(is.unsorted(cycle_id) || any(vapply(slots,function(ii)any(diff(grid_index[ii])<=0),logical(1))))
    stop("Calendar sampler slots must follow chronological cycle and phase order")
  step_size<-.season2_positive(step_size,1L,"step_size")
  # Normalize constraint rows for conditioning; their nullspace remains the
  # exact per-cycle zero-sum subspace. No dense reduced-coordinate basis.
  C<-Matrix::sparseMatrix(i=cycle_id,j=seq_len(N),
    x=1/sqrt(lengths(slots)[cycle_id]),dims=c(M,N))
  B<-Matrix::sparseMatrix(i=seq_len(N),j=grid_index,x=1,dims=c(N,J))
  A<-B-Matrix::t(C)%*%(C%*%B)
  h<-diff(grid);rows<-3:J
  if(distance_scale){
    w<-sqrt((h[seq_len(J-2L)]+h[2:(J-1L)])/2)
    left<-1/h[seq_len(J-2L)]/w;right<-1/h[2:(J-1L)]/w
    coefficients<-as.numeric(t(cbind(left,-left-right,right)))
  } else coefficients<-rep(c(1,-2,1),J-2L)
  L<-Matrix::sparseMatrix(i=c(1:2,rep(rows,each=3L)),
    j=c(1:2,as.numeric(t(cbind(rows-2L,rows-1L,rows)))),
    x=c(1,1,coefficients),dims=c(J,J))
  # Complete cycles sharing the same selected phase nodes share one centering
  # direction. Cache their reduced sufficient-statistic design once.
  signature<-vapply(slots,function(ii)paste(grid_index[ii],collapse=":"),character(1))
  pattern<-match(signature,unique(signature))
  groups<-lapply(seq_len(max(pattern)),function(g){
    members<-which(pattern==g);ii<-grid_index[slots[[members[1]]]]
    ng<-length(ii);mg<-length(members)
    Bg<-Matrix::sparseMatrix(i=seq_len(ng),j=ii,x=1,dims=c(ng,J))
    centered<-Bg-Matrix::Matrix(1,ng,1,sparse=TRUE)%*%Matrix::Matrix(colSums(Bg)/ng,nrow=1,sparse=TRUE)
    list(indices=ii,multiplicity=mg,slots=matrix(unlist(slots[members]),nrow=ng),
      reduced_design=sqrt(mg)*centered)
  })
  reduced<-do.call(rbind,lapply(groups,`[[`,"reduced_design"))
  centers<-Matrix::sparseMatrix(i=unlist(lapply(groups,`[[`,"indices")),
    j=rep(seq_along(groups),lengths(lapply(groups,`[[`,"indices"))),
    x=unlist(lapply(groups,function(g)rep(sqrt(g$multiplicity/length(g$indices)),length(g$indices)))),
    dims=c(J,length(groups)))
  Dtime<-if(deviation_model%in%c("none","seasonality2","seasonality3"))Matrix::Matrix(0,0L,N,sparse=TRUE) else if(deviation_model=="iid")Matrix::Diagonal(N) else Matrix::sparseMatrix(i=rep(seq_len(N-2L),3L),
    j=c(seq_len(N-2L),2:(N-1L),3:N),
    x=c(rep(1,N-2L),rep(-2,N-2L),rep(1,N-2L))/step_size^(3/2),
    dims=c(N-2L,N))
  result<-list(grid=as.numeric(grid),cycle_id=cycle_id,grid_index=grid_index,
    distance_scale=distance_scale,M=M,J=J,Td=N,cycle_lengths=lengths(slots),
    template_design=A,template_operator=L,
    template_innovation_scale=c(1,1,rep(if(distance_scale)J^(3/2) else 1,J-2L)),
    delta_operator=Dtime,
    delta_precision=Matrix::crossprod(Dtime),delta_constraints=C,
    deviation_model=deviation_model,
    delta_rank=if(deviation_model=="none")0L else if(deviation_model=="RW2" && M==1L)N-2L else N-M,
    delta_flat_modes=if(deviation_model=="RW2" && M==1L)1L else 0L,
    step_size=step_size,template_groups=groups,
    template_counts=tabulate(grid_index,nbins=J),template_centers=centers,
    template_reduced_design=reduced)
  if(deviation_model%in%c("seasonality2","seasonality3")) {
    if(!is.numeric(deviation_period) || length(deviation_period)!=1L ||
       !is.finite(deviation_period) || deviation_period<3 ||
       deviation_period!=as.integer(deviation_period))
      stop("Supply an integer deviation_period >= 3 in collection slots")
    ii<-if(is.null(deviation_slot_index))seq_len(N) else deviation_slot_index
    if(!is.numeric(ii) || !length(ii) || any(!is.finite(ii)) ||
       any(ii!=as.integer(ii)) || any(ii<1 | ii>N) || any(diff(ii)!=1))
      stop("deviation_slot_index must identify consecutive chronological fitting slots")
    dims<-.season2_dimensions(length(ii),deviation_period)
    result$delta_period<-as.integer(deviation_period)
    result$delta_slot_index<-as.integer(ii)
    result$delta_input_length<-length(ii)
    result$delta_geometry<-if(deviation_model=="seasonality3")
      list(M=dims[["M"]],k=dims[["k"]]) else .season2_geometry(dims[["M"]],dims[["k"]])
    result$delta_rank<-(dims[["M"]]-1L)*(dims[["k"]]-1L)
    # These time-RW2 operators/constraints do not describe the fixed-block prior.
    result$delta_operator<-NULL
    result$delta_precision<-NULL
    result$delta_constraints<-NULL
  }
  result
}


# L remains the selected raw distance-weighted operator. Its prior SD receives
# the matching coordinate-unit factor. On a regular J-point phase grid, both
# L's interior rows and their SDs grow by J^(3/2), recovering ordinary D2 exactly.
# HS hierarchy parameters are for the standardized innovations Lh / scale.
.seasoncal_template_variances <- function(design,sigma,relative_sd) {
  sigma<-.season2_positive(sigma,1L,"template observation scale")
  relative_sd<-.season2_positive(relative_sd,design$J,"template relative SD")
  (sigma*relative_sd*design$template_innovation_scale)^2
}
.seasoncal_template_innovations <- function(design,raw) {
  as.numeric(design$template_operator%*%raw)/design$template_innovation_scale
}

.seasoncal_data <- function(data,design,obserror) {
  if(!is.numeric(data) || length(data)!=design$Td || any(is.infinite(data)))
    stop("data must contain one finite value or NA per supplied cycle slot")
  sigma<-.season2_positive(obserror$sigma_e,1L,"obserror$sigma_e")
  variance<-.season2_positive(obserror$sigma_et^2,design$Td,"obserror$sigma_et^2")
  list(data=as.numeric(data),observed=!is.na(data),sigma=sigma,variance=variance)
}

# U' U=Q. This function samples the actual Gaussian conditional directly.
.seasoncal_gaussian <- function(Q,linear) {
  n<-length(linear)
  if(n==0L)return(numeric())
  if(any(!is.finite(linear)))stop("Nonfinite Gaussian linear term")
  U<-Matrix::chol(Matrix::forceSymmetric(Q),pivot=FALSE)
  z<-as.numeric(Matrix::solve(U,Matrix::solve(Matrix::t(U),as.numeric(linear))+stats::rnorm(n)))
  if(any(!is.finite(z)))stop("Nonfinite calendar Gaussian draw; no repair applied")
  z
}

.seasoncal_template_conditional <- function(input,design,delta,prior_variance) {
  o<-input$observed;A<-design$template_design[o,,drop=FALSE]
  residual<-input$data[o]-delta[o]
  v<-input$variance[o]
  L<-design$template_operator
  vp<-.season2_positive(prior_variance,design$J,"template prior variance")
  Q<-Matrix::crossprod(A,Matrix::Diagonal(x=1/v)%*%A)+
    Matrix::crossprod(L,Matrix::Diagonal(x=1/vp)%*%L)
  list(Q=Q,linear=as.numeric(Matrix::crossprod(A,residual/v)))
}


# Apply automatic precision/square-root Gaussian sampling to the
# template's actual centered design and prior operator. Same Gaussian target.
.seasoncal_draw_template_direct <- function(input,design,delta,prior_variance) {
  .gaussian_sampler_select(
    function(){
      q<-.seasoncal_template_conditional(input,design,delta,prior_variance)
      .seasoncal_gaussian(q$Q,q$linear)
    },
    function(){
      o<-input$observed
      as.numeric(.gaussian_sqrt_draw(design$template_design[o,,drop=FALSE],
        input$data[o]-delta[o],design$template_operator,
        input$variance[o],prior_variance))
    },sampler="auto",block="Calendar recurring template")
}


# Auxiliary observations use chronological collection slots.
.seasoncal_template_auxiliary <- function(input,design,delta,current_raw) {
  signal<-as.numeric(design$template_design%*%current_raw)
  z<-.season2_auxiliary(input$data-delta,input$variance,
    matrix(signal,nrow=1L),1L,design$Td)
  z<-as.numeric(z$data);c0<-min(input$variance)
  means<-lapply(design$template_groups,function(g)
    rowMeans(matrix(z[g$slots],nrow=nrow(g$slots))))
  rhs<-unlist(Map(function(g,m)sqrt(g$multiplicity)*m,design$template_groups,means))
  linear<-numeric(design$J)
  for(g in seq_along(means)) {
    info<-design$template_groups[[g]];m<-means[[g]]
    linear[info$indices]<-linear[info$indices]+info$multiplicity*(m-mean(m))/c0
  }
  list(variance=c0,response=rhs,linear=linear)
}

# Banded Gaussian draw with low-rank corrections for distinct
# calendar cycle layouts. No T-by-J cross-product is formed per Gibbs sweep.
.seasoncal_template_structured_precision <- function(design,aux,prior_variance) {
  L<-design$template_operator
  K<-Matrix::crossprod(L,Matrix::Diagonal(x=1/prior_variance)%*%L)+
    Matrix::Diagonal(x=design$template_counts/aux$variance)
  U<-Matrix::chol(Matrix::forceSymmetric(K),pivot=FALSE)
  F<-as.matrix(design$template_centers)/sqrt(aux$variance)
  solved<-as.matrix(Matrix::solve(U,Matrix::solve(Matrix::t(U),cbind(aux$linear,F))))
  m0<-solved[,1L];Z<-solved[,-1L,drop=FALSE]
  # S=I-F'K^-1 F. Use S=S0'S0+Z'QZ, a sum of Gram matrices,
  # to avoid cancellation in the nearly constant template direction.
  S0<-diag(ncol(F))-crossprod(F,Z)
  Rprior<-as.matrix((L%*%Z)/sqrt(prior_variance))
  Robs<-as.matrix(design$template_reduced_design%*%Z)/sqrt(aux$variance)
  S<-crossprod(S0)+crossprod(Rprior)+crossprod(Robs)
  R<-base::chol(S)
  correction<-backsolve(R,forwardsolve(t(R),as.numeric(crossprod(F,m0))))
  mu<-as.numeric(m0+Z%*%correction)
  # Exact mean identity from a centered likelihood and the shared boundary
  # variance: E(h[1]+h[2] | data)=0. This is not a constraint on random draws.
  mu<-mu-mean(mu[1:2])
  mu+as.numeric(Matrix::solve(U,stats::rnorm(design$J)))+
    as.numeric(Z%*%backsolve(R,stats::rnorm(ncol(F))))
}

.seasoncal_draw_template <- function(input,design,delta,prior_variance,current_raw) {
  aux<-.seasoncal_template_auxiliary(input,design,delta,current_raw)
  .gaussian_sampler_select(
    function().seasoncal_template_structured_precision(design,aux,prior_variance),
    function()as.numeric(.gaussian_sqrt_draw(design$template_reduced_design,
      aux$response,design$template_operator,aux$variance,prior_variance)),
    sampler="auto",block="Calendar recurring template (structured auxiliary)")
}

.seasoncal_delta_conditional <- function(input,design,raw,lambda) {
  o<-input$observed
  residual<-input$data[o]-as.numeric(design$template_design[o,,drop=FALSE]%*%raw)
  lambda<-.season2_positive(lambda,1L,"delta variance")
  weight<-linear<-numeric(design$Td)
  weight[o]<-1/input$variance[o];linear[o]<-residual*weight[o]
  Q<-Matrix::Diagonal(x=weight)+design$delta_precision/lambda
  if(design$deviation_model=="RW2" && sum(o)<2L) {
    if(design$M==1L && (!any(o) ||
        all((seq_len(design$Td)-(design$Td+1)/2)[o]==0)))
      stop("Single-cycle intrinsic RW2 has an unidentified linear mode: observations must identify its slope; no boundary prior was added")
    # Two cycle sums identify the affine RW2 nullspace. A penalty only in
    # already constrained directions is zero on C delta=0 and is removed by
    # exact conditioning. It changes no prior or boundary assumption.
    c0<-design$delta_constraints[seq_len(min(2L,design$M)),,drop=FALSE]
    Q<-Q+Matrix::crossprod(c0)/lambda
  }
  list(Q=Q,linear=linear)
}


# iid N(0,lambda) prior, conditioned on per-cycle zero sums.
# Independent Gaussian observation/prior product draw.
# Under SV, covariance-weighted conditioning is required; ordinary unweighted
# mean subtraction would give a different posterior. All operations are O(N).
.seasoncal_iid_delta_draw <- function(input,design,raw,lambda) {
  lambda<-.season2_positive(lambda,1L,"delta variance")
  o<-input$observed
  residual<-numeric(design$Td)
  residual[o]<-input$data[o]-as.numeric(design$template_design[o,,drop=FALSE]%*%raw)
  v<-input$variance;v[!o]<-Inf
  post_var<-1/(1/v+1/lambda)
  z<-sampleOutlier(residual,v,rep(lambda,design$Td),design$Td)
  sums<-as.numeric(rowsum(z,design$cycle_id,reorder=FALSE))
  variance_sums<-as.numeric(rowsum(post_var,design$cycle_id,reorder=FALSE))
  delta<-z-post_var*(sums/variance_sums)[design$cycle_id]
  if(any(!is.finite(delta)))stop("Nonfinite iid calendar deviation draw; no repair applied")
  delta
}

.seasoncal_constrained_draw <- function(Q,linear,C) {
  U<-Matrix::chol(Matrix::forceSymmetric(Q),pivot=FALSE)
  solve_Q<-function(b)Matrix::solve(U,Matrix::solve(Matrix::t(U),b))
  x<-as.numeric(solve_Q(linear)+Matrix::solve(U,stats::rnorm(length(linear))))
  Z<-as.matrix(solve_Q(Matrix::t(C)))
  S<-as.matrix(C%*%Z)
  R<-base::chol(S)
  multiplier<-backsolve(R,forwardsolve(t(R),as.numeric(C%*%x)))
  delta<-as.numeric(x-Z%*%multiplier)
  if(any(!is.finite(delta)))stop("Nonfinite constrained calendar deviation draw; no repair applied")
  delta
}

# Auxiliary observation construction and spectral Gaussian draw.
# Only the recurring mean is evaluated on the calendar grid before subtraction.
.seasoncal_fixed_delta_auxiliary <- function(input,design,raw,delta_grid) {
  ii<-design$delta_slot_index
  recurring<-as.numeric(design$template_design[ii,,drop=FALSE]%*%raw)
  .season2_auxiliary(input$data[ii]-recurring,input$variance[ii],delta_grid,
    design$delta_geometry$M,design$delta_geometry$k)
}

.seasoncal_delta_energy <- function(state) {
  if(state$design$deviation_model=="seasonality3")
    sum(.season3_innovations(state$delta_grid)^2) else
  if(state$design$deviation_model=="seasonality2")
    .season2_energy(state$delta_grid) else
    sum(as.numeric(state$design$delta_operator%*%state$delta)^2)
}

.seasoncal_finish <- function(state) {
  d<-state$design
  if(d$deviation_model%in%c("seasonality2","seasonality3")) {
    # Calendar-only pads carry no deviation likelihood. Fixed-block end pads
    # remain latent in delta_grid and contribute to its prior and constraints.
    state$delta<-numeric(d$Td)
    state$delta[d$delta_slot_index]<-as.numeric(t(state$delta_grid))[seq_len(d$delta_input_length)]
  }
  state$template<-state$template_raw-mean(state$template_raw)
  state$common_at_observations<-as.numeric(d$template_design%*%state$template_raw)
  state$deviation_at_observations<-state$delta
  state$s_full<-state$common_at_observations+state$delta
  state$s_mu<-matrix(state$common_at_observations+state$deviation_at_observations,ncol=1L)
  innovations<-.seasoncal_template_innovations(d,state$template_raw)
  scales<-c(state$s_evolParams0$sigma_w0,state$s_evolParamsShape$sigma_wt)
  state$n_squared_sum<-sum((innovations/as.numeric(scales))^2)
  state$delta_energy<-.seasoncal_delta_energy(state)
  state
}

# Public-to-the-source-file entry point, separate from the main fitter.
# Grid inputs are explicit; distance-weighted differences use matched SD units.
init_Sbeta_cal <- function(data,obserror,grid,cycle_id,grid_index,distance_scale=TRUE,
    evol_error="HS",delta_shape=0.001,delta_rate=0.001,
    delta_variance_init=NULL,init_evol_variance=0.01,step_size=1,design=NULL,
    deviation_model=c("RW2","iid","none","seasonality2","seasonality3"),
    deviation_period=NULL,deviation_slot_index=NULL) {
  deviation_model<-match.arg(deviation_model)
  if(!identical(evol_error,"HS"))stop("Calendar template currently supports evol_error='HS' only")
  d<-if(is.null(design))build_seasonality_cal(grid,cycle_id,grid_index,distance_scale,step_size,deviation_model,
    deviation_period,deviation_slot_index) else design
  if(!identical(d$deviation_model,deviation_model))stop("Prebuilt design and requested deviation model differ")
  input<-.seasoncal_data(data,d,obserror)
  shape<-if(deviation_model%in%c("none","seasonality3"))NA_real_ else .season2_positive(delta_shape,1L,"delta_shape")
  rate<-if(deviation_model%in%c("none","seasonality3"))NA_real_ else .season2_positive(delta_rate,1L,"delta_rate")
  initial<-.season2_positive(init_evol_variance,1L,"init_evol_variance")
  lambda<-if(deviation_model%in%c("none","seasonality3"))NA_real_ else if(is.null(delta_variance_init))input$sigma^2 else
    .season2_positive(delta_variance_init,1L,"delta_variance_init")
  delta<-numeric(d$Td)
  raw<-.seasoncal_draw_template(input,d,delta,
    .seasoncal_template_variances(d,input$sigma,sqrt(initial)),numeric(d$J))
  e<-.seasoncal_template_innovations(d,raw)
  p0<-dsp_initEvol0(raw[1:2]/input$sigma,commonSD=TRUE)
  ps<-dsp_initEvolParams(e[-(1:2)]/input$sigma,evol_error="HS")
  scales<-c(p0$sigma_w0,p0$px_sigma_w0,ps$tauLambdaj,ps$xiLambdaj,ps$tauLambda,ps$xiLambda)
  if(any(!is.finite(scales)) || any(scales<=0))stop("Horseshoe initialization produced invalid scales")
  state<-list(model=if(deviation_model=="seasonality3")"calendar_shared_template_group_HS_seasonal_differences" else if(deviation_model=="seasonality2")"calendar_shared_template_fixed_blocks" else if(deviation_model=="none")"calendar_shared_template_only" else if(deviation_model=="iid")"calendar_shared_template_iid_deviation" else "calendar_shared_template_time_RW2",
    deviation_model=deviation_model,Td=d$Td,J=d$J,M=d$M,
    design=d,distance_scale=distance_scale,template_raw=raw,delta=delta,
    delta_variance=lambda,delta_shape=shape,delta_rate=rate,delta_rank=d$delta_rank,
    prior_rank=d$J,s_evolParams0=p0,s_evolParamsShape=ps,evol_error=evol_error,
    colname="SeasonalCalendar",sampler_backend="precision",step_size=d$step_size,
    template_sampler="auto",template_sampler_backend=attr(raw,"sampler_backend"),
    template_update="homoskedastic_auxiliary_structured",
    template_innovation_scale=d$template_innovation_scale)
  if(deviation_model%in%c("seasonality2","seasonality3"))
    state$delta_grid<-matrix(0,d$delta_geometry$M,d$delta_geometry$k)
  if(deviation_model=="seasonality3") {
    state$delta_hs<-list(tau2=1,lambda2=rep(1,d$delta_geometry$M-1L),
      nu=rep(1,d$delta_geometry$M-1L),xi=1)
    state$delta_sampler_backend<-"initial_zero"
  }
  .seasoncal_finish(state)
}

# One exact block Gibbs sweep; observation and trend updates belong outside.
fit_Sbeta_cal <- function(data,sParam,obserror,evol_error=sParam$evol_error) {
  if(!(sParam$model%in%c("calendar_shared_template_time_RW2","calendar_shared_template_iid_deviation","calendar_shared_template_only","calendar_shared_template_fixed_blocks","calendar_shared_template_group_HS_seasonal_differences")))stop("Initialize a calendar seasonal state first")
  if(!identical(evol_error,"HS"))stop("Calendar template currently supports evol_error='HS' only")
  d<-sParam$design;input<-.seasoncal_data(data,d,obserror)
  relative_sd<-as.numeric(c(sParam$s_evolParams0$sigma_w0,sParam$s_evolParamsShape$sigma_wt))
  raw<-.seasoncal_draw_template(input,d,sParam$delta,
    .seasoncal_template_variances(d,input$sigma,relative_sd),sParam$template_raw)
  sParam$template_sampler_backend<-attr(raw,"sampler_backend")
  sParam$template_raw<-raw
  sParam$s_evolParams0<-.season2_update_initial(raw,input$sigma,sParam$s_evolParams0)
  e<-.seasoncal_template_innovations(d,raw)
  sParam$s_evolParamsShape<-.season2_update_hs(e[-(1:2)]/input$sigma,sParam$s_evolParamsShape)
  if(d$delta_rank){
    if(d$deviation_model%in%c("seasonality2","seasonality3")) {
      aux<-.seasoncal_fixed_delta_auxiliary(input,d,raw,sParam$delta_grid)
      if(d$deviation_model=="seasonality3") {
        sParam$delta_grid<-.season3_draw_delta(aux$data,aux$variance,
          sParam$delta_hs$tau2*sParam$delta_hs$lambda2,sampler="auto")
        sParam$delta_sampler_backend<-attr(sParam$delta_grid,"sampler_backend")
      } else sParam$delta_grid<-.season2_draw_delta(aux$data,aux$variance,
        sParam$delta_variance,d$delta_geometry)
    } else if(d$deviation_model=="iid") {
      sParam$delta<-.seasoncal_iid_delta_draw(input,d,raw,sParam$delta_variance)
    } else {
      q<-.seasoncal_delta_conditional(input,d,raw,sParam$delta_variance)
      sParam$delta<-.seasoncal_constrained_draw(q$Q,q$linear,d$delta_constraints)
    }
  }
  if(d$deviation_model=="seasonality3") {
    sParam$delta_hs<-.season3_update_scales(sParam$delta_grid,sParam$delta_hs)
    sParam$delta_innovation_variance<-sParam$delta_hs$tau2*sParam$delta_hs$lambda2
    return(.seasoncal_finish(sParam))
  }
  if(d$deviation_model=="none")return(.seasoncal_finish(sParam))
  energy<-.seasoncal_delta_energy(sParam)
  precision<-stats::rgamma(1,shape=sParam$delta_shape+d$delta_rank/2,
    rate=sParam$delta_rate+energy/2)
  sParam$delta_variance<-1/precision
  if(!is.finite(sParam$delta_variance) || sParam$delta_variance<=0)
    stop("Invalid calendar deviation variance draw; no clipping applied")
  .seasoncal_finish(sParam)
}


# ---- Adaptive recurring template (periodic for long cycles) ----
# Adaptive calendar shapes. Long-cycle templates default to a periodic cubic
# basis with q=ceiling(k/3) coefficients/intervals. Periodicity is built into
# evaluation.
# Ordinary coefficient D2 and the first-two-coefficient shared boundary scale
# retain conjugate updates. The coefficient prior is anchored, not circular.
# spline_boundary="clamped" selects a clamped basis.

.seasoncal2_periodic_basis <- function(x,intervals,derivative=0L) {
  q<-.calendar_integer(intervals,"periodic knot intervals",4L)
  derivative<-.calendar_integer(derivative,"derivative",0L)
  if(derivative>2L || !is.numeric(x) || any(!is.finite(x)) || any(x<0 | x>1))
    stop("Periodic basis needs x in [0,1] and derivative 0, 1 or 2")
  # Extend uniform knots three intervals past each boundary; identify columns
  # separated by one full period. At both endpoints value/derivatives agree.
  extended_knots<-seq.int(-3L,q+3L)/q
  B<-splines::splineDesign(extended_knots,x,ord=4L,
    derivs=rep(derivative,length(x)),outer.ok=FALSE,sparse=TRUE)
  wrap<-Matrix::sparseMatrix(i=seq_len(q+3L),
    j=(seq_len(q+3L)-1L)%%q+1L,x=1,dims=c(q+3L,q))
  B%*%wrap
}

build_seasonality_cal_adaptive <- function(design,cycle_length,spline_boundary=c("periodic","clamped")) {
  spline_boundary<-match.arg(spline_boundary)
  if(!identical(design$deviation_model,"seasonality3"))
    stop("Adaptive calendar templates require group-horseshoe seasonal deviations")
  k<-.calendar_integer(cycle_length,"cycle_length",3L)
  design$template_kind<-if(k<100L)"grid" else "bspline"
  design$template_cycle_length<-k
  if(k<100L)return(design)
  q<-as.integer(ceiling(k/3));knots<-seq(0,1,length.out=q+1L)[-c(1L,q+1L)]
  B<-if(spline_boundary=="periodic") .seasoncal2_periodic_basis(design$grid,q) else
    Matrix::Matrix(splines::bs(design$grid,knots=knots,degree=3L,
      intercept=TRUE,Boundary.knots=c(0,1)),sparse=TRUE)
  P<-ncol(B)
  R<-design$template_reduced_design%*%B
  design$spline<-list(basis=B,degree=3L,intervals=q,internal_knots=knots,
    boundary_knots=c(0,1),coefficients=P,boundary=spline_boundary,
    coefficient_penalty="ordinary_second_differences_with_shared_first_two_scale",
    operator=.gaussian_difference_operator(P,2L),reduced_design=R,
    observation_gram=Matrix::crossprod(R))
  design
}

.seasoncal2_coefficient_draw <- function(input,design,delta,variance,current_theta,
                                        sampler="auto") {
  b<-design$spline
  raw<-as.numeric(b$basis%*%current_theta)
  aux<-.seasoncal_template_auxiliary(input,design,delta,raw)
  linear<-as.numeric(Matrix::crossprod(b$basis,aux$linear))
  .gaussian_sampler_select(function(){
    Q<-b$observation_gram/aux$variance+
      Matrix::crossprod(b$operator,Matrix::Diagonal(x=1/variance)%*%b$operator)
    .seasoncal_gaussian(Q,linear)
  },function()as.numeric(.gaussian_sqrt_draw(b$reduced_design,aux$response,
    b$operator,aux$variance,variance)),sampler=sampler,
    block="Calendar adaptive B-spline recurring template")
}

.seasoncal2_variances <- function(state,sigma) {
  (sigma*as.numeric(c(state$s_evolParams0$sigma_w0,
    state$s_evolParamsShape$sigma_wt)))^2
}

.seasoncal2_finish <- function(state) {
  # Evaluate the same calendar-centred function and restore fixed-block pads.
  # Prior energy must be computed on coefficients, never on grid heights.
  d<-state$design
  state$template_raw<-as.numeric(d$spline$basis%*%state$theta)
  state$template<-state$template_raw-mean(state$template_raw)
  state$delta<-numeric(d$Td)
  state$delta[d$delta_slot_index]<-as.numeric(t(state$delta_grid))[seq_len(d$delta_input_length)]
  state$common_at_observations<-as.numeric(d$template_design%*%state$template_raw)
  state$deviation_at_observations<-state$delta
  state$s_full<-state$common_at_observations+state$delta
  state$s_mu<-matrix(state$s_full,ncol=1L)
  state$delta_energy<-sum(.season3_innovations(state$delta_grid)^2)
  e<-as.numeric(state$design$spline$operator%*%state$theta)
  sd<-as.numeric(c(state$s_evolParams0$sigma_w0,state$s_evolParamsShape$sigma_wt))
  state$n_squared_sum<-sum((e/sd)^2)
  state
}

init_Sbeta_cal_adaptive <- function(data,obserror,grid,cycle_id,grid_index,
    distance_scale=TRUE,evol_error="HS",init_evol_variance=.01,step_size=1,
    design,deviation_model="seasonality3") {
  if(!identical(deviation_model,"seasonality3") || !identical(design$deviation_model,"seasonality3"))
    stop("Adaptive calendar templates require group-horseshoe seasonal deviations")
  if(identical(design$template_kind,"grid"))return(init_Sbeta_cal(data,obserror,
    grid,cycle_id,grid_index,distance_scale=distance_scale,evol_error=evol_error,
    init_evol_variance=init_evol_variance,step_size=step_size,design=design,
    deviation_model=deviation_model))
  if(!identical(design$template_kind,"bspline"))stop("Build an adaptive calendar template design first")
  if(!identical(evol_error,"HS"))stop("Adaptive calendar templates support horseshoe priors only")
  input<-.seasoncal_data(data,design,obserror)
  initial<-.season2_positive(init_evol_variance,1L,"init_evol_variance")
  P<-design$spline$coefficients
  theta<-.seasoncal2_coefficient_draw(input,design,numeric(design$Td),
    rep(input$sigma^2*initial,P),numeric(P))
  backend<-attr(theta,"sampler_backend")
  p0<-dsp_initEvol0(theta[1:2]/input$sigma,commonSD=TRUE)
  ps<-dsp_initEvolParams(diff(theta,differences=2)/input$sigma,evol_error="HS")
  if(any(!is.finite(c(p0$sigma_w0,p0$px_sigma_w0,ps$tauLambdaj,
      ps$xiLambdaj,ps$tauLambda,ps$xiLambda))) ||
     any(c(p0$sigma_w0,p0$px_sigma_w0,ps$tauLambdaj,ps$xiLambdaj,ps$tauLambda,ps$xiLambda)<=0))
    stop("Horseshoe initialization produced invalid scales")
  state<-list(model="calendar_adaptive_bspline_group_HS_deviation",
    deviation_model=deviation_model,Td=design$Td,J=design$J,M=design$M,
    design=design,distance_scale=distance_scale,theta=as.numeric(theta),
    template_raw=as.numeric(design$spline$basis%*%theta),delta=numeric(design$Td),
    delta_grid=matrix(0,design$delta_geometry$M,design$delta_geometry$k),
    delta_hs=list(tau2=1,lambda2=rep(1,design$delta_geometry$M-1L),
      nu=rep(1,design$delta_geometry$M-1L),xi=1),
    delta_variance=NA_real_,delta_shape=NA_real_,delta_rate=NA_real_,
    delta_rank=design$delta_rank,prior_rank=P,
    s_evolParams0=p0,s_evolParamsShape=ps,evol_error=evol_error,
    colname="SeasonalCalendar",step_size=step_size,
    template_sampler="auto",template_sampler_backend=backend,
    template_update="homoskedastic_auxiliary_bspline",delta_sampler_backend="initial_zero")
  .seasoncal2_finish(state)
}

fit_Sbeta_cal_adaptive <- function(data,sParam,obserror,evol_error=sParam$evol_error) {
  if(!identical(sParam$design$template_kind,"bspline"))
    return(fit_Sbeta_cal(data,sParam,obserror,evol_error))
  if(!identical(evol_error,"HS"))stop("Adaptive calendar templates support horseshoe priors only")
  d<-sParam$design;input<-.seasoncal_data(data,d,obserror)
  theta<-.seasoncal2_coefficient_draw(input,d,sParam$delta,
    .seasoncal2_variances(sParam,input$sigma),sParam$theta,
    sampler=sParam$template_sampler)
  sParam$template_sampler_backend<-attr(theta,"sampler_backend")
  sParam$theta<-as.numeric(theta)
  sParam$s_evolParams0<-.season2_update_initial(theta,input$sigma,sParam$s_evolParams0)
  sParam$s_evolParamsShape<-.season2_update_hs(diff(theta,differences=2)/input$sigma,
    sParam$s_evolParamsShape)
  raw<-as.numeric(d$spline$basis%*%theta)
  if(d$delta_rank) {
    aux<-.seasoncal_fixed_delta_auxiliary(input,d,raw,sParam$delta_grid)
    sParam$delta_grid<-.season3_draw_delta(aux$data,aux$variance,
      sParam$delta_hs$tau2*sParam$delta_hs$lambda2,sampler="auto")
    sParam$delta_sampler_backend<-attr(sParam$delta_grid,"sampler_backend")
  }
  sParam$delta_hs<-.season3_update_scales(sParam$delta_grid,sParam$delta_hs)
  sParam$delta_innovation_variance<-sParam$delta_hs$tau2*sParam$delta_hs$lambda2
  .seasoncal2_finish(sParam)
}
