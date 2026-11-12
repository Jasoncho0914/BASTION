# Infer fixed deviation periods from complete cycles on the declared schedule.
# Missing values/dates inside the fitting span still count as scheduled slots.
.calendar_deviation_periods <- function(layouts,labels,overrides=NULL) {
  if(!is.null(overrides)) {
    if(is.null(names(overrides)) || anyDuplicated(names(overrides)) ||
       any(!names(overrides)%in%labels))
      stop("deviation_periods overrides must be uniquely named seasonalities")
    overrides<-vapply(overrides,function(x).calendar_integer(x,"deviation_periods entry",3L),integer(1))
  }
  details<-lapply(seq_along(layouts),function(j) {
    d<-layouts[[j]];span<-range(d$original_to_slot)
    slots<-split(seq_along(d$cycle_id),d$cycle_id)
    complete<-vapply(slots,function(ii)min(ii)>=span[1] && max(ii)<=span[2],logical(1))
    lengths<-lengths(slots)[complete];counts<-table(lengths)
    label<-labels[j]
    override<-!is.null(overrides) && label%in%names(overrides)
    modes<-if(length(counts))as.integer(names(counts)[counts==max(counts)]) else integer()
    if(!override && length(modes)!=1L)
      stop("Cannot infer a unique deviation period for '",label,
        "': no complete cycles or tied modal lengths. Supply a named deviation_periods override.")
    period<-if(override)overrides[[label]] else .calendar_integer(modes,"inferred deviation period",3L)
    list(period=period,source=if(override)"explicit_override" else "modal_complete_cycle_length",
      complete_cycle_lengths=as.integer(lengths),
      cycle_length_counts=data.frame(length=as.integer(names(counts)),count=as.integer(counts)),
      complete_cycle_ids=as.integer(names(slots)[complete]))
  })
  names(details)<-labels
  list(periods=setNames(vapply(details,`[[`,integer(1),"period"),labels),details=details)
}

# Calendar fitter helpers.
.calendar_integer <- function(x,name,minimum) {
  if(!is.numeric(x) || length(x)!=1L || !is.finite(x) || x!=as.integer(x) || x<minimum)
    stop(name," must be an integer >= ",minimum)
  as.integer(x)
}

# Weighted version of the same Gaussian regression conditional/prior.
.calendar_regression_draw <- function(residual,X,relative_variance,obs) {
  p<-ncol(X);vp<-obs$sigma_e^2*as.numeric(relative_variance)
  .gaussian_sampler_select(
    function(){
      w<-1/obs$sigma_et^2
      .seasoncal_gaussian(Matrix::crossprod(X,Matrix::Diagonal(x=w)%*%X)+
        Matrix::Diagonal(x=1/vp),as.numeric(Matrix::crossprod(X,residual*w)))
    },
    function()as.numeric(.gaussian_sqrt_draw(X,residual,Matrix::Diagonal(p),
      obs$sigma_et^2,vp)),sampler="auto",block="Calendar regression")
}

.calendar_regression_init <- function(data,X,obs) {
  beta<-as.numeric(.calendar_regression_draw(data,X,rep(1,ncol(X)),obs))
  params<-dsp_initEvol0(beta/obs$sigma_e)
  list(beta=beta,beta_params=params,Td=ncol(X),s_mu=as.numeric(X%*%beta),
    n_squared_sum=sum((beta/params$sigma_w0)^2))
}
.calendar_regression_update <- function(data,X,state,obs) {
  state$beta<-as.numeric(.calendar_regression_draw(data,X,state$beta_params$sigma_w0^2,obs))
  state$beta_params<-dsp_sampleEvol0(state$beta/obs$sigma_e,state$beta_params,commonSD=FALSE)
  state$s_mu<-as.numeric(X%*%state$beta)
  state$n_squared_sum<-sum((state$beta/state$beta_params$sigma_w0)^2)
  state
}

# The only missing-value likelihood augmentation. Its observation count is N,
# because these sampled latent observations are conditioned on during the sweep.
# Integrating them out restores the likelihood using observed data alone.
.calendar_complete_data <- function(data,mean,sd) {
  missing<-is.na(data);completed<-data
  if(any(missing))completed[missing]<-stats::rnorm(sum(missing),mean[missing],sd[missing])
  completed
}

# Calendar pads lie outside the global fitting range. They remain NA, rather
# than contributing artificial observations or observation-variance dimensions.
# Their positive SD placeholder is ignored by the calendar sampler likelihood.
.calendar_season_inputs <- function(residual,obs,layout) {
  sd<-rep(obs$sigma_e,layout$n_slots)
  sd[layout$original_to_slot]<-obs$sigma_et
  list(data=calendar_expand(residual,layout),
       obserror=list(sigma_e=obs$sigma_e,sigma_et=sd))
}

#' BASTION decomposition with calendar seasonality
#'
#' Bayesian trend and seasonal decomposition using observation dates and
#' calendar seasonal specifications.
#'
#' @details BASTION decomposes a time series into a locally adaptive trend, one or more
#' seasonal components, and optional regression and outlier components. The
#' calendar fitter uses observation dates to align seasonal patterns with
#' calendar cycles. All components are estimated jointly using Bayesian
#' inference.
#'
#' Supply observations through `y` and matching `POSIXct` timestamps through
#' `dates`. Use `data_unit` to describe the observation interval, such as
#' `list(unit = "day")` for daily data or `list(unit = "hour")` for hourly data.
#'
#' Specify seasonal patterns as a named list through `seasonalities`. For
#' example, `list(weekly = list(unit = "week"), yearly = list(unit = "year"))`
#' fits weekly and yearly seasonality together. Calendar alignment accommodates
#' different cycle lengths, such as months with different numbers of days and
#' leap years. Each seasonal component combines a smooth recurring pattern with
#' deviations that allow the pattern to evolve across cycles.
#'
#' The trend uses second differences by default (`D_trend = 2`). Use `D_trend =
#' 1` for first differences. Setting `sparse = TRUE` applies stronger shrinkage
#' to trend differences while retaining local adaptivity.
#'
#' Use `obsSV = "const"` for constant observation variance or `obsSV = "SV"` to
#' estimate observation variability over time. Set `Outlier = TRUE` to estimate a
#' separate outlier component, and supply `X` to include regression effects.
#'
#' Missing observations are estimated during fitting. The returned timeline
#' includes all scheduled timestamps, and `imputation` identifies those with
#' missing observations. Component estimates and observation standard deviations
#' are returned in the units of the supplied data. Posterior means and credible
#' intervals at the level specified by `cl` are returned whether or not posterior
#' draws are retained. Set `save_samples = TRUE` to retain the draws.
#'
#' Use `print()` for an overview of the fit, `summary()` for the component
#' estimates, `plot()` for the decomposition, and `plot_component()` to inspect
#' an individual component with credible intervals. Seasonal component names
#' follow the names supplied in `seasonalities`, such as `seasonal_weekly` and
#' `seasonal_yearly`.
#'
#' @section Supplying observation dates and intervals:
#' `y` must be a numeric vector, and `dates` must be a `POSIXct` vector of the
#' same length, with a valid time zone. Each timestamp identifies the observation
#' in the corresponding position of `y`. Convert character dates before calling
#' the fitter, for example with `as.POSIXct(date_text, tz = "UTC")`.
#'
#' `data_unit` describes how often observations are collected. Supply a list
#' containing `unit` and, optionally, `every` and `anchor`. For example,
#' `list(unit = "day")` means one observation per day, `list(unit = "day", every
#' = 2)` means one observation every two days, and `list(unit = "minute", every =
#' 15)` means one observation every 15 minutes. The default value of `every` is
#' 1.
#'
#' Use the intended observation interval even when some observations are missing.
#' For daily data with missing days, keep `data_unit = list(unit = "day")`.
#' Missing responses may be supplied as `NA`; omitted timestamps are also
#' recognized as missing positions in that interval. If `X` is supplied, include
#' every timestamp and its covariates, using `NA` in `y` for missing responses.
#'
#' Allowed units are `"millisecond"`, `"second"`, `"minute"`, `"hour"`, `"day"`,
#' `"week"`, `"month"`, `"quarter"`, and `"year"`. Day, week, month, quarter, and
#' year use calendar intervals and require a positive integer `every`. Subdaily
#' units use elapsed time and allow positive fractional values, such as
#' `list(unit = "second", every = 0.01)` for observations every 10 milliseconds.
#'
#' @section Specifying seasonal patterns:
#' `seasonalities` is a named list with one list per seasonal component. Each
#' component has `unit`, optional `every` (default 1), optional `anchor`, and
#' optional `distance_scale` (default `TRUE`). The allowed units are the same as
#' for `data_unit`. The list names label the results; for example, `biweekly =
#' list(unit = "week", every = 2)` creates a component named `seasonal_biweekly`.
#'
#' Here `unit` and `every` describe how often a seasonal pattern repeats. With
#' daily observations, `list(unit = "week")` repeats every week, `list(unit =
#' "week", every = 2)` repeats every two weeks, and `list(unit = "month")`
#' repeats every calendar month. With hourly observations, `list(unit = "day")`
#' repeats every day. Supplying multiple entries estimates those seasonal
#' components together.
#'
#' Keep the observation interval and the seasonal period separate: daily
#' observations with a biweekly pattern use `data_unit = list(unit = "day")` and
#' `seasonalities = list(biweekly = list(unit = "week", every = 2))`. The
#' seasonal cycle must contain enough observations to estimate its shape; a daily
#' pattern cannot be estimated from one observation per day.
#'
#' `anchor` sets a reference boundary for an interval or seasonal cycle and must
#' be a single `POSIXct` timestamp. The observation interval defaults to the
#' earliest supplied timestamp. Seasonal boundaries default to midnight on
#' January 1, 1970, except weekly cycles use Monday, January 5, 1970, in the
#' selected time zone. For a two-week cycle, set a custom anchor to specify which
#' week starts the cycle.
#'
#' `distance_scale = TRUE` accounts for spacing between seasonal positions when
#' smoothing a grid-based recurring pattern. Set it to `FALSE` to use unscaled
#' second differences for that pattern. This option is supplied separately within
#' each seasonal specification; it does not change the observation interval or
#' the number of days in a calendar cycle.
#'
#' @param y Numeric observations, allowing NA.
#' @param dates POSIXct dates, one per supplied data value.
#' @param data_unit List with unit, optional every (default 1), optional anchor.
#' @param seasonalities Named list of lists with unit/every/anchor and optional
#'   distance_scale (default TRUE, with matched innovation SDs). Shorthand list("daily") works.
#' @param X Optional numeric covariate matrix, one row per supplied date. Supply
#'   every collection date (y may be NA) when using X, so covariates are known.
#' @param Outlier Include horseshoe+ outliers at all observations (default FALSE).
#' @param cl Pointwise credible interval probability (default .95).
#' @param sparse Fix the global trend scale to 1/(1000*(T-D_trend))
#'   when TRUE (default FALSE).
#' @param obsSV Observation variance model, "const" or "SV" (default "const").
#' @param nchains Number of chains (default 2).
#' @param nsave Saved draws per chain (default 1000).
#' @param nburn Burn-in sweeps per chain (default 1000).
#' @param nskip Discarded sweeps between saved draws (default 4).
#' @param verbose Display progress (default TRUE).
#' @param save_samples Return posterior draws (default FALSE).
#' @param D_trend Trend difference order, 1 or 2 (default 2).
#' @param deviation_periods Optional named fixed-lag overrides; otherwise infer
#'   the unique modal complete-cycle length in collection slots. Ties or no
#'   complete cycles require an explicit override.
#' @param tz Optional explicit calendar time zone, otherwise attached dates zone.
#' @param tolerance_seconds Optional preprocessing timestamp matching tolerance.
#' @return A `bastion_fit` with the same components/intervals/seasonality/noise
#'   structure as `fit_BASTION`, with named calendar seasonalities. Main rows
#'   cover all scheduled timestamps, including missing responses. Recurring
#'   summaries have J rows and mean/lower/upper columns; deviations have T rows.
#'   Native phase positions, cycle-centering mappings, input specifications,
#'   and fitting arguments are stored under `metadata`. `imputation` identifies
#'   timestamps with missing responses. Residuals retain latent completed-data
#'   draws at missing responses. Intervals are returned with or without draws.
#' @seealso [simulate_BASTION_calendar()]
#' @examples
#' \dontrun{
#' # Daily observations with weekly and yearly patterns
#' sim <- simulate_BASTION_calendar("daily", seed = 20)
#' fit_daily <- fit_BASTION_calendar(
#'   y = sim$y,
#'   dates = sim$dates,
#'   data_unit = list(unit = "day"),
#'   seasonalities = list(
#'     weekly = list(unit = "week"),
#'     yearly = list(unit = "year")
#'   ),
#'   obsSV = "const", nchains = 1
#' )
#' plot(fit_daily)
#'
#' # Daily observations; the pattern repeats every two weeks
#' sim <- simulate_BASTION_calendar("biweekly", seed = 21)
#' fit_biweekly <- fit_BASTION_calendar(
#'   y = sim$y,
#'   dates = sim$dates,
#'   data_unit = list(unit = "day"),
#'   seasonalities = list(
#'     biweekly = list(unit = "week", every = 2, anchor = sim$dates[1])
#'   ),
#'   obsSV = "const", nchains = 1
#' )
#' plot_component(fit_biweekly, "seasonal_biweekly")
#'
#' # Hourly observations with daily and weekly patterns
#' sim <- simulate_BASTION_calendar("hourly", seed = 22)
#' fit_hourly <- fit_BASTION_calendar(
#'   y = sim$y,
#'   dates = sim$dates,
#'   data_unit = list(unit = "hour"),
#'   seasonalities = list(
#'     daily = list(unit = "day"),
#'     weekly = list(unit = "week")
#'   ),
#'   obsSV = "const", nchains = 1
#' )
#' plot(fit_hourly)
#' }
#' @export
fit_BASTION_calendar <- function(y,dates,data_unit,seasonalities,X=NULL,
    Outlier=FALSE,cl=.95,sparse=FALSE,obsSV="const",nchains=2L,
    nsave=1000L,nburn=1000L,nskip=4L,verbose=TRUE,save_samples=FALSE,
    D_trend=2L,tz=NULL,tolerance_seconds=NULL,
    deviation_periods=NULL) {
  if(missing(dates) || is.null(dates))
    stop("dates must be supplied as a nonempty POSIXct vector")
  # Shared recurring templates and group-horseshoe seasonal differences.
  deviation_model<-"seasonality3"
  template_model<-"cal2"
  spline_boundary<-"periodic"
  D_trend<-.validate_D_trend(D_trend)
  nchains<-.calendar_integer(nchains,"nchains",1L)
  nsave<-.calendar_integer(nsave,"nsave",1L)
  nburn<-.calendar_integer(nburn,"nburn",0L)
  nskip<-.calendar_integer(nskip,"nskip",0L)
  for(name in c("Outlier","sparse","verbose","save_samples")) {
    z<-get(name);if(!is.logical(z) || length(z)!=1L || is.na(z))stop(name," must be TRUE or FALSE")
  }
  if(length(cl)!=1L || !is.numeric(cl) || !is.finite(cl) || cl<=0 || cl>=1)stop("cl must lie strictly between 0 and 1")
  if(!is.character(obsSV) || length(obsSV)!=1L || is.na(obsSV) || !(obsSV%in%c("const","SV")))stop("obsSV must be 'const' or 'SV'")
  supplied_y<-y
  initial<-prepare_calendar(supplied_y,dates,data_unit,seasonalities,tz,tolerance_seconds)
  # All component layouts use the same declared collection times. Use the
  # first layout only to recover the common grid between the endpoint data.
  first<-initial$seasonalities[[1]]
  span<-seq.int(min(first$original_to_slot),max(first$original_to_slot))
  grid_dates<-first$slot_time[span]
  original_to_grid<-match(first$original_to_slot,span)
  N<-length(span);n<-length(supplied_y)
  if(N<=D_trend)stop("Too few collection points for D_trend")
  if(Outlier && N<5L)stop("Outlier model requires at least five collection points")
  y<-rep(NA_real_,N);y[original_to_grid]<-supplied_y
  observed<-!is.na(y);missing<-which(!observed)
  if(sum(observed)<2L)stop("At least two observed data values are required")
  offset<-mean(y[observed]);scale<-stats::sd(y[observed])
  if(!is.finite(scale) || scale<=0)stop("Observed data must have positive finite sample SD")
  z<-(y-offset)/scale
  Xgrid<-NULL
  if(!is.null(X)) {
    if(!is.matrix(X) || !is.numeric(X) || nrow(X)!=n || ncol(X)<1L || any(!is.finite(X)))
      stop("X must be a finite numeric matrix with one row per supplied date")
    if(N!=n)stop("X is needed at every collection time: add omitted dates with y=NA and supply their covariates")
    Xgrid<-matrix(NA_real_,N,ncol(X),dimnames=list(NULL,colnames(X)))
    Xgrid[original_to_grid,]<-X
  }
  prepared<-prepare_calendar(z,grid_dates,data_unit,seasonalities,tz,tolerance_seconds)
  layouts<-prepared$seasonalities
  labels<-names(layouts);season_names<-paste0("Seasonal_",labels);p<-length(layouts)
  names(layouts)<-season_names
  period_inference<-NULL
  if(deviation_model=="seasonality3") {
    inferred<-.calendar_deviation_periods(layouts,labels,deviation_periods)
    deviation_periods<-inferred$periods;period_inference<-inferred$details
  } else if(deviation_model=="seasonality2") {
    if(is.null(deviation_periods) || length(deviation_periods)!=p ||
       is.null(names(deviation_periods)) || anyDuplicated(names(deviation_periods)) ||
       !setequal(names(deviation_periods),labels))
      stop("Supply deviation_periods named exactly as seasonalities, in collection slots")
    deviation_periods<-vapply(deviation_periods[labels],function(x)
      .calendar_integer(x,"deviation_periods entry",3L),integer(1))
  } else if(!is.null(deviation_periods)) {
    stop("deviation_periods applies only to fixed-period deviation models")
  }
  component_names<-c(season_names,"Trend",if(!is.null(Xgrid))"Regression",if(Outlier)"Outlier")
  error_names<-c(season_names,"Trend",if(Outlier)"Outlier")
  designs<-lapply(seq_len(p),function(j){
    d<-layouts[[j]]
    build_seasonality_cal(d$grid,d$cycle_id,d$grid_index,d$distance_scale,
      step_size=initial$data_unit$every,deviation_model=deviation_model,
      deviation_period=if(deviation_model%in%c("seasonality2","seasonality3"))deviation_periods[[j]] else NULL,
      deviation_slot_index=if(deviation_model%in%c("seasonality2","seasonality3"))d$original_to_slot else NULL)
  })
  if(template_model=="cal2")designs<-lapply(seq_len(p),function(j)
    build_seasonality_cal_adaptive(designs[[j]],deviation_periods[[j]],spline_boundary))
  names(designs)<-season_names
  # Calendar variance-prior hyperparameters and initialization settings.
  settings<-list(data_unit=initial$data_unit,seasonalities=lapply(layouts,`[[`,"specification"),
    D_trend=D_trend,Outlier=Outlier,sparse=sparse,obsSV=obsSV,cl=cl,
    nchains=nchains,nsave=nsave,nburn=nburn,nskip=nskip,
    delta_shape=if(deviation_model%in%c("none","seasonality3"))NA_real_ else .001,
    delta_rate=if(deviation_model%in%c("none","seasonality3"))NA_real_ else .001,init_evol_variance=.01,
    template_sampler="auto",template_update="homoskedastic_auxiliary_structured",
    template_cycle_patterns=vapply(designs,function(d)length(d$template_groups),integer(1)),
    template_distance_scale=vapply(layouts,`[[`,logical(1),"distance_scale"),
    template_prior_scale_match="ordinary_D2_at_same_grid_size",
    template_innovation_scale=lapply(designs,`[[`,"template_innovation_scale"),
    deviation_model=if(deviation_model=="seasonality3")"seasonality3_fixed_blocks" else if(deviation_model=="seasonality2")"seasonality2_fixed_blocks" else if(deviation_model=="none")"none" else if(deviation_model=="iid")"iid" else "time_RW2",
    deviation_difference_order=if(deviation_model%in%c("seasonality2","seasonality3"))1L else if(deviation_model=="none")NA_integer_ else if(deviation_model=="iid")0L else 2L,
    deviation_spacing_power=if(deviation_model%in%c("seasonality2","seasonality3"))0 else if(deviation_model=="none")NA_real_ else if(deviation_model=="iid")0 else 3,deviation_prior_rank=vapply(designs,`[[`,numeric(1),"delta_rank"),
    deviation_step_size=initial$data_unit$every,
    deviation_step_unit=initial$data_unit$unit)
  if(deviation_model%in%c("seasonality2","seasonality3")) {
    settings$deviation_periods<-deviation_periods
    settings$deviation_anchor<-grid_dates[1L]
    settings$deviation_blocks<-vapply(designs,function(d)d$delta_geometry$M,integer(1))
    settings$deviation_constraints<-"row_and_column_sums_on_fixed_blocks"
    settings$deviation_sampler<-if(deviation_model=="seasonality3")"Seasonality3_tridiagonal_with_QR_fallback" else "Seasonality2_exact_spectral_auxiliary"
    if(deviation_model=="seasonality3") {
      settings$deviation_period_inference<-period_inference
      settings$deviation_global_half_Cauchy_scale<-1
      settings$deviation_local_half_Cauchy_scale<-1
      settings$deviation_scale_units<-"standardized_data"
      settings$deviation_scaled_by_observation_variance<-FALSE
      settings$deviation_local_scale<-"one_per_cycle_transition"
    }
  }
  spline_components<-which(vapply(designs,function(d)identical(d$template_kind,"bspline"),logical(1)))
  if(template_model=="cal2") {
    settings$template_model<-"cal2"
    settings$spline_boundary<-spline_boundary
    settings$template_cutoff<-100L
    settings$template_slots_per_knot_interval<-3L
    settings$template_kind<-setNames(vapply(designs,`[[`,character(1),"template_kind"),season_names)
    settings$template_basis<-setNames(lapply(designs,function(d)if(is.null(d$spline))NULL else
      d$spline[c("degree","intervals","internal_knots","boundary_knots","coefficients","boundary","coefficient_penalty")]),season_names)
    settings$template_coefficient_prior<-"ordinary_D2_HS_observation_scaled_shared_first_two_boundary_scale"
    for(j in spline_components)settings$template_innovation_scale[[j]]<-rep(1,designs[[j]]$spline$coefficients)
    settings$template_update<-"grid_structured_or_bspline_cached_Gram_auxiliary"
    settings$template_distance_scale_applies_to<-"grid_templates_only"
  }
  total<-nchains*nsave;C<-length(component_names)
  samples<-list(beta=array(NA_real_,c(total,N,C),dimnames=list(NULL,NULL,component_names)),
    beta_combined=matrix(NA_real_,total,N),obs_sigma_t2=matrix(NA_real_,total,N),
    evol_sigma_t2=array(NA_real_,c(total,N,length(error_names)),dimnames=list(NULL,NULL,error_names)),
    remainder=matrix(NA_real_,total,N),yhat=matrix(NA_real_,total,N),
    common_seasonality=setNames(lapply(layouts,function(d)matrix(NA_real_,total,d$grid_size)),season_names),
    seasonal_template_variance=setNames(lapply(layouts,function(d)matrix(NA_real_,total,d$grid_size)),season_names),
    seasonal_delta_variance=matrix(NA_real_,total,p,dimnames=list(NULL,season_names)),
    template_sampler_backend=matrix(NA_character_,total,p,dimnames=list(NULL,season_names)),
    imputed_data=matrix(NA_real_,total,length(missing)),
    chain=rep(seq_len(nchains),each=nsave))
  if(deviation_model=="seasonality3") {
    samples$seasonal_delta_global_variance<-matrix(NA_real_,total,p,dimnames=list(NULL,season_names))
    samples$seasonal_delta_local_variance<-setNames(lapply(designs,function(d)
      matrix(NA_real_,total,d$delta_geometry$M-1L)),season_names)
    samples$seasonal_delta_innovation_variance<-samples$seasonal_delta_local_variance
    samples$seasonal_delta_sampler_backend<-matrix(NA_character_,total,p,dimnames=list(NULL,season_names))
  }
  if(!is.null(Xgrid))samples$reg_coef<-matrix(NA_real_,total,ncol(Xgrid),dimnames=list(NULL,colnames(Xgrid)))
  if(obsSV=="SV")samples$sv_parameters<-matrix(NA_real_,total,3L,
    dimnames=list(NULL,c("mu","phi","innovation_sd")))
  if(template_model=="cal2") {
    samples$seasonal_spline_coefficients<-setNames(lapply(designs[spline_components],function(d)
      matrix(NA_real_,total,d$spline$coefficients)),season_names[spline_components])
    samples$seasonal_spline_coefficient_variance<-samples$seasonal_spline_coefficients
  }
  nsteps<-nburn+(nskip+1L)*nsave
  started<-proc.time()[["elapsed"]]
  for(chain in seq_len(nchains)) {
    .bastion_progress_chain(chain,verbose)
    # Initialization only: missing values at standardized mean. They are sampled
    # before every parameter sweep; no interpolated values enter as fixed data.
    completed<-z;completed[missing]<-0
    obs<-init_sigmaE_0(completed)
    comp<-matrix(0,N,C,dimnames=list(NULL,component_names))
    sp<-vector("list",p)
    for(j in seq_len(p)) {
      d<-layouts[[j]];input<-.calendar_season_inputs(completed,obs,d)
      initialize_season<-if(template_model=="cal2")init_Sbeta_cal_adaptive else init_Sbeta_cal
      sp[[j]]<-initialize_season(input$data,input$obserror,d$grid,d$cycle_id,d$grid_index,
        distance_scale=d$distance_scale,step_size=initial$data_unit$every,design=designs[[j]],
        deviation_model=deviation_model)
      comp[,j]<-as.numeric(calendar_restore(sp[[j]]$s_mu,d))
    }
    tp<-init_Tbeta(completed,obs,"HS",D_trend,sparse)
    comp[,"Trend"]<-as.numeric(tp$s_mu)
    bp<-NULL;zp<-NULL;sv<-NULL
    if(!is.null(Xgrid)) {bp<-.calendar_regression_init(completed,Xgrid,obs);comp[,"Regression"]<-bp$s_mu}
    if(Outlier) {zp<-init_Outlier(completed,obs);comp[,"Outlier"]<-as.numeric(zp$s_mu)}
    if(obsSV=="SV") {sv<-dsp_initSV(completed-rowSums(comp));obs$sigma_e<-1;obs$sigma_et<-as.numeric(sv$sigma_wt)}
    saved<-0L
    if(verbose)pb<-.bastion_progress_bar(nsteps)
    for(iteration in seq_len(nsteps)) {
      tryCatch({
        completed<-.calendar_complete_data(z,rowSums(comp),obs$sigma_et)
        residual<-completed-rowSums(comp)
        if(obsSV=="SV") {
          sv<-dsp_sampleSVparams(residual,sv);obs$sigma_et<-as.numeric(sv$sigma_wt)
        } else {
          rank<-N+sum(vapply(sp,`[[`,numeric(1),"prior_rank"))+
            if(is.null(bp))0 else bp$Td
          if(Outlier)rank<-rank+zp$Td
          energy<-tp$n_squared_sum+sum(vapply(sp,`[[`,numeric(1),"n_squared_sum"))
          if(!is.null(bp))energy<-energy+bp$n_squared_sum
          if(Outlier)energy<-energy+zp$n_squared_sum
          obs$sigma_e<-sqrt(extraDistr::rinvgamma(1,(N+rank)/2,(sum(residual^2)+energy)/2))
          if(!is.finite(obs$sigma_e) || obs$sigma_e<=0)stop("Invalid observation variance draw")
          obs$sigma_et<-rep(obs$sigma_e,N)
        }
        if(!is.null(bp)) {
          partial<-residual+comp[,"Regression"]
          bp<-.calendar_regression_update(partial,Xgrid,bp,obs)
          comp[,"Regression"]<-bp$s_mu;residual<-partial-bp$s_mu
        }
        for(j in seq_len(p)) {
          partial<-residual+comp[,j]
          input<-.calendar_season_inputs(partial,obs,layouts[[j]])
          update_season<-if(template_model=="cal2")fit_Sbeta_cal_adaptive else fit_Sbeta_cal
          sp[[j]]<-update_season(input$data,sp[[j]],input$obserror)
          comp[,j]<-as.numeric(calendar_restore(sp[[j]]$s_mu,layouts[[j]]))
          residual<-partial-comp[,j]
        }
        partial<-residual+comp[,"Trend"]
        tp<-fit_Tbeta(partial,tp,obs,"HS",D_trend,sparse)
        comp[,"Trend"]<-as.numeric(tp$s_mu);residual<-partial-comp[,"Trend"]
        if(Outlier) {
          partial<-residual+comp[,"Outlier"]
          zp<-fit_Outlier(partial,zp,obs)
          comp[,"Outlier"]<-as.numeric(zp$s_mu);residual<-partial-comp[,"Outlier"]
        }
      },error=function(e)stop("Calendar fit: chain ",chain,", sweep ",iteration,": ",conditionMessage(e),call.=FALSE))
      if(iteration>nburn && (iteration-nburn)%%(nskip+1L)==0L) {
        saved<-saved+1L;i<-(chain-1L)*nsave+saved
        physical<-comp*scale;physical[,"Trend"]<-physical[,"Trend"]+offset
        samples$beta[i,,]<-physical
        mu<-rowSums(physical)
        samples$beta_combined[i,]<-mu
        samples$obs_sigma_t2[i,]<-(obs$sigma_et*scale)^2
        if(obsSV=="SV") {
          samples$sv_parameters[i,]<-as.numeric(sv$svParams[,1L])
          samples$sv_parameters[i,"mu"]<-samples$sv_parameters[i,"mu"]+2*log(scale)
        }
        samples$remainder[i,]<-(completed*scale+offset-mu)
        samples$yhat[i,]<-(mu+stats::rnorm(N,sd=obs$sigma_et*scale))
        samples$imputed_data[i,]<-completed[missing]*scale+offset
        samples$evol_sigma_t2[i,,"Trend"]<-(obs$sigma_e*scale)^2*
          as.numeric(c(tp$s_evolParams0$sigma_w0,tp$s_evolParams$sigma_wt))^2
        if(Outlier) {
          vv<-(obs$sigma_e*scale*zp$s_evolParams$sigma_wt)^2
          samples$evol_sigma_t2[i,,"Outlier"]<-vv
        }
        for(j in seq_len(p)) {
          samples$template_sampler_backend[i,j]<-sp[[j]]$template_sampler_backend
          samples$common_seasonality[[j]][i,]<-sp[[j]]$template*scale
          if(j%in%spline_components) {
            samples$seasonal_spline_coefficients[[season_names[j]]][i,]<-sp[[j]]$theta*scale
            samples$seasonal_spline_coefficient_variance[[season_names[j]]][i,]<-.seasoncal2_variances(sp[[j]],obs$sigma_e*scale)
          } else samples$seasonal_template_variance[[j]][i,]<-.seasoncal_template_variances(
            sp[[j]]$design,obs$sigma_e*scale,
            as.numeric(c(sp[[j]]$s_evolParams0$sigma_w0,sp[[j]]$s_evolParamsShape$sigma_wt)))
          samples$seasonal_delta_variance[i,j]<-sp[[j]]$delta_variance*scale^2
          if(deviation_model=="seasonality3") {
            samples$seasonal_delta_global_variance[i,j]<-sp[[j]]$delta_hs$tau2*scale^2
            samples$seasonal_delta_local_variance[[j]][i,]<-sp[[j]]$delta_hs$lambda2
            samples$seasonal_delta_innovation_variance[[j]][i,]<-sp[[j]]$delta_hs$tau2*sp[[j]]$delta_hs$lambda2*scale^2
            samples$seasonal_delta_sampler_backend[i,j]<-sp[[j]]$delta_sampler_backend
          }
        }
        if(!is.null(bp))samples$reg_coef[i,]<-bp$beta*scale
      }
      if(verbose).bastion_progress_tick(pb,iteration,nsteps)
    }
  }
  mapping<-.bastion_calendar_mapping(initial,layouts,designs,seq_len(N))
  original_observed<-rep(NA_real_,N);original_observed[original_to_grid]<-supplied_y
  .bastion_fitted_object(samples,original_observed,grid_dates,
    lapply(layouts,function(d)length(d$grid)),X,Outlier,obsSV,cl,
    nchains,nsave,nburn,nskip,save_samples,scale=1,
    calendar=mapping,recurring_positions=lapply(layouts,`[[`,"grid"),
    imputation_time=grid_dates[missing],
    input_metadata=list(seasonalities=seasonalities),sparse=sparse,D_trend=D_trend)
}
