# Fitted result formatting. Sampling and priors live in the component files;
# this layer summarizes saved states.
.bastion_draw_summary <- function(x,cl) {
  x<-as.matrix(x)
  alpha<-(1-cl)/2
  data.frame(mean=colMeans(x),
    lower=apply(x,2L,stats::quantile,probs=alpha,names=FALSE),
    upper=apply(x,2L,stats::quantile,probs=1-alpha,names=FALSE),
    row.names=NULL,check.names=FALSE)
}

.bastion_fitted_object <- function(samples,observed,time,Ks,X,Outlier,obsSV,
    cl,nchains,nsave,nburn,nskip,save_samples,scale,
    calendar=NULL,recurring_positions=NULL,imputation_time=NULL,
    input_metadata=list(),sparse=FALSE,D_trend=2L) {
  T<-length(observed);D<-dim(samples$beta)[1L]
  old_season_names<-paste0("Seasonal",unlist(Ks,use.names=FALSE))
  season_names<-paste0("seasonal_",unlist(Ks,use.names=FALSE))
  if(!is.null(calendar)) {
    season_names<-names(calendar$seasonality)
    old_season_names<-sub("^seasonal_","Seasonal_",season_names)
  }
  old_names<-c("Trend",old_season_names,if(!is.null(X))"Regression",if(Outlier)"Outlier")
  component_names<-c("trend",season_names,if(!is.null(X))"regression",if(Outlier)"outlier",
    "residual")
  M<-length(component_names)
  # Preserve all array dimensions when there is only one saved draw.
  components<-array(NA_real_,c(D,T,M),
    dimnames=list(NULL,as.character(time),component_names))
  components[,,seq_len(M-1L)]<-samples$beta[,,old_names,drop=FALSE]
  components[,,"residual"]<-matrix(samples$remainder,D,T)*scale
  component_summaries<-lapply(seq_len(M),function(j)
    .bastion_draw_summary(matrix(components[,,j],D,T),cl))
  table_for<-function(field) {
    value<-as.data.frame(setNames(lapply(component_summaries,`[[`,field),
      component_names),check.names=FALSE)
    rownames(value)<-NULL
    value
  }
  result<-list(time=time,observed=as.numeric(observed),
    components=table_for("mean"),
    intervals=list(lower=table_for("lower"),upper=table_for("upper")),
    seasonality=setNames(vector("list",length(Ks)),season_names))
  seasonal_draws<-setNames(vector("list",length(Ks)),season_names)
  for(j in seq_along(Ks)) {
    k<-if(is.null(calendar))as.integer(Ks[[j]]) else length(recurring_positions[[j]])
    recurring<-matrix(samples$common_seasonality[[old_season_names[j]]],D,k)
    colnames(recurring)<-as.character(seq_len(k))
    phase<-if(is.null(calendar))(seq_len(T)-1L)%%k+1L else
      calendar$seasonality[[j]]$phase_index
    deviation<-matrix(components[,,season_names[j]],D,T)-recurring[,phase,drop=FALSE]
    if(!is.null(calendar)) {
      mapping<-calendar$seasonality[[j]]
      for(g in seq_along(mapping$patterns)) {
        ii<-which(mapping$pattern_id==g)
        if(length(ii)) {
          center<-rowMeans(recurring[,mapping$patterns[[g]],drop=FALSE])
          deviation[,ii]<-deviation[,ii,drop=FALSE]+center
        }
      }
    }
    colnames(deviation)<-as.character(time)
    result$seasonality[[j]]<-list(
      recurring=.bastion_draw_summary(recurring,cl),
      deviation=.bastion_draw_summary(deviation,cl))
    if(save_samples)seasonal_draws[[j]]<-list(recurring=recurring,deviation=deviation)
  }
  regression_draws<-NULL
  if(!is.null(X)) {
    p<-ncol(X)
    coefficients<-matrix(samples$reg_coef,D,p)*scale
    terms<-colnames(X)
    if(is.null(terms))terms<-paste0("beta_",seq_len(p))
    colnames(coefficients)<-terms
    result$regression<-list(coefficients=data.frame(term=terms,
      .bastion_draw_summary(coefficients,cl),row.names=NULL,check.names=FALSE))
    if(save_samples)regression_draws<-list(coefficients=coefficients)
  }
  sd<-if(obsSV=="SV")sqrt(matrix(samples$obs_sigma_t2,D,T)) else
    sqrt(matrix(samples$obs_sigma_t2[,1L],D,1L))
  colnames(sd)<-if(obsSV=="SV")as.character(time) else "sd"
  result$noise<-list(sd=.bastion_draw_summary(sd,cl))
  noise_draws<-list(sd=sd)
  if(obsSV=="SV") {
    parameters<-matrix(samples$sv_parameters,D,3L,
      dimnames=list(NULL,c("mu","phi","innovation_sd")))
    result$noise$parameters<-data.frame(parameter=colnames(parameters),
      .bastion_draw_summary(parameters,cl),row.names=NULL,check.names=FALSE)
    noise_draws$parameters<-parameters
  }
  positions<-if(is.null(calendar))lapply(Ks,function(k)seq_len(k)) else recurring_positions
  names(positions)<-season_names
  result$metadata<-list(input=input_metadata,
    fitting=list(obsSV=obsSV,Outlier=Outlier,sparse=sparse,D_trend=D_trend,
      cl=cl,alpha=1-cl,nchains=nchains,nburn=nburn,nsave=nsave,nskip=nskip),
    mapping=list(recurring_positions=positions,
      recurring_prior=setNames(rep("second differences",length(season_names)),season_names)))
  if(!is.null(calendar)) {
    result$metadata$input$data_unit<-calendar$data_unit
    result$metadata$input$data_anchor<-calendar$data_anchor
    result$metadata$input$tz<-calendar$tz
    result$metadata$mapping$seasonality<-calendar$seasonality
    result$metadata$mapping$recurring_prior<-vapply(calendar$seasonality,
      `[[`,character(1),"recurring_prior")
    result$imputation<-imputation_time
  }
  if(save_samples) {
    iterations<-nburn+(nskip+1L)*seq_len(nsave)
    result$draws<-list(index=data.frame(draw=seq_len(D),
      chain=rep(seq_len(nchains),each=nsave),iteration=rep(iterations,nchains)),
      components=components,seasonality=seasonal_draws)
    if(!is.null(regression_draws))result$draws$regression<-regression_draws
    result$draws$noise<-noise_draws
  }
  class(result)<-c("bastion_fit","list")
  result
}

# Retain only the phase selection and full-cycle centering patterns, rather
# than dense evaluation matrices or repeated copies of the observation data.
.bastion_calendar_mapping <- function(initial,layouts,designs,original_to_grid) {
  mapping<-lapply(seq_along(layouts),function(j) {
    layout<-layouts[[j]];groups<-designs[[j]]$template_groups
    slot_pattern<-integer(layout$n_slots)
    for(g in seq_along(groups))slot_pattern[as.numeric(groups[[g]]$slots)]<-g
    slots<-layout$original_to_slot[original_to_grid]
    list(unit=layout$specification$unit,every=layout$specification$every,
      anchor=layout$resolved_anchor,
      recurring_prior=if(identical(designs[[j]]$template_kind,"bspline"))
        "periodic cubic spline" else "second differences",
      phase_index=layout$grid_index[slots],pattern_id=slot_pattern[slots],
      patterns=lapply(groups,`[[`,"indices"))
  })
  names(mapping)<-sub("^Seasonal_","seasonal_",names(layouts))
  list(data_unit=initial$data_unit,data_anchor=initial$data_anchor,
    tz=initial$tz,seasonality=mapping)
}

#' Print a BASTION fit
#'
#' Display the model priors, sampling settings, credible interval level,
#' and the first five rows of component estimates with their time values.
#' @param x A fitted object of class `bastion_fit`.
#' @param ... Additional arguments passed to printing the component preview.
#' @return The fitted object, invisibly.
#' @export
print.bastion_fit <- function(x,...) {
  n<-nrow(x$components)
  seasonal_names<-names(x$seasonality)
  f<-x$metadata$fitting
  calendar<-x$metadata$mapping$seasonality
  detail<-if(is.null(calendar)) {
    paste0("period = ",paste(unlist(x$metadata$input$Ks),collapse=", "))
  } else {
    labels<-vapply(seq_along(calendar),function(j) {
      spec<-calendar[[j]]
      paste0(sub("^seasonal_","",seasonal_names[j]),": ",
        if(spec$every==1)"" else paste0(spec$every," "),spec$unit)
    },character(1))
    paste0("calendar: ",paste(labels,collapse=", "))
  }
  cat(n,"observations ·",length(seasonal_names),"seasonal components",
    if(length(seasonal_names))paste0(" (",detail,")") else "", "\n\n",sep=" ")
  cat("Prior:\n")
  if(!is.null(f)) {
    trend<-if(f$D_trend==1L)"First Differences" else "Second Differences"
    cat("- Trend: Horseshoe on ",trend,
      if(isTRUE(f$sparse))" (fixed global scale)" else "","\n",sep="")
  }
  priors<-x$metadata$mapping$recurring_prior
  if(is.null(priors) && !is.null(calendar))priors<-vapply(calendar,function(spec)
    if(is.null(spec$recurring_prior))"unspecified" else spec$recurring_prior,character(1))
  if(is.null(priors))priors<-setNames(rep("second differences",length(seasonal_names)),seasonal_names)
  for(name in seasonal_names) {
    prior<-priors[[name]]
    label<-if(identical(prior,"second differences"))"Horseshoe on Second Differences" else
      if(identical(prior,"periodic cubic spline"))
        "Periodic Cubic Spline (Horseshoe on Coefficient Second Differences)" else prior
    cat("- ",sub("^seasonal_","Seasonal_",name),":\n",
      "    Recurring: ",label,"\n",
      "    Deviation: Group Horseshoe on Seasonal Differences\n",sep="")
  }
  if(!is.null(f)) {
    cat("- Noise:",if(f$obsSV=="SV")"Stochastic Volatility" else "Constant Variance","\n")
    cat("- Outlier:",if(isTRUE(f$Outlier))"Horseshoe+" else "False","\n")
  }
  cat("\n")
  if(!is.null(x$regression))cat("Regression coefficients:",
    nrow(x$regression$coefficients),"\n")
  if(length(x$imputation))cat("Imputed observations:",length(x$imputation),"\n")
  if(!is.null(f))cat("Sampling:",f$nchains,"chains ·",f$nburn,
    "burn-in ·",f$nsave,"saved per chain · skip",f$nskip,"\n")
  if(!is.null(f$cl))cat("Credible intervals:",
    paste0(format(100*f$cl,trim=TRUE),"%"),
    paste0("(alpha = ",format(1-f$cl,trim=TRUE),")"),"\n")
  cat("Posterior draws retained:",if(is.null(x$draws))"no" else "yes","\n\n")
  rows<-seq_len(min(5L,n))
  preview<-data.frame(time=x$time[rows],x$components[rows,,drop=FALSE],
    check.names=FALSE)
  print(preview,...)
  if(n>5L)cat("...",n-5L,"more rows\n")
  invisible(x)
}

#' Component summaries of a BASTION fit
#' @param object A fitted object of class `bastion_fit`.
#' @param ... Additional arguments, currently unused.
#' @return The full posterior-mean component data frame, identical to
#'   `object$components`.
#' @export
summary.bastion_fit <- function(object,...) {
  object$components
}
