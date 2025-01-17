#' MCMC Sampler for Dynamic Shrinkage with Changepoint
#'
#' Run the MCMC for Bayesian trend filtering with a penalty on zeroth (D = 0),
#' first (D = 1), or second (D = 2) differences of the conditional expectation.
#' The user has the option to utilize threshold shrinkage to identify changepoints
#' for D = 1 or 2.
#' The penalty is determined by the prior on the evolution errors, which include:
#' \itemize{
#' \item the dynamic horseshoe prior ('DHS');
#' \item the static horseshoe prior ('HS');
#' \item the Bayesian lasso ('BL');
#' \item the normal stochastic volatility model ('SV');
#' \item the normal-inverse-gamma prior ('NIG').
#' }
#' In each case, the evolution error is a scale mixture of Gaussians.
#' Note that only 'DHS' works for threshold shrinkage with changepoints.
#' Sampling is accomplished with a (parameter-expanded) Gibbs sampler,
#' mostly relying on a dynamic linear model representation.

#' @param y a numeric vector of the \code{T x 1} vector of time series observations
#' @param family a string specifying expontial family for observation likelihood.
#' Defaults to "gaussian" and implementation also available for "negbinomial"
#' @param cp a logical flag (default is FALSE) indicating to determine whether to use threshold shrinkage with changepoints.
#' @param evol_error the evolution error distribution; must be one of
#' \itemize{
#' \item the dynamic horseshoe prior ('DHS');
#' \item the static horseshoe prior ('HS');
#' \item the Bayesian lasso ('BL');
#' \item the normal stochastic volatility model ('SV');
#' \item the normal-inverse-gamma prior ('NIG').
#' }
#' @param D integer scalar indicating degree of differencing defaults to 1; implementation is available D = 0, D = 1, or D = 2
#' @param obsSV Options for modeling the error variance. It must be one of the following:
#' \itemize{
#' \item const: Constant error variance for all time points.
#' \item SV: Stochastic Volatility model.
#' \item ASV: Adaptive Stochastic Volatility model.
#' }
#' @param useAnom logical (default FALSE); if TRUE, include an anomaly component in the observation equation
#' (only for threshold shrinkage with changepoints.)
#' @param nsave integer scalar (default = 1000); number of MCMC iterations to record
#' @param nburn integer scalar (default = 1000); number of MCMC iterations to discard (burn-in)
#' @param nskip integer scalar (default = 4); number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param mcmc_params list of character scalars naming parameters for which we store the MCMC output;
#' must be one or more of:
#' \itemize{
#' \item "mu" (conditional mean)
#' \item "yhat" (posterior predictive distribution)
#' \item "evol_sigma_t2" (evolution error variance)
#' \item "obs_sigma_t2" (observation error variance)
#' \item "dhs_phi" (DHS AR(1) coefficient)
#' \item "dhs_mean" (DHS AR(1) unconditional mean)
#' }
#' defaults to list("mu", "omega", "r")
#' @param computeDIC logical; if TRUE (default), compute the deviance information criterion \code{DIC}
#' and the effective number of parameters \code{p_d}
#' @param verbose logical; should R report extra information on progress? Defaults to FALSE
#' @param cp_thres Proportion of posterior samples of latent indicator being 1 needed to declare a changepoint; defaults to 0.4
#' @param ... optional additional arguments to pass to \code{\link{btf_nb}} when family = "negbinomial"
#'
#' @return \code{dsp_fit} returns an object of class "\code{dsp}".
#'
#' An object of class "\code{dsp}" is defined as a list containing at least the following components:
#'    \item{pars}{a list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}} ## TODO change so matches
#'    \item{cp}{if threshold shrinkage with changepoints is used, also return detected changepoint locations; otherwise FALSE}
#'    \item{DIC}{}
#'    \item{family}{value supplied for family argument}
#'    \item{evol_error}{value supplied for evol_error argument}
#'    \item{D}{value supplied for D argument}
#'    \item{obsSV}{value supplied for obsSV argument}
#'    \item{mcpar}{named vector of supplied nsave, nburn, and nskip}
#'    \item{cp_thres}{value supplied for cp_thres argument}
#'
#'
#' @note The data \code{y} may contain NAs, which will be treated with a simple imputation scheme
#' via an additional Gibbs sampling step. In general, rescaling \code{y} to have unit standard
#' deviation is recommended to avoid numerical issues when family is "gaussian".
#'
#' @examples
#'
#' # Example 1: Change in mean with stochastic volatility
#' signal = c(rep(0, 50), rep(10, 50))
#' noise = rep(1, 100)
#' noise_var = rep(1, 100)
#' for (k in 2:100){
#'   noise_var[k] = exp(0.9*log(noise_var[k-1]) + rnorm(1, 0, 0.5))
#'   noise[k] = rnorm(1, 0, sqrt(noise_var[k])) }
#'
#' y = signal + noise
#' mcmc_output = dsp_fit(y, cp=TRUE, mcmc_params = list('yhat', 'mu', "omega", "r"))
#' cp = mcmc_output$cp
#' print(paste0('Changepoint Locations: ', cp))
#' plot_fitted(y, mu = colMeans(mcmc_output$mu), postY = mcmc_output$yhat, y_true = signal)
#'
#' # Example 2: Change in mean with Outliers
#' signal = c(rep(0, 100), rep(5, 100))
#' y = c(rep(0, 100), rep(5, 100)) + rnorm(200)
#' y[50] = 10
#' y[150] = -10
#'
#' mcmc_output = dsp_fit(y, cp=TRUE, useAnom = TRUE, mcmc_params = list('yhat', 'mu', "omega", "r"))
#' cp = mcmc_output$cp
#' plot_fitted(y, mu = colMeans(mcmc_output$mu), postY = mcmc_output$yhat, y_true = signal)
#'
#' # Example 3: Change in linear trend
#' signal = c(seq(1, 50), seq(51, 2))
#' y = c(seq(1, 50), seq(51, 2)) + rnorm(100)
#'
#' mcmc_output = dsp_fit(y, cp=TRUE, D=2, mcmc_params = list('yhat', 'mu', "omega", "r"))
#' cp = mcmc_output$cp
#' plot_fitted(y, mu = colMeans(mcmc_output$mu), postY = mcmc_output$yhat, y_true = signal)
#'
#'
#' @export
dsp_fit = function(y, family = "gaussian",
                   cp = FALSE, evol_error = 'DHS',
                   D = 1, obsSV = "const", useAnom = FALSE,
                   nsave = 1000, nburn = 1000, nskip = 4,
                   mcmc_params = list("mu", "omega", "r"),
                   computeDIC = TRUE,
                   verbose = FALSE,
                   cp_thres = 0.4, ...){

  if(!((evol_error == "DHS") || (evol_error == "HS") || (evol_error == "BL") || (evol_error == "SV") || (evol_error == "NIG"))) stop('Error type must be one of DHS, HS, BL, SV, or NIG')
  if(!((D == 0) || (D == 1) || (D == 2))) stop('D must be 0, 1 or 2')
  if(!family %in% c("gaussian", "negbinomial")) stop('family must be gaussian or negbinomial')

  if(family == "gaussian"){
    if (cp == TRUE) {
      if (evol_error == 'DHS' & (D == 1 || D == 2)) {
        # Add in omega and r for finding changepoints
        if (is.na(match('omega', mcmc_params))) {
          mcmc_params = append(mcmc_params, list('omega'))
        }
        if (is.na(match('r', mcmc_params))) {
          mcmc_params = append(mcmc_params, list('r'))
        }
        mcmc_output = abco(y, D = D, obsSV = obsSV, useAnom = useAnom, nsave = nsave,
                           nburn = nburn, nskip = nskip, mcmc_params = mcmc_params, verbose = verbose, cp_thres = cp_thres)
      }
    } else {
      mcmc_output = btf(y, evol_error = evol_error, D = D, obsSV = obsSV, nsave = nsave, nburn = nburn, nskip = nskip,
                        mcmc_params = mcmc_params, computeDIC = computeDIC, verbose = verbose)
    }
  }else{
    if(family == "negbinomial"){
      if(!all(y >= 0)) stop("Negative Binomial likelihood only appropriate for positive data")

      mcmc_output = btf_nb(y = y,
                           evol_error = evol_error,
                           D = D,
                           nsave = nsave, nburn = nburn, nskip = nskip,
                           mcmc_params = mcmc_params,
                           computeDIC = computeDIC,
                           verbose = verbose, ...)
    }
  }


  structure(c(mcmc_output, list(cp = cp_thres, #TODO change documentation or this so matches
                 DIC = mcmc_output[c("DIC", "p_d")],
                 D = D,
                 obsSV = obsSV,
                 mcpar = c(nsave = nsave, nburn = nburn, nskip = nskip),
                 cp_thres = cp_thres)),
            class = c(mcmc_output$class, c("dsp")))

}
#' MCMC Sampler for Bayesian Trend Filtering
#'
#' Run the MCMC for Bayesian trend filtering with a penalty on zeroth (D = 0),
#' first (D = 1), or second (D = 2) differences of the conditional expectation.
#' The penalty is determined by the prior on the evolution errors, which include:
#' \itemize{
#' \item the dynamic horseshoe prior ('DHS');
#' \item the static horseshoe prior ('HS');
#' \item the Bayesian lasso ('BL');
#' \item the normal stochastic volatility model ('SV');
#' \item the normal-inverse-gamma prior ('NIG').
#' }
#' In each case, the evolution error is a scale mixture of Gaussians.
#' Sampling is accomplished with a (parameter-expanded) Gibbs sampler,
#' mostly relying on a dynamic linear model representation.

#' @param y the \code{T x 1} vector of time series observations
#' @param evol_error the evolution error distribution; must be one of
#' 'DHS' (dynamic horseshoe prior; the default), 'HS' (horseshoe prior), 'BL' (Bayesian lasso), or 'NIG' (normal-inverse-gamma prior)
#' @param D degree of differencing (D = 0, D = 1, or D = 2)
#' @param obsSV Options for modeling the error variance. It must be one of the following:
#' \itemize{
#' \item const: Constant error variance for all time points.
#' \item SV: Stochastic Volatility model.
#' \item ASV: Adaptive Stochastic Volatility model.
#' }
#' for the observation error variance
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burnin)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param mcmc_params named list of parameters for which we store the MCMC output;
#' must be one or more of:
#' \itemize{
#' \item "mu" (conditional mean)
#' \item "yhat" (posterior predictive distribution)
#' \item "evol_sigma_t2" (evolution error variance)
#' \item "obs_sigma_t2" (observation error variance)
#' \item "dhs_phi" (DHS AR(1) coefficient)
#' \item "dhs_mean" (DHS AR(1) unconditional mean)
#' }
#' @param computeDIC logical; if TRUE, compute the deviance information criterion \code{DIC}
#' and the effective number of parameters \code{p_d}
#' @param verbose logical; should R report extra information on progress?
#'
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#'
#' @note The data \code{y} may contain NAs, which will be treated with a simple imputation scheme
#' via an additional Gibbs sampling step. In general, rescaling \code{y} to have unit standard
#' deviation is recommended to avoid numerical issues.
#'
#' @examples
#' \dontrun{
#' # Example 1: Bumps Data
#' simdata = simUnivariate(signalName = "bumps", T = 128, RSNR = 7, include_plot = TRUE)
#' y = simdata$y
#'
#' out = btf(y, D=1)
#' plot_fitted(y, mu = colMeans(out$mu), postY = out$yhat, y_true = simdata$y_true)
#'
#' # Example 2: Doppler Data; longer series, more noise
#' simdata = simUnivariate(signalName = "doppler", T = 500, RSNR = 5, include_plot = TRUE)
#' y = simdata$y
#'
#' out = btf(y)
#' plot_fitted(y, mu = colMeans(out$mu), postY = out$yhat, y_true = simdata$y_true)
#'
#' # And examine the AR(1) parameters for the log-volatility w/ traceplots:
#' plot(as.ts(out$dhs_phi)) # AR(1) coefficient
#' plot(as.ts(out$dhs_mean)) # Unconditional mean
#'
#'# Example 3: Blocks data (locally constant)
#' simdata = simUnivariate(signalName = "blocks", T = 1000, RSNR = 3, include_plot = TRUE)
#' y = simdata$y
#'
#' out = btf(y, D = 1) # try D = 1 to approximate the locally constant behavior
#' plot_fitted(y, mu = colMeans(out$mu), postY = out$yhat, y_true = simdata$y_true)
#' }
#'
#' @import progress
#' @export
btf = function(y, evol_error = 'DHS', D = 2, obsSV = "const",
               nsave = 1000, nburn = 1000, nskip = 4,
               mcmc_params = list("mu", "yhat","evol_sigma_t2", "obs_sigma_t2", "dhs_phi", "dhs_mean"),
               computeDIC = TRUE,
               verbose = TRUE){

  # Convert to upper case:
  evol_error = toupper(evol_error)

  # For D = 0, return special case:
  if(D == 0){
    return(btf0(y = y, evol_error = evol_error, obsSV = obsSV,
                nsave = nsave, nburn = nburn, nskip = nskip,
                mcmc_params = mcmc_params,
                computeDIC = computeDIC,
                verbose = verbose))
  }

  # Time points (in [0,1])
  T = length(y); t01 = seq(0, 1, length.out=T);

  # Initialize bandsparse locations
  loc = t_create_loc(length(y)-D, 1)
  loc_obs = t_create_loc(length(y), D)

  # Begin by checking for missing values, then imputing (for initialization)
  is.missing = which(is.na(y)); any.missing = (length(is.missing) > 0)

  # Impute the active "data"
  y = approxfun(t01, y, rule = 2)(t01)

  # Initial SD (implicitly assumes a constant mean)
  sigma_e = sd(y, na.rm=TRUE); sigma_et = rep(sigma_e, T)

  # Compute the Cholesky term (uses random variances for a more conservative sparsity pattern)
  chol0 = NULL

  # Initialize the conditional mean, mu, via sampling:
  mu = sampleBTF(y, obs_sigma_t2 = sigma_et^2, evol_sigma_t2 = 0.01*sigma_et^2, D = D, loc_obs = loc_obs, chol0)

  # Compute the evolution errors:
  omega = diff(mu, differences = D)

  # And the initial states:
  mu0 = as.matrix(mu[1:D,])

  # Initialize the evolution error variance paramters:
  evolParams = initEvolParams(omega, evol_error = evol_error)

  # Initial variance parameters:
  evolParams0 = initEvol0(mu0)

  # SV parameters, if necessary:
  if(obsSV == "SV") {svParams = initSV(y - mu); sigma_et = svParams$sigma_wt}
  else if(obsSV == "ASV"){sParams = init_paramsASV(y-mu, evol_error = "HS", D = 2); sigma_et = exp(sParams$s_mu/2)}

  # For HS MCMC comparisons:
  # evolParams$dhs_phi = 0

  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('mu', mcmc_params)) || computeDIC) post_mu = array(NA, c(nsave, T))
  if(!is.na(match('yhat', mcmc_params))) post_yhat = array(NA, c(nsave, T))
  if(!is.na(match('obs_sigma_t2', mcmc_params)) || computeDIC) post_obs_sigma_t2 = array(NA, c(nsave, T))
  if(!is.na(match('evol_sigma_t2', mcmc_params))) post_evol_sigma_t2 = array(NA, c(nsave, T))
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi = numeric(nsave)
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean = numeric(nsave)
  post_loglike = numeric(nsave)

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  if(verbose){
    pb <- progress::progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                     total = nstot,
                                     complete = "=",   # Completion bar character
                                     incomplete = "-", # Incomplete bar character
                                     current = ">",    # Current bar character
                                     clear = FALSE,    # If TRUE, clears the bar when finish
                                     width = 100)      # Width of the progress bar
  }
  for(nsi in 1:nstot){
    if(verbose){
      if(nsi < 10){
        pb$tick()
      }
      else if(((nsi%%100) == 0)){
        pb$tick(100)
      }
    }
    # Impute missing values, if any:
    if(any.missing) y[is.missing] = mu[is.missing] + sigma_et[is.missing]*rnorm(length(is.missing))

    # Sample the states:
    mu = sampleBTF(y, obs_sigma_t2 = sigma_et^2,
                   evol_sigma_t2 = c(evolParams0$sigma_w0^2,
                                     evolParams$sigma_wt^2),
                   D = D, loc_obs, chol0)

    # Compute the evolution errors:
    omega = diff(mu, differences = D)

    # And the initial states:
    mu0 = as.matrix(mu[1:D,])

    # Sample the initial variance parameters:
    evolParams0 = sampleEvol0(mu0, evolParams0, A = 1)

    # Sample the (observation and evolution) variances and associated parameters:
    if(obsSV == "SV"){
      # Evolution error variance + params:
      evolParams = sampleEvolParams(omega, evolParams, 1/sqrt(T), evol_error, loc)

      # Observation error variance + params:
      svParams = sampleSVparams(omega = y - mu, svParams = svParams)
      sigma_et = svParams$sigma_wt

    }else if(obsSV == "const"){
      # Evolution error variance + params:
      evolParams = sampleEvolParams(omega, evolParams, sigma_e/sqrt(T), evol_error, loc)

      if(evol_error == 'DHS') {
        sigma_e = uni.slice(sigma_e, g = function(x){
          -(T+2)*log(x) - 0.5*sum((y - mu)^2, na.rm=TRUE)/x^2 - log(1 + (sqrt(T)*exp(evolParams$dhs_mean0/2)/x)^2)
        }, lower = 0, upper = Inf)[1]
      }
      if(evol_error == 'HS') sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2 + length(evolParams$xiLambda), rate = sum((y - mu)^2, na.rm=TRUE)/2 + T*sum(evolParams$xiLambda)))
      if(evol_error == 'BL') sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2 + length(evolParams$tau_j)/2, rate = sum((y - mu)^2, na.rm=TRUE)/2 + T*sum((omega/evolParams$tau_j)^2)/2))
      if((evol_error == 'NIG') || (evol_error == 'SV')) sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2, rate = sum((y - mu)^2, na.rm=TRUE)/2))

      # Replicate for coding convenience:
      sigma_et = rep(sigma_e, T)
    }else if(obsSV == "ASV"){
      evolParams = sampleEvolParams(omega, evolParams, 1/sqrt(T), evol_error, loc)

      # Observation error variance + params:
      sParams = fit_paramsASV(y-mu,sParams,evol_error = "HS", D = 2)
      sigma_et = exp(sParams$s_mu/2)
    }else{
      stop('obsSV has to be one of const, SV, or ASV')
    }

    # Store the MCMC output:
    if(nsi > nburn){
      # Increment the skip counter:
      skipcount = skipcount + 1

      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1

        # Save the MCMC samples:
        if(!is.na(match('mu', mcmc_params)) || computeDIC) post_mu[isave,] = mu
        if(!is.na(match('yhat', mcmc_params))) post_yhat[isave,] = mu + sigma_et*rnorm(T)
        if(!is.na(match('obs_sigma_t2', mcmc_params)) || computeDIC) post_obs_sigma_t2[isave,] = sigma_et^2
        if(!is.na(match('evol_sigma_t2', mcmc_params))) post_evol_sigma_t2[isave,] = c(evolParams0$sigma_w0^2, evolParams$sigma_wt^2)
        if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi[isave] = evolParams$dhs_phi
        if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean[isave] = evolParams$dhs_mean
        post_loglike[isave] = sum(dnorm(y, mean = mu, sd = sigma_et, log = TRUE))

        # And reset the skip counter:
        skipcount = 0
      }
    }
  }

  if(!is.na(match('mu', mcmc_params))) mcmc_output$mu = post_mu
  if(!is.na(match('yhat', mcmc_params))) mcmc_output$yhat = post_yhat
  if(!is.na(match('obs_sigma_t2', mcmc_params))) mcmc_output$obs_sigma_t2 = post_obs_sigma_t2
  if(!is.na(match('evol_sigma_t2', mcmc_params))) mcmc_output$evol_sigma_t2 = post_evol_sigma_t2
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_phi = post_dhs_phi
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_mean = post_dhs_mean

  # Also include the log-likelihood:
  mcmc_output$loglike = post_loglike

  if(computeDIC){
    # Log-likelihood evaluated at posterior means:
    loglike_hat = sum(dnorm(y,
                            mean = colMeans(post_mu),
                            sd = colMeans(sqrt(post_obs_sigma_t2)),
                            log = TRUE))

    # Effective number of parameters (Note: two options)
    p_d = c(2*(loglike_hat - mean(post_loglike)),
            2*var(post_loglike))
    # DIC:
    DIC = -2*loglike_hat + 2*p_d

    # Store the DIC and the effective number of parameters (p_d)
    mcmc_output$DIC = DIC; mcmc_output$p_d = p_d
  }

  return (mcmc_output);
}
#' MCMC Sampler for Bayesian Trend Filtering: D = 0
#'
#' Run the MCMC for Bayesian trend filtering with a penalty on the conditional expectation.
#' The penalty is determined by the prior on the evolution errors, which include:
#' \itemize{
#' \item the dynamic horseshoe prior ('DHS');
#' \item the static horseshoe prior ('HS');
#' \item the Bayesian lasso ('BL');
#' \item the normal stochastic volatility model ('SV');
#' \item the normal-inverse-gamma prior ('NIG').
#' }
#' In each case, the evolution error is a scale mixture of Gaussians.
#' Sampling is accomplished with a (parameter-expanded) Gibbs sampler,
#' mostly relying on a dynamic linear model representation.

#' @param y the \code{T x 1} vector of time series observations
#' @param evol_error the evolution error distribution; must be one of
#' 'DHS' (dynamic horseshoe prior), 'HS' (horseshoe prior), 'BL' (Bayesian lasso), or 'NIG' (normal-inverse-gamma prior)
#' @param obsSV Options for modeling the error variance. It must be one of the following:
#' \itemize{
#' \item const: Constant error variance for all time points.
#' \item SV: Stochastic Volatility model.
#' \item ASV: Adaptive Stochastic Volatility model.
#' }
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burnin)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param mcmc_params named list of parameters for which we store the MCMC output;
#' must be one or more of:
#' \itemize{
#' \item "mu" (conditional mean)
#' \item "yhat" (posterior predictive distribution)
#' \item "evol_sigma_t2" (evolution error variance)
#' \item "obs_sigma_t2" (observation error variance)
#' \item "dhs_phi" (DHS AR(1) coefficient)
#' \item "dhs_mean" (DHS AR(1) unconditional mean)
#' }
#' @param computeDIC logical; if TRUE, compute the deviance information criterion \code{DIC}
#' and the effective number of parameters \code{p_d}
#' @param verbose logical; should R report extra information on progress?
#'
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#'
#' @note The data \code{y} may contain NAs, which will be treated with a simple imputation scheme
#' via an additional Gibbs sampling step. In general, rescaling \code{y} to have unit standard
#' deviation is recommended to avoid numerical issues.
btf0 = function(y, evol_error = 'DHS', obsSV = "const",
                nsave = 1000, nburn = 1000, nskip = 4,
                mcmc_params = list("mu", "yhat","evol_sigma_t2", "obs_sigma_t2", "dhs_phi", "dhs_mean"),
                computeDIC = TRUE,
                verbose = TRUE){

  # Time points (in [0,1])
  T = length(y); t01 = seq(0, 1, length.out=T);
  loc = t_create_loc(length(y), 0)

  # Begin by checking for missing values, then imputing (for initialization)
  is.missing = which(is.na(y)); any.missing = (length(is.missing) > 0)

  # Impute the active "data"
  y = approxfun(t01, y, rule = 2)(t01)

  # Initial SD (implicitly assumes a constant mean)
  sigma_e = sd(y, na.rm=TRUE); sigma_et = rep(sigma_e, T)

  # Initialize the conditional mean, mu, via sampling:
  mu = sampleBTF(y, obs_sigma_t2 = sigma_et^2, evol_sigma_t2 = 0.01*sigma_et^2, D = 0)

  # Initialize the evolution error variance paramters:
  evolParams = initEvolParams(omega = mu, evol_error = evol_error)

  # SV parameters, if necessary:
  if(obsSV == "SV") {svParams = initSV(y - mu); sigma_et = svParams$sigma_wt}
  else if(obsSV == "ASV"){sParams = init_paramsASV(y-mu, evol_error = "HS", D = 2); sigma_et = exp(sParams$s_mu/2)}

  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('mu', mcmc_params)) || computeDIC) post_mu = array(NA, c(nsave, T))
  if(!is.na(match('yhat', mcmc_params))) post_yhat = array(NA, c(nsave, T))
  if(!is.na(match('obs_sigma_t2', mcmc_params)) || computeDIC) post_obs_sigma_t2 = array(NA, c(nsave, T))
  if(!is.na(match('evol_sigma_t2', mcmc_params))) post_evol_sigma_t2 = array(NA, c(nsave, T))
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi = numeric(nsave)
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean = numeric(nsave)
  post_loglike = numeric(nsave)

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  if(verbose){
    pb <- progress::progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                     total = nstot,
                                     complete = "=",   # Completion bar character
                                     incomplete = "-", # Incomplete bar character
                                     current = ">",    # Current bar character
                                     clear = FALSE,    # If TRUE, clears the bar when finish
                                     width = 100)      # Width of the progress bar
  }
  for(nsi in 1:nstot){
    if(verbose){
      if(nsi < 10){
        pb$tick()
      }
      else if(((nsi%%100) == 0)){
        pb$tick(100)
      }
    }
    # Impute missing values, if any:
    if(any.missing) y[is.missing] = mu[is.missing] + sigma_et[is.missing]*rnorm(length(is.missing))

    # Sample the states:
    mu = sampleBTF(y, obs_sigma_t2 = sigma_et^2, evol_sigma_t2 = evolParams$sigma_wt^2, D = 0)

    # Sample the (observation and evolution) variances and associated parameters:
    if(obsSV == "SV"){
      # Evolution error variance + params:
      evolParams = sampleEvolParams(mu, evolParams, 1/sqrt(T), evol_error, loc)

      # Observation error variance + params:
      svParams = sampleSVparams(y - mu, svParams = svParams)
      sigma_et = svParams$sigma_wt

    }else if(obsSV == "const"){
      # Evolution error variance + params:
      evolParams = sampleEvolParams(mu, evolParams, sigma_e/sqrt(T), evol_error, loc)

      # Sample the observation error SD:
      if(evol_error == 'DHS') {
        sigma_e = uni.slice(sigma_e, g = function(x){
          -(T+2)*log(x) - 0.5*sum((y - mu)^2, na.rm=TRUE)/x^2 - log(1 + (sqrt(T)*exp(evolParams$dhs_mean0/2)/x)^2)
        }, lower = 0, upper = Inf)[1]
      }
      if(evol_error == 'HS') sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2 + length(evolParams$xiLambda),
                                                     rate = sum((y - mu)^2, na.rm=TRUE)/2 + T*sum(evolParams$xiLambda)))
      if(evol_error == 'BL') sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2 + length(evolParams$tau_j)/2,
                                                     rate = sum((y - mu)^2, na.rm=TRUE)/2 + T*sum((mu/evolParams$tau_j)^2)/2))
      if((evol_error == 'NIG') || (evol_error == 'SV')) sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2,
                                                                                rate = sum((y - mu)^2, na.rm=TRUE)/2))


      # Replicate for coding convenience:
      sigma_et = rep(sigma_e, T)
    }else if(obsSV == "ASV"){
      evolParams = sampleEvolParams(mu, evolParams, 1/sqrt(T), evol_error, loc)

      # Observation error variance + params:
      sParams = fit_paramsASV(y-mu,sParams,evol_error = "HS", D = 2)
      sigma_et = exp(sParams$s_mu/2)
    }else{
      stop('obsSV has to be one of const, SV, or ASV')
    }

    # Store the MCMC output:
    if(nsi > nburn){
      # Increment the skip counter:
      skipcount = skipcount + 1

      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1

        # Save the MCMC samples:
        if(!is.na(match('mu', mcmc_params)) || computeDIC) post_mu[isave,] = mu
        if(!is.na(match('yhat', mcmc_params))) post_yhat[isave,] = mu + sigma_et*rnorm(T)
        if(!is.na(match('obs_sigma_t2', mcmc_params)) || computeDIC) post_obs_sigma_t2[isave,] = sigma_et^2
        if(!is.na(match('evol_sigma_t2', mcmc_params))) post_evol_sigma_t2[isave,] =  evolParams$sigma_wt^2
        if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi[isave] = evolParams$dhs_phi
        if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean[isave] = evolParams$dhs_mean
        post_loglike[isave] = sum(dnorm(y, mean = mu, sd = sigma_et, log = TRUE))

        # And reset the skip counter:
        skipcount = 0
      }
    }
  }

  if(!is.na(match('mu', mcmc_params))) mcmc_output$mu = post_mu
  if(!is.na(match('yhat', mcmc_params))) mcmc_output$yhat = post_yhat
  if(!is.na(match('obs_sigma_t2', mcmc_params))) mcmc_output$obs_sigma_t2 = post_obs_sigma_t2
  if(!is.na(match('evol_sigma_t2', mcmc_params))) mcmc_output$evol_sigma_t2 = post_evol_sigma_t2
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_phi = post_dhs_phi
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_mean = post_dhs_mean

  # Also include the log-likelihood:
  mcmc_output$loglike = post_loglike

  if(computeDIC){
    # Log-likelihood evaluated at posterior means:
    loglike_hat = sum(dnorm(y,
                            mean = colMeans(post_mu),
                            sd = colMeans(sqrt(post_obs_sigma_t2)),
                            log = TRUE))

    # Effective number of parameters (Note: two options)
    p_d = c(2*(loglike_hat - mean(post_loglike)),
            2*var(post_loglike))
    # DIC:
    DIC = -2*loglike_hat + 2*p_d

    # Store the DIC and the effective number of parameters (p_d)
    mcmc_output$DIC = DIC; mcmc_output$p_d = p_d
  }

  return (mcmc_output);
}

#' Run the MCMC for sparse Bayesian trend filtering
#'
#' Sparse Bayesian trend filtering has two penalties:
#' (1) a penalty on the first (D = 1) or second (D = 2) differences of the conditional expectation and
#' (2) a penalty on the conditional expectation, i.e., shrinkage to zero.
#'
#' Each penalty is determined by a prior, which include:
#' \itemize{
#' \item the dynamic horseshoe prior ('DHS');
#' \item the static horseshoe prior ('HS');
#' \item the Bayesian lasso ('BL');
#' \item the normal stochastic volatility model ('SV');
#' \item the normal-inverse-gamma prior ('NIG').
#' }
#' In each case, the prior is a scale mixture of Gaussians.
#' Sampling is accomplished with a (parameter-expanded) Gibbs sampler,
#' mostly relying on a dynamic linear model representation.

#' @param y the \code{T x 1} vector of time series observations
#' @param evol_error the evolution error distribution; must be one of
#' 'DHS' (dynamic horseshoe prior), 'HS' (horseshoe prior), 'BL' (Bayesian lasso), or 'NIG' (normal-inverse-gamma prior)
#' @param zero_error the shrinkage-to-zero distribution; must be one of
#' 'DHS' (dynamic horseshoe prior), 'HS' (horseshoe prior), 'BL' (Bayesian lasso), or 'NIG' (normal-inverse-gamma prior)
#' @param D degree of differencing (D = 1, or D = 2)
#' @param obsSV Options for modeling the error variance. It must be one of the following:
#' \itemize{
#' \item const: Constant error variance for all time points.
#' \item SV: Stochastic Volatility model.
#' \item ASV: Adaptive Stochastic Volatility model.
#' }
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burin-in)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param mcmc_params named list of parameters for which we store the MCMC output;
#' must be one or more of:
#' \itemize{
#' \item "mu" (conditional mean)
#' \item "yhat" (posterior predictive distribution)
#' \item "evol_sigma_t2" (evolution error variance)
#' \item "zero_sigma_t2" (shrink-to-zero error variance)
#' \item "obs_sigma_t2" (observation error variance)
#' \item "dhs_phi" (DHS AR(1) coefficient for evolution error)
#' \item "dhs_mean" (DHS AR(1) unconditional mean for evolution error)
#' \item "dhs_phi_zero" (DHS AR(1) coefficient for shrink-to-zero error)
#' \item "dhs_mean_zero" (DHS AR(1) unconditional mean for shrink-to-zero error)
#' }
#' @param computeDIC logical; if TRUE, compute the deviance information criterion \code{DIC}
#' and the effective number of parameters \code{p_d}
#' @param verbose logical; should R report extra information on progress?
#'
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#'
#' @note The data \code{y} may contain NAs, which will be treated with a simple imputation scheme
#' via an additional Gibbs sampling step. In general, rescaling \code{y} to have unit standard
#' deviation is recommended to avoid numerical issues.
#'
#' @examples
#' # Example 1: Bumps Data
#' y = make.signal(name = "bumps", n = 128, snr = 7)
#'
#' out = btf_sparse(y)
#' plot_fitted(y, mu = colMeans(out$mu), postY = out$yhat)
#'
#' # Example 2: Doppler Data; longer series, more noise
#' y = make.signal(name = "doppler", n = 500, snr = 7)
#'
#' out = btf_sparse(y)
#' plot_fitted(y, mu = colMeans(out$mu), postY = out$yhat)
#'
#' # And examine the AR(1) parameters for the log-volatility w/ traceplots:
#' plot(as.ts(out$dhs_phi)) # AR(1) coefficient
#' plot(as.ts(out$dhs_mean)) # Unconditional mean
#'
#' # Example 3: Blocks data (locally constant)
#' y = make.signal(name = "blocks", n = 1000, snr = 3)
#'
#' out = btf_sparse(y, D = 1) # try D = 1 to approximate the locally constant behavior
#' plot_fitted(y, mu = colMeans(out$mu), postY = out$yhat)
#'
#'
#' @import spam progress
#' @export
btf_sparse = function(y, evol_error = 'DHS', zero_error = 'DHS', D = 2, obsSV = "const",
                      nsave = 1000, nburn = 1000, nskip = 4,
                      mcmc_params = list("mu", "yhat","evol_sigma_t2", "obs_sigma_t2", "zero_sigma_t2", "dhs_phi", "dhs_mean","dhs_phi_zero", "dhs_mean_zero"),
                      computeDIC = TRUE,
                      verbose = TRUE){

  # Check differencing:
  if((D != 1) && (D != 2)) stop('btf_sparse() requires D = 1 or D = 2')

  # Convert to upper case:
  evol_error = toupper(evol_error)
  zero_error = toupper(zero_error)

  # Time points (in [0,1])
  T = length(y); t01 = seq(0, 1, length.out=T);

  # Begin by checking for missing values, then imputing (for initialization)
  is.missing = which(is.na(y)); any.missing = (length(is.missing) > 0)

  # Impute the active "data"
  y = approxfun(t01, y, rule = 2)(t01)

  # Initial SD (implicitly assumes a constant mean)
  sigma_e = sd(y, na.rm=TRUE); sigma_et = rep(sigma_e, T)

  # Compute the Cholesky term (uses random variances for a more conservative sparsity pattern)
  chol0 = initChol.spam(T = T, D = D)

  # Initialize the conditional mean, mu, via sampling:
  mu = sampleBTF(y, obs_sigma_t2 = sigma_et^2, evol_sigma_t2 = 0.01*sigma_et^2, D = D, chol0 = chol0)

  # Compute the evolution errors:
  omega = diff(mu, differences = D)

  # And the initial states:
  mu0 = as.matrix(mu[1:D,])

  # Initialize the evolution error variance paramters:
  evolParams = initEvolParams(omega, evol_error = evol_error)

  # Initial variance parameters:
  evolParams0 = initEvol0(mu0)

  # Initialize the shrink-to-zero parameters:
  zeroParams = initEvolParams(omega = mu, evol_error = zero_error)

  # For HS MCMC comparisons:
  # evolParams$dhs_phi = 0
  if(obsSV == "SV") {svParams = initSV(y - mu); sigma_et = svParams$sigma_wt}
  else if(obsSV == "ASV"){sParams = init_paramsASV(y-mu, evol_error = "HS", D = 2); sigma_et = exp(sParams$s_mu/2)}

  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('mu', mcmc_params)) || computeDIC) post_mu = array(NA, c(nsave, T))
  if(!is.na(match('yhat', mcmc_params))) post_yhat = array(NA, c(nsave, T))
  if(!is.na(match('obs_sigma_t2', mcmc_params)) || computeDIC) post_obs_sigma_t2 = array(NA, c(nsave, T))
  if(!is.na(match('evol_sigma_t2', mcmc_params))) post_evol_sigma_t2 = array(NA, c(nsave, T))
  if(!is.na(match('zero_sigma_t2', mcmc_params))) post_zero_sigma_t2 = array(NA, c(nsave, T))
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi = numeric(nsave)
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean = numeric(nsave)
  if(!is.na(match('dhs_phi_zero', mcmc_params)) && zero_error == "DHS") post_dhs_phi_zero = numeric(nsave)
  if(!is.na(match('dhs_mean_zero', mcmc_params)) && zero_error == "DHS") post_dhs_mean_zero = numeric(nsave)
  post_loglike = numeric(nsave)

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  if(verbose){
    pb <- progress::progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                     total = nstot,
                                     complete = "=",   # Completion bar character
                                     incomplete = "-", # Incomplete bar character
                                     current = ">",    # Current bar character
                                     clear = FALSE,    # If TRUE, clears the bar when finish
                                     width = 100)      # Width of the progress bar
  }
  for(nsi in 1:nstot){
    if(verbose){
      if(nsi < 10){
        pb$tick()
      }
      else if(((nsi%%100) == 0)){
        pb$tick(100)
      }
    }
    # Impute missing values, if any:
    if(any.missing) y[is.missing] = mu[is.missing] + sigma_et[is.missing]*rnorm(length(is.missing))


    # Sample the observation SD:
    if(obsSV == "SV"){
      # Observation error variance + params:
      svParams = sampleSVparams(omega = y - mu, svParams = svParams)
      sigma_et = svParams$sigma_wt
    }else if(obsSV == "const"){
      sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2, rate = sum((y - mu)^2, na.rm=TRUE)/2))
      sigma_et = rep(sigma_e, T) # For coding convenience
    }else if(obsSV == "ASV"){
      # Observation error variance + params:
      sParams = fit_paramsASV(y-mu,sParams,evol_error = "HS", D = 2)
      sigma_et = exp(sParams$s_mu/2)
    }else{
      stop('obsSV has to be one of const, SV, or ASV')
    }
    # Sample the states:
    #mu = sampleBTF(y, obs_sigma_t2 = sigma_et^2, evol_sigma_t2 = c(evolParams0$sigma_w0^2, evolParams$sigma_wt^2), D = D, chol0 = chol0)
    mu = sampleBTF_sparse(y,
                          obs_sigma_t2 = sigma_et^2,
                          evol_sigma_t2 = c(evolParams0$sigma_w0^2, evolParams$sigma_wt^2),
                          zero_sigma_t2 = zeroParams$sigma_wt^2,
                          D = D, chol0 = chol0)

    # Compute the evolution errors:
    omega = diff(mu, differences = D)

    # And the initial states:
    mu0 = as.matrix(mu[1:D,])

    # Sample the initial variance parameters:
    evolParams0 = sampleEvol0(mu0, evolParams0, A = 1)

    # Sample the evolution variance and associated parameters:
    evolParams = sampleEvolParams(omega = omega, evolParams = evolParams, sigma_e = 1/sqrt(T), evol_error = evol_error)

    # Sample the shrink-to-zero variance and associated parameters:
    zeroParams = sampleEvolParams(omega = mu, evolParams = zeroParams, sigma_e = 1/sqrt(T), evol_error = zero_error)

    # Store the MCMC output:
    if(nsi > nburn){
      # Increment the skip counter:
      skipcount = skipcount + 1

      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1

        # Save the MCMC samples:
        if(!is.na(match('mu', mcmc_params)) || computeDIC) post_mu[isave,] = mu
        if(!is.na(match('yhat', mcmc_params))) post_yhat[isave,] = mu + sigma_et*rnorm(T)
        if(!is.na(match('obs_sigma_t2', mcmc_params)) || computeDIC) post_obs_sigma_t2[isave,] = sigma_et^2
        if(!is.na(match('evol_sigma_t2', mcmc_params))) post_evol_sigma_t2[isave,] = c(evolParams0$sigma_w0^2, evolParams$sigma_wt^2)
        if(!is.na(match('zero_sigma_t2', mcmc_params))) post_zero_sigma_t2[isave,] = zeroParams$sigma_wt^2
        if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi[isave] = evolParams$dhs_phi
        if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean[isave] = evolParams$dhs_mean
        if(!is.na(match('dhs_phi_zero', mcmc_params)) && zero_error == "DHS") post_dhs_phi_zero[isave] = zeroParams$dhs_phi
        if(!is.na(match('dhs_mean_zero', mcmc_params)) && zero_error == "DHS") post_dhs_mean_zero[isave] = zeroParams$dhs_mean
        post_loglike[isave] = sum(dnorm(y, mean = mu, sd = sigma_et, log = TRUE))

        # And reset the skip counter:
        skipcount = 0
      }
    }
  }

  if(!is.na(match('mu', mcmc_params))) mcmc_output$mu = post_mu
  if(!is.na(match('yhat', mcmc_params))) mcmc_output$yhat = post_yhat
  if(!is.na(match('obs_sigma_t2', mcmc_params))) mcmc_output$obs_sigma_t2 = post_obs_sigma_t2
  if(!is.na(match('evol_sigma_t2', mcmc_params))) mcmc_output$evol_sigma_t2 = post_evol_sigma_t2
  if(!is.na(match('zero_sigma_t2', mcmc_params))) mcmc_output$zero_sigma_t2 = post_zero_sigma_t2
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_phi = post_dhs_phi
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_mean = post_dhs_mean
  if(!is.na(match('dhs_phi_zero', mcmc_params)) && zero_error == "DHS") mcmc_output$dhs_phi_zero = post_dhs_phi_zero
  if(!is.na(match('dhs_mean_zero', mcmc_params)) && zero_error == "DHS") mcmc_output$dhs_mean_zero = post_dhs_mean_zero

  # Also include the log-likelihood:
  mcmc_output$loglike = post_loglike

  if(computeDIC){
    # Log-likelihood evaluated at posterior means:
    loglike_hat = sum(dnorm(y,
                            mean = colMeans(post_mu),
                            sd = colMeans(sqrt(post_obs_sigma_t2)),
                            log = TRUE))

    # Effective number of parameters (Note: two options)
    p_d = c(2*(loglike_hat - mean(post_loglike)),
            2*var(post_loglike))
    # DIC:
    DIC = -2*loglike_hat + 2*p_d

    # Store the DIC and the effective number of parameters (p_d)
    mcmc_output$DIC = DIC; mcmc_output$p_d = p_d
  }

  return (mcmc_output);
}
#----------------------------------------------------------------------------
#' MCMC Sampler for Bayesian Trend Filtering: Regression
#'
#' Run the MCMC for Bayesian trend filtering regression with a penalty on
#' first (D=1) or second (D=2) differences of each dynamic regression coefficient.
#' The penalty is determined by the prior on the evolution errors, which include:
#' \itemize{
#' \item the dynamic horseshoe prior ('DHS');
#' \item the static horseshoe prior ('HS');
#' \item the Bayesian lasso ('BL');
#' \item the normal stochastic volatility model ('SV');
#' \item the normal-inverse-gamma prior ('NIG').
#' }
#' In each case, the evolution error is a scale mixture of Gaussians.
#' Sampling is accomplished with a (parameter-expanded) Gibbs sampler,
#' mostly relying on a dynamic linear model representation.
#'
#' @param y the \code{T x 1} vector of time series observations
#' @param X the \code{T x p} matrix of time series predictors
#' @param evol_error the evolution error distribution; must be one of
#' 'DHS' (dynamic horseshoe prior), 'HS' (horseshoe prior), 'BL' (Bayesian lasso), or 'NIG' (normal-inverse-gamma prior)
#' @param D degree of differencing (D = 1 or D = 2)
#' @param obsSV Options for modeling the error variance. It must be one of the following:
#' \itemize{
#' \item const: Constant error variance for all time points.
#' \item SV: Stochastic Volatility model.
#' \item ASV: Adaptive Stochastic Volatility model.
#' }
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burin-in)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param mcmc_params named list of parameters for which we store the MCMC output;
#' must be one or more of:
#' \itemize{
#' \item "mu" (conditional mean)
#' \item "yhat" (posterior predictive distribution)
#' \item "beta" (dynamic regression coefficients)
#' \item "evol_sigma_t2" (evolution error variance)
#' \item "obs_sigma_t2" (observation error variance)
#' \item "dhs_phi" (DHS AR(1) coefficient)
#' \item "dhs_mean" (DHS AR(1) unconditional mean)
#' }
#' @param use_backfitting logical; if TRUE, use backfitting to sample the predictors j=1,...,p
#' (faster, but usually less MCMC efficient)
#' @param computeDIC logical; if TRUE, compute the deviance information criterion \code{DIC}
#' and the effective number of parameters \code{p_d}
#' @param verbose logical; should R report extra information on progress?
#'
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#'
#' @note The data \code{y} may contain NAs, which will be treated with a simple imputation scheme
#' via an additional Gibbs sampling step. In general, rescaling \code{y} to have unit standard
#' deviation is recommended to avoid numerical issues.
#'
#' @examples
#'
#' # Example 1: all signals
#' simdata = simRegression(T = 200, p = 5, p_0 = 0)
#' y = simdata$y; X = simdata$X
#' out = btf_reg(y, X)
#' for(j in 1:ncol(X))
#'  plot_fitted(rep(0, length(y)),
#'              mu = colMeans(out$beta[,,j]),
#'              postY = out$beta[,,j],
#'              y_true = simdata$beta_true[,j])
#'
#'
#' # Example 2: some noise, longer series
#' simdata = simRegression(T = 500, p = 10, p_0 = 5)
#' y = simdata$y; X = simdata$X
#' out = btf_reg(y, X, nsave = 1000, nskip = 0) # Short MCMC run for a quick example
#' for(j in 1:ncol(X))
#'   plot_fitted(rep(0, length(y)),
#'               mu = colMeans(out$beta[,,j]),
#'               postY = out$beta[,,j],
#'               y_true = simdata$beta_true[,j])
#'
#'
#' @import spam progress
#' @export
btf_reg = function(y, X = NULL, evol_error = 'DHS', D = 1, obsSV = "const",
                   nsave = 1000, nburn = 1000, nskip = 4,
                   mcmc_params = list("mu", "yhat","beta","evol_sigma_t2", "obs_sigma_t2", "dhs_phi", "dhs_mean"),
                   use_backfitting = FALSE,
                   computeDIC = TRUE,
                   verbose = TRUE){

  # Convert to upper case:
  evol_error = toupper(evol_error)

  # If no predictors are specified, just call the univariate trend filtering model: btf()
  if(is.null(X)) {
    mcmc_params[match("beta", mcmc_params)] = NULL # Remove "beta", since it does not apply for btf()
    return(btf(y = y, evol_error = evol_error, D = D, nsave = nsave, nburn = nburn,
               nskip = nskip, mcmc_params = mcmc_params, computeDIC = computeDIC, verbose = verbose))
  }

  # Time points (in [0,1])
  T = length(y); t01 = seq(0, 1, length.out=T);

  # Begin by checking for missing values, then imputing (for initialization)
  is.missing = which(is.na(y)); any.missing = (length(is.missing) > 0)

  # Impute the active "data"
  if(any.missing) y[is.missing] = mean(y, na.rm=TRUE)

  # Obtain XtX
  XtX = build_XtX(X);

  # Number of predictors:
  p = ncol(X)

  # Initial SD (implicitly assumes a constant mean)
  sigma_e = sd(y, na.rm=TRUE); sigma_et = rep(sigma_e, T)

  # Compute the Cholesky term: use random variances for a more conservative sparsity pattern
  chol0 = initCholReg.spam(obs_sigma_t2 = abs(rnorm(T)),
                           evol_sigma_t2 = matrix(abs(rnorm(T*p)), nrow = T),
                           XtX = XtX, D = D)

  # Initialize the dynamic regression coefficients, beta, via sampling:
  beta = sampleBTF_reg(y, X, obs_sigma_t2 = sigma_et^2, evol_sigma_t2 = matrix(0.01*sigma_et^2, nrow = T, ncol = p), XtX = XtX, D = D, chol0 = chol0)

  # Conditional mean:
  mu = rowSums(X*beta)

  # Compute the evolution errors:
  omega = diff(beta, differences = D)

  # And the initial states:
  beta0 = matrix(beta[1:D,], nrow = D)

  # Initialize the evolution error variance paramters:
  evolParams = initEvolParams(omega, evol_error = evol_error)

  # Initial variance parameters:
  evolParams0 = initEvol0(beta0, commonSD = FALSE)

  # SV parameters, if necessary:
  if(obsSV == "SV") {svParams = initSV(y - mu); sigma_et = svParams$sigma_wt}
  else if(obsSV == "ASV"){sParams = init_paramsASV(y-mu, evol_error = "HS", D = 2); sigma_et = exp(sParams$s_mu/2)}

  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('mu', mcmc_params)) || computeDIC) post_mu = array(NA, c(nsave, T))
  if(!is.na(match('yhat', mcmc_params))) post_yhat = array(NA, c(nsave, T))
  if(!is.na(match('beta', mcmc_params))) post_beta = array(NA, c(nsave, T, p))
  if(!is.na(match('obs_sigma_t2', mcmc_params)) || computeDIC) post_obs_sigma_t2 = array(NA, c(nsave, T))
  if(!is.na(match('evol_sigma_t2', mcmc_params))) post_evol_sigma_t2 = array(NA, c(nsave, T, p))
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi = array(NA, c(nsave, p))
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean = array(NA, c(nsave, p))
  post_loglike = numeric(nsave)

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  if(verbose){
    pb <- progress::progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                     total = nstot,
                                     complete = "=",   # Completion bar character
                                     incomplete = "-", # Incomplete bar character
                                     current = ">",    # Current bar character
                                     clear = FALSE,    # If TRUE, clears the bar when finish
                                     width = 100)      # Width of the progress bar
  }
  for(nsi in 1:nstot){
    if(verbose){
      if(nsi < 10){
        pb$tick()
      }
      else if(((nsi%%100) == 0)){
        pb$tick(100)
      }
    }
    # Impute missing values, if any:
    if(any.missing) y[is.missing] = mu[is.missing] + sigma_et[is.missing]*rnorm(length(is.missing))

    # Sample the dynamic regression coefficients, beta:
    # Backfitting is faster, but likely less efficient
    if(use_backfitting){
      beta = sampleBTF_reg_backfit(y, X, beta, obs_sigma_t2 = sigma_et^2, evol_sigma_t2 = rbind(matrix(evolParams0$sigma_w0^2, nrow = D), evolParams$sigma_wt^2), D = D)
    } else beta = sampleBTF_reg(y, X, obs_sigma_t2 = sigma_et^2, evol_sigma_t2 = rbind(matrix(evolParams0$sigma_w0^2, nrow = D), evolParams$sigma_wt^2), XtX = XtX, D = D, chol0 = chol0)

    # Conditional mean:
    mu = rowSums(X*beta)

    # Compute the evolution errors:
    omega = diff(beta, differences = D)

    # And the initial states:
    beta0 = matrix(beta[1:D,], nrow = D)

    # Sample the initial variance parameters:
    evolParams0 = sampleEvol0(beta0, evolParams0, A = 1, commonSD = FALSE)

    # Sample the (observation and evolution) variances and associated parameters:
    if(obsSV == "SV"){
      # Evolution error variance + params:
      evolParams = sampleEvolParams(omega, evolParams, 1/sqrt(T*p), evol_error)

      # Observation error variance + params:
      svParams = sampleSVparams(omega = y - mu, svParams = svParams)
      sigma_et = svParams$sigma_wt

    }else if(obsSV == "const"){
      # Evolution error variance + params:
      evolParams = sampleEvolParams(omega, evolParams, sigma_e/sqrt(T*p), evol_error)

      # Sample the observation error SD:
      if(evol_error == 'DHS') {
        sigma_e = uni.slice(sigma_e, g = function(x){
          -(T+2)*log(x) - 0.5*sum((y - mu)^2, na.rm=TRUE)/x^2 - log(1 + (sqrt(T*p)*exp(evolParams$dhs_mean0/2)/x)^2)
        }, lower = 0, upper = Inf)[1]
      }
      #if(evol_error == 'HS') sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2 + length(evolParams$xiLambda), rate = sum((y - mu)^2, na.rm=TRUE)/2 + T*p*sum(evolParams$xiLambda)))
      if(evol_error == 'HS') sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2, rate = sum((y - mu)^2, na.rm=TRUE)/2))
      if(evol_error == 'BL') sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2 + length(evolParams$tau_j)/2, rate = sum((y - mu)^2, na.rm=TRUE)/2 + T*p*sum((omega/evolParams$tau_j)^2)/2))
      if((evol_error == 'NIG') || (evol_error == 'SV')) sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2, rate = sum((y - mu)^2, na.rm=TRUE)/2))

      # Replicate for coding convenience:
      sigma_et = rep(sigma_e, T)
    }else if(obsSV == "ASV"){
      evolParams = sampleEvolParams(omega, evolParams, 1/sqrt(T*p), evol_error)

      # Observation error variance + params:
      sParams = fit_paramsASV(y-mu,sParams,evol_error = "HS", D = 2)
      sigma_et = exp(sParams$s_mu/2)
    }else{
      stop('obsSV has to be one of const, SV, or ASV')
    }

    # Store the MCMC output:
    if(nsi > nburn){
      # Increment the skip counter:
      skipcount = skipcount + 1

      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1

        # Save the MCMC samples:
        if(!is.na(match('mu', mcmc_params)) || computeDIC) post_mu[isave,] = mu
        if(!is.na(match('yhat', mcmc_params))) post_yhat[isave,] = mu + sigma_et*rnorm(T)
        if(!is.na(match('beta', mcmc_params))) post_beta[isave,,] = beta
        if(!is.na(match('obs_sigma_t2', mcmc_params)) || computeDIC) post_obs_sigma_t2[isave,] = sigma_et^2
        if(!is.na(match('evol_sigma_t2', mcmc_params))) {
          if(D == 0){
            post_evol_sigma_t2[isave,,] = evolParams$sigma_wt^2
          } else post_evol_sigma_t2[isave,,] = rbind(matrix(evolParams0$sigma_w0^2, nrow = D), evolParams$sigma_wt^2)
        }
        if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi[isave,] = evolParams$dhs_phi
        if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean[isave,] = evolParams$dhs_mean
        post_loglike[isave] = sum(dnorm(y, mean = mu, sd = sigma_et, log = TRUE))

        # And reset the skip counter:
        skipcount = 0
      }
    }
  }

  if(!is.na(match('mu', mcmc_params))) mcmc_output$mu = post_mu
  if(!is.na(match('yhat', mcmc_params))) mcmc_output$yhat = post_yhat
  if(!is.na(match('beta', mcmc_params))) mcmc_output$beta = post_beta
  if(!is.na(match('obs_sigma_t2', mcmc_params))) mcmc_output$obs_sigma_t2 = post_obs_sigma_t2
  if(!is.na(match('evol_sigma_t2', mcmc_params))) mcmc_output$evol_sigma_t2 = post_evol_sigma_t2
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_phi = post_dhs_phi
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_mean = post_dhs_mean

  # Also include the log-likelihood:
  mcmc_output$loglike = post_loglike

  if(computeDIC){
    # Log-likelihood evaluated at posterior means:
    loglike_hat = sum(dnorm(y,
                            mean = colMeans(post_mu),
                            sd = colMeans(sqrt(post_obs_sigma_t2)),
                            log = TRUE))
    # Effective number of parameters (Note: two options)
    p_d = c(2*(loglike_hat - mean(post_loglike)),
            2*var(post_loglike))
    # DIC:
    DIC = -2*loglike_hat + 2*p_d

    # Store the DIC and the effective number of parameters (p_d)
    mcmc_output$DIC = DIC; mcmc_output$p_d = p_d
  }

  return (mcmc_output);
}
#----------------------------------------------------------------------------
#' MCMC Sampler for B-spline Bayesian Trend Filtering
#'
#' Run the MCMC for B-spline fitting with a Bayesian trend filtering model on the
#' coefficients, i.e., a penalty on zeroth (D=0), first (D=1), or second (D=2)
#' differences of the B-spline basis coefficients.
#' The penalty is determined by the prior on the evolution errors, which include:
#' \itemize{
#' \item the dynamic horseshoe prior ('DHS');
#' \item the static horseshoe prior ('HS');
#' \item the Bayesian lasso ('BL');
#' \item the normal stochastic volatility model ('SV');
#' \item the normal-inverse-gamma prior ('NIG').
#' }
#' In each case, the evolution error is a scale mixture of Gaussians.
#' Sampling is accomplished with a (parameter-expanded) Gibbs sampler,
#' mostly relying on a dynamic linear model representation.

#' @param y the \code{T x 1} vector of time series observations
#' @param x the \code{T x 1} vector of observation points; if NULL, assume equally spaced
#' @param num_knots the number of knots; if NULL, use the default of \code{max(20, min(ceiling(T/4), 150))}
#' @param evol_error the evolution error distribution; must be one of
#' 'DHS' (dynamic horseshoe prior), 'HS' (horseshoe prior), 'BL' (Bayesian lasso), or 'NIG' (normal-inverse-gamma prior)
#' @param D degree of differencing (D = 0, D = 1, or D = 2)
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burin-in)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param mcmc_params named list of parameters for which we store the MCMC output;
#' must be one or more of:
#' \itemize{
#' \item "mu" (conditional mean)
#' \item "beta" (B-spline basis coefficients)
#' \item "yhat" (posterior predictive distribution)
#' \item "evol_sigma_t2" (evolution error variance)
#' \item "obs_sigma_t2" (observation error variance)
#' \item "dhs_phi" (DHS AR(1) coefficient)
#' \item "dhs_mean" (DHS AR(1) unconditional mean)
#' }
#' @param computeDIC logical; if TRUE, compute the deviance information criterion \code{DIC}
#' and the effective number of parameters \code{p_d}
#' @param verbose logical; should R report extra information on progress?
#'
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#'
#' @note The data \code{y} may contain NAs, which will be treated with a simple imputation scheme
#' via an additional Gibbs sampling step. In general, rescaling \code{y} to have unit standard
#' deviation is recommended to avoid numerical issues.
#'
#' @note The primary advantages of \code{btf_bspline} over \code{btf} are
#' \enumerate{
#' \item Unequally-spaced points are handled automatically and
#' \item Computations are linear in the number of basis coefficients, which may be
#' substantially fewer than the number of time points.
#' }
#'
#' @examples
#'
#' # Example 1: Blocks data
#' simdata <- simUnivariate(signalName = "blocks", T = 1000, RSNR = 3, include_plot = FALSE)
#' y <- simdata$y
#' out <- btf_bspline(y, D = 1)
#' plot_fitted(y, mu = colMeans(out$mu), postY = out$yhat, y_true = simdata$y_true)
#'
#'
#' # Example 2: motorcycle data (unequally-spaced points)
#' data("mcycle", package = "MASS")
#' y <- scale(mcycle$accel) # Center and Scale for numerical stability
#' x <- mcycle$times
#' plot(x, y, xlab = 'Time (ms)', ylab='Acceleration (g)', main = 'Motorcycle Crash Data')
#' out <- btf_bspline(y = y, x = x)
#' plot_fitted(y, mu = colMeans(out$mu), postY = out$yhat, t01 = x)
#'
#'
#' @import fda progress
#' @export
btf_bspline = function(y, x = NULL, num_knots = NULL, evol_error = 'DHS', D = 2,
                       nsave = 1000, nburn = 1000, nskip = 4,
                       mcmc_params = list("mu", "yhat", "beta", "evol_sigma_t2", "obs_sigma_t2", "dhs_phi", "dhs_mean"),
                       computeDIC = TRUE,
                       verbose = TRUE){

  # Convert to upper case:
  evol_error = toupper(evol_error)

  # For D = 0, return special case:
  if(D == 0){
    return(btf_bspline0(y = y, x = x, num_knots = num_knots, evol_error = evol_error,
                        nsave = nsave, nburn = nburn, nskip = nskip,
                        mcmc_params = mcmc_params,
                        computeDIC = computeDIC, verbose = verbose))
  }

  # Length of time series
  T = length(y);

  # Observation points
  if(is.null(x)) x = seq(0, 1, length.out=T);

  # Rescale to (0,1):
  t01 = (x - min(x))/diff(range(x))

  # Compute B-spline basis matrix:
  if(is.null(num_knots)) num_knots = max(20, min(ceiling(T/4), 150))
  X = eval.basis(t01, create.bspline.basis(c(0,1), nbasis = num_knots))
  p = ncol(X)

  # In place of XtX, store the 4-bands of XtX:
  XtX = crossprod(X) #   XtX = Matrix(crossprod(X));
  XtX_bands = list(XtX_0 = diag(XtX),
                   XtX_1 = diag(XtX[-1,]),
                   XtX_2 = diag(XtX[-(1:2),]),
                   XtX_3 = diag(XtX[-(1:3),]))
  rm(XtX)

  # Begin by checking for missing values, then imputing (for initialization)
  is.missing = which(is.na(y)); any.missing = (length(is.missing) > 0)

  # Impute the active "data"
  y = approxfun(t01, y, rule = 2)(t01)

  # Initial SD (implicitly assumes a constant mean)
  sigma_e = sd(y, na.rm=TRUE)

  # Initialize the B-spline coefficients via sampling
  Xty = crossprod(X,y)
  beta = sampleBTF_bspline(y, X, obs_sigma2 = sigma_e^2, evol_sigma_t2 = rep(0.01*sigma_e^2, p), XtX_bands = XtX_bands, Xty = Xty, D = D)

  # And the conditional expectation:
  mu = as.numeric(X%*%beta)

  # Compute the evolution errors:
  omega = diff(beta, differences = D)

  # And the initial states:
  beta0 = matrix(beta[1:D,], nrow = D)

  # Initialize the evolution error variance paramters:
  evolParams = initEvolParams(omega, evol_error = evol_error)

  # Initial variance parameters:
  evolParams0 = initEvol0(beta0)

  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('mu', mcmc_params)) || computeDIC) post_mu = array(NA, c(nsave, T))
  if(!is.na(match('yhat', mcmc_params))) post_yhat = array(NA, c(nsave, T))
  if(!is.na(match('beta', mcmc_params))) post_beta = array(NA, c(nsave, p))
  if(!is.na(match('obs_sigma_t2', mcmc_params)) || computeDIC) post_obs_sigma_t2 = array(NA, c(nsave, T))
  if(!is.na(match('evol_sigma_t2', mcmc_params))) post_evol_sigma_t2 = array(NA, c(nsave, p))
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi = numeric(nsave)
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean = numeric(nsave)
  post_loglike = numeric(nsave)

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  if(verbose){
    pb <- progress::progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                     total = nstot,
                                     complete = "=",   # Completion bar character
                                     incomplete = "-", # Incomplete bar character
                                     current = ">",    # Current bar character
                                     clear = FALSE,    # If TRUE, clears the bar when finish
                                     width = 100)      # Width of the progress bar
  }
  for(nsi in 1:nstot){
    if(verbose){
      if(nsi < 10){
        pb$tick()
      }
      else if(((nsi%%100) == 0)){
        pb$tick(100)
      }
    }
    # Impute missing values, if any:
    if(any.missing) {
      y[is.missing] = mu[is.missing] + sigma_e*rnorm(length(is.missing))
      Xty = crossprod(X,y) # Recompute when missing values present
    }

    # Sample the dynamic regression coefficients, beta:
    beta = sampleBTF_bspline(y, X, obs_sigma2 = sigma_e^2, evol_sigma_t2 = rbind(matrix(evolParams0$sigma_w0^2, nrow = D), evolParams$sigma_wt^2), XtX_bands = XtX_bands, Xty = Xty, D = D)

    # And the conditional expectation:
    mu = as.numeric(X%*%beta)

    # Compute the evolution errors:
    omega = diff(beta, differences = D)

    # And the initial states:
    beta0 = matrix(beta[1:D,], nrow = D)

    # Sample the evolution error variance (and associated parameters):
    evolParams = sampleEvolParams(omega, evolParams, sigma_e/sqrt(p), evol_error)

    # Sample the observation error SD:
    if(evol_error == 'DHS') {
      sigma_e = uni.slice(sigma_e, g = function(x){
        -(T+2)*log(x) - 0.5*sum((y - mu)^2, na.rm=TRUE)/x^2 - log(1 + (sqrt(p)*exp(evolParams$dhs_mean0/2)/x)^2)
      }, lower = 0, upper = Inf)[1]
    }
    if(evol_error == 'HS') sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2 + length(evolParams$xiLambda), rate = sum((y - mu)^2, na.rm=TRUE)/2 + p*sum(evolParams$xiLambda)))
    if(evol_error == 'BL') sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2 + length(evolParams$tau_j)/2, rate = sum((y - mu)^2, na.rm=TRUE)/2 + p*sum((omega/evolParams$tau_j)^2)/2))
    if((evol_error == 'NIG') || (evol_error == 'SV')) sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2, rate = sum((y - mu)^2, na.rm=TRUE)/2))

    # Sample the initial variance parameters:
    evolParams0 = sampleEvol0(beta0, evolParams0, A = 1)

    # Store the MCMC output:
    if(nsi > nburn){
      # Increment the skip counter:
      skipcount = skipcount + 1

      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1

        # Save the MCMC samples:
        if(!is.na(match('mu', mcmc_params)) || computeDIC) post_mu[isave,] = mu
        if(!is.na(match('yhat', mcmc_params))) post_yhat[isave,] = mu + sigma_e*rnorm(T)
        if(!is.na(match('beta', mcmc_params))) post_beta[isave,] = matrix(beta)
        if(!is.na(match('obs_sigma_t2', mcmc_params)) || computeDIC) post_obs_sigma_t2[isave,] = sigma_e^2
        if(!is.na(match('evol_sigma_t2', mcmc_params))) post_evol_sigma_t2[isave,] =  rbind(matrix(evolParams0$sigma_w0^2, nrow = D), evolParams$sigma_wt^2)
        if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi[isave] = evolParams$dhs_phi
        if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean[isave] = evolParams$dhs_mean
        post_loglike[isave] = sum(dnorm(y, mean = mu, sd = sigma_e, log = TRUE))

        # And reset the skip counter:
        skipcount = 0
      }
    }
  }

  if(!is.na(match('mu', mcmc_params))) mcmc_output$mu = post_mu
  if(!is.na(match('yhat', mcmc_params))) mcmc_output$yhat = post_yhat
  if(!is.na(match('beta', mcmc_params))) mcmc_output$beta = post_beta
  if(!is.na(match('obs_sigma_t2', mcmc_params))) mcmc_output$obs_sigma_t2 = post_obs_sigma_t2
  if(!is.na(match('evol_sigma_t2', mcmc_params))) mcmc_output$evol_sigma_t2 = post_evol_sigma_t2
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_phi = post_dhs_phi
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_mean = post_dhs_mean

  # Also include the log-likelihood:
  mcmc_output$loglike = post_loglike

  if(computeDIC){
    # Log-likelihood evaluated at posterior means:
    loglike_hat = sum(dnorm(y,
                            mean = colMeans(post_mu),
                            sd = colMeans(sqrt(post_obs_sigma_t2)),
                            log = TRUE))
    # Effective number of parameters (Note: two options)
    p_d = c(2*(loglike_hat - mean(post_loglike)),
            2*var(post_loglike))
    # DIC:
    DIC = -2*loglike_hat + 2*p_d

    # Store the DIC and the effective number of parameters (p_d)
    mcmc_output$DIC = DIC; mcmc_output$p_d = p_d
  }

  return (mcmc_output);
}

#----------------------------------------------------------------------------
#' MCMC Sampler for B-spline Bayesian Trend Filtering: D = 0
#'
#' Run the MCMC for B-spline fitting with a penalty the B-spline basis coefficients.
#' The penalty is determined by the prior on the evolution errors, which include:
#' \itemize{
#' \item the dynamic horseshoe prior ('DHS');
#' \item the static horseshoe prior ('HS');
#' \item the Bayesian lasso ('BL');
#' \item the normal stochastic volatility model ('SV');
#' \item the normal-inverse-gamma prior ('NIG').
#' }
#' In each case, the evolution error is a scale mixture of Gaussians.
#' Sampling is accomplished with a (parameter-expanded) Gibbs sampler,
#' mostly relying on a dynamic linear model representation.

#' @param y the \code{T x 1} vector of time series observations
#' @param x the \code{T x 1} vector of observation points; if NULL, assume equally spaced
#' @param num_knots the number of knots; if NULL, use the default of \code{max(20, min(ceiling(T/4), 150))}
#' @param evol_error the evolution error distribution; must be one of
#' 'DHS' (dynamic horseshoe prior), 'HS' (horseshoe prior), 'BL' (Bayesian lasso), or 'NIG' (normal-inverse-gamma prior)
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burin-in)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param mcmc_params named list of parameters for which we store the MCMC output;
#' must be one or more of:
#' \itemize{
#' \item "mu" (conditional mean)
#' \item "beta" (B-spline basis coefficients)
#' \item "yhat" (posterior predictive distribution)
#' \item "evol_sigma_t2" (evolution error variance)
#' \item "obs_sigma_t2" (observation error variance)
#' \item "dhs_phi" (DHS AR(1) coefficient)
#' \item "dhs_mean" (DHS AR(1) unconditional mean)
#' }
#' @param computeDIC logical; if TRUE, compute the deviance information criterion \code{DIC}
#' and the effective number of parameters \code{p_d}
#' @param verbose logical; should R report extra information on progress?
#'
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#'
#' @note The data \code{y} may contain NAs, which will be treated with a simple imputation scheme
#' via an additional Gibbs sampling step. In general, rescaling \code{y} to have unit standard
#' deviation is recommended to avoid numerical issues.
#'
#' @note The primary advantages of \code{btf_bspline} over \code{btf} are
#' \enumerate{
#' \item Unequally-spaced points are handled automatically and
#' \item Computations are linear in the number of basis coefficients, which may be
#' substantially fewer than the number of time points.
#' }
#'
#'
#' @import fda progress
btf_bspline0 = function(y, x = NULL, num_knots = NULL, evol_error = 'DHS',
                        nsave = 1000, nburn = 1000, nskip = 4,
                        mcmc_params = list("mu", "yhat", "beta", "evol_sigma_t2", "obs_sigma_t2", "dhs_phi", "dhs_mean"),
                        computeDIC = TRUE,
                        verbose = TRUE){

  # Length of time series
  T = length(y);

  # Observation points
  if(is.null(x)) x = seq(0, 1, length.out=T);

  # Rescale to (0,1):
  t01 = (x - min(x))/diff(range(x))

  # Compute B-spline basis matrix:
  if(is.null(num_knots)) num_knots = max(20, min(ceiling(T/4), 150))
  X = eval.basis(t01, create.bspline.basis(c(0,1), nbasis = num_knots))
  p = ncol(X)

  # In place of XtX, store the 4-bands of XtX:
  XtX = crossprod(X) #   XtX = Matrix(crossprod(X));
  XtX_bands = list(XtX_0 = diag(XtX),
                   XtX_1 = diag(XtX[-1,]),
                   XtX_2 = diag(XtX[-(1:2),]),
                   XtX_3 = diag(XtX[-(1:3),]))
  rm(XtX)

  # Begin by checking for missing values, then imputing (for initialization)
  is.missing = which(is.na(y)); any.missing = (length(is.missing) > 0)

  # Impute the active "data"
  y = approxfun(t01, y, rule = 2)(t01)

  # Initial SD (implicitly assumes a constant mean)
  sigma_e = sd(y, na.rm=TRUE)

  # Initialize the B-spline coefficients via sampling
  Xty = crossprod(X,y)
  beta = sampleBTF_bspline(y, X, obs_sigma2 = sigma_e^2, evol_sigma_t2 = rep(0.01*sigma_e^2, p), XtX_bands = XtX_bands, Xty = Xty, D = 0)

  # And the conditional expectation:
  mu = as.numeric(X%*%beta)

  # Initialize the evolution error variance paramters:
  evolParams = initEvolParams(omega = beta, evol_error = evol_error)

  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('mu', mcmc_params)) || computeDIC) post_mu = array(NA, c(nsave, T))
  if(!is.na(match('yhat', mcmc_params))) post_yhat = array(NA, c(nsave, T))
  if(!is.na(match('beta', mcmc_params))) post_beta = array(NA, c(nsave, p))
  if(!is.na(match('obs_sigma_t2', mcmc_params)) || computeDIC) post_obs_sigma_t2 = array(NA, c(nsave, T))
  if(!is.na(match('evol_sigma_t2', mcmc_params))) post_evol_sigma_t2 = array(NA, c(nsave, p))
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi = numeric(nsave)
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean = numeric(nsave)
  post_loglike = numeric(nsave)

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  if(verbose){
    pb <- progress::progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                     total = nstot,
                                     complete = "=",   # Completion bar character
                                     incomplete = "-", # Incomplete bar character
                                     current = ">",    # Current bar character
                                     clear = FALSE,    # If TRUE, clears the bar when finish
                                     width = 100)      # Width of the progress bar
  }
  for(nsi in 1:nstot){
    if(verbose){
      if(nsi < 10){
        pb$tick()
      }
      else if(((nsi%%100) == 0)){
        pb$tick(100)
      }
    }
    # Impute missing values, if any:
    if(any.missing) {
      y[is.missing] = mu[is.missing] + sigma_e*rnorm(length(is.missing))
      Xty = crossprod(X,y) # Recompute when missing values present
    }

    # Sample the dynamic regression coefficients, beta:
    beta = sampleBTF_bspline(y, X, obs_sigma2 = sigma_e^2, evol_sigma_t2 = evolParams$sigma_wt^2, XtX_bands = XtX_bands, Xty = Xty, D = 0)

    # And the conditional expectation:
    mu = as.numeric(X%*%beta)

    # Sample the evolution error variance (and associated parameters):
    evolParams = sampleEvolParams(omega = beta, evolParams, sigma_e/sqrt(p), evol_error)

    # Sample the observation error SD:
    if(evol_error == 'DHS') {
      sigma_e = uni.slice(sigma_e, g = function(x){
        -(T+2)*log(x) - 0.5*sum((y - mu)^2, na.rm=TRUE)/x^2 - log(1 + (sqrt(p)*exp(evolParams$dhs_mean0/2)/x)^2)
      }, lower = 0, upper = Inf)[1]
    }
    if(evol_error == 'HS') sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2 + length(evolParams$xiLambda), rate = sum((y - mu)^2, na.rm=TRUE)/2 + p*sum(evolParams$xiLambda)))
    if(evol_error == 'BL') sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2 + length(evolParams$tau_j)/2, rate = sum((y - mu)^2, na.rm=TRUE)/2 + p*sum((beta/evolParams$tau_j)^2)/2))
    if((evol_error == 'NIG') || (evol_error == 'SV'))  sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2, rate = sum((y - mu)^2, na.rm=TRUE)/2))

    # Store the MCMC output:
    if(nsi > nburn){
      # Increment the skip counter:
      skipcount = skipcount + 1

      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1

        # Save the MCMC samples:
        if(!is.na(match('mu', mcmc_params)) || computeDIC) post_mu[isave,] = mu
        if(!is.na(match('yhat', mcmc_params))) post_yhat[isave,] = mu + sigma_e*rnorm(T)
        if(!is.na(match('beta', mcmc_params))) post_beta[isave,] = matrix(beta)
        if(!is.na(match('obs_sigma_t2', mcmc_params)) || computeDIC) post_obs_sigma_t2[isave,] = sigma_e^2
        if(!is.na(match('evol_sigma_t2', mcmc_params))) post_evol_sigma_t2[isave,] = evolParams$sigma_wt^2
        if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi[isave] = evolParams$dhs_phi
        if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean[isave] = evolParams$dhs_mean
        post_loglike[isave] = sum(dnorm(y, mean = mu, sd = sigma_e, log = TRUE))

        # And reset the skip counter:
        skipcount = 0
      }
    }
  }

  if(!is.na(match('mu', mcmc_params))) mcmc_output$mu = post_mu
  if(!is.na(match('yhat', mcmc_params))) mcmc_output$yhat = post_yhat
  if(!is.na(match('beta', mcmc_params))) mcmc_output$beta = post_beta
  if(!is.na(match('obs_sigma_t2', mcmc_params))) mcmc_output$obs_sigma_t2 = post_obs_sigma_t2
  if(!is.na(match('evol_sigma_t2', mcmc_params))) mcmc_output$evol_sigma_t2 = post_evol_sigma_t2
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_phi = post_dhs_phi
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_mean = post_dhs_mean

  # Also include the log-likelihood:
  mcmc_output$loglike = post_loglike

  if(computeDIC){
    # Log-likelihood evaluated at posterior means:
    loglike_hat = sum(dnorm(y,
                            mean = colMeans(post_mu),
                            sd = colMeans(sqrt(post_obs_sigma_t2)),
                            log = TRUE))
    # Effective number of parameters (Note: two options)
    p_d = c(2*(loglike_hat - mean(post_loglike)),
            2*var(post_loglike))
    # DIC:
    DIC = -2*loglike_hat + 2*p_d

    # Store the DIC and the effective number of parameters (p_d)
    mcmc_output$DIC = DIC; mcmc_output$p_d = p_d
  }

  return (mcmc_output);
}
