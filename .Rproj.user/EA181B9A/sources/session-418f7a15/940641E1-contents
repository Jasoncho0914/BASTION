## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @importFrom grDevices dev.new
#' @importFrom graphics lines par polygon
#' @importFrom methods as
#' @importFrom stats approxfun arima dbeta dnorm mad quantile rbinom rexp rgamma rnorm runif sd var
#' @useDynLib dsp
## usethis namespace: end
NULL
#----------------------------------------------------------------------------
#' Adaptive Bayesian Changepoint with Outliers
#'
#' Run the MCMC sampler for ABCO with a penalty on
#' first (D = 1), or second (D = 2) differences of the conditional expectation.
#' The penalty utilizes the dynamic horseshoe prior on the evolution errors.
#' Sampling is accomplished with a (parameter-expanded) Gibbs sampler,
#' mostly relying on a dynamic linear model representation.

#' @param y the \code{T} vector of time series observations
#' @param D degree of differencing (D = 1, or D = 2)
#' @param obsSV Options for modeling the error variance. It must be one of the following:
#' \itemize{
#' \item const: Constant error variance for all time points.
#' \item SV: Stochastic Volatility model.
#' \item ASV: Adaptive Stochastic Volatility model.
#' }
#' @param useAnom logical; if TRUE, include an anomaly component in the observation equation
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
#' \item "zeta_sigma_t2" (outlier error variance)
#' \item "dhs_phi" (DHS AR(1) coefficient)
#' \item "dhs_mean" (DHS AR(1) unconditional mean)
#' }
#' @param verbose logical; should R report extra information on progress?
#' @param cp_thres Percentage of posterior samples needed to declare a changepoint
#'
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#'
#' @import Matrix progress
#' @export
abco = function(y, D = 1, useAnom=TRUE, obsSV = "const",
                nsave = 1000, nburn = 1000, nskip = 4,
                mcmc_params = list('mu', "omega", "r"),
                verbose = TRUE, cp_thres = 0.5){
  # Time points (in [0,1])
  T = length(y); t01 = seq(0, 1, length.out=T)
  evol_error = 'DHS'

  # Store the original (w/ missing) data, then impute the active "data"
  yna = y; y = approxfun(t01, y, rule = 2)(t01)

  # Compute the Cholesky term (uses random variances for a more conservative sparsity pattern)
  chol0 = NULL
  loc = t_create_loc(length(y)-D, 1)
  loc_obs = t_create_loc(length(y), D)

  # Initial SD (implicitly assumes a constant mean)
  sigma_e = rep(sd(y, na.rm=TRUE), T)
  sigma_et = sd(y)

  evolParams0 = list(sigma_w0 = rep(1000,D))
  mu = t_sampleBTF(y, obs_sigma_t2 = sigma_e^2, evol_sigma_t2 = 0.01*sigma_e^2, D = D, loc_obs)
  mu[1:D] = y[1:D]

  if (useAnom){
    zeta = t_sampleBTF(y[-(1:D)]-mu[-(1:D)], obs_sigma_t2 = sigma_e[-(1:D)]^2, evol_sigma_t2 = 0.01*sigma_e[-(1:D)]^2, D = 0, loc_obs)
    zeta_params = t_initEvolZeta_ps(zeta)
  }

  # Compute the evolution errors:
  omega = diff(mu, differences = D)

  evolParams = t_initEvolParams_no(y, D, omega)
  lower_b = evolParams$lower_b
  upper_b = evolParams$upper_b

  if (useAnom) {
    if(obsSV == "SV"){
      svParams = t_initSV(y-mu-c(rep(0,D),zeta))
      sigma_e = svParams$sigma_wt
    }else if(obsSV == "ASV"){
      sParams = init_paramsASV(y-mu - c(rep(0,D),zeta), evol_error = "HS", D = 2)
      sigma_e = exp(sParams$s_mu/2)
    }
  }else{
    if(obsSV == "SV"){
      svParams = t_initSV(y-mu)
      sigma_e = svParams$sigma_wt
    }else if(obsSV == "ASV"){
      sParams = init_paramsASV(y-mu, evol_error = "HS", D = 2)
      sigma_e = exp(sParams$s_mu/2)
    }
  }


  #Array to store omega
  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('mu', mcmc_params))) post_mu = array(NA, c(nsave, T))
  if(!is.na(match('zeta', mcmc_params))) post_zeta = array(NA, c(nsave, T-D))
  if(!is.na(match('yhat', mcmc_params))) post_yhat = array(NA, c(nsave, T))
  if(!is.na(match('h', mcmc_params))) post_h = array(NA, c(nsave, T-D))
  if(!is.na(match('r', mcmc_params)) && evol_error == "DHS") post_r = numeric(nsave)
  if(!is.na(match('omega', mcmc_params))) post_omega = array(NA, c(nsave, T-D))
  if(!is.na(match('zeta_sigma_t2', mcmc_params))) post_zeta_sigma_t2 = array(NA, c(nsave, T-D))
  if(!is.na(match('obs_sigma_t2', mcmc_params))) post_obs_sigma_t2 = array(NA, c(nsave, T))
  if(!is.na(match('evol_sigma_t2', mcmc_params))) post_evol_sigma_t2 = array(NA, c(nsave, T))
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi = numeric(nsave)
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean = numeric(nsave)

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
    # Sample the states:
    if (useAnom){
      mu = t_sampleBTF(y-c(rep(0, D), zeta), obs_sigma_t2 = sigma_e^2, evol_sigma_t2 = c(evolParams0$sigma_w0^2, evolParams$sigma_wt^2), D = D, loc_obs)
      zeta = t_sampleBTF(y[-(1:D)]-mu[-(1:D)], obs_sigma_t2 = sigma_e[-(1:D)]^2, evol_sigma_t2 = zeta_params$sigma_wt^2, D = 0)
      zeta_params = t_sampleEvolZeta_ps(zeta, zeta_params)
    } else{
      mu = t_sampleBTF(y, obs_sigma_t2 = sigma_e^2, evol_sigma_t2 = c(evolParams0$sigma_w0^2, evolParams$sigma_wt^2), D = D, loc_obs)
    }
    omega = diff(mu, differences = D)

    evolParams = t_sampleEvolParams(omega, evolParams, D, 1/T, lower_b, upper_b, loc)
    #Observational Params
    if (useAnom){
      if (obsSV == "SV"){
        svParams = t_sampleSVparams(omega = y - mu- c(rep(0, D), zeta), svParams = svParams)
        sigma_e = svParams$sigma_wt
      } else if(obsSV == "const"){
        sigma_et = uni.slice(sigma_et, g = function(x){
          -(T+2)*log(x) - 0.5*sum((y - mu - c(rep(0,D), zeta))^2, na.rm=TRUE)/x^2
        }, lower = 0, upper = Inf)[1]
        sigma_e = rep(sigma_et, T)
      } else if(obsSV == "ASV"){
        sParams = fit_paramsASV(y-mu - c(rep(0, D), zeta),sParams,evol_error = "HS", D = 2)
        sigma_e = exp(sParams$s_mu/2)
      } else{
        stop('obsSV has to be one of const, SV, or ASV')
      }
    }else{
      if(obsSV == "SV"){
        svParams = t_sampleSVparams(omega = y - mu, svParams = svParams)
        sigma_e = svParams$sigma_wt
      }else if(obsSV == "const"){
        sigma_et = uni.slice(sigma_et, g = function(x){
          -(T+2)*log(x) - 0.5*sum((y - mu)^2, na.rm=TRUE)/x^2 - log(1 + (sqrt(T)*exp(evolParams$dhs_mean/2)/x)^2)
        }, lower = 0, upper = Inf)[1]
        sigma_e = rep(sigma_et, T)
      }else if(obsSV == "ASV"){
        sParams = fit_paramsASV(y-mu,sParams,evol_error = "HS", D = 2)
        sigma_e = exp(sParams$s_mu/2)
      }else{
        stop('obsSV has to be one of const, SV, or ASV')
      }
    }

    if(nsi > nburn){
      # Increment the skip counter:
      skipcount = skipcount + 1

      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1

        if(!is.na(match('mu', mcmc_params))) post_mu[isave,] = mu
        if(!is.na(match('yhat', mcmc_params))) post_yhat[isave,] = mu + rnorm(T, 0, sigma_e)
        if (useAnom) {
          if(!is.na(match('zeta', mcmc_params))) post_zeta[isave,] = zeta
          if(!is.na(match('zeta_sigma_t2', mcmc_params))) post_zeta_sigma_t2[isave,] = zeta_params$sigma_wt^2
        }
        if(!is.na(match('h', mcmc_params))) post_h[isave,] = evolParams$ht
        if(!is.na(match('omega', mcmc_params))) post_omega[isave,] = omega
        if(!is.na(match('r', mcmc_params)) && evol_error == "DHS") post_r[isave] = evolParams$r
        if(!is.na(match('obs_sigma_t2', mcmc_params))) post_obs_sigma_t2[isave,] = sigma_e^2
        if(!is.na(match('evol_sigma_t2', mcmc_params))) post_evol_sigma_t2[isave,] = c(evolParams0$sigma_w0^2, evolParams$sigma_wt^2)
        if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi[isave] = evolParams$dhs_phi
        if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean[isave] = evolParams$dhs_mean

        skipcount = 0
      }
    }
  }
  if(!is.na(match('zeta', mcmc_params))) mcmc_output$zeta = post_zeta
  if(!is.na(match('mu', mcmc_params))) mcmc_output$mu = post_mu
  if(!is.na(match('yhat', mcmc_params))) mcmc_output$yhat = post_yhat
  if(!is.na(match('r', mcmc_params)) && evol_error == "DHS") mcmc_output$r = post_r
  if(!is.na(match('h', mcmc_params))) mcmc_output$h = post_h
  if(!is.na(match('omega', mcmc_params))) mcmc_output$omega = post_omega
  if(!is.na(match('zeta_sigma_t2', mcmc_params))) mcmc_output$zeta_sigma_t2 = post_zeta_sigma_t2
  if(!is.na(match('obs_sigma_t2', mcmc_params))) mcmc_output$obs_sigma_t2 = post_obs_sigma_t2
  if(!is.na(match('evol_sigma_t2', mcmc_params))) mcmc_output$evol_sigma_t2 = post_evol_sigma_t2
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_phi = post_dhs_phi
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_mean = post_dhs_mean

  mcmc_output$cp = identify_cp(D, mcmc_output, cp_thres)
  return (mcmc_output)
}
#----------------------------------------------------------------------------
#' Initialize the anomaly component parameters
#'
#' Compute initial values for either a horseshoe prior or horseshoe+ prior
#' for the anomaly component.
#'
#' @param zeta \code{T} vector of initial estimates.
#' @return List of relevant components: \code{sigma_wt}, the \code{T} vector of standard deviations,
#' and additional parameters for inverse gamma priors (shape and scale).
#'
#' @import MCMCpack
#' @export
t_initEvolZeta_ps = function(zeta){
  zeta = as.matrix(zeta)
  n = nrow(zeta)

  xi = rinvgamma(1, 1/2, 1)
  tau_t2 = rinvgamma(1, 1/2, 1/xi)

  #eta_v = rinvgamma(n, 1/2, 1)
  #eta_t2 = rinvgamma(n, 1/2, 1/eta_v)

  #v = rinvgamma(n, 1/2, 1/eta_t2/tau_t2)
  v = rinvgamma(n, 1/2, 1/tau_t2)
  lambda_t2 = rinvgamma(n, 1/2, 1/v)

  #list(xi=xi, v=v, eta_v = eta_v, tau_t2 = tau_t2, lambda_t2 = lambda_t2, eta_t2 = eta_t2,
  #     sigma_wt = sqrt(lambda_t2))
  list(xi=xi, v=v, tau_t2 = tau_t2, lambda_t2 = lambda_t2, sigma_wt = sqrt(lambda_t2))
}
#----------------------------------------------------------------------------
#' Sampler for the anomaly component parameters
#'
#' Compute one draw of the anomaly component parameters.
#'
#' @param omega \code{T} vector of errors
#' @param evolParams list of parameters to be updated
#' @return List of relevant components in \code{evolParams}:  \code{sigma_wt},
#' the \code{T} vector of standard deviations, and additional parameters for inverse gamma priors (shape and scale).
#'
#' @import MCMCpack
#' @export
t_sampleEvolZeta_ps = function(omega, evolParams){
  #omega = as.matrix(omega)
  n = length(omega)

  hsInput2 = omega^2

  evolParams$lambda_t2 = rinvgamma(n, 1, 1/evolParams$v + hsInput2/2)
  #evolParams$v = rinvgamma(n, 1, 1/evolParams$lambda_t2 + 1/evolParams$tau_t2/evolParams$eta_t2)
  evolParams$v = rinvgamma(n, 1, 1/evolParams$lambda_t2 + 1/evolParams$tau_t2)

  #evolParams$tau_t2 = rinvgamma(1, (n+1)/2, 1/evolParams$xi + sum(1/evolParams$v/evolParams$eta_t2))
  evolParams$tau_t2 = rinvgamma(1, (n+1)/2, 1/evolParams$xi + sum(1/evolParams$v))
  evolParams$xi = rinvgamma(1, 1, 1+1/evolParams$tau_t2)

  #evolParams$eta_t2 = rinvgamma(n, 1, 1/evolParams$eta_v+1/evolParams$v/evolParams$tau_t2)
  #evolParams$eta_v = rinvgamma(n, 1, 1/16+1/evolParams$eta_t2)

  evolParams$sigma_wt = sqrt(evolParams$lambda_t2)

  evolParams
}
#----------------------------------------------------------------------------
#' Sampler for first or second order random walk (RW) Gaussian dynamic linear model (DLM)
#'
#' Compute one draw of the \code{T} state variable \code{mu} in a DLM using back-band substitution methods.
#' This model is equivalent to the Bayesian trend filtering (BTF) model, assuming appropriate
#' (shrinkage/sparsity) priors for the evolution errors.
#'
#' @param y the \code{T x 1} vector of time series observations
#' @param obs_sigma_t2 the \code{T x 1} vector of observation error variances
#' @param evol_sigma_t2 the \code{T x 1} vector of evolution error variances
#' @param D the degree of differencing (one or two)
#' @param loc_obs list of the row and column indices to fill in a band-sparse matrix
#' @return \code{T x 1} vector of simulated states
#'
#' @note Missing entries (NAs) are not permitted in \code{y}. Imputation schemes are available.
#'
#' @import Matrix Rcpp RcppZiggurat
#' @export
t_sampleBTF = function(y, obs_sigma_t2, evol_sigma_t2, D = 1, loc_obs){

  # Some quick checks:
  if((D < 0) || (D != round(D)))  stop('D must be a positive integer')

  if(any(is.na(y))) stop('y cannot contain NAs')

  T = length(y)

  # Linear term:
  linht = y/obs_sigma_t2

  # Quadratic terms and solutions are computed differently, depending on D:

  if(D == 0){
    # Special case: no differencing

    # Posterior SDs and posterior means:
    postSD = 1/sqrt(1/obs_sigma_t2 + 1/evol_sigma_t2)
    postMean = (linht)*postSD^2

    # Sample the states:
    mu = rnorm(n = T, mean = postMean, sd = postSD)

  } else {
    # Original sampler, based on Matrix package:
    if (D == 1) {
      diag1 = 1/obs_sigma_t2 + 1/evol_sigma_t2 + c(1/evol_sigma_t2[-1], 0)
      diag2 = -1/evol_sigma_t2[-1]
      rd = RcppZiggurat::zrnorm(T)
      mu = sample_mat_c(loc_obs$r, loc_obs$c, c(diag1, diag2, diag2), length(diag1), length(loc_obs$r), c(linht), rd, D)
    } else {
      diag1 = 1/obs_sigma_t2 + 1/evol_sigma_t2 + c(0, 4/evol_sigma_t2[-(1:2)], 0) + c(1/evol_sigma_t2[-(1:2)], 0, 0)
      diag2 = c(-2/evol_sigma_t2[3], -2*(1/evol_sigma_t2[-(1:2)] + c(1/evol_sigma_t2[-(1:3)],0)))
      diag3 = 1/evol_sigma_t2[-(1:2)]
      rd = RcppZiggurat::zrnorm(T)
      mu = sample_mat_c(loc_obs$r, loc_obs$c, c(diag1, diag2, diag2, diag3, diag3), length(diag1), length(loc_obs$r), c(linht), rd, D)
    }
  }

  # And return the states:
  mu
}
#----------------------------------------------------------------------------
#' Sample the thresholded dynamic shrinkage process parameters
#'
#' Compute one draw for each of the parameters in the thresholded dynamic shrinkage process
#' for the special case in which the shrinkage parameter \code{kappa ~ Beta(alpha, beta)}
#' with \code{alpha = beta = 1/2}.
#'
#' @param omega \code{T} vector of evolution errors
#' @param evolParams list of parameters to be updated (see Value below)
#' @param D the degree of differencing (one or two)
#' @param sigma_e the observation error standard deviation; for (optional) scaling purposes
#' @param upper_b the upper bound in the uniform prior of the threshold variable
#' @param lower_b the lower bound in the uniform prior of the threshold variable
#' @param loc list of the row and column indices to fill in a band-sparse matrix
#' @param prior_dhs_phi the parameters of the prior for the log-volatility AR(1) coefficient \code{dhs_phi};
#' either \code{NULL} for uniform on [-1,1] or a 2-dimensional vector of (shape1, shape2) for a Beta prior
#' on \code{[(dhs_phi + 1)/2]}
#' @param alphaPlusBeta For the symmetric prior kappa ~ Beta(alpha, beta) with alpha=beta,
#' specify the sum [alpha + beta]
#' @return List of relevant components:
#' \itemize{
#' \item the \code{T} evolution error standard deviations \code{sigma_wt},
#' \item the \code{T} log-volatility \code{ht},
#' \item the \code{1} log-vol unconditional mean(s) \code{dhs_mean},
#' \item the \code{1} log-vol AR(1) coefficient(s) \code{dhs_phi},
#' \item the \code{1} log-vol correction coefficient(s) \code{dhs_phi2},
#' \item the \code{T} log-vol innovation standard deviations \code{sigma_eta_t} from the Polya-Gamma priors,
#' \item the \code{1} initial log-vol SD \code{sigma_eta_0},
#' \item the \code{1} threshold parameter r
#' }
#'
#' @note The priors induced by \code{prior_dhs_phi} all imply a stationary (log-) volatility process.
#'
#' @import pgdraw
t_sampleEvolParams = function(omega, evolParams, D = 1, sigma_e = 1, lower_b, upper_b, loc, prior_dhs_phi = c(20,1), alphaPlusBeta = 1){
  # Store the DSP parameters locally:
  ht = evolParams$ht; dhs_mean = evolParams$dhs_mean; dhs_phi = evolParams$dhs_phi; dhs_phi2 = evolParams$dhs_phi2
  sigma_eta_t = evolParams$sigma_eta_t; sigma_eta_0 = evolParams$sigma_eta_0; r=evolParams$r

  # "Local" number of time points
  n = length(ht); p = 1
  st = (log(omega^2 + 0.0001) >= r)

  # Sample the log-volatilities using AWOL sampler
  ht = t_sampleLogVols(h_y = omega, h_prev = ht, h_mu = dhs_mean, h_phi=dhs_phi, h_phi2 = dhs_phi2, h_sigma_eta_t = sigma_eta_t, h_sigma_eta_0 = sigma_eta_0, h_st = st, loc = loc)

  # Compute centered log-vols for the samplers below:
  ht_tilde = ht - dhs_mean

  # Sample AR(1) parameters
  # Note: dhs_phi = 0 means non-dynamic HS, while dhs_phi = 1 means RW, in which case we don't sample either
  coef = t_sampleAR1(h_yc = ht_tilde, h_phi = dhs_phi, h_phi2 = dhs_phi2, h_sigma_eta_t = sigma_eta_t, h_st = st, prior_dhs_phi = prior_dhs_phi)
  dhs_phi = coef[1]
  dhs_phi2 = coef[2]

  # Sample the unconditional mean(s), unless dhs_phi = 1 (not defined)
  dhs_mean = t_sampleLogVolMu(h = ht, h_mu = dhs_mean, h_phi = dhs_phi, h_phi2 = dhs_phi2, h_sigma_eta_t = sigma_eta_t, h_sigma_eta_0 = sigma_eta_0, h_st = st, h_log_scale = log(sigma_e^2));

  # Sample the evolution error SD of log-vol (i.e., Polya-Gamma mixing weights)
  eta_t = ht_tilde[-1] - (dhs_phi+ dhs_phi2*st[-n])*ht_tilde[-n]       # Residuals
  #sigma_eta_t = matrix(1/sqrt(rpg(num = (n-1), h = alphaPlusBeta, z = eta_t)), nc = 1) # Sample
  #sigma_eta_0 = 1/sqrt(rpg(num = 1, h = 1, z = ht_tilde[1]))                # Sample the inital
  sigma_eta_t = 1/sqrt(pgdraw(alphaPlusBeta, eta_t))
  sigma_eta_0 = 1/sqrt(pgdraw(1, ht_tilde[1]))

  r = t_sampleR_mh(h_yc = ht_tilde, h_phi = dhs_phi, h_phi2 = dhs_phi2, h_sigma_eta_t = sigma_eta_t,
                   h_sigma_eta_0 = sigma_eta_0, h_st = st, h_r = r, lower_b = lower_b, upper_b = upper_b, omega = omega, D = D)

  # Evolution error SD:
  sigma_wt = exp(ht/2)

  # Return the same list, but with the new values
  list(sigma_wt = sigma_wt, ht = ht, dhs_mean = dhs_mean, dhs_phi = dhs_phi, dhs_phi2 = dhs_phi2, sigma_eta_t = sigma_eta_t, sigma_eta_0 = sigma_eta_0, r=r)
}
#----------------------------------------------------------------------------
#' Sample the latent log-volatilities
#'
#' Compute one draw of the log-volatilities using a discrete mixture of Gaussians
#' approximation to the likelihood (see Omori, Chib, Shephard, and Nakajima, 2007)
#' where the log-vols are assumed to follow an TAR(1) model with time-dependent
#' innovation variances. More generally, the code operates for \code{p} independent
#' TAR(1) log-vol processes to produce an efficient joint sampler in \code{O(Tp)} time.
#'
#' @param h_y the \code{T } vector of data, which follow independent SV models
#' @param h_prev the \code{T} vector of the previous log-vols
#' @param h_mu the \code{1} vector of log-vol unconditional means
#' @param h_phi the \code{1} vector of log-vol AR(1) coefficients
#' @param h_phi2 the \code{1} vector of previous penalty coefficient(s)
#' @param h_sigma_eta_t the \code{T} vector of log-vol innovation standard deviations
#' @param h_sigma_eta_0 the \code{1} vector of initial log-vol innovation standard deviations
#' @param h_st the \code{T} vector of indicators on whether each time-step exceed the estimated threshold
#' @param loc list of the row and column indices to fill in the band-sparse matrix in the sampler
#'
#' @return \code{T x p} vector of simulated log-vols
#'
#' @note For Bayesian trend filtering, \code{p = 1}. More generally, the sampler allows for
#' \code{p > 1} but assumes (contemporaneous) independence across the log-vols for \code{j = 1,...,p}.
#'
#' @import Matrix RcppZiggurat
t_sampleLogVols = function(h_y, h_prev, h_mu, h_phi, h_phi2, h_sigma_eta_t, h_sigma_eta_0, h_st, loc){

  # Compute dimensions:
  n = length(h_prev); p = 1

  yoffset = 10^-8

  # This is the response in our DLM, log(y^2)
  ystar = log(h_y^2 + yoffset)

  # Sample the mixture components
  #z = sapply(ystar-h_prev, ncind, m_st, sqrt(v_st2), q)
  #z = sapply(ystar-h_prev, ncind, m_st, sqrt(v_st2), q)
  # z = c(draw_indicators_generic(ystar-h_prev, rep(0, length(ystar)), length(ystar),
  #                             q, m_st, sqrt(v_st2), length(q)))

  # Subset mean and variances to the sampled mixture components; (n x p) matrices
  # m_st_all = m_st[z]
  # v_st2_all = v_st2[z]
  mixture_component = sample_j_wrap(n,ystar-hprev)

  # Subset mean and variances to the sampled mixture components; (n x p) matrices
  m_st_all = mixture_component$mean
  v_st2_all = mixture_component$var

  # Joint AWOL sampler for j=1,...,p:

  # Linear term:
  linht = (ystar - m_st_all - h_mu)/v_st2_all

  # Evolution precision matrix (n x p)
  evol_prec_mat = rep(0, n);
  evol_prec_mat[1] = 1/h_sigma_eta_0^2;
  evol_prec_mat[-1] = 1/h_sigma_eta_t^2;

  # Lagged version, with zeros as appropriate (needed below)
  evol_prec_lag_mat = rep(0, n)
  evol_prec_lag_mat[1:(n-1)] = evol_prec_mat[-1]

  # Diagonal of quadratic term:
  Q_diag = 1/v_st2_all +  evol_prec_mat + (h_phi+h_phi2*h_st)^2*evol_prec_lag_mat

  # Off-diagonal of quadratic term:
  Q_off = (-(h_phi+h_phi2*h_st)*evol_prec_lag_mat)[-n]

  # Sample the log-vols:
  rd = RcppZiggurat::zrnorm(length(linht))
  hsamp = h_mu + sample_mat_c(loc$r, loc$c, c(Q_diag, Q_off, Q_off), length(Q_diag), length(loc$r), c(linht), rd, 1)

  # Return the (uncentered) log-vols
  hsamp
}
#----------------------------------------------------------------------------
#' Sample the TAR(1) coefficients
#'
#' Compute one draw of the TAR(1) coefficients in a model with Gaussian innovations
#' and time-dependent innovation variances. In particular, we use the sampler for the
#' log-volatility TAR(1) process with the parameter-expanded Polya-Gamma sampler. The sampler also applies
#' to a multivariate case with independent components.
#'
#' @param h_yc the \code{T} vector of centered log-volatilities
#' (i.e., the log-vols minus the unconditional means \code{dhs_mean})
#' @param h_phi the \code{1} vector of previous AR(1) coefficient(s)
#' @param h_phi2 the \code{1} vector of previous penalty coefficient(s)
#' @param h_sigma_eta_t the \code{T} vector of log-vol innovation standard deviations
#' @param h_st the \code{T} vector of indicators on whether each time-step exceed the estimated threshold
#' @param prior_dhs_phi the parameters of the prior for the log-volatility AR(1) coefficient \code{dhs_phi};
#' either \code{NULL} for uniform on [-1,1] or a 2-dimensional vector of (shape1, shape2) for a Beta prior
#' on \code{[(dhs_phi + 1)/2]}
#'
#' @return \code{2} vector of sampled TAR(1) coefficient(s)
#'
#' @import truncdist
#' @export
t_sampleAR1 = function(h_yc, h_phi, h_phi2, h_sigma_eta_t, h_st, prior_dhs_phi = NULL){

  # Compute dimensions:
  n = length(h_yc);

  # Compute "regression" terms for dhs_phi_j:
  y_ar = h_yc[-1]/h_sigma_eta_t # Standardized "response"
  x_ar = h_yc[-n]/h_sigma_eta_t # Standardized "predictor"

  dhs_phi01 = (h_phi + 1)/2 # ~ Beta(prior_dhs_phi[1], prior_dhs_phi[2])

  # Slice sampler when using Beta prior:
  dhs_phi01 = uni.slice(dhs_phi01, g = function(x){
    -0.5*sum((y_ar - (2*x - 1 + h_phi2*h_st[-n])*x_ar)^2) +
      dbeta(x, shape1 = prior_dhs_phi[1], shape2 = prior_dhs_phi[2], log = TRUE)
  }, lower = 0, upper = 1)[1]#}, lower = 0.005, upper = 0.995)[1] #

  h_phi = 2*dhs_phi01 - 1

  index = which(h_st[-n] == 1)
  h_phi2 = uni.slice(h_phi2, g = function(x) {
    result = dnorm(x, mean = -1, sd=.5, log = TRUE)
    for (i in index){
      result = result - 0.5*(y_ar[i] - (h_phi + x)*x_ar[i])^2
    }
    result
  }, lower = -5, upper = 0)[1]


  c(h_phi, h_phi2)
}
#----------------------------------------------------------------------------
#' Sample the TAR(1) unconditional means
#'
#' Compute one draw of the unconditional means in an TAR(1) model with Gaussian innovations
#' and time-dependent innovation variances. In particular, we use the sampler for the
#' log-volatility TAR(1) process with the parameter-expanded Polya-Gamma sampler. The sampler also applies
#' to a multivariate case with independent components.
#'
#' @param h the \code{T} vector of log-volatilities
#' @param h_mu the \code{1} vector of previous means
#' @param h_phi the \code{1} vector of AR(1) coefficient(s)
#' @param h_phi2 the \code{1} vector of previous penalty coefficient(s)
#' @param h_sigma_eta_t the \code{T} vector of log-vol innovation standard deviations
#' @param h_sigma_eta_0 the standard deviations of initial log-vols
#' @param h_st the \code{T} vector of indicators on whether each time-step exceed the estimated threshold
#' @param h_log_scale prior mean from scale mixture of Gaussian (Polya-Gamma) prior, e.g. log(sigma_e^2) or dhs_mean0
#'
#' @return the sampled mean(s) \code{dhs_mean}
#'
#' @import pgdraw
#' @export
t_sampleLogVolMu = function(h, h_mu, h_phi, h_phi2, h_sigma_eta_t, h_sigma_eta_0, h_st, h_log_scale = 0){

  # Compute "local" dimensions:
  n = length(h)

  # Sample the precision term(s)
  dhs_mean_prec_j = pgdraw(1, h_mu - h_log_scale)
  #dhs_mean_prec_j = pgdraw(1, h_mu - h_log_scale)

  # Now, form the "y" and "x" terms in the (auto)regression
  y_mu = (h[-1] - (h_phi+h_phi2*h_st[-n])*h[-n])/h_sigma_eta_t;
  x_mu = (1 - h_phi-h_phi2*h_st[-n])/h_sigma_eta_t

  # Include the initial sd?
  #if(!is.null(h_sigma_eta_0)){y_mu = rbind(h[1,]/h_sigma_eta_0, y_mu); x_mu = rbind(1/h_sigma_eta_0, x_mu)}
  y_mu = c(h[1]/h_sigma_eta_0, y_mu);
  x_mu = c(1/h_sigma_eta_0, x_mu)

  # Posterior SD and mean:
  postSD = 1/sqrt(sum(x_mu^2) + dhs_mean_prec_j)
  postMean = (sum(x_mu*y_mu) + h_log_scale*dhs_mean_prec_j)*postSD^2
  dhs_mean = rnorm(n = 1, mean = postMean, sd = postSD)

  dhs_mean
}
#----------------------------------------------------------------------------
#' Sample the threshold parameter
#'
#' Compute one draw of thee threshold parameter in th TAR(1) model with Gaussian innovations
#' and time-dependent innovation variances. The sampler utilizes metropolis hasting to draw
#' from uniform prior.
#'
#' @param h_yc the \code{T} vector of centered log-volatilities
#' (i.e., the log-vols minus the unconditional means \code{dhs_mean})
#' @param h_phi the \code{1} vector of previous AR(1) coefficient(s)
#' @param h_phi2 the \code{1} vector of previous penalty coefficient(s)
#' @param h_sigma_eta_t the \code{T} vector of log-vol innovation standard deviations
#' @param h_sigma_eta_0 the \code{1} vector of initial log-vol innovation standard deviations
#' @param h_st the \code{T} vector of indicators on whether each time-step exceed the estimated threshold
#' @param h_r \code{1} the previous draw of the threshold parameter
#' @param upper_b the upper bound in the uniform prior of the threshold variable
#' @param lower_b the lower bound in the uniform prior of the threshold variable
#' @param omega \code{T} vector of evolution errors
#' @param D the degree of differencing (one or two)
#'
#' @return the sampled threshold value \code{r}
t_sampleR_mh = function(h_yc, h_phi, h_phi2, h_sigma_eta_t, h_sigma_eta_0, h_st, h_r, lower_b, upper_b, omega, D){
  n = length(h_yc);

  y_ar = h_yc[-1]/h_sigma_eta_t
  x_ar = h_yc[-n]/h_sigma_eta_t

  new_cand = runif(1, lower_b, upper_b)
  num = sum(dnorm(y_ar, (h_phi + h_phi2*(log(omega^2+0.0001)[-n] >= new_cand))*x_ar, 1, log=TRUE))
  den = sum(dnorm(y_ar, (h_phi + h_phi2*(log(omega^2+0.0001)[-n] >= h_r))*x_ar, 1, log=TRUE))

  if (num > den){
    r = new_cand
  }else{
    prob = runif(1, 0, 1)
    if (prob < exp(num-den)){
      r = new_cand
    } else{
      r = h_r
    }
  }

  r
}
#----------------------------------------------------------------------------
#' Initialize the evolution error variance parameters
#'
#' Compute initial values for evolution error variance parameters under the
#' dynamic horseshoe prior
#'
#' @param y the \code{T} vector of time series observations
#' @param D degree of differencing (D = 1, or D = 2)
#' @param omega \code{T} vector of evolution errors
#'
#' @return List of relevant components: \code{sigma_wt}, the \code{T} vector of evolution standard deviations,
#' and additional parameters associated with the DHS priors.
#'
#' @import msm
#' @export
t_initEvolParams_no = function(y, D, omega){

  # "Local" number of time points
  yoffset = tcrossprod(rep(1,length(omega)),
                       apply(as.matrix(omega), 2,
                             function(x) any(x^2 < 10^-16)*max(10^-8, mad(x)/10^6)))
  ht = log(omega^2+0.0001)
  n = length(omega)

  # Initialize the AR(1) model to obtain unconditional mean and AR(1) coefficient
  arCoefs = arima(ht, c(1,0,0), method="ML")$coef
  dhs_mean = arCoefs[2]; dhs_phi = arCoefs[1]
  dhs_phi2 = rtnorm(1, mean = -0.5, sd = .5, upper = 0, lower = -1)

  # Initialize the SD of log-vol innovations simply using the expectation:
  sigma_eta_t = matrix(pi, nrow = n-1, ncol = 1)
  sigma_eta_0 = rep(pi) # Initial value

  # Evolution error SD:
  sigma_wt = exp(ht/2)

  diff_y = log(diff(y, differences=D)^2+0.0001)
  lower_b = log(0.0001)
  upper_b = max(diff_y)
  r = runif(1, min=lower_b, max = upper_b)

  list(sigma_wt = sigma_wt, ht = ht, dhs_mean = dhs_mean, dhs_phi = dhs_phi, dhs_phi2 = dhs_phi2, sigma_eta_t = sigma_eta_t, sigma_eta_0 = sigma_eta_0, r=r, lower_b = lower_b, upper_b = upper_b)
}
#----------------------------------------------------------------------------
#' Initialize the stochastic volatility parameters
#'
#' Compute initial values for normal stochastic volatility parameters.
#' The model assumes an AR(1) for the log-volatility.
#'
#' @param omega \code{T} vector of errors
#' @return List of relevant components: \code{sigma_wt}, the \code{T} vector of standard deviations,
#' and additional parameters (unconditional mean, AR(1) coefficient, and standard deviation).
#' @importFrom methods is
t_initSV = function(omega){

  # Make sure omega is (n x p) matrix
  omega = omega; n = length(omega)

  # log-volatility:
  ht = log(omega^2 + 0.0001)

  # AR(1) pararmeters: check for error in initialization too
  ar_fit = try(arima(ht, c(1,0,0)), silent = TRUE)
  if(is(ar_fit, "try-error")) {
    params = c(ar_fit$coef[2], ar_fit$coef[1], sqrt(ar_fit$sigma2))
  } else params = c(mean(ht)/(1 - 0.8),0.8, 1)
  svParams = params

  # SDs, log-vols, and other parameters:
  return(list(sigma_wt = exp(ht/2), ht = ht, svParams = svParams))
}
#----------------------------------------------------------------------------
#' Sampler for the stochastic volatility parameters
#'
#' Compute one draw of the normal stochastic volatility parameters.
#' The model assumes an AR(1) for the log-volatility.
#'
#' @param omega \code{T} vector of errors
#' @param svParams list of parameters to be updated
#' @return List of relevant components in \code{svParams}: \code{sigma_wt}, the \code{T} vector of standard deviations,
#' and additional parameters associated with SV model.
#'
#' @import stochvol
t_sampleSVparams = function(omega, svParams){

  # Make sure omega is (n x p) matrix
  n = length(omega)

  # First, check for numerical issues:
  svInput = omega; #if(all(svInput==0)) {svInput = 10^-8} else svInput = svInput + sd(svInput)/10^8
  #hsOffset = tcrossprod(rep(1,n), apply(omega, 2, function(x) any(x^2 < 10^-16)*max(10^-8, mad(x)/10^6)))

  # Sample the SV parameters:
  svsamp = stochvol::svsample_fast_cpp(svInput,
                                       startpara = list(
                                         mu = svParams$svParams[1],
                                         phi = svParams$svParams[2],
                                         sigma = svParams$svParams[3]),
                                       startlatent = svParams$ht)# ,priorphi = c(10^4, 10^4));
  # Update the parameters:
  svParams$svParams = svsamp$para[1:3];
  svParams$ht = svsamp$latent[1,]

  # Finally, up the evolution error SD:
  svParams$sigma_wt = exp(svParams$ht/2)

  # Check for numerically large values:
  svParams$sigma_wt[which(svParams$sigma_wt > 10^3, arr.ind = TRUE)] = 10^3

  return(svParams)
}
#----------------------------------------------------------------------------
#' Initializer for location indices for filling in band-sparse matrix
#'
#' Create row and column indices for locations of symmetric band-sparse matrix.
#' Starts with the locations of the diagonal, proceed with upper-diagonals,
#' followed by lower-diagonals.
#'
#' @param len length of the diagonal of the band-sparse matrix
#' @param D number of super-diagonals to include for the band-sparse
#'
#' @return a list containing
#' \itemize{
#' \item the row indices \code{r} and
#' \item the column indices \code{c}
#' }
#'
#' @export
t_create_loc <- function(len, D){
  if (D == 0 || D == 1){
    row_ind = c()
    col_ind = c()
    for (i in 0:(len-1)){
      row_ind = c(row_ind, i)
      col_ind = c(col_ind, i)
    }
    for (i in 0:(len-2)){
      row_ind = c(row_ind, i)
      col_ind = c(col_ind, i+1)
    }
    for (i in 0:(len-2)){
      row_ind = c(row_ind, i+1)
      col_ind = c(col_ind, i)
    }
  }
  if (D == 2) {
    row_ind = c()
    col_ind = c()
    for (i in 0:(len-1)){
      row_ind = c(row_ind, i)
      col_ind = c(col_ind, i)
    }
    for (i in 0:(len-2)){
      row_ind = c(row_ind, i)
      col_ind = c(col_ind, i+1)
    }
    for (i in 0:(len-2)){
      row_ind = c(row_ind, i+1)
      col_ind = c(col_ind, i)
    }
    for (i in 0:(len-3)){
      row_ind = c(row_ind, i)
      col_ind = c(col_ind, i+2)
    }
    for (i in 0:(len-3)){
      row_ind = c(row_ind, i+2)
      col_ind = c(col_ind, i)
    }
  }
  list(r = row_ind, c = col_ind)
}
