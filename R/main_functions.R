#' @keywords internal
#' @importFrom stats approxfun rnorm
fit_ASD = function(y,
                   Ks,
                   X = NULL,
                   Outlier = FALSE,
                   sparse = FALSE,
                   obsSV = "const",
                   nsave = 1000,
                   nburn = 1000,
                   nskip = 4,
                   verbose = TRUE) {
  D = 2;evol_error = "HS";
  SVm = (obsSV == "SV")
  ASVm = (obsSV == "ASV")
  reg = !is.null(X) # is it regression
  if (reg & !is.matrix(X)) {
    stop("X needs to be a matrix")
  }
  nKs = length(Ks) # number of seasonality
  TT = length(y)
  t01 = seq(0, 1, length.out = TT)
  is.missing = which(is.na(y))
  any.missing = (length(is.missing) > 0)
  y = stats::approxfun(t01, y, rule = 2)(t01)
  offset_y = mean(y)
  y = y - offset_y
  # list of parameters
  params_list = list()

  # parameter matrices
  dims_er = c(TT, nKs + 1)
  dims_b =  c(TT, nKs + 1)
  dimnames_er = c(paste0("Seasonal", Ks), "Trend")
  dimnames_b = c(paste0("Seasonal", Ks), "Trend")
  col_trend = nKs + 1
  if (reg & Outlier) {
    col_reg = nKs + 2
    col_out_b = nKs + 3
    col_out_e = nKs + 2
    dims_b[2] = dims_b[2] + 2
    dims_er[2] = dims_er[2] + 1
    dimnames_b = c(dimnames_b, "Regression", "Outlier")
    dimnames_er = c(dimnames_er, "Regression")
  } else if (Outlier) {
    col_out_b = nKs + 2
    col_out_e = nKs + 2
    dims_b[2] = dims_b[2] + 1
    dims_er[2] = dims_er[2] + 1
    dimnames_b = c(dimnames_b, "Outlier")
    dimnames_er = c(dimnames_er, "Outlier")
  } else if (reg) {
    col_reg = nKs + 2
    dims_b[2] = dims_b[2] + 1
    dimnames_b = c(dimnames_b, "Regression")
  }
  error_mat = array(
    data = NA,
    dim = dims_er,
    dimnames = list(NULL, dimnames_er)
  )
  beta_mat  = array(
    data = NA,
    dim = dims_b,
    dimnames = list(NULL, dimnames_b)
  )


  # Initializing parameters
  ## observation error
  obserror = init_sigmaE_0(y)

  ## Seasonality
  for (ik in 1:nKs) {
    sParam =  init_Sbeta(y, obserror, evol_error = "HS", Ks[[ik]])
    cn = sParam$colname
    params_list[[cn]] = sParam
    beta_mat[, cn] = sParam$s_mu
    error_mat[, cn] = c(
      sParam$s_evolParams23$sigma_w0,
      sParam$s_evolParams3k$sigma_wt,
      sParam$s_evolParamskT$sigma_wt
    ) ^ 2
    # error_mat[-1, cn] = c(
    #   sParam$s_evolParams23$sigma_w0,
    #   sParam$s_evolParams3k$sigma_wt,
    #   sParam$s_evolParamskT$sigma_wt
    # ) ^ 2
  }
  ## trend Parameter Sv
  tParam = init_Tbeta(y, obserror, evol_error, D, sparse)
  params_list[["Trend"]] = tParam
  beta_mat[, "Trend"] = tParam$s_mu
  error_mat[, "Trend"] = c(tParam$s_evolParams0$sigma_w0 ^ 2,
                           tParam$s_evolParams$sigma_wt ^ 2)

  ## Regression
  if (reg) {
    bParam = init_Regression(y, X, obserror)
    params_list[["Regression"]] = bParam
    beta_mat[, "Regression"] = bParam$s_mu
  }

  if (Outlier) {
    zParam = init_Outlier(y, obserror)
    params_list[["Outlier"]] = zParam
    beta_mat[, "Outlier"] = zParam$s_mu
    #error_mat[,"Outlier"] = zParam$s_evolParams$sigma_wt^2
    error_mat[-c(1:4), "Outlier"] = zParam$s_evolParams$sigma_wt ^ 2
  }

  if (SVm) {
    svParam = dsp_initSV((y - rowSums(beta_mat))/obserror$sigma_e)
    obserror$sigma_et = obserror$sigma_e * svParam$sigma_wt
  }else if(ASVm){
    svParam = init_paramsASV(log(((y - rowSums(beta_mat))/obserror$sigma_e)^2))
    obserror$sigma_et = obserror$sigma_e * svParam$sigma_wt
  }

  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = list()
  post_obs_sigma_t2 = array(NA, c(nsave, TT))

  #Mean Matrix
  post_s_beta = array(NA, c(nsave, dims_b), dimnames = list(NULL, NULL, dimnames_b))
  #Error Matrix
  post_s_evol_sigma_t2 = array(NA, c(nsave, dims_er), dimnames = list(NULL, NULL, dimnames_er))

  #remainder
  post_remainder = array(NA, c(nsave, TT))

  # Regression
  if (reg) {
    p = ncol(X)
    post_reg = array(NA, c(nsave, p))
  }

  #combined beta
  post_beta_combined = array(NA, c(nsave, TT))
  #yhat
  post_yhat = array(NA, c(nsave, TT))

  # Total number of MCMC simulations:
  nstot = nburn + (nskip + 1) * (nsave)
  skipcount = 0
  isave = 0 # For counting

  # Run the MCMC:
  #if(verbose) timer0 = proc.time()[3] # For timing the sampler
  if (verbose) {
    pb <- progress::progress_bar$new(
      format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
      total = nstot,
      complete = "=",
      incomplete = "-",
      current = ">",
      clear = FALSE,
      width = 100
    )
  }

  for (nsi in 1:nstot) {
    if (verbose) {
      if (nsi < 10) {
        pb$tick()
      }
      else if (((nsi %% 100) == 0)) {
        pb$tick(100)
      }
    }

    if(SVm){
      obserror = fit_sigmaE_0_m_SV(y,params_list,TT,svParam = svParam)
      svParam = dsp_sampleSVparams((y - rowSums(beta_mat))/obserror$sigma_e,svParam)
      obserror$sigma_et = obserror$sigma_e *svParam$sigma_wt

    }else if(ASVm){
      obserror = fit_sigmaE_0_m_SV(y,params_list,TT,svParam = svParam)
      svParam = fit_paramsASV(log(((y - rowSums(beta_mat))/obserror$sigma_e)^2),
                              svParam)
      obserror$sigma_et = obserror$sigma_e *svParam$sigma_wt
    }else{
      obserror = fit_sigmaE_0_m(y, params_list, TT)
    }

    if (reg) {
      bParam = fit_Regression(y - rowSums(beta_mat[, -col_reg, drop = FALSE]), X, params_list[["Regression"]], obserror)
      params_list[["Regression"]] = bParam
      beta_mat[, "Regression"] = bParam$s_mu
    }

    for (ik in 1:nKs) {
      # sParam = fit_Sbeta(y - rowSums(beta_mat[, -ik, drop = FALSE]),
      #                    params_list[[ik]],
      #                    obserror,
      #                    evol_error,
      #                    Ks[[ik]])
      sParam = fit_Sbeta(y - rowSums(beta_mat[, -ik, drop = FALSE]),
                           params_list[[ik]],
                           obserror,
                           evol_error,
                           Ks[[ik]])
      # sParam = fit_Sbeta0(y - rowSums(beta_mat[, -ik, drop = FALSE]),
      #                    params_list[[ik]],
      #                    obserror,
      #                    Ks[[ik]])
      cn = sParam$colname
      params_list[[cn]] = sParam
      beta_mat[, cn] = sParam$s_mu
      error_mat[, cn] = c(
        sParam$s_evolParams23$sigma_w0,
        sParam$s_evolParams3k$sigma_wt,
        sParam$s_evolParamskT$sigma_wt
      ) ^ 2
      # error_mat[-1, cn] = c(
      #   sParam$s_evolParams23$sigma_w0,
      #   sParam$s_evolParams3k$sigma_wt,
      #   sParam$s_evolParamskT$sigma_wt
      # ) ^ 2
    }
    tParam_data = y - rowSums(beta_mat[, -(col_trend), drop = FALSE])
    tParam = fit_Tbeta(tParam_data, tParam, obserror, evol_error, D, sparse)
    params_list[["Trend"]] = tParam
    beta_mat[, "Trend"] = tParam$s_mu
    error_mat[, "Trend"] = c(tParam$s_evolParams0$sigma_w0 ^ 2,
                             tParam$s_evolParams$sigma_wt ^ 2)

    if (Outlier) {
      zParam = fit_Outlier(y - rowSums(beta_mat[, -(col_out_b)]), zParam, obserror)
      params_list[["Outlier"]] = zParam
      beta_mat[, "Outlier"] = zParam$s_mu
      #error_mat[,"Outlier"] = zParam$s_evolParams$sigma_wt^2
      error_mat[-c(1:4), "Outlier"] = zParam$s_evolParams$sigma_wt ^ 2
    }
    #plot(log(((y - rowSums(beta_mat))/obserror$sigma_e)^2))
    # Stor the MCMC output:
    if (nsi > nburn) {
      # Increment the skip counter:
      skipcount = skipcount + 1

      # Save the iteration:
      if (skipcount > nskip) {
        # Increment the save index
        isave = isave + 1

        # observation error
        post_obs_sigma_t2[isave, ] = obserror$sigma_et ^ 2

        #Mean Matrix
        post_s_beta[isave, , ] = beta_mat

        #Error Matrix
        post_s_evol_sigma_t2[isave, , ] = error_mat

        #combined beta
        beta_combined = rowSums(beta_mat)
        post_beta_combined[isave, ] = beta_combined

        #Regression
        if (reg) {
          post_reg[isave, ] = bParam$beta
        }

        #remainder
        post_remainder[isave, ] = y - beta_combined

        #yhat
        post_yhat[isave, ] = beta_combined + obserror$sigma_et * stats::rnorm(TT)
        # And reset the skip counter:
        skipcount = 0
      }
    }
  }
  posterior_samples = list()
  posterior_samples$beta_combined = post_beta_combined
  post_s_beta[, , "Trend"] = post_s_beta[, , "Trend"] + offset_y
  posterior_samples$beta = post_s_beta

  posterior_samples$evol_sigma_t2 = post_s_evol_sigma_t2
  posterior_samples$obs_sigma_t2 = post_obs_sigma_t2
  posterior_samples$remainder = post_remainder
  if (reg) {
    posterior_samples$reg_coef = post_reg
  }
  posterior_samples$yhat = post_yhat

  mcmc_output$samples = posterior_samples
  return(mcmc_output)

}
#' @keywords internal
#' @importFrom Matrix chol diag t
robust_cholesky <- function(matrix) {
  tryCatch({
        # Attempt Cholesky decomposition
        return(suppressWarnings(Matrix::chol(matrix)))  # Return decomposition if successful
        },
      error = function(e) {
        saveRDS(matrix,file = "matrix_pert.rds")
        sus = which(rowSums(matrix)<0)
        for(i in sus){
          rowsum = sum(matrix[i,])
          while(rowsum <0){
            matrix[i,i]  = matrix[i,i] - rowsum*10
            rowsum = sum(matrix[i,])
          }
        }
        return(suppressWarnings(Matrix::chol(matrix)))
      }
  )
}
#' @keywords internal
#' @importFrom stats quantile
summarize_output <- function(mcmc_output,y,Ks,cl,reg,Outlier){
  nKs = length(Ks)
  beta = mcmc_output$beta
  remainder = mcmc_output$remainder
  #summary statistics
  posterior_summary = list()
  #
  posterior_mean = data.frame(
    y = y,
    apply(beta, c(2, 3), mean),
    Remainder = colMeans(remainder)
  )
  posterior_summary$p_means = posterior_mean
  alpha_d2 = (1 - cl) / 2
  lower = apply(beta, c(2, 3), stats::quantile, alpha_d2)
  upper = apply(beta, c(2, 3), stats::quantile, 1 - alpha_d2)
  posterior_summary$Trend_sum = data.frame(Mean = posterior_mean[, "Trend"],
                                           CR_lower = lower[, "Trend"],
                                           CR_upper = upper[, "Trend"])
  for (ik in 1:nKs) {
    nm = paste0("Seasonal", Ks[[ik]])
    posterior_summary[[paste0(nm, "_sum")]] = data.frame(Mean = posterior_mean[, nm],
                                                         CR_lower = lower[, nm],
                                                         CR_upper = upper[, nm])
  }
  if (reg) {
    posterior_summary$Regression_sum = data.frame(Mean = posterior_mean[, "Regression"],
                                                  CR_lower = lower[, "Regression"],
                                                  CR_upper = upper[, "Regression"])
  }
  if (Outlier) {
    posterior_summary$Outlier_sum = data.frame(Mean = posterior_mean[, "Outlier"],
                                               CR_lower = lower[, "Outlier"],
                                               CR_upper = upper[, "Outlier"])
  }
  # summarizing volatility
  posterior_summary$Volatility = data.frame(
    Mean = apply(sqrt(mcmc_output$obs_sigma_t2),2,mean),
    CR_lower =  apply(sqrt(mcmc_output$obs_sigma_t2),2,quantile,alpha_d2),
    CR_upper =  apply(sqrt(mcmc_output$obs_sigma_t2),2,quantile,1 - alpha_d2))

  return(posterior_summary)
}
#' @keywords internal
debug_precision <- function(QHt_Matrix){
  sus = which(rowSums(QHt_Matrix)<0)
  for(i in sus){
    rowsum = sum(QHt_Matrix[i,])
    while(rowsum <0){
      QHt_Matrix[i,i]  = QHt_Matrix[i,i] - rowsum*10
      rowsum = sum(QHt_Matrix[i,])
    }
  }
  return(QHt_Matrix)
}
#' @keywords internal
coverage <- function(upper,lower,a){
  n = length(a)
  return(sum((upper > a) & (lower < a))/n)
}
#' @keywords internal
run_with_retries <- function(func, retries = 10, delay = 1, ...) {
  for (attempt in 1:retries) {
    tryCatch(
      {
        # Attempt to execute the function with additional arguments
        result <- func(...)
        return(result)  # If successful, return the result
      },
      error = function(e) {
        if (attempt == retries) {
          stop("Function failed after ", retries, " attempts: ", e$message)
        } else {
          message("Error on attempt ", attempt, ": ", e$message)
          Sys.sleep(delay)  # Wait before retrying
        }
      }
    )
  }
}
#' Decomposition of time series data with Bayesian Adaptive Seasonality Trend decomposition Incorporating Outlier and Noise (BASTION)
#'
#' decompose time series data using BASTION. The time series should be a vector.
#' @param y numeric vector of the \code{T x 1} vector of time series observations
#' @param Ks list of values containing the seasonal periods
#' @param X matrix for additional covariate for regression (default is NULL)
#' @param Outlier logical flag (default is FALSE) to model outlier
#' @param cl scalar between (0,1) for confidence leve (default is 0.95)
#' @param sparse logical flag (default is FALSE) to induce additional shrinkage to the trend estimate
#' @param obsSV options for modeling the error variance. It must be one of the following:
#' \itemize{
#' \item const: Constant error variance for all time points.
#' \item SV: Stochastic Volatility model.
#' }
#' @param nchains integer scalar for the number of chains for MCMC sampling (default is 2)
#' @param nsave integer scalar (default = 1000); number of MCMC iterations to record
#' @param nburn integer scalar (default = 1000); number of MCMC iterations to discard (burn-in)
#' @param nskip integer scalar (default = 4); number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param verbose logical; report extra information on progress if true. Defaults to FALSE
#' @param save_samples logical; save and return posterior samples of each components if true. Default is FALSE
#' @return \code{fit_BASTION} returns an object class list.
#'
#' @section `summary`:
#' A list providing summarized posterior estimates:
#' - `p_means`: A matrix containing the posterior means of the observed data (`y`) and each decomposed component:
#'   **Trend**, **Seasonality**, and **Outlier**, if applicable.
#' - `Trend_sum`: A matrix containing the posterior mean and 95% credible interval of the trend estimate.
#' - `Seasonal"k"_sum`: A matrix containing the posterior mean and 95% credible interval of the seasonal estimate,
#'   where `"k"` is determined by the input.
#' - `Outlier_sum`: A matrix containing the posterior mean and 95% credible interval of the outlier component.
#'
#' @section `samples`: (returned only when `save_samples` = TRUE)
#' A list containing MCMC samples of the relevant parameters:
#' - `beta_combined`: Posterior samples of **Trend + Seasonality**.
#' - `beta`: Posterior samples of each component excluding the remainder.
#' - `obs_sigma_t2`: Posterior samples of the variance of the observation equation.
#' - `evol_sigma_t2`: Posterior samples of the evolution error term.
#' - `remainder`: Posterior samples of the remainder term.
#' - `yhat`: Posterior samples of the signal + error term.
#'
#' @import progress spam abind
#' @export
fit_BASTION = function(y,Ks,X=NULL,Outlier=FALSE,cl=0.95,sparse = FALSE,obsSV = "const",
                       nchains = 2,nsave = 1000, nburn= 1000, nskip = 4,
                       verbose = TRUE,save_samples = FALSE){
  #v2
  y = as.vector(y)
  if (!is.list(Ks)) {
    stop("Ks needs to be a list")
  }
  if (!(obsSV %in% c("const", "SV"))) {
    stop("obsSV needs to be either SV or const")
  }
  reg = !is.null(X) # is it regression
  if (reg & !is.matrix(X)) {
    stop("X needs to be a matrix")
  }
  elements_to_combine <- list(
    beta_combined = rbind,
    beta = function(x, y) abind::abind(x, y, along = 1),
    evol_sigma_t2 = function(x, y) abind::abind(x, y, along = 1),
    obs_sigma_t2 = rbind,
    remainder = rbind,
    Yhat = rbind
  )

  for(i in 1:nchains){
    print(paste("Chain",i))
    # model = run_with_retries(fit_ASD,
    #                          retries = retries,
    #                          delay = 1,
    #                          y = y,
    #                          X = X,
    #                          Ks = Ks,
    #                          Outlier = Outlier,
    #                          obsSV = obsSV,
    #                          ...)
    model = fit_ASD(y = y,
                    X = X,
                    Ks = Ks,
                    Outlier = Outlier,
                    obsSV = obsSV,
                    nsave = nsave,
                    nburn = nburn,
                    nskip = nskip,
                    verbose = verbose)
    if(i ==1){
      combined_samples = model
    }else{
      for (element in names(elements_to_combine)) {
        combined_samples$samples[[element]] <- elements_to_combine[[element]](
          combined_samples$samples[[element]],
          model$samples[[element]]
        )
      }
      if(reg){
        combined_samples$samples$reg_coef <- rbind(combined_samples$samples$reg_coef,
                                                   model$samples$reg_coef)
      }
    }
  }
  summary = summarize_output(mcmc_output = combined_samples$samples,
                             y = y,
                             Ks = Ks,
                             cl = cl,
                             reg = reg,
                             Outlier = Outlier)
  if(save_samples){
    return(list(summary = summary,
                samples = combined_samples$samples))
  }else{
    list(summary = summary)
  }
}


