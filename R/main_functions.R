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
    error_mat[, "Outlier"] = zParam$s_evolParams$sigma_wt ^ 2
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
    pb <- .bastion_progress_bar(nstot)
  }

  for (nsi in 1:nstot) {
    if (verbose) .bastion_progress_tick(pb,nsi,nstot)

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
      error_mat[, "Outlier"] = zParam$s_evolParams$sigma_wt ^ 2
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
#' @importFrom stats approxfun rnorm
fit_ASD_SV = function(y,
                      Ks,
                      X = NULL,
                      Outlier = FALSE,
                      sparse = FALSE,
                      obsSV = "const",
                      nsave = 1000,
                      nburn = 1000,
                      nskip = 4,
                      verbose = TRUE, D_trend = 2L,
                      seasonality_model = c("seasonality3", "seasonality")) {
  seasonality_model = match.arg(seasonality_model)
  D = .validate_D_trend(D_trend);evol_error = "HS";
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
  norm_y = sd(y)
  y = (y - offset_y)/norm_y

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
    dimnames_er = c(dimnames_er, "Outlier")
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

  #initializing estimates
  # yts = forecast::msts(y,seasonal.periods = unlist(Ks))
  # ie = forecast::mstl(yts)
  # Initializing parameters
  ## observation error
  obserror = init_sigmaE_0(y)

  ## Seasonality
  for (ik in 1:nKs) {
    # Differenced-state sampler: sParam = init_Sbeta(y, obserror, evol_error="HS", Ks[[ik]])
    sParam = if(seasonality_model == "seasonality") {
      init_Sbeta(y, obserror, evol_error="HS", k=Ks[[ik]])
    } else {
      init_Sbeta3(y, obserror, evol_error="HS", k=Ks[[ik]])
    }
    cn = sParam$colname
    params_list[[cn]] = sParam
    beta_mat[, cn] = sParam$s_mu
    # The former T innovation variances do not exist in this model.
    error_mat[, cn] = if(seasonality_model == "seasonality") c(
        sParam$s_evolParams23$sigma_w0, sParam$s_evolParams3k$sigma_wt,
        sParam$s_evolParamskT$sigma_wt)^2 else NA_real_
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
    error_mat[, "Outlier"] = zParam$s_evolParams$sigma_wt ^ 2
  }

  if (SVm) {
    svParam = dsp_initSV(y - rowSums(beta_mat))
    obserror$sigma_et = svParam$sigma_wt
    obserror$sigma_e = 1
  }else if(ASVm){
    svParam = init_paramsASV(log((y - rowSums(beta_mat))^2))
    obserror$sigma_et = svParam$sigma_wt
    obserror$sigma_e = 1
  }

  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = list()
  post_obs_sigma_t2 = array(NA, c(nsave, TT))

  #Mean Matrix
  post_s_beta = array(NA, c(nsave, dims_b), dimnames = list(NULL, NULL, dimnames_b))
  #Error Matrix
  post_s_evol_sigma_t2 = array(NA, c(nsave, dims_er), dimnames = list(NULL, NULL, dimnames_er))

  post_template_variance = lapply(Ks, function(k) matrix(NA_real_, nsave, k))
  names(post_template_variance) = paste0("Seasonal", Ks)
  post_common_seasonality = lapply(Ks, function(k) matrix(NA_real_, nsave, k))
  names(post_common_seasonality) = paste0("Seasonal", Ks)
  post_sv_parameters = if(SVm) matrix(NA_real_, nsave, 3L,
    dimnames=list(NULL,c("mu","phi","innovation_sd"))) else NULL
  post_delta_variance = matrix(NA_real_, nsave, nKs,
    dimnames=list(NULL,paste0("Seasonal", Ks)))
  if(seasonality_model == "seasonality3") {
    post_delta_global_variance = matrix(NA_real_, nsave, nKs,
      dimnames=list(NULL,paste0("Seasonal", Ks)))
    post_delta_local_variance = lapply(Ks,function(k)
      matrix(NA_real_,nsave,ceiling(TT/k)-1L))
    names(post_delta_local_variance) = paste0("Seasonal", Ks)
    post_delta_innovation_variance = post_delta_local_variance
    post_delta_sampler_backend = matrix(NA_character_,nsave,nKs,
      dimnames=list(NULL,paste0("Seasonal", Ks)))
  }

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
    pb <- .bastion_progress_bar(nstot)
  }

  for (nsi in 1:nstot) {
    if (verbose) .bastion_progress_tick(pb,nsi,nstot)

    if(SVm){
      #obserror = fit_sigmaE_0_m_SV(y,params_list,TT,svParam = svParam)
      svParam = dsp_sampleSVparams(y - rowSums(beta_mat),svParam)
      obserror$sigma_et = svParam$sigma_wt

    }else if(ASVm){
      #obserror = fit_sigmaE_0_m_SV(y,params_list,TT,svParam = svParam)
      svParam = fit_paramsASV(log(((y - rowSums(beta_mat)))^2),svParam)
      obserror$sigma_et = svParam$sigma_wt

    }else{
      # The raw seasonal template has rank k; its independent delta
      # variance does not contribute to the sigma_e conditional.
      sigma_params = lapply(params_list, function(param) {
        if (!is.null(param$prior_rank)) param$Td = param$prior_rank
        param
      })
      obserror = fit_sigmaE_0_m(y, sigma_params, TT)
    }

    if (reg) {
      bParam = fit_Regression(y - rowSums(beta_mat[, -col_reg, drop = FALSE]),
                              X,
                              params_list[["Regression"]],
                              obserror)
      params_list[["Regression"]] = bParam
      beta_mat[, "Regression"] = bParam$s_mu
    }

    for (ik in 1:nKs) {
      # Differenced-state sampler: sParam = fit_Sbeta(residual, params_list[[ik]], obserror, evol_error, Ks[[ik]])
      sParam = if(seasonality_model == "seasonality") {
        fit_Sbeta(y - rowSums(beta_mat[, -ik, drop=FALSE]),
          params_list[[ik]], obserror, evol_error, k=Ks[[ik]])
      } else {
        fit_Sbeta3(y - rowSums(beta_mat[, -ik, drop=FALSE]),
          params_list[[ik]], obserror, evol_error, k=Ks[[ik]])
      }
      cn = sParam$colname
      params_list[[cn]] = sParam
      beta_mat[, cn] = sParam$s_mu
      error_mat[, cn] = if(seasonality_model == "seasonality") c(
        sParam$s_evolParams23$sigma_w0, sParam$s_evolParams3k$sigma_wt,
        sParam$s_evolParamskT$sigma_wt)^2 else NA_real_
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
      error_mat[, "Outlier"] = zParam$s_evolParams$sigma_wt ^ 2
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
        if(SVm)post_sv_parameters[isave, ] = as.numeric(svParam$svParams[,1L])

        #Mean Matrix
        post_s_beta[isave, , ] = beta_mat

        #Error Matrix
        post_s_evol_sigma_t2[isave, , ] = error_mat

        for (ik in seq_len(nKs)) {
          seasonal_state = params_list[[ik]]
          if(seasonality_model == "seasonality") next
          post_common_seasonality[[ik]][isave, ] = seasonal_state$template * norm_y
          post_template_variance[[ik]][isave, ] = obserror$sigma_e^2 * c(
            rep(seasonal_state$s_evolParams0$sigma_w0[1]^2, 2L),
            seasonal_state$s_evolParamsShape$sigma_wt^2) * norm_y^2
          post_delta_variance[isave, ik] = seasonal_state$delta_variance * norm_y^2
          if(seasonality_model == "seasonality3") {
            post_delta_global_variance[isave,ik] = seasonal_state$delta_hs$tau2 * norm_y^2
            post_delta_local_variance[[ik]][isave,] = seasonal_state$delta_hs$lambda2
            post_delta_innovation_variance[[ik]][isave,] = seasonal_state$delta_innovation_variance * norm_y^2
            post_delta_sampler_backend[isave,ik] = seasonal_state$delta_sampler_backend
          }
        }

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
  posterior_samples$beta_combined = post_beta_combined*norm_y
  post_s_beta = post_s_beta*norm_y
  post_s_beta[, , "Trend"] = post_s_beta[, , "Trend"] + offset_y
  posterior_samples$beta = post_s_beta

  posterior_samples$common_seasonality = post_common_seasonality
  if(SVm) {
    post_sv_parameters[,"mu"] = post_sv_parameters[,"mu"] + 2*log(norm_y)
    posterior_samples$sv_parameters = post_sv_parameters
  }
  posterior_samples$seasonal_template_variance = post_template_variance
  posterior_samples$seasonal_delta_variance = post_delta_variance
  if(seasonality_model == "seasonality3") {
    posterior_samples$seasonal_delta_global_variance = post_delta_global_variance
    posterior_samples$seasonal_delta_local_variance = post_delta_local_variance
    posterior_samples$seasonal_delta_innovation_variance = post_delta_innovation_variance
    posterior_samples$seasonal_delta_sampler_backend = post_delta_sampler_backend
  }
  posterior_samples$evol_sigma_t2 = post_s_evol_sigma_t2
  posterior_samples$obs_sigma_t2 = robust_prod(post_obs_sigma_t2,norm_y^2)
  posterior_samples$remainder = post_remainder
  if (reg) {
    posterior_samples$reg_coef = post_reg
  }
  posterior_samples$yhat = post_yhat

  mcmc_output$samples = posterior_samples
  mcmc_output$normalization = list(mean=offset_y,sd=norm_y)
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
        message("Numerically Unstable...")
        #saveRDS(matrix,file = "matrix_pert.rds")
        sus = which(rowSums(matrix)<0)
        for(i in sus){
          rowsum = sum(matrix[i,])
          while(rowsum <0){
            matrix[i,i]  = matrix[i,i] - rowsum*100
            rowsum = sum(matrix[i,])
          }
        }
        message("stabilized matrix")
        return(suppressWarnings(Matrix::chol(matrix)))
      }
  )
}
#' @keywords internal
robust_prod <- function(a,b){
  return(exp(log(a) + log(b)))
}
#' @keywords internal
robust_div <- function(a,b){
  return(exp(log(a) - log(b)))
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
#' @return An object of class `bastion_fit` with `time`, `observed`,
#'   `components`, `intervals`, `seasonality`, `noise`, optional `regression`,
#'   and optional `draws`. `components` is a T by M data frame of posterior
#'   means, with matching lower/upper data frames under `intervals`.
#' @details BASTION decomposes a time series into a locally adaptive trend, one or more
#' seasonal components, and optional regression and outlier components. The
#' components are estimated jointly using Bayesian inference.
#'
#' Supply seasonal periods in numbers of observations through `Ks`. For example,
#' use `list(12)` for yearly seasonality in monthly data, or `list(7, 365)` for
#' weekly and yearly seasonality in daily data. Each seasonal component combines
#' a common recurring pattern with deviations that allow the pattern to evolve
#' across cycles.
#'
#' The trend uses second differences by default (`D_trend = 2`). Use `D_trend =
#' 1` for first differences. Setting `sparse = TRUE` applies stronger shrinkage
#' to trend differences while retaining local adaptivity.
#'
#' Use `obsSV = "const"` for constant observation variance or `obsSV = "SV"` to
#' estimate observation variability over time. Set `Outlier = TRUE` to estimate a
#' separate outlier component, and supply `X` to include regression effects.
#'
#' Component estimates and observation standard deviations are returned in the
#' units of the supplied data. Posterior means and credible intervals at the
#' level specified by `cl` are returned whether or not posterior draws are
#' retained. Set `save_samples = TRUE` to retain the draws.
#'
#' Use `print()` for an overview of the fit, `summary()` for the component
#' estimates, `plot()` for the decomposition, and `plot_component()` to inspect
#' an individual component with credible intervals.
#'
#' @param D_trend Trend difference order: 1 for first differences or 2 for
#'   second differences (default). Does not change seasonal differencing.
#' @import progress spam abind
#' @seealso [simulate_BASTION()]
#' @examples
#' \dontrun{
#' # Daily observations with fixed weekly and yearly periods
#' sim <- simulate_BASTION("daily", seed = 20)
#' fit_daily <- fit_BASTION(
#'   y = sim$y,
#'   Ks = list(7, 365),
#'   obsSV = "const", nchains = 1
#' )
#' plot(fit_daily)
#' plot_component(fit_daily, "seasonal_365", part = "recurring")
#'
#' # A pattern repeating every 14 observations
#' sim <- simulate_BASTION("biweekly", seed = 21)
#' fit_biweekly <- fit_BASTION(
#'   y = sim$y,
#'   Ks = list(14),
#'   obsSV = "const", nchains = 1
#' )
#' plot_component(fit_biweekly, "seasonal_14")
#'
#' # Hourly observations with daily and weekly periods
#' sim <- simulate_BASTION("hourly", seed = 22)
#' fit_hourly <- fit_BASTION(
#'   y = sim$y,
#'   Ks = list(24, 168),
#'   obsSV = "const", nchains = 1
#' )
#' plot(fit_hourly)
#' }
#' @export
fit_BASTION = function(y,Ks,X=NULL,Outlier=FALSE,cl=0.95,sparse = FALSE,obsSV = "const",
                       nchains = 2,nsave = 1000, nburn= 1000, nskip = 4,
                       verbose = TRUE,save_samples = FALSE,D_trend = 2L){
  # Shared recurring templates and group-horseshoe seasonal differences.
  seasonality_model = "seasonality3"
  # seasonality_model = "seasonality" # Differenced-state seasonal model.
  D_trend = .validate_D_trend(D_trend)
  observation_time = if(stats::is.ts(y)) as.numeric(stats::time(y)) else seq_along(y)
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
    common_seasonality = function(x,y) Map(rbind,x,y),
    seasonal_template_variance = function(x,y) Map(rbind,x,y),
    seasonal_delta_variance = rbind,
    sv_parameters = rbind,
    obs_sigma_t2 = rbind,
    remainder = rbind,
    yhat = rbind
  )
  if(seasonality_model == "seasonality3") {
    elements_to_combine = c(elements_to_combine,list(
      seasonal_delta_global_variance = rbind,
      seasonal_delta_local_variance = function(x,y)Map(rbind,x,y),
      seasonal_delta_innovation_variance = function(x,y)Map(rbind,x,y),
      seasonal_delta_sampler_backend = rbind))
  }

  for(i in 1:nchains){
    .bastion_progress_chain(i,verbose)
    model = fit_ASD_SV(y = y,
                       X = X,
                       Ks = Ks,
                       Outlier = Outlier,
                       sparse = sparse,
                       obsSV = obsSV,
                       nsave = nsave,
                       nburn = nburn,
                       nskip = nskip,
                       verbose = verbose,
                       D_trend = D_trend,
                       seasonality_model = seasonality_model)
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
  .bastion_fitted_object(combined_samples$samples,y,observation_time,Ks,X,
    Outlier,obsSV,cl,nchains,nsave,nburn,nskip,save_samples,
    combined_samples$normalization$sd,
    input_metadata=list(Ks=Ks),sparse=sparse,D_trend=D_trend)
}



# Validate before initialization or random-number generation.
.validate_D_trend <- function(D_trend) {
  if(!is.numeric(D_trend) || length(D_trend)!=1L ||
     !is.finite(D_trend) || !(D_trend %in% c(1,2)))
    stop("D_trend must be 1 or 2")
  as.integer(D_trend)
}
