#' @keywords internal
dsp_initEvol0 = function(mu0, commonSD = TRUE){

  p = length(mu0)

  # Common or distinct:
  if(commonSD) {
    sigma_w0 = rep(mean(abs(mu0)), p)
  } else  sigma_w0 = abs(mu0)

  # Initialize at 1 for simplicity:
  px_sigma_w0 = rep(1, p)

  sigma_00 = px_sigma_00 = 1

  list(sigma_w0 = sigma_w0, px_sigma_w0 = px_sigma_w0, sigma_00 = sigma_00, px_sigma_00 = px_sigma_00)
}
#' @keywords internal
#' @importFrom stats mad rgamma
dsp_sampleEvol0 = function(mu0, evolParams0, commonSD = FALSE, A = 1){
  # Store length locally:
  p = length(mu0)

  # For numerical stability:
  mu02offset = any(mu0^2 < 10^-16)*max(10^-8, mad(mu0)/10^6)
  mu02 = mu0^2 + mu02offset

  if(commonSD){
    # (Common) standard deviations:
    evolParams0$sigma_w0 = rep(1/sqrt(rgamma(n = 1, shape = p/2 + 1/2, rate = sum(mu02)/2 + evolParams0$px_sigma_w0[1])), p)

    # (Common) paramater expansion:
    evolParams0$px_sigma_w0 = rep(rgamma(n = 1, shape = 1/2 + 1/2, rate = 1/evolParams0$sigma_w0[1]^2 + 1/A^2), p)

  } else {
    # (Distinct) standard deviations:
    evolParams0$sigma_w0 = 1/sqrt(rgamma(n = p, shape = 1/2 + 1/2, rate = mu02/2 + evolParams0$px_sigma_w0))

    # (Distinct) paramater expansion:
    #evolParams0$px_sigma_w0 = rgamma(n = p, shape = 1/2 + 1/2, rate = 1/evolParams0$sigma_w0^2 + 1/A^2)
    evolParams0$px_sigma_w0 = rgamma(n = p, shape = 1/2 + 1/2, rate = 1/evolParams0$sigma_w0^2 + 1/evolParams0$sigma_00^2)

    # Global standard deviations:
    evolParams0$sigma_00 = 1/sqrt(rgamma(n = 1, shape = p/2 + 1/2, rate = sum(evolParams0$px_sigma_w0) + evolParams0$px_sigma_00))

    # (Global) parameter expansion:
    evolParams0$px_sigma_00 = rgamma(n = 1, shape = 1/2 + 1/2, rate = 1/evolParams0$sigma_00^2 + 1/A^2)
  }

  # And return the list:
  evolParams0
}
#' @keywords internal
#' @importFrom stats sd
#' @importFrom mgcv rig
dsp_initEvolParams = function(omega, evol_error = "HS"){
  # Check:
  if(!((evol_error == "HS") || (evol_error == "BL") ||(evol_error == "NIG"))){
    stop('Error type must be one of DHS, HS, BL, SV, or NIG')
  }

  # Make sure omega is (n x p) matrix
  omega = as.matrix(omega); n = nrow(omega); p = ncol(omega)

  if(evol_error == "HS"){
    tauLambdaj = 1/omega^2;
    xiLambdaj = 1/(2*tauLambdaj); tauLambda = 1/(2*colMeans(xiLambdaj)); xiLambda = 1/(tauLambda + 1)

    # Parameters to store/return:
    return(list(sigma_wt = 1/sqrt(tauLambdaj), tauLambdaj = tauLambdaj, xiLambdaj = xiLambdaj, tauLambda = tauLambda, xiLambda = xiLambda))
  }
  if(evol_error == "BL"){
    tau_j = abs(omega); lambda2 = mean(tau_j)
    return(list(sigma_wt = tau_j, tau_j = tau_j, lambda2 = lambda2))
  }
  if(evol_error == "NIG")
    return(list(sigma_wt = tcrossprod(rep(1,n), apply(omega, 2, function(x) sd(x, na.rm=TRUE)))))
}
#' @keywords internal
#' @importFrom stats mad rgamma
dsp_sampleEvolParams = function(omega, evolParams,  sigma_e = 1, evol_error = "HS", loc = NULL){

  # Check:
  if(!((evol_error == "HS") || (evol_error == "BL") || (evol_error == "NIG"))) stop('Error type must be one of HS, BL, or NIG')

  # Make sure omega is (n x p) matrix
  omega = as.matrix(omega); n = nrow(omega); p = ncol(omega)

  if(evol_error == "HS"){

    # For numerical reasons, keep from getting too small
    hsOffset = tcrossprod(rep(1,n), apply(omega, 2, function(x) any(x^2 < 10^-16)*max(10^-8, mad(x)/10^6)))
    hsInput2 = omega^2 + hsOffset

    # Local scale params:
    evolParams$tauLambdaj = matrix(rgamma(n = n*p, shape = 1, rate = evolParams$xiLambdaj + hsInput2/2), nrow = n)
    evolParams$xiLambdaj = matrix(rgamma(n = n*p, shape = 1, rate = evolParams$tauLambdaj + tcrossprod(rep(1,n), evolParams$tauLambda)), nrow = n)

    # Global scale params:
    evolParams$tauLambda = rgamma(n = p, shape = 0.5 + n/2, colSums(evolParams$xiLambdaj) + evolParams$xiLambda)
    #evolParams$xiLambda = rgamma(n = p, shape = 1, rate = evolParams$tauLambda + 1/sigma_e^2)
    evolParams$xiLambda = rgamma(n = p, shape = 1, rate = evolParams$tauLambda + 1)

    evolParams$sigma_wt = 1/sqrt(evolParams$tauLambdaj)

    return(evolParams)
  }
  if(evol_error == "BL"){

    # For numerical reasons, keep from getting too small
    hsOffset = tcrossprod(rep(1,n), apply(omega, 2, function(x) any(x^2 < 10^-16)*max(10^-8, mad(x)/10^6)))
    hsInput2 = omega^2 + hsOffset

    # 1/tau_j^2 is inverse-gaussian (NOTE: this is very slow!)
    evolParams$tau_j = matrix(sapply(matrix(hsInput2), function(x){1/sqrt(rig(n = 1,
                                                                                    mean = sqrt(evolParams$lambda2*sigma_e^2/x), # already square the input
                                                                                    scale = 1/evolParams$lambda2))}), nrow = n)
    # Note: should be better priors for lambda2
    evolParams$lambda2 = rgamma(n = 1,
                                shape = 1 + n*p,
                                rate = 2 + sum(evolParams$tau_j^2)/2)

    # For Bayesian lasso, scale by sigma_e:
    evolParams$sigma_wt = sigma_e*evolParams$tau_j

    return(evolParams)
  }
  if(evol_error == "SV") return(dsp_sampleSVparams(omega = omega, svParams = evolParams))
  if(evol_error == "NIG") {
    evolParams = list(sigma_wt = tcrossprod(rep(1,n),
                                            apply(omega, 2,
                                                  function(x) 1/sqrt(rgamma(n = 1, shape = n/2 + 0.01, rate = sum(x^2)/2 + 0.01)))))
    return(evolParams)
  }
}
#' @keywords internal
#' @importFrom stats arima
#' @import methods
dsp_initSV = function(omega){

  # Make sure omega is (n x p) matrix
  omega = as.matrix(omega); n = nrow(omega); p = ncol(omega)

  # log-volatility:
  ht = log(omega^2 + 0.0001)

  # AR(1) pararmeters: check for error in initialization too
  svParams = apply(ht, 2, function(x){
    ar_fit = try(arima(x, c(1,0,0)), silent = TRUE)
    if(methods::is(ar_fit, "try-error")) {
      params = c(ar_fit$coef[2], ar_fit$coef[1], sqrt(ar_fit$sigma2))
    } else params = c(mean(x)/(1 - 0.8),0.8, 1)
    params
  }); rownames(svParams) = c("intercept", "ar1", "sig")

  # SDs, log-vols, and other parameters:
  return(list(sigma_wt = exp(ht/2), ht = ht, svParams = svParams))
}
#' @keywords internal
#' @import stochvol
dsp_sampleSVparams = function(omega, svParams){

  # Make sure omega is (n x p) matrix
  omega = as.matrix(omega); n = nrow(omega); p = ncol(omega)

  for(j in 1:p){
    # First, check for numerical issues:
    svInput = omega[,j]; #if(all(svInput==0)) {svInput = 10^-8} else svInput = svInput + sd(svInput)/10^8
    #hsOffset = tcrossprod(rep(1,n), apply(omega, 2, function(x) any(x^2 < 10^-16)*max(10^-8, mad(x)/10^6)))

    # Sample the SV parameters:
    svsamp = stochvol::svsample_fast_cpp(svInput,
                                         startpara = list(
                                           mu = svParams$svParams[1,j],
                                           phi = svParams$svParams[2,j],
                                           sigma = svParams$svParams[3,j]),
                                         startlatent = svParams$ht[,j])# ,priorphi = c(10^4, 10^4));
    # Update the parameters:
    svParams$svParams[,j] = svsamp$para[1:3];
    svParams$ht[,j] = svsamp$latent
  }
  # Finally, up the evolution error SD:
  svParams$sigma_wt = exp(svParams$ht/2)

  # Check for numerically large values:
  svParams$sigma_wt[which(svParams$sigma_wt > 10^3, arr.ind = TRUE)] = 10^3

  return(svParams)
}
#' @keywords internal
sample_jfast <- function(T,obs= NULL){
  m_st  = c(1.92677, 1.34744, 0.73504, 0.02266, -0.85173, -1.97278, -3.46788, -5.55246, -8.68384, -14.65000)
  v_st2 = c(0.11265, 0.17788, 0.26768, 0.40611,  0.62699,  0.98583,  1.57469,  2.54498,  4.16591,   7.33342)
  q     = c(0.00609, 0.04775, 0.13057, 0.20674,  0.22715,  0.18842,  0.12047,  0.05591,  0.01575,   0.00115)
  #z = draw.indicators(res = ystar-h_prev, nmix = list(m = m_st, v = v_st2, p = q))
  if(is.null(obs)){
    z = sample.int(10,T,replace = T, prob = q)
  }else{
    z = sapply(obs, ncind, m_st, sqrt(v_st2), q)
  }
  # Subset mean and variances to the sampled mixture components; (n x p) matrices
  return(
    data.frame(mean = m_st[z],
               var = v_st2[z]))
}
#' @keywords internal
#' @importFrom stats dnorm
ncind = function(y,mu,sig,q){
  tryCatch({
    sample(1:length(q),
           size = 1,
           prob = q*dnorm(y,mu,sig))
  },
  error = function(e){
    which.max(q*dnorm(y,mu,sig))
  }
  )

}
#' @keywords internal
sample_j_wrap <- function(Td,obs=NULL){
  # Omori, Chib, Shephard, Nakajima (2007) 10-component mixture:
  m_st  = c(1.92677, 1.34744, 0.73504, 0.02266, -0.85173, -1.97278, -3.46788, -5.55246, -8.68384, -14.65000)
  v_st2 = c(0.11265, 0.17788, 0.26768, 0.40611,  0.62699,  0.98583,  1.57469,  2.54498,  4.16591,   7.33342)
  q     = c(0.00609, 0.04775, 0.13057, 0.20674,  0.22715,  0.18842,  0.12047,  0.05591,  0.01575,   0.00115)

  # if(is.null(obs)){
  #   z = sample.int(10,Td,replace = Td, prob = q)
  # }else{
  #   z = c(draw_indicators_generic(obs, rep(0, Td), Td, q, m_st, sqrt(v_st2), 10))
  # }
  z = sample.int(10,Td,replace = Td, prob = q)
  return(
    data.frame(mean = m_st[z],
               var = v_st2[z]))

}
