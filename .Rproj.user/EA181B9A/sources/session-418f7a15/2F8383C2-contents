#' @keywords internal
#' @importFrom stats sd
init_sigmaE_0 <- function(data) {
  T = length(data)
  sigma_e = sd(data, na.rm = TRUE)
  sigma_et = rep(sigma_e, T)
  return(list(sigma_e = sigma_e, sigma_et = sigma_et))
}
#' @keywords internal
#' @import extraDistr
fit_sigmaE_0 <- function(data, tParam) {
  T = length(data)
  s_mu = tParam$s_mu
  s_omega = tParam$s_omega
  s_mu0 = tParam$s_mu0

  lambda0_2 = tParam$s_evolParams0$sigma_w0 ^ 2
  omega_sigma_wt2 = tParam$s_evolParams$sigma_wt ^ 2
  a = T + 0.001
  b = 0.001 + 1 / 2 * (sum((data - s_mu) ^ 2) + sum(s_omega ^ 2 / omega_sigma_wt2) + sum(s_mu0 ^
                                                                                           2 / lambda0_2))
  sigma_e = sqrt(extraDistr::rinvgamma(1, a, b))
  sigma_et = rep(sigma_e, T)
  return(list(sigma_e = sigma_e, sigma_et = sigma_et))
}
#' @keywords internal
#' @import extraDistr
fit_sigmaE_0_m <- function(data,
                           params_list,
                           TT,
                           a = 0,
                           b = 0) {
  a = a
  b = b
  offset = data
  for (param in params_list) {
    offset = offset - param$s_mu
    a = a + param$Td / 2
    b = b + param$n_squared_sum / 2
  }
  a = a + TT / 2
  b = b + sum((offset) ^ 2) / 2
  sigma_e = sqrt(extraDistr::rinvgamma(1, a, b))
  sigma_et = rep(sigma_e, TT)
  return(list(sigma_e = sigma_e, sigma_et = sigma_et))
}
#' @keywords internal
#' @import extraDistr
fit_sigmaE_0_m_SV <- function(data,
                              params_list,
                              TT,
                              svParam,
                              a = 0,
                              b = 0) {
  a = a
  b = b
  offset = data
  for (param in params_list) {
    offset = offset - param$s_mu
    a = a + param$Td / 2
    b = b + param$n_squared_sum / 2
  }
  a = a + TT / 2
  b = b + sum((offset/svParam$sigma_wt)^2)/2
  sigma_e = sqrt(extraDistr::rinvgamma(1, a, b))
  sigma_et = rep(sigma_e, TT)
  return(list(sigma_e = sigma_e, sigma_et = sigma_et))
}
#' @keywords internal
init_paramsASV <- function(data, D = 2){
  Td = length(data)
  s_p_error_term = sample_jfast(Td)
  s_mu = sampleTrend(data- s_p_error_term$mean,
                     obs_sigma_t2 = s_p_error_term$var,
                     evol_sigma_t2 = 0.01*rep(1,Td),
                     D = D,
                     Td = Td)
  s_omega = diff(s_mu,differences = D)
  s_mu0 = as.matrix(s_mu[1:D,])
  s_evolParams0 = dsp_initEvol0(s_mu0)
  s_evolParams = dsp_initEvolParams(s_omega,"HS")
  return(list(
    s_p_error_term = s_p_error_term,
    s_mu = s_mu,
    sigma_wt = exp(s_mu/2),
    s_evolParams0 = s_evolParams0,
    s_evolParams = s_evolParams,
    Td = Td
  ))
}
#' @keywords internal
fit_paramsASV <- function(data, sParams , D=2){
  Td = sParams$Td
  s_p_error_term = sample_jfast(Td,data-sParams$s_mu)
  s_mu = sampleTrend(
    data - s_p_error_term$mean,
    obs_sigma_t2 = s_p_error_term$var,
    evol_sigma_t2 = c(sParams$s_evolParams0$sigma_w0^2,
                      sParams$s_evolParams$sigma_wt^2),
    D = D,
    Td = Td)
  s_omega = diff(s_mu, differences = D)
  s_mu0 = as.matrix(s_mu[1:D,])
  s_evolParams0 = dsp_sampleEvol0(s_mu0, sParams$s_evolParams0)
  s_evolParams = dsp_sampleEvolParams(omega = s_omega,
                                      evolParams = sParams$s_evolParams,
                                      sigma_e = 1,
                                      evol_error = "HS")
  sParams$s_p_error_term = s_p_error_term
  sParams$s_mu = s_mu
  sParams$s_evolParams0 = s_evolParams0
  sParams$s_evolParams = s_evolParams
  sParams$sigma_wt =exp(s_mu/2)
  return(sParams)
  #
  # list(s_p_error_term = s_p_error_term,
  #      s_mu = s_mu,
  #      s_evolParams0 = s_evolParams0,
  #      s_evolParams = s_evolParams)
}
