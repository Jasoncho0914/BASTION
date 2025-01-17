#' @keywords internal
#' @importFrom stats rnorm
sampleOutlier <- function(data, obs_sigma_t2, evol_sigma_t2, Td) {
  linht = data / obs_sigma_t2
  postSD = 1 / sqrt(1 / obs_sigma_t2 + 1 / evol_sigma_t2)
  postMean = (linht) * postSD ^ 2
  return(stats::rnorm(n = Td, mean = postMean, sd = postSD))
}
#' @keywords internal
#' @import extraDistr
t_initEvolZeta_ps <- function(zeta, Td) {
  zeta = as.matrix(zeta)
  xi = extraDistr::rinvgamma(1, 1 / 2, 1)
  tau_t2 = extraDistr::rinvgamma(1, 1 / 2, 1 / xi)
  v = extraDistr::rinvgamma(Td, 1 / 2, 1 / tau_t2)
  lambda_t2 = extraDistr::rinvgamma(Td, 1 / 2, 1 / v)
  return(list(
    xi = xi,
    v = v,
    tau_t2 = tau_t2,
    lambda_t2 = lambda_t2,
    sigma_wt = sqrt(lambda_t2)
  ))
}
#' @keywords internal
#' @import extraDistr
t_sampleEvolZeta_ps <- function(zeta, Td, evolParams) {
  hsInput2 = zeta ^ 2
  evolParams$lambda_t2 = extraDistr::rinvgamma(Td, 1, 1 / evolParams$v +
                                                 hsInput2 / 2)
  evolParams$v = extraDistr::rinvgamma(Td, 1, 1 / evolParams$lambda_t2 +
                                         1 / evolParams$tau_t2)
  evolParams$tau_t2 = extraDistr::rinvgamma(1, (Td + 1) / 2, 1 / evolParams$xi +
                                              sum(1 / evolParams$v))
  evolParams$xi = extraDistr::rinvgamma(1, 1, 1 + 1 / evolParams$tau_t2)
  evolParams$sigma_wt = pmax(sqrt(evolParams$lambda_t2), 1e-6)
  evolParams
}
#' @keywords internal
init_Outlier <- function(data, obserror) {
  # Td = length(data)
  # s_mu  = sampleOutlier(data,
  #                       obs_sigma_t2 = obserror$sigma_et^2,
  #                       evol_sigma_t2 = 0.01*rep(1,Td),
  #                       Td)
  # s_evolParams = initEvolParams_HS_sparse(s_mu/obserror$sigma_et,Td)
  # n_squared_sum = sum((s_mu/s_evolParams$sigma_wt)^2)
  Td = length(data) - 4
  idx14 = c(1:4)
  data_5T = data[-idx14]
  zeta_5T  = sampleOutlier(
    data_5T,
    obs_sigma_t2 = obserror$sigma_et[-idx14] ^ 2,
    evol_sigma_t2 = 0.01 * rep(1, Td),
    Td
  )
  #s_evolParams = initEvolParams_HS_sparse(zeta_5T/obserror$sigma_et[-idx14],Td)
  s_evolParams = t_initEvolZeta_ps(zeta_5T / obserror$sigma_et[-idx14], Td)
  n_squared_sum = sum((zeta_5T / s_evolParams$sigma_wt) ^ 2)
  s_mu = c(0, 0, 0, 0, zeta_5T)
  list(
    s_mu = s_mu,
    s_evolParams = s_evolParams,
    Td = Td,
    n_squared_sum = n_squared_sum,
    colname = "Outlier"
  )
}
#' @keywords internal
fit_Outlier <- function(data, zParam, obserror) {
  Td = zParam$Td
  idx14 = c(1:4)
  obs_sigma_e = obserror$sigma_e
  zeta_5T  = sampleOutlier(
    data[-idx14],
    obs_sigma_t2 = obserror$sigma_et[-idx14] ^ 2,
    evol_sigma_t2 = obs_sigma_e ^ 2 *
      zParam$s_evolParams$sigma_wt ^ 2,
    Td
  )
  # s_evolParams = sampleEvolParams_HS_sparse(zeta_5T/obserror$sigma_et[-idx14],
  #                                           Td = Td,
  #                                           zParam$s_evolParams)
  s_evolParams = t_sampleEvolZeta_ps(zeta_5T/obs_sigma_e, Td = Td, zParam$s_evolParams)
  n_squared_sum = sum((zeta_5T/s_evolParams$sigma_wt) ^ 2)
  s_mu = c(0, 0, 0, 0, zeta_5T)
  # Td = zParam$Td
  # s_mu  = sampleOutlier(data,
  #                       obs_sigma_t2 = obserror$sigma_et^2,
  #                       evol_sigma_t2 = obserror$sigma_et^2*zParam$s_evolParams$sigma_wt^2,
  #                       Td)
  # s_evolParams = sampleEvolParams_HS_sparse(s_mu/obserror$sigma_et,
  #                                           Td = Td,
  #                                           zParam$s_evolParams)
  # n_squared_sum = sum((s_mu/s_evolParams$sigma_wt)^2)
  zParam$s_mu = s_mu
  zParam$s_evolParams = s_evolParams
  zParam$n_squared_sum = n_squared_sum
  return(zParam)
}
