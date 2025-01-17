#' @keywords internal
#' @importFrom Matrix bandSparse
build_Q_trend = function (obs_sigma_t2, evol_sigma_t2, D = 1, Td) {
  # T = length(evol_sigma_t2)
  if (D == 1) {
    Q = Matrix::bandSparse(
      Td,
      k = c(0, 1),
      diagonals = list(
        1 / obs_sigma_t2 +
          1 / evol_sigma_t2 + c(1 /
                                  evol_sigma_t2[-1], 0),
        -1 / evol_sigma_t2[-1]
      ),
      symmetric = TRUE
    )
  } else if (D == 2) {
    Q = Matrix::bandSparse(
      Td,
      k = c(0, 1, 2),
      diagonals = list(
        1 / obs_sigma_t2 +
          1 / evol_sigma_t2 + c(0, 4 / evol_sigma_t2[-(1:2)], 0) + c(1 /evol_sigma_t2[-(1:2)], 0, 0),
        c(-2 / evol_sigma_t2[3], -2 * (1 / evol_sigma_t2[-(1:2)] + c(1 /evol_sigma_t2[-(1:3)], 0))),
        1 / evol_sigma_t2[-(1:2)]
      ),
      symmetric = TRUE
    )
  } else if (D == 3) {
    d0 = c(1 / evol_sigma_t2) + c(0, 0, 9 / evol_sigma_t2[4:Td], 0) + c(0, 9/evol_sigma_t2[4:Td], 0, 0) +
      c(1 / evol_sigma_t2[4:Td], 0, 0, 0) + 1 / obs_sigma_t2
    d1 = -c(3 / evol_sigma_t2[4:Td], 0, 0) - c(0, 9 / evol_sigma_t2[4:Td], 0) -
      c(0, 0, 3 / evol_sigma_t2[4:Td])
    d2 = c(3 / evol_sigma_t2[4:Td], 0) + c(0, 3 / evol_sigma_t2[4:Td])
    d3 = -1 / evol_sigma_t2[4:Td]
    Q = Matrix::bandSparse(
      Td,
      k = c(0, 1, 2, 3),
      diagonals = list(d0, d1, d2, d3),
      symmetric = TRUE
    )
  } else{
    stop("build_Q requires D = 1 or D = 2")
  }
  return(Q)
}
#' @keywords internal
#' @importFrom Matrix t chol  solve
#' @importFrom stats rnorm
sampleTrend <- function(data, obs_sigma_t2, evol_sigma_t2, D = 1, Td) {
  if ((D < 0) || (D != round(D)))
    stop("D must be a positive integer")
  if (any(is.na(data))) {
    stop("y cannot contain NAs")
  }
  linht = data / obs_sigma_t2
  QHt_Matrix = build_Q_trend(obs_sigma_t2 = obs_sigma_t2,
                             evol_sigma_t2 = pmax(evol_sigma_t2,1e-16),
                             D = D,
                             Td)
  chQht_Matrix <- Matrix::chol(QHt_Matrix)
  mu = as.matrix(Matrix::solve(chQht_Matrix, Matrix::solve(Matrix::t(chQht_Matrix),
                                                           linht) + stats::rnorm(Td)))
  return(mu)
}
#' @keywords internal
#' @import extraDistr
initEvolParams_HS <- function(omega, Td) {
  tau = 1
  tau_x = extraDistr::rinvgamma(1, 1, 1 + 1 / tau ^ 2)
  lambda_x = rep(1, Td)
  lambda = sqrt(extraDistr::rinvgamma(Td, 1, 1 / lambda_x + sum(omega ^
                                                                  2) / (2 * tau)))
  return(list(
    sigma_wt = pmax(tau * lambda, 1e-7),
    tau = tau,
    tau_x = tau_x,
    lambda = lambda,
    lambda_x = lambda_x
  ))
}
#' @keywords internal
#' @import extraDistr
sampleEvolParams_HS <- function(omega, evolParams, Td) {
  tau = evolParams$tau
  tau_x = evolParams$tau_x
  lambda = evolParams$lambda
  lambda_x = evolParams$lambda_x

  tau = sqrt(extraDistr::rinvgamma(1, (Td + 1) / 2, 1 / tau_x + sum((omega /
                                                                       lambda) ^ 2) / 2))
  tau_x = extraDistr::rinvgamma(1, 1, 1 + 1 / tau ^ 2)
  lambda = sqrt(extraDistr::rinvgamma(Td, 1, 1 / lambda_x + sum(omega ^
                                                                  2) / (2 * tau)))
  lambda_x = extraDistr::rinvgamma(Td, 1, 1 + 1 / lambda ^ 2)

  evolParams$sigma_wt = pmax(tau * lambda, 1e-7)
  evolParams$tau = tau
  evolParams$tau_x = tau_x
  evolParams$lambda = lambda
  evolParams$lambda_x = lambda_x
  return(evolParams)
}
# Regularized Horseshoe for outliers
#' @keywords internal
#' @import extraDistr
initEvolParams_HS_sparse <- function(omega,
                                     Td,
                                     tau = 1 / (100 * Td)) {
  #nmin1 = 0
  omega_norm = omega / tau
  x_lambda_t = rep(100, Td)
  lambda_2 = extraDistr::rinvgamma(Td, 1, 1 / x_lambda_t + omega_norm ^
                                     2 / 2)
  #sigma_wt = pmax(sqrt(lambda_2)*tau,1e-7)
  sigma_wt = pmax(sqrt(lambda_2) * tau, 1e-6)
  return(list(
    sigma_wt = sigma_wt,
    tau = tau,
    lambda_2 = lambda_2
  ))
}
#' @keywords internal
#' @import extraDistr
sampleEvolParams_HS_sparse <- function(omega,
                                       Td,
                                       evolParams,
                                       tau = 1 / (100 * Td)){
  tau = evolParams$tau
  lambda_2 = evolParams$lambda_2
  omega_norm = omega / tau
  x_lambda_t = extraDistr::rinvgamma(Td, 1, 1 + 1 / lambda_2)
  lambda_2 = extraDistr::rinvgamma(Td, 1, 1 / x_lambda_t + omega_norm ^
                                     2 / 2)
  sigma_wt = pmax(sqrt(lambda_2) * tau, 1e-6)
  #sigma_wt = sqrt(lambda_2)*tau
  return(list(
    sigma_wt = sigma_wt,
    tau = tau,
    lambda_2 = lambda_2
  ))
}
#' @keywords internal
init_Tbeta <- function(data, obserror, evol_error, D, sparse) {
  Td = length(data)
  s_mu = sampleTrend(
    data,
    obs_sigma_t2 = obserror$sigma_et ^ 2,
    evol_sigma_t2 = 0.01 * rep(1, Td),
    D = D,
    Td
  )
  s_omega = diff(s_mu, differences = D)
  s_mu0 = as.matrix(s_mu[1:D, ])
  s_evolParams0 = dsp_initEvol0(s_mu0 / obserror$sigma_et[c(1:D)])
  if (sparse) {
    # s_evolParams = t_initEvolZeta_ps(s_omega/obserror$sigma_et[-c(1:D)],
    #                                         Td-D)
    s_evolParams = initEvolParams_HS_sparse(s_omega / obserror$sigma_et[-c(1:D)], Td -
                                              D)
  } else{
    s_evolParams = dsp_initEvolParams(s_omega / obserror$sigma_et[-c(1:D)], evol_error = "HS")
  }
  n_squared_sum = sum(c(s_mu[c(1:D)], s_omega) ^ 2 /
                        c(s_evolParams0$sigma_w0 ^ 2, s_evolParams$sigma_wt ^
                            2))
  return(
    list(
      s_mu = s_mu,
      s_evolParams0 = s_evolParams0,
      s_evolParams = s_evolParams,
      Td = Td,
      n_squared_sum = n_squared_sum,
      colname = "Trend"
    )
  )
}
#' @keywords internal
fit_Tbeta <- function(data,
                      tParam,
                      obserror,
                      evol_error,
                      D,
                      sparse) {
  obs_sigma_e  = obserror$sigma_e

  s_mu = sampleTrend(
    data,
    obs_sigma_t2 = obserror$sigma_et ^ 2,
    evol_sigma_t2 =(obs_sigma_e *
                      c(tParam$s_evolParams0$sigma_w0, tParam$s_evolParams$sigma_wt)
    ) ^2,
    D = D,
    Td = tParam$Td
  )

  s_omega = diff(s_mu, differences = D)
  idx_d = 1:D
  s_mu0 = as.matrix(s_mu[idx_d, ])
  s_evolParams0 = dsp_sampleEvol0(s_mu0/obs_sigma_e, tParam$s_evolParams0)
  if (sparse) {
    # s_evolParams = t_sampleEvolZeta_ps(s_omega/obserror$sigma_et[-c(1:D)],
    #                                    tParam$Td-D,
    #                                    tParam$s_evolParams)
    s_evolParams = sampleEvolParams_HS_sparse(s_omega/obs_sigma_e,
                                              tParam$Td - D,
                                              tParam$s_evolParams)
  } else{
    s_evolParams = dsp_sampleEvolParams(
      omega = s_omega/obs_sigma_e,
      evolParams = tParam$s_evolParams,
      sigma_e = 1,
      evol_error = evol_error
    )
  }

  n_squared_sum = sum(c(s_mu[idx_d], s_omega) ^ 2 /
                        c(s_evolParams0$sigma_w0 ^ 2, s_evolParams$sigma_wt ^
                            2))
  # storing
  tParam$s_mu = s_mu
  tParam$s_evolParams0 = s_evolParams0
  tParam$s_evolParams = s_evolParams
  tParam$n_squared_sum = n_squared_sum
  return(tParam)
}
