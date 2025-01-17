#' @keywords internal
#' @importFrom Matrix bandSparse
build_Q_season <- function(obs_sigma_t2 = NULL, evol_sigma_t2, Td, k) {
  prec_prior = 1 / evol_sigma_t2
  # initializing indices to reduce runtime
  i0_1 = 1:(k-3)
  i0_2 = 1:(Td-k)
  i0_3 = 2:(k-2)
  diag0 = prec_prior
  diag0[i0_1] = diag0[i0_1] + prec_prior[i0_1 + 2]
  diag0[i0_2] = diag0[i0_2] + prec_prior[i0_2 + k]
  diag0[i0_3] = diag0[i0_3] + 4*prec_prior[i0_3 + 1]

  diag1 = numeric(Td-1)
  diag1[i0_1] = diag1[i0_1] -2*prec_prior[i0_1+2]
  diag1[i0_3] = diag1[i0_3] -2*prec_prior[i0_3+1]

  diag2 = numeric(Td-2)
  diag2[i0_1] = diag2[i0_1] + prec_prior[i0_1 + 2]

  diagk = -prec_prior[-c(1:k)]

  # Add observation variance to main diagonal if provided
  if (!is.null(obs_sigma_t2)) {
    diag0 = diag0 + 1 / obs_sigma_t2
  }

  # Create sparse matrix using the bandSparse function
  Q = Matrix::bandSparse(
    Td,
    k = c(0, 1, 2, k),
    diagonals = list(diag0, diag1, diag2, diagk),
    symmetric = TRUE
  )
  row_indices <- c(rep(1:(k-1), each = (k-1)), 1:(k-1), rep(k, k-1))
  col_indices <- c(rep(1:(k-1), (k-1)), rep(k, k-1), 1:(k-1))
  values <- rep(prec_prior[k], length(row_indices))
  update <- Matrix::sparseMatrix(i = row_indices, j = col_indices, x = values,
                         dims = c(Td, Td))
  Q <- Q + update
  return(Q)
}
#' @keywords internal
#' @importFrom Matrix t chol  solve
#' @importFrom stats rnorm
sampleBeta_season <- function(data, obs_sigma_t2, evol_sigma_t2, Td, k) {
  Q = build_Q_season(obs_sigma_t2, evol_sigma_t2, Td, k)
  linht = data / obs_sigma_t2
  chQht_Matrix <- Matrix::chol(Q)
  mu = as.matrix(Matrix::solve(chQht_Matrix, Matrix::solve(Matrix::t(chQht_Matrix), linht) + stats::rnorm(Td)))
  return(mu)
}
#' @keywords internal
init_Sbeta <- function(data, obserror, evol_error, k) {
  Td = length(data)
  obs_sigma_2T = obserror$sigma_et
  sigma_e = obserror$sigma_e
  s_mu = sampleBeta_season(
    data,
    obs_sigma_t2 = obs_sigma_2T ^ 2,
    evol_sigma_t2 =  sigma_e ^ 2 * rep(0.01, Td),
    Td = Td,
    k = k
  )

  # S 12
  s_evolParams23 = dsp_initEvol0(s_mu[1:2]/sigma_e,
                                 commonSD = TRUE)
  # S(3:km1)
  i1_km1 = c(1:(k-1))
  omega_2 = diff(s_mu[i1_km1], differences = 2)
  s_evolParams3k = dsp_initEvolParams(omega_2/sigma_e,
                                      evol_error = "HS")
  # S k+1,..,T
  omega_kT = matrix(c((s_mu[k] + sum(s_mu[i1_km1])),
                      diff(s_mu, lag = k))
  )
  s_evolParamskT = dsp_initEvolParams(omega_kT/sigma_e,
                                      evol_error = "HS")
  # s_evolParamskT = initEvolParams_HS_sparse(omega_kT/sigma_e,Td-k+1,
  #                                           1/(Td%%k))
  #
  #normalized squared sum
  n_squared_sum = sum((
    c(s_mu[1:2], omega_2, omega_kT) /
      c(
        s_evolParams23$sigma_w0,
        s_evolParams3k$sigma_wt,
        s_evolParamskT$sigma_wt
      )
  ) ^ 2)
  return(
    list(
      s_mu = s_mu,
      s_evolParams23 = s_evolParams23,
      s_evolParams3k = s_evolParams3k,
      s_evolParamskT = s_evolParamskT,
      Td = Td,
      n_squared_sum = n_squared_sum,
      colname = paste0("Seasonal", k)
    )
  )
}
#' @keywords internal
fit_Sbeta <- function(data, sParam, obserror, evol_error, k) {
  Td = sParam$Td
  obs_sigma_2T = obserror$sigma_et
  obs_sigma_e  = obserror$sigma_e
  s_mu = sampleBeta_season(
    data,
    obs_sigma_t2 = obs_sigma_2T ^ 2,
    evol_sigma_t2 = (
      obs_sigma_e *
        c(sParam$s_evolParams23$sigma_w0,
          sParam$s_evolParams3k$sigma_wt,
          sParam$s_evolParamskT$sigma_wt
        )
    ) ^ 2,
    Td = Td,
    k = k
  )
  s_evolParams23 = dsp_sampleEvol0(s_mu[1:2]/obs_sigma_e,
                                   sParam$s_evolParams23,
                                   commonSD = TRUE)

  # S(3:km1)
  i1_km1 = c(1:(k-1))
  omega_2 = diff(s_mu[i1_km1], differences = 2)
  s_evolParams3k = dsp_sampleEvolParams(
    omega_2 / obs_sigma_e,
    evolParams = sParam$s_evolParams3k,
    evol_error = "HS"
  )

  # S k+1,..,T
  omega_kT = matrix(c((s_mu[k] + sum(s_mu[i1_km1])),
                      diff(s_mu, lag = k)))
  s_evolParamskT = dsp_sampleEvolParams(
    omega_kT /obs_sigma_e,
    evolParams = sParam$s_evolParamskT,
    evol_error = "HS"
  )

  #normalized squared sum
  n_squared_sum = sum((
    c(s_mu[1:2], omega_2, omega_kT) /
      c(
        s_evolParams23$sigma_w0,
        s_evolParams3k$sigma_wt,
        s_evolParamskT$sigma_wt
      )
  ) ^ 2)

  sParam$s_mu = s_mu
  sParam$s_evolParams23 = s_evolParams23
  sParam$s_evolParams3k = s_evolParams3k
  sParam$s_evolParamskT = s_evolParamskT
  sParam$n_squared_sum = n_squared_sum
  return(sParam)
}
