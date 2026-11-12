#' @keywords internal
#' @importFrom Matrix chol solve t
#' @importFrom stats rnorm
sample_RC_precision <- function(y, X, beta_sigma_2, sigma_e, Td) {
  A = crossprod(X)
  diag(A) = diag(A) + 1 / beta_sigma_2
  linht = t(X) %*% matrix(y / sigma_e)
  chQht_Matrix <- Matrix::chol(A)
  return(as.matrix(Matrix::solve(
    chQht_Matrix, Matrix::solve(Matrix::t(chQht_Matrix), linht) + stats::rnorm(Td)
  )) * sigma_e)
}
#' @keywords internal
init_Regression <- function(data, X, obserror,sampler=c("auto","precision","sqrt")) {
  sampler<-match.arg(sampler)
  Td = ncol(X)
  beta = sample_RC(data, X, rep(1, Td), obserror$sigma_e, Td,sampler=sampler)
  beta_params = dsp_initEvol0(beta / obserror$sigma_e)
  n_squared_sum = crossprod(beta * 1 / beta_params$sigma_w0)[[1]]

  return(
    list(
      gaussian_sampler = sampler,
      s_mu = X %*% beta,
      beta = beta,
      Td = Td,
      beta_params = beta_params,
      n_squared_sum = n_squared_sum
    )
  )
}
#' @keywords internal
fit_Regression <- function(data, X, bParam, obserror,
    sampler=if(is.null(bParam$gaussian_sampler)) "auto" else bParam$gaussian_sampler) {
  sampler<-match.arg(sampler,c("auto","precision","sqrt"))
  bParam$gaussian_sampler<-sampler
  beta = sample_RC(data,
                   X,
                   bParam$beta_params$sigma_w0 ^ 2,
                   obserror$sigma_e,
                   bParam$Td,sampler=sampler)
  beta_params = dsp_sampleEvol0(beta / obserror$sigma_e, bParam$beta_params, commonSD = FALSE)
  n_squared_sum = crossprod(beta * 1 / beta_params$sigma_w0)[[1]]

  bParam$s_mu = X %*% beta
  bParam$beta = beta
  bParam$beta_params = beta_params
  bParam$n_squared_sum = n_squared_sum
  return(bParam)

}

sample_RC_sqrt <- function(y,X,beta_sigma_2,sigma_e,Td,
                           column_scale=TRUE,qr_order=0L) {
  if(Td!=ncol(X))stop("Td must equal ncol(X)")
  sigma_e<-.gaussian_positive(sigma_e,1L,"sigma_e")
  .gaussian_sqrt_draw(X,as.numeric(y)/sigma_e,Matrix::Diagonal(Td),
    1,beta_sigma_2,column_scale,qr_order)*sigma_e
}
sample_RC <- function(y,X,beta_sigma_2,sigma_e,Td,
                       sampler=c("auto","precision","sqrt"),notify=TRUE,
                       column_scale=TRUE,qr_order=0L) {
  .gaussian_sampler_select(
    function()sample_RC_precision(y,X,beta_sigma_2,sigma_e,Td),
    function()sample_RC_sqrt(y,X,beta_sigma_2,sigma_e,Td,column_scale,qr_order),
    sampler,notify,"Regression coefficients")
}
