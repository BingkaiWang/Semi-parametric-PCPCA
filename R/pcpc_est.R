source("R/pcpc_test.R")
#' Consistent estimating the common principle components (CPCs) with or without knowing the number of CPCs
#' in partial common principle component analysis
#'
#' @param cov_data a list sample covariance matrices with dimension p.
#' @param ns the number of samples to estimate each covariance matrix.
#' @param k the number of CPCs (an integer between 1 and p). If not specified, a seqential testing procedure is used to estimate k.
#'
#' @return a p-by-k (or p-by-\hat{k} if k is not given) orthogonal matrix, which is an estimate of covariacne matrices.
#' @export
#'
#' @examples
#' set.seed(1234)
#' p <- 4 # dimension of covariance matrices
#' k <- 2 # number of CPCs
#' n <- 20 # numbder of covariance matrices
#' ns <- 20 # number of samples the estimate each matrix
#' lambda <- exp(4:1) # eigenvalues
#' gamma_true <- svd(matrix(rnorm(p^2), p, p))$u[,1:k] # true common eigenvectors
#' cov <- map(1:n, function(j){
#'   lambda_common <- sample(lambda, 2)
#'   lambda_individual <- setdiff(lambda,  lambda_common)
#'   temp <- gramSchmidt(cbind(gamma_true, matrix(rnorm(p*(p-k)), p, p-k)))$Q[,(k+1):p]
#'   gamma_true %*% diag(lambda_common, nrow = k) %*% t(gamma_true) +
#'   temp %*% diag(lambda_individual, nrow = p-k) %*% t(temp)
#' }) # generate covariance matrices
#' cov_data <- map(cov, ~cov(MASS::mvrnorm(ns, rep(0,p), .))) # generate sample covariance matrices from a wishart distribution
#' gamma_hat <- pcpc_est(cov_data, ns, k)
#' t(gamma_hat) %*% gamma_true
#' gamma_hat <- pcpc_est(cov_data, ns)
#' t(gamma_hat) %*% gamma_true
pcpc_est <- function(cov_data, ns, k = NULL){
  n <- length(cov_data)
  p <- nrow(cov_data[[1]])
  avg_cov <- avg_mat(cov_data)
  pc_candi <- svd(avg_cov)$u
  
  reconstruct_cov_data <- map(cov_data, ~t(pc_candi) %*% . %*% pc_candi)
  eig_avg_cov <- eig(avg_cov)
  cpc_error <- map_dbl(1:p, function(i){
    mean(map_dbl(reconstruct_cov_data, ~sum(.[i,-i]^2/(eig_avg_cov[i]*eig_avg_cov[-i])))) * ns / (p-1)
  })
  pc_candi <- pc_candi[,order(cpc_error)]
  
  if(!is.null(k)){
    return(pc_candi[,1:k])
  } else {
    return(pcpc_test(pc_candi, cov_data, ns)$cpc)
  }
}
