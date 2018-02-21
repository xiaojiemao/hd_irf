source("../fun_collection.R")
library(MTS)
library(glmnet)
library(assertthat)
library(purrr)
library(expm)
library(testthat)

n = 20
p = 6
d = 1
alpha = 0.3
sparsity = 0.2
snr = 0.8
rho = 0.5
A = make_A(p, alpha, sparsity)
sim_series = simulateVAR(n, A, snr, rho)

test_that("Test VAR Simulations", {
  expect_equal(dim(sim_series$series), c(n, p))
  expect_equal(dim(sim_series$noises), c(n, p))
  expect_equal(sim_series$A, A)
  expect_equal(sim_series$d, 1)
})

test_that("Design Matrix for regression equation by equation", {
  # n = 20
  # p = 6
  # alpha = 0.3
  # sparsity = 0.2
  # snr = 0.8
  # rho = 0.5
  # A = make_A(p, alpha, sparsity)
  # sim_series = simulateVAR(n, A, snr, rho)
  
  design = design_matrix_ebye(sim_series$series)
  expect_equal(design$Y, sim_series$series[-1, ])
  expect_equal(design$X, sim_series$series[-n, ])
})

test_that("Basic estimation suites", {
  ## lasso 
  design = design_matrix_ebye(sim_series$series, d = 1)
  est = estimate_lasso_ebye(design$X, design$Y, nfolds = 3)
  expect_equal(dim(est$coefs), c(6, 6))
  expect_equal(dim(est$residual), c(n - 1, 6))
  
  design = design_matrix_ebye(sim_series$series, d = 2)
  est = estimate_lasso_ebye(design$X, design$Y, nfolds = 3)
  expect_equal(dim(est$coefs), c(6, 12))
  expect_equal(dim(est$residual), c(n - 2, 6))
  
  ## ols
  design = design_matrix_ebye(sim_series$series, d = 1)
  est = estimate_ols_ebye(design$X, design$Y)
  expect_equal(dim(est$coefs), c(6, 6))
  expect_equal(dim(est$residual), c(n - 1, 6))
  
  design = design_matrix_ebye(sim_series$series, d = 2)
  est = estimate_ols_ebye(design$X, design$Y)
  expect_equal(dim(est$coefs), c(6, 12))
  expect_equal(dim(est$residual), c(n - 2, 6))
  
  ## ridge 
  design = design_matrix_ebye(sim_series$series, d = 1)
  est = estimate_ridge_ebye(design$X, design$Y)
  expect_equal(dim(est$coefs), c(6, 6))
  expect_equal(dim(est$residual), c(n - 1, 6))
  
  design = design_matrix_ebye(sim_series$series, d = 2)
  est = estimate_ridge_ebye(design$X, design$Y)
  expect_equal(dim(est$coefs), c(6, 12))
  expect_equal(dim(est$residual), c(n - 2, 6))
  
})


test_that("Test the design matrix for eye-by-eye regression", {
  # n = 20
  # p = 6
  # alpha = 0.3
  # sparsity = 0.2
  # snr = 0.8
  # rho = 0.5
  # A = make_A(p, alpha, sparsity)
  # sim_series = simulateVAR(n, A, snr, rho)
  
  ### First Stage
  design = design_matrix_ebye(sim_series$series, d = 2)
  expect_equal(design$Y, sim_series$series[c(-1, -2), ])
  expect_equal(design$X, cbind(sim_series$series[c(-1, -n), ], sim_series$series[c(-(n-1), -n), ]))
  
  design = design_matrix_ebye(sim_series$series, d = 4)
  expect_equal(design$Y, sim_series$series[-c(1:4), ])
  expect_equal(design$X, cbind(sim_series$series[-c(1:3, n), ], sim_series$series[-c(1:2, n-1, n), ], 
                               sim_series$series[-c(1, n-2, n-1, n), ], sim_series$series[-c(n-3, n-2, n-1, n), ]))
  ### Second Stage
  design = design_matrix_ebye(sim_series$series, d = 1)
  est = estimate_lasso_ebye(design$X, design$Y, nfolds = 3)
  design2 = design_matrix_ebye_2nd(design$Y, est$residual, lag = 3)
  expect_equal(design2$Yhat, design$Y[-(1:3), ] - est$residual[-c(1:3), ])
  expect_equal(design2$Edesign, cbind(est$residual[-c(1:2, n-1), ], est$residual[-c(1, n-2, n-1), ], 
                                      est$residual[-c(n-3, n-2, n-1), ]))
  
  design = design_matrix_ebye(sim_series$series, d = 2)
  est = estimate_lasso_ebye(design$X, design$Y, nfolds = 3)
  design2 = design_matrix_ebye_2nd(design$Y, est$residual, lag = 2)
  expect_equal(design2$Yhat, design$Y[-(1:2), ] - est$residual[-c(1:2), ])
  expect_equal(design2$Edesign, cbind(est$residual[-c(1, n-2), ], est$residual[-c(n-3, n-2), ]))
})

test_that("IRF estimation", {
  design = design_matrix_ebye(sim_series$series, d = 1)
  est = estimate_lasso_ebye(design$X, design$Y, nfolds = 3)
  design2 = design_matrix_ebye_2nd(design$Y, est$residual, lag = 3)
  irf = irf_2nd(design2$Yhat, design2$Edesign, estimate_lasso_ebye)
  expect_equal(dim(irf), c(p, p*3))
})

test_that("VARMAirf", {
  irf = VARMAirf(Phi = A, lag = 3)
  expect_equal(dim(irf), c(p, p*3))
  expect_equal(irf[, 1:p], A)
  expect_equal(irf[, (p+1):(2*p)], A %^% 2)
  expect_equal(irf[, (2*p+1):(3*p)], A %^% 3)

  
  B = matrix(rnorm(p*p), p, p)
  irf = VARMAirf(Phi = cbind(A, B), lag = 3)
  expect_equal(dim(irf), c(p, p*3))
  expect_equal(irf[, 1:p], A)
  expect_equal(irf[, (p+1):(2*p)], irf[, 1:p] %*% A + B)
  expect_equal(irf[, (2*p+1):(3*p)], irf[, (p+1):(2*p)] %*% A + irf[, 1:p] %*% B)
  
  irf = VARMAirf(Phi = A, Theta = B, lag = 3)
  expect_equal(dim(irf), c(p, p*3))
  expect_equal(irf[, 1:p],  A - B)
  expect_equal(irf[, (p+1):(2*p)], A %*% irf[, 1:p])
  expect_equal(irf[, (2*p+1):(3*p)], A %*% irf[, (p+1):(2*p)])
})

test_that("Test the design matrix construction for local projection", {
  Y <- sim_series$series
  h = 1; d = 1
  design1 = design_matrix_lp_ebye(Y, h = h, d = d)
  expect_equal(design1$Y, Y[-c(1:(h + d - 1)), ])
  expect_equal(design1$X, Y[-n, ])
  h = 1; d = 2
  design2 = design_matrix_lp_ebye(Y, h = h, d = d)
  expect_equal(design2$Y, Y[-c(1:(h + d - 1)), ])
  expect_equal(design2$X, cbind(Y[-c(1, n),], Y[-c(n-1, n), ]))
  h = 2; d = 2
  design3 = design_matrix_lp_ebye(Y, h = h, d = d)
  expect_equal(design3$Y, Y[-c(1:(h + d - 1)), ])
  expect_equal(design3$X, cbind(Y[-c(1, n-1, n),], Y[-c((n - 2):n), ]))
  h = 2; d = 3
  design4 = design_matrix_lp_ebye(Y, h = h, d = d)
  expect_equal(design4$Y, Y[-c(1:(h + d - 1)), ])
  expect_equal(design4$X, cbind(Y[-c(1, 2, n-1, n),], Y[-c(1, (n - 2):n), ], Y[-c((n - 3):n), ]))
})

test_that("test local projection regression", {

  Y <- sim_series$series
  
  seed = Sys.time()
  set.seed(seed)
  irf = irf_lp(Y, d = 1, estimator = estimate_ridge_ebye, h = 1)
  expect_equal(dim(irf), c(p, p))
  
  design = design_matrix_lp_ebye(Y, h = 1, d = 1)
  set.seed(seed)
  temp = estimate_ridge_ebye(design$X, design$Y)
  expect_equal(irf, temp$coefs[, 1:p])
  
  seed = Sys.time()
  set.seed(seed)
  irf = irf_lp(Y, d = 1, estimator = estimate_ridge_ebye, h = 2)
  expect_equal(dim(irf), c(p, p*2))
  design1 = design_matrix_lp_ebye(Y, h = 1, d = 1)
  design2 = design_matrix_lp_ebye(Y, h = 2, d = 1)
  set.seed(seed)
  temp1 = estimate_ridge_ebye(design1$X, design1$Y)
  temp2 = estimate_ridge_ebye(design2$X, design2$Y)
  expect_equal(irf, cbind(temp1$coefs[, 1:p], temp2$coefs[, 1:p]))
  
})

test_that("intermediate and helper function", {
  estimator = estimate_lasso_ebye
  int_res = intermediate(A, n, d, estimator, lag = 3, rep = 2)
  
  est_A = list(matrix(rnorm(p*p), p, p), matrix(rnorm(p*p), p, p))
  true_A = matrix(rnorm(p*p), p, p)
  A_error = compute_A_error(true_A, est_A)
  err1 = sqrt(sum((est_A[[1]] - true_A)^2)/sum(true_A^2))
  err2 = sqrt(sum((est_A[[2]] - true_A)^2)/sum(true_A^2))
  expect_equal(A_error, mean(c(err1, err2)))
  
  res_error = compute_res_error(int_res$est_res, int_res$true_res)
  err1 = compute_rel_error(int_res$est_res[[1]], int_res$true_res[[1]][-1, ])
  err2 = compute_rel_error(int_res$est_res[[2]], int_res$true_res[[2]][-1, ])
  expect_equal(res_error, mean(c(err1, err2)))
  
  irf1_error = compute_irf_total_error(A = A, est_irf = int_res$est_irf1, lag = 3)
  true_irf = VARMAirf(Phi = A, lag = 3)
  expect_equal(irf1_error$mean_relative_error, mean(c(compute_rel_error(int_res$est_irf1[[1]], true_irf), 
                                                     compute_rel_error(int_res$est_irf1[[2]], true_irf))))
  expect_equal(irf1_error$sd_relative_error,  sd(c(compute_rel_error(int_res$est_irf1[[1]], true_irf), 
                                                    compute_rel_error(int_res$est_irf1[[2]], true_irf))))
})

test_that("Main function", {
  nseq = 50
  pseq = 40
  d1 = 1
  lag = 3
  estimator = estimate_lasso_ebye
  res = main(nseq, pseq, d1, lag, estimator, alpha = 0.15, sparsity = 0.1, rep = 3, snr = 3, rho = 0.8, d = 1)
})
