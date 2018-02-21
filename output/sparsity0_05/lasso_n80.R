source("../fun_collection.R")
library(MTS)
library(glmnet)
library(assertthat)
library(purrr)
library(expm)
library(testthat)



nseq = 80
pseq = seq(5, 60, 5)
d1 = 1
lag = 10
rep_num = 30
estimator = estimate_lasso_ebye
res = main(nseq, pseq, d1, lag, estimator, alpha = 0.12, sparsity = 0.05, rep = rep_num, snr = 3, rho = 0.8, d = 1)
saveRDS(res, "output/lasso_n80.rds")