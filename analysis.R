#######################
# Plotting function
#######################
# max(lasso_result_n$irf_lp_error$abs_mean)
# max(lasso_result_n$irf_lp_error$abs_sd)
plot_mean_irf <- function(lasso_result_n, ridge_result_n, n, pseq){
  # this function plot the mean absolute error of irf
  plot(pseq, rep(0, length(pseq)), ylim = c(0, 7), xlab = "number of series", ylab = "absolute error",
       type = 'n', main = c("sample size = ", n))
  lines(pseq, lasso_result_n$irf1_error$abs_mean, col = 1)
  lines(pseq, lasso_result_n$irf2_error$abs_mean, col = 2)
  lines(pseq, lasso_result_n$irf_lp_error$abs_mean, col = 3)
  lines(pseq, ridge_result_n$irf1_error$abs_mean, col = 1, lty = 2)
  lines(pseq, ridge_result_n$irf2_error$abs_mean, col = 2, lty = 2)
  lines(pseq, ridge_result_n$irf_lp_error$abs_mean, col = 3, lty = 2)
}
plot_sd_irf <- function(lasso_result_n, ridge_result_n, n, pseq){
  # this function plot the mean absolute error of irf
  plot(pseq, rep(0, length(pseq)), ylim = c(0, 1), xlab = "number of series", ylab = "sd",
       type = 'n', main = c("sample size = ", n))
  lines(pseq, lasso_result_n$irf1_error$abs_sd, col = 1)
  lines(pseq, lasso_result_n$irf2_error$abs_sd, col = 2)
  lines(pseq, lasso_result_n$irf_lp_error$abs_sd, col = 3)
  lines(pseq, ridge_result_n$irf1_error$abs_sd, col = 1, lty = 2)
  lines(pseq, ridge_result_n$irf2_error$abs_sd, col = 2, lty = 2)
  lines(pseq, ridge_result_n$irf_lp_error$abs_sd, col = 3, lty = 2)
}
plot_rel_mean_irf <- function(lasso_result_n, ridge_result_n, n, pseq){
  # this function plot the mean absolute error of irf
  plot(pseq, rep(0, length(pseq)), ylim = c(0, 7), xlab = "number of series", ylab = "relative error",
       type = 'n', main = c("sample size = ", n))
  lines(pseq, lasso_result_n$irf1_error$rel_mean, col = 1)
  lines(pseq, lasso_result_n$irf2_error$rel_mean, col = 2)
  lines(pseq, lasso_result_n$irf_lp_error$rel_mean, col = 3)
  lines(pseq, ridge_result_n$irf1_error$rel_mean, col = 1, lty = 2)
  lines(pseq, ridge_result_n$irf2_error$rel_mean, col = 2, lty = 2)
  lines(pseq, ridge_result_n$irf_lp_error$rel_mean, col = 3, lty = 2)
}
plot_rel_sd_irf <- function(lasso_result_n, ridge_result_n, n, pseq){
  # this function plot the mean absolute error of irf
  plot(pseq, rep(0, length(pseq)), ylim = c(0, 1), xlab = "number of series", ylab = "sd",
       type = 'n', main = c("sample size = ", n))
  lines(pseq, lasso_result_n$irf1_error$rel_sd, col = 1)
  lines(pseq, lasso_result_n$irf2_error$rel_sd, col = 2)
  lines(pseq, lasso_result_n$irf_lp_error$rel_sd, col = 3)
  lines(pseq, ridge_result_n$irf1_error$rel_sd, col = 1, lty = 2)
  lines(pseq, ridge_result_n$irf2_error$rel_sd, col = 2, lty = 2)
  lines(pseq, ridge_result_n$irf_lp_error$rel_sd, col = 3, lty = 2)
}
plot_A_res <- function(lasso_result_n, ridge_result_n, n, pseq){
  plot(pseq, rep(0, length(pseq)), ylim = c(0, 5), xlab = "number of series", ylab = "relative error",
       type = 'n', main = c("sample size = ", n))
  lines(pseq, lasso_result_n$A_error, col = 1)
  lines(pseq, lasso_result_n$res_error, col = 2)
  lines(pseq, ridge_result_n$A_error, col = 1, lty = 2)
  lines(pseq, ridge_result_n$res_error, col = 2, lty = 2)
}

op = par()



#####################
#  Sparsity: 0.05   #
#####################

lasso_n40 = readRDS("output/sparsity0_05/output/lasso_n40.rds")
lasso_n60 = readRDS("output/sparsity0_05/output/lasso_n60.rds")
lasso_n80 = readRDS("output/sparsity0_05/output/lasso_n80.rds")

ridge_n40 = readRDS("output/sparsity0_05/output/ridge_n40.rds")
ridge_n60 = readRDS("output/sparsity0_05/output/ridge_n60.rds")
ridge_n80 = readRDS("output/sparsity0_05/output/ridge_n80.rds")

pseq = seq(5, 60, 5)
par(mfrow = c(3, 2))
plot_mean_irf(lasso_n40, ridge_n40, 40, pseq)
plot_sd_irf(lasso_n40, ridge_n40, 40, pseq)
plot_mean_irf(lasso_n60, ridge_n60, 60, pseq)
plot_sd_irf(lasso_n60, ridge_n60, 60, pseq)
plot_mean_irf(lasso_n80, ridge_n80, 80, pseq)
plot_sd_irf(lasso_n80, ridge_n80, 80, pseq)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = T)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("1s-lasso", "1s-ridge", "2s-lasso", "2s-ridge", "lp-lasso", "lp-ridge"), 
       col = c(1, 1, 2, 2, 3, 3), lty = c(1, 2, 1, 2, 1, 2), xpd = T, horiz = T, inset = c(0, 0),
       bty = "n", cex = 0.5)
par(op)

pseq = seq(5, 60, 5)
par(mfrow = c(3, 2))
plot_rel_mean_irf(lasso_n40, ridge_n40, 40, pseq)
plot_rel_sd_irf(lasso_n40, ridge_n40, 40, pseq)
plot_rel_mean_irf(lasso_n60, ridge_n60, 60, pseq)
plot_rel_sd_irf(lasso_n60, ridge_n60, 60, pseq)
plot_rel_mean_irf(lasso_n80, ridge_n80, 80, pseq)
plot_rel_sd_irf(lasso_n80, ridge_n80, 80, pseq)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = T)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("1s-lasso", "1s-ridge", "2s-lasso", "2s-ridge", "lp-lasso", "lp-ridge"), 
       col = c(1, 1, 2, 2, 3, 3), lty = c(1, 2, 1, 2, 1, 2), xpd = T, horiz = T, inset = c(0, 0),
       bty = "n", cex = 0.5)
par(op)

pseq = seq(5, 60, 5)
par(mfrow = c(1, 3))
plot_A_res(lasso_n40, ridge_n40, 40, pseq)
plot_A_res(lasso_n60, ridge_n60, 60, pseq)
plot_A_res(lasso_n80, ridge_n80, 80, pseq)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = T)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("A-lasso", "A-ridge", "res-lasso", "res-ridge"), 
       col = c(1, 1, 2, 2), lty = c(1, 2, 1, 2), xpd = T, horiz = T, inset = c(0, 0),
       bty = "n", cex = 0.5)
par(op)

