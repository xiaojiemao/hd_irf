# Collections of functions

#################
#  Simulations  #
#################
make_A <- function(p, alpha, sparsity, maxIter = 1000){
  # make a VAR coefficient matrix for VAR(1) system with dimension p where 
  # we randomly sample some entries to put nonzero values alpha; it's hard 
  # to generate VAR(d) coefficient matrix corresponding to a stable process
  # alpha is the values for nonzero elements in A: it can be either a constant
  #   or a vector (in later simulations I usually use constants)
  # sparsity is the proportion of nonzero elements
  # the function gives error when the associated VAR process is not stable and
  # it will repeat process of creating A until either A is stable or it's repeated
  # up to 1000 times
  
  A = matrix(0, p, p)
  
  A[sample(length(A), round(length(A)*sparsity), replace =F)] <- alpha
  max(abs(eigen(A)$values))
  iter = 0
  
  while (max(abs(eigen(A)$values)) >= 1 & iter < maxIter){
    A[sample(length(A), round(length(A)*sparsity), replace =F)] <- alpha
    iter = iter + 1
  }
  
  if (iter < maxIter) {
    return(A)
  } else {
    stop("The VAR is unstable")
  }
  
  return(A)
}
simulateVAR <- function(n, A, snr, rho){
  # this function simulates a length-n series from VAR process with VA coeffcient matrix A and toplitz 
  # error covariance matrix
  # snr: signal noise ratio, which is the ratio of the norm of the AR coefficient matrix and the norm
  # of the contemporaneous error covariance matrix E_cov
  # rho: the constant in the toplitz error covariance matrix E_cov
  p = nrow(A)
  d = ncol(A)/p
  E_cov = toeplitz(rho^(0:(p - 1)))
  E_cov = E_cov/(max(svd(E_cov)$d))*(max(svd(A)$d)/snr)  # scale the error covariance matrix such that 
  # the ratio of the norm of the AR coefficient matrix and the norm
  # of the contemporaneous error covariance matrix is equal to the given snr
  sim = VARMAsim(n, arlags = d, phi = A, sigma = E_cov)
  
  return(list(series = sim$series, noises = sim$noises, A = A, E_cov = E_cov, d = d))
}
simulateVARMA <- function(n, A, B, snr, rho){
  # this function simulates a length-n series from VARMA process with VA coeffcient matrix A, MA coefficient matirx B
  # and toplitz error covariance matrix
  # snr: signal noise ratio, which is the ratio of the (||A|| + ||B||) and the norm
  # of the contemporaneous error covariance matrix E_cov
  # rho: the constant in the toplitz error covariance matrix
  
  p = nrow(A)
  d = ncol(A)/p
  q = ncol(B)/p
  
  sigma1 = max(svd(A)$d)
  sigma2 = max(svd(B)$d)
  
  E_cov = toeplitz(rho^(0:(p - 1)))
  E_cov = E_cov/(max(svd(E_cov)$d))*((sigma1 + sigma2)/snr)
  
  sim = VARMAsim(n, arlags = d, malags = q, phi = A, theta = B, sigma = E_cov)
  return(list(series = sim$series, noises = sim$noises, A = A, B = B, E_cov = E_cov, d = d, q = q))
}


###############
#  Estimation #
###############
design_matrix_ebye <- function(y, d = 1) {
  # y: the matrix form of y with rows representing observation for different
  # time and cols representing observation for different series; it can be the
  # series component from the VARMAsim()
  # d: the order of VAR model
  # return: regress the ith column of Y on the whole X gives the ith row of 
  #    [A_1, A_2, ..., A_d]
  n = dim(y)[1]
  p = dim(y)[2]
  
  Y = y[-(1:d), ]
  X = matrix(NA, nrow = n - d, ncol = d*p)
  if (d == 1) {
    X = y[-n, ] 
  } else if (d == 2) {
    X[, 1:p] = y[-c(1, n), ]
    X[, (p + 1):(2 * p)] = y[-c(n - 1, n), ]
  } else if (d > 2) {
    # here the first column blocp and the last column blocp have to be dealt with separately
    # because they only exclude either the head or tail; other column blocps exclude both head
    # and tail
    X[, 1:p] = y[-c(1:(d - 1), n), ]
    for (i in 2:(d - 1)) {
      X[, ((i - 1)*p + 1):(i*p)] = y[-c(1:(d - i), (n - i + 1):n), ]
    }
    X[, ((d - 1)*p + 1):(d*p)] = y[-c((n - d + 1):n), ]
  }
  
  assert_that(are_equal(sum(is.na(X)), 0))
  
  return(list(Y = Y, X = X))
}
estimate_ols_ebye  <- function(X, Y) {
  n = nrow(Y)
  p = ncol(Y)
  d = ncol(X)/p
  
  est_residual = matrix(0, n, p)
  est_coef = matrix(0, p, p*d)
  
  for (i in 1:p) {
    regression = lm(Y[, i] ~ X + 0)
    est_coef[i, ] = regression$coefficients
    est_residual[, i] = residuals(regression)
  }
  if (sum(is.na(est_coef)) > 0){
    message("Multicollinearity!")
  }
  return(list(coefs = est_coef, residual = est_residual))
}
estimate_lasso_ebye <- function(X, Y, nfolds = 5){
  # estimate the coefficient matrix 
  n = nrow(Y)
  p = ncol(Y)
  d = ncol(X)/p
  
  est_residual = matrix(0, n, p)
  est_coef = matrix(0, p, p*d)
  
  for (i in 1:p) {
    ysd = sd(Y[, i])
    if (ysd > 0){
      regression = cv.glmnet(x = X, y = Y[, i], alpha = 1, intercept = F, nfolds = nfolds, 
                             lambda = exp(seq(log(0.001), log(10), length.out=100)))
      est_coef[i, ] =  coef(regression, s = 'lambda.min')[-1]
      est_residual[, i] = Y[, i] - X %*% est_coef[i, ]
    } else {
      # if the response variable is constant across all observations, then glmnet will give error;
      # this can happen if the first stage estimation gives all 0 coefficients to the regression
      # for a certain series then the residuals for the regression would be the same as the response
      # variable values; then because in the design_matrix_ebye_2nd() we enforce that the contemporaneous
      # residuals have identity matrix, the resulting response variable y-e in the second stage would be 
      # all 0
      est_coef[i, ] = rep(0, length(est_coef[i, ]))
      est_residual[, i] = Y[, i]
    }
  }
  
  return(list(coefs = est_coef, residual = est_residual))
}
estimate_ridge_ebye <- function(X, Y, nfolds = 5){
  n = nrow(Y)
  p = ncol(Y)
  d = ncol(X)/p
  
  est_residual = matrix(0, n, p)
  est_coef = matrix(0, p, p*d)
  
  for (i in 1:p) {
    ysd = sd(Y[, i])
    if (ysd > 0){
      regression = cv.glmnet(x = X, y = Y[, i], alpha = 0, intercept = F, nfolds = nfolds,
                             lambda = exp(seq(log(0.001), log(5), length.out=100)))
      est_coef[i, ] =  coef(regression, s = 'lambda.min')[-1]
      est_residual[, i] = Y[, i] - X %*% est_coef[i, ]
    } else {
      est_coef[i, ] = rep(0, length(est_coef[i, ]))
      est_residual[, i] = Y[, i]
    }
  }
  
  return(list(coefs = est_coef, residual = est_residual))
}
VARMAirf <- function(Phi = NULL, Theta = NULL, Sigma = NULL, lag = 6, 
                      orth = FALSE) {
  # this function is modified from the VARMAirf() from the MTS package
  #   by muting the automatic plotting side-effect
  # it returns the irf based on the AR and MA coefficient matrix
  #  VARMA process: (note that both sides use minus by default)
  #   (I - Phi_1L - Ph_2L^2 - ... - Phi_dL^d)y_t = (I - Theta_1L - ... - Theta_dL^d)e_t
  q = 0
  p = 0
  k = 0
  if (length(Theta) > 0) {
    k = dim(Theta)[1]
    k1 = dim(Theta)[2]
    q = floor(k1/k)
  }
  if (length(Phi) > 0) {
    k = dim(Phi)[1]
    k1 = dim(Phi)[2]
    p = floor(k1/k)
  }
  if (is.null(Sigma)) {
    Sigma = diag(rep(1, k))
  }
  if (orth) {
    m1 = eigen(Sigma)
    v1 = sqrt(m1$values)
    vv = diag(v1)
    Pmtx = m1$vectors
    Sh = Pmtx %*% vv %*% t(Pmtx)
  }
  if (k < 1) 
    k = 1
  PSI = diag(rep(1, k))
  if (orth) {
    WGT = c(PSI %*% Sh)
  }
  else {
    WGT = c(PSI)
  }
  for (il in 1:lag) {
    ilk = il * k
    tmp = matrix(0, k, k)
    if ((q > 0) && (il <= q)) 
      tmp = -Theta[, (ilk - k + 1):ilk]
    if (p > 0) {
      iend = min(il, p)
      for (j in 1:iend) {
        jdx = (il - j)
        kdx = j * k
        tmp = tmp + Phi[, (kdx - k + 1):kdx] %*% PSI[, 
                                                     (jdx * k + 1):(jdx * k + k)]
      }
    }
    PSI = cbind(PSI, tmp)
    if (orth) {
      WGT = cbind(WGT, c(tmp %*% Sh))
    }
    else {
      WGT = cbind(WGT, c(tmp))
    }
  }
  wk1 = WGT
  for (i in 1:k^2) {
    wk1[i, ] = cumsum(WGT[i, ])
  }
  tdx = c(1:(lag + 1)) - 1
  # par(mfcol = c(k, k), mai = c(0.3, 0.3, 0.3, 0.3))
  if (orth) {
    gmax = max(WGT)
    gmin = min(WGT)
    cx = (gmax - gmin)/10
    gmax = gmax + cx
    gmin = gmin - cx
    # for (j in 1:k^2) {
    #   plot(tdx, WGT[j, ], type = "l", xlab = "lag", ylab = "IRF", 
    #        ylim = c(gmin, gmax), cex.axis = 0.8)
    #   points(tdx, WGT[j, ], pch = "*", cex = 0.8)
    #   title(main = "Orth. innovations")
    #}
    # cat("Press return to continue ", "\n")
    # readline()
    gmax = max(wk1)
    gmin = min(wk1)
    cx = (gmax - gmin)/10
    gmax = gmax + cx
    gmin = gmin - cx
    # for (j in 1:k^2) {
    #   plot(tdx, wk1[j, ], type = "l", xlab = "lag", ylab = "Acu-IRF", 
    #        ylim = c(gmin, gmax), cex.axis = 0.8)
    #   points(tdx, wk1[j, ], pch = "*", cex = 0.8)
    #   title(main = "Orth. innovations")
    #}
  }
  else {
    gmax = max(WGT)
    gmin = min(WGT)
    cx = (gmax - gmin)/10
    gmax = gmax + cx
    gmin = gmin - cx
    # for (j in 1:k^2) {
    #   plot(tdx, WGT[j, ], type = "l", xlab = "lag", ylab = "IRF", 
    #        ylim = c(gmin, gmax), cex.axis = 0.8)
    #   points(tdx, WGT[j, ], pch = "*", cex = 0.8)
    #   title(main = "Orig. innovations")
    # }
    # cat("Press return to continue ", "\n")
    # readline()
    gmax = max(wk1)
    gmin = min(wk1)
    cx = (gmax - gmin)/10
    gmax = gmax + cx
    gmin = gmin - cx
    # for (j in 1:k^2) {
    #   plot(tdx, wk1[j, ], type = "l", xlab = "lag", ylab = "Acu-IRF", 
    #        ylim = c(gmin, gmax), cex.axis = 0.8)
    #   points(tdx, wk1[j, ], pch = "*", cex = 0.8)
    #   title(main = "Orig. innovations")
    # }
  }
  # par(mfcol = c(1, 1))
  return(PSI[, -c(1:nrow(Phi))])
}
irf_2nd <- function(Yhat, Edesign, estimator){
  # this function uses different methods for estimating the second stage coefficient:
  # method: estimate_ols_ebye, estimate_ridge_ebye, estimate_lasso_ebye
  # Yhat, Edesign:: the response matrix and design matrixin the 2nd stage regression
  #  they should be output from design_matrix_ebye_2nd()
  # return: the estimated impulse response function, which is exactly the coef matrix in the second
  #   stage regression
  
  est = estimator(Edesign, Yhat)
  return(est$coefs)
}
irf_lp <- function(y, d, estimator, h = 6){
  # use the estimator to estimate impulse response function of order 1:h
  # y: the original series
  # d: the lagged order used for local projection; for VAR(d) process this should be d
  # estimator: the estimator used for equation by equation regression
  # h: 1:h impulse response functions are estimated
  p = ncol(y)
  
  irf = matrix(NA, p, p*h)
  
  for (i in 1:h){
    design = design_matrix_lp_ebye(y, i, d)
    estimate <- estimator(design$X, design$Y)
    irf[, ((i - 1)*p + 1):(i*p)] <- estimate$coefs[, 1:p]
  }
  
  assert_that(are_equal(sum(is.na(irf)), 0))
  
  return(irf)
}



###################
#  design matrix #
##################
design_matrix_ebye <- function(y, d = 1) {
  # y: the matrix form of y with rows representing observation for different
  # time and cols representing observation for different series; it can be the
  # series component from the VARMAsim(), e.g., the series component from the 
  # output of simulateVAR()
  # d: the order of VAR model
  # return: response matrix Y and design matrix X such that regress the ith c
  #   olumn of Y on the whole X gives the ith row of 
  #    [A_1, A_2, ..., A_d]
  #   the design matrix is used in the first stage to get the residuals  
  n = dim(y)[1]
  p = dim(y)[2]
  
  Y = y[-(1:d), ]
  X = matrix(NA, nrow = n - d, ncol = d*p)
  if (d == 1) {
    X = y[-n, ] 
  } else if (d == 2) {
    X[, 1:p] = y[-c(1, n), ]
    X[, (p + 1):(2 * p)] = y[-c(n - 1, n), ]
  } else if (d > 2) {
    # here the first column blocp and the last column blocp have to be dealt with separately
    # because they only exclude either the head or tail; other column blocps exclude both head
    # and tail
    X[, 1:p] = y[-c(1:(d - 1), n), ]
    for (i in 2:(d - 1)) {
      X[, ((i - 1)*p + 1):(i*p)] = y[-c(1:(d - i), (n - i + 1):n), ]
    }
    X[, ((d - 1)*p + 1):(d*p)] = y[-c((n - d + 1):n), ]
  }
  
  assert_that(are_equal(sum(is.na(X)), 0))
  
  return(list(Y = Y, X = X))
}
design_matrix_ebye_2nd <- function(Y, E, lag = 6){
  # this function makes the deisgn matrix for the second stage regression
  #   y_t = e_{t} + B_1e_{t-1} + ... + B_qe_{t-q}
  # Y(n-d * p): this Y should come from the response matrix in the first stage, which is returned by
  #   design_matrix_ebye()
  # residuals E: it should come from the residuals estimated by different equation-by-equation methods
  #   in the first stage
  # lag: the maximum order of impulse response function that we compute, which is q in the above regression
  # return:
  # y_hat: this is essentially (Y - E)[{q+1}:n, ]; it enforces the constraint that the contemporaneous 
  #   residual have identity coefficient matrix; note that compared to the original simulated series
  #   we drop the first d terms in the first stage regression and then we further drop the first q terms for
  #   the second stage regression
  # E_design: this is the design matrix for the second stage regression
  
  n = nrow(Y)
  p = ncol(Y)
  Yhat = Y[-c(1:lag), ] - E[-c(1:lag), ]
  Edesign = matrix(NA, n - lag, p*lag)
  
  if (lag == 1) {
    Edesign = E[-n, ] 
  } else if (lag == 2) {
    Edesign[, 1:p] = E[-c(1, n), ]
    Edesign[, (p + 1):(2 * p)] = E[-c(n - 1, n), ]
  } else if (lag > 2) {
    Edesign[, 1:p] = E[-c(1:(lag - 1), n), ]
    for (i in 2:(lag - 1)) {
      Edesign[, ((i - 1)*p + 1):(i*p)] = E[-c(1:(lag - i), (n - i + 1):n), ]
    }
    Edesign[, ((lag - 1)*p + 1):(lag*p)] = E[-c((n - lag + 1):n), ]
  }
  
  assert_that(are_equal(sum(is.na(Edesign)), 0))
  
  return(list(Yhat = Yhat, Edesign = Edesign))
}
design_matrix_lp_ebye <- function(y, h, d = 1){
  # construct design matrix for local projection regression 
  # y_{t+h} on y_t, y_{t-1}, ..., y_{t-d} where h is given by lag and d is given by d
  # y: the input multivariate series
  # h: the order for impulse response function to be estimated
  # d: the order of lag term included in the local projection
  
  n = dim(y)[1]
  p = dim(y)[2]
  
  Y = y[-(1:(d + h - 1)), ]
  X = matrix(NA, nrow = n - d - h + 1, ncol = d*p)
  
  if (d == 1) {
    X = y[-((n - h + 1):n), ] 
  } else if (d >= 2) {
    # here the first column blocp and the last column blocp have to be dealt with separately
    # because they only exclude either the head or tail; other column blocps exclude both head
    # and tail
    for (i in 1:(d - 1)) {
      X[, ((i - 1)*p + 1):(i*p)] = y[-c(1:(d - i), (n - h - i + 2):n), ]
    }
    X[, ((d - 1)*p + 1):(d*p)] = y[-c((n - h - d + 2):n), ]
  }
  
  assert_that(are_equal(sum(is.na(X)), 0))
  
  return(list(Y = Y, X = X))
}


#################
#  Analysis 
#################
### two stage vs local projection under var process 
intermediate <- function(A, n, d,  estimator, lag = 6, rep = 5, snr = 3, rho = 0.8){
  # main: simulate, estimate and store the result for given n, p(coefficient matrices)
  
  # INPUT:
  # A: the true VAR transition matrix which is used for generating the series and computing the
  #    true IRF
  # n: sample size 
  # p: dimension of time series
  # d: the order of VAR model used in the first stage to get the residuals and the order used in local projection
  # estimator
  # lag: the order up to which the irf is computed
  # rep: the number of times the whole experiment is conducted
  # snr & rho: the parameters for simulating the var systems
  
  # OUTPUT:
  # the list containing all estimated residuals, transition matrix and irf
  #   each list has rep many components with each component storing the result for one repetition
  # est_A: estimated coefficient matrix for each repetition
  # est_res: the estimated residual for each repetition
  # est_irf1: the estimated irf based on stage 1 estimated VAR transition matrix
  # est_irf2: the estimated irf based on stage 2 regression w.r.t. the residuals estimated from
  #   stage 1
  # true_res: the residual in the simulation process; note that it has length n (est_res has length
  #   n - d instead)
  
  est_A = vector('list', rep)
  est_res = vector('list', rep)
  est_irf1 = vector('list', rep)
  est_irf2 = vector('list', rep)
  est_irf_lp = vector('list', rep)  
  true_res = vector('list', rep)
  p = nrow(A)
  
  for (i in 1:rep){
    cat("repetition", i, "\n")
    
    simulation <- simulateVAR(n, A, snr, rho)
    true_res[[i]] = simulation$noises
    # compute first-stage irf estimator 
    design1 <- design_matrix_ebye(simulation$series, d)
    estimate <- estimator(design1$X, design1$Y)
    est_A[[i]] <- estimate$coefs
    est_res[[i]] <- estimate$residual
    est_irf1[[i]] <- VARMAirf(Phi = estimate$coefs, lag = lag)
    # compute the second-stage irf estimator
    design2 <- design_matrix_ebye_2nd(Y = design1$Y, E = estimate$residual, lag = lag)
    est_irf2[[i]] <- irf_2nd(Yhat = design2$Yhat, Edesign = design2$Edesign, estimator = estimator)
    # compute the local projection irf estimator 
    est_irf_lp[[i]] <- irf_lp(simulation$series, d, estimator, h = lag)
  }
  
  return(list(est_A = est_A, est_res = est_res, true_res = true_res, est_irf1 = est_irf1,
              est_irf2 = est_irf2, est_irf_lp = est_irf_lp))
}
main <- function(nseq, pseq, d1, lag, estimator, alpha = 0.3, sparsity = 0.2, rep = 5, snr = 3, rho = 0.8, d = 1){
  # This function return the error measures for different n and p
  
  # nseq, pseq: the grid of sample size and time series dimension 
  # d1: order for fitting the first stage regression
  # lag: the order up to which the irf is computed
  # estimator: the estimator used in both the first and second stage; lasso, ridge or ols
  # d: order of the true VAR process; note that now we can only simulate VAR(1) because we only implement that in make_A
  # alpha, sparsity: parameters in generating A matrix
  # snr, rho: parameters in generating the time series
  # rep: number of repetitions of the experiments
  
  A_error = matrix(NA, length(nseq), length(pseq))
  res_error = matrix(NA, length(nseq), length(pseq))
  irf1_error = init_irf_error(nseq, pseq)
  irf2_error = init_irf_error(nseq, pseq)
  irf_lp_error = init_irf_error(nseq, pseq)
  
  for (j in 1:length(pseq)){
    cat("p = ", pseq[j], "\n")
    
    A = make_A(p = pseq[j], alpha = alpha, sparsity = sparsity)
    
    for (i in 1:length(nseq)){
      cat("n = ", nseq[i], "\n")  
      
      est_result <- intermediate(A = A, n = nseq[i], d = d1, lag = lag, estimator = estimator, rep = rep, snr = snr, rho = rho) 
      A_error[i, j] <- compute_A_error(true_A = A, est_A = est_result$est_A)
      res_error[i, j] <- compute_res_error(est_res = est_result$est_res, true_res = est_result$true_res)
      # store result for irf1
      result = compute_irf_total_error(A = A, est_irf = est_result$est_irf1, lag = lag)
      irf1_error$rel_mean[i, j] <- result$mean_relative_error
      irf1_error$rel_sd[i , j] <- result$sd_relative_error
      irf1_error$abs_mean[i , j] <- result$mean_abs_error
      irf1_error$abs_sd[i, j] <- result$sd_abs_error
      # store result for irf2
      result = compute_irf_total_error(A = A, est_irf = est_result$est_irf2, lag = lag)
      irf2_error$rel_mean[i, j] <- result$mean_relative_error
      irf2_error$rel_sd[i , j] <- result$sd_relative_error
      irf2_error$abs_mean[i , j] <- result$mean_abs_error
      irf2_error$abs_sd[i, j] <- result$sd_abs_error
      # store result for irf_lp
      result = compute_irf_total_error(A = A, est_irf = est_result$est_irf_lp, lag = lag)
      irf_lp_error$rel_mean[i, j] <- result$mean_relative_error
      irf_lp_error$rel_sd[i , j] <- result$sd_relative_error
      irf_lp_error$abs_mean[i , j] <- result$mean_abs_error
      irf_lp_error$abs_sd[i, j] <- result$sd_abs_error
    }
  }
  
  return(list(A_error = A_error, res_error = res_error, irf1_error = irf1_error, 
              irf2_error = irf2_error, irf_lp_error = irf_lp_error))
}


##############
# Helper 
##############
init_irf_error <- function(nseq, pseq){
  # this function initialize the list associated with each type irf estimator that records the mean & sd
  return(list(rel_mean = matrix(NA, length(nseq), length(pseq)),
              rel_sd = matrix(NA, length(nseq), length(pseq)),
              abs_mean = matrix(NA, length(nseq), length(pseq)),
              abs_sd = matrix(NA, length(nseq), length(pseq))))
}
compute_rel_error <- function(estA, trueA) {
  # return the relative error between two matrices 
  #    ||estA - trueA||_F/||trueA||_F
  norm(estA - trueA, "F")/norm(trueA, "F")
}

compute_A_error <- function(true_A, est_A){
  # compute average relative errors for the estimated coefficient matrix
  # est_A: the list of estimated coefficient matrix for each repetition of the experiment 
  
  if (sum(map_dbl(est_A, ~sum(is.na(.)))) > 0){
    message("NA entries in the estimated coefficient matrix")
  }
  mean(map_dbl(est_A, function(x) compute_rel_error(x, true_A)))
  
} 
compute_res_error <- function(est_res, true_res){
  # compute the average error across repetitions where the error is 
  #    ||est_res - true_res||_F/||true_res||_F
  if (sum(map_dbl(est_res, ~sum(is.na(.)))) > 0){
    message("NA entries in the estimated residuals")
  }
  d1 = nrow(true_res[[1]]) - nrow(est_res[[1]])
  err_vector = map2_dbl(est_res, true_res, 
                        function(x, y) compute_rel_error(x, y[-c(1:d1), ]))
  mean(err_vector)
}

compute_irf_total_error <- function(A, est_irf, lag, B = NULL){
  # compute the mean and the standard deviation of the relative/asbsolute error of irf across different repetitions
  if (sum(map_dbl(est_irf, ~sum(is.na(.)))) > 0){
    message("NA entries in the estimated IRF")
  }
  
  p = nrow(est_irf[[1]])
  lag = ncol(est_irf[[1]])/p
  true_irf = VARMAirf(Phi = A, Theta = B, lag = lag)
  
  mean_relative_error = mean(map_dbl(est_irf, function(x) compute_rel_error(x, true_irf)))
  sd_relative_error = sd(map_dbl(est_irf, function(x) compute_rel_error(x, true_irf)))
  mean_abs_error = mean(map_dbl(est_irf, function(x) norm(x - true_irf, "F")))
  sd_abs_error = sd(map_dbl(est_irf, function(x) norm(x - true_irf, "F")))
  
  return(list(mean_relative_error =  mean_relative_error, sd_relative_error = sd_relative_error,
              mean_abs_error = mean_abs_error, sd_abs_error = sd_abs_error))
}


##################
## Plotting 
##################
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

