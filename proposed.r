library(nnet)

## We first define some auxiliary functions "DRM.est", "cal.KS.stat", "DRM.gof", and then proceed to the main function "DRM.multiclass.Youden".

# This function estimates the DRM parameters using three-class multinomial logistic regression.

DRM.est = function(obs, qt = 1) {
  
  # Preliminary validation
  m = length(obs) - 1
  if (m != 2) { stop("The number of classes must equal 3.") }
  if (!(qt %in% c(1, 2 ,3 ,4))) { stop("The basis function is not defined.") }
  
  n0 = length(unlist(obs[1]))
  n1 = length(unlist(obs[2])) 
  n2 = length(unlist(obs[3])) 
  n = n0 + n1 + n2
  nn = c(n0, n1, n2) # sample sizes
  d = length(qt)
  group = c(rep(1, n0), rep(2, n1), rep(3, n2))  # labels
  pool = unlist(obs) # pooled sample
  
  # Perform multinomial logistic regression with selected basis functions
  if (qt == 1){ out = summary(multinom(factor(group) ~ pool, trace = F))
  Qt_matrix = cbind(1, pool) 
  }
  if (qt == 2){ out = summary(multinom(factor(group) ~ log(pool), trace = F))
  Qt_matrix = cbind(1, log(pool)) 
  }
  if (qt == 3){ out = summary(multinom(factor(group) ~ pool + I(pool^2), trace = F))
  Qt_matrix = cbind(1, pool, pool**2) 
  }
  if (qt == 4){ out = summary(multinom(factor(group) ~ log(pool)+I((log(pool))^2), trace = F))
  Qt_matrix = cbind(1, log(pool), log(pool)**2) 
  }
  
  coef = matrix(coef(out), nrow = m)
  coef[1, 1] = coef[1, 1] - log(n1 / n0)
  coef[2, 1] = coef[2, 1] - log(n2 / n0)
  rownames(coef) = c('θ1', 'θ2')
  
  # obtain estimated distributions based on pooled sample
  P0 = 1 / (n0 + n1 * exp(Qt_matrix %*% coef[1, ]) + n2 * exp(Qt_matrix %*% coef[2, ]))
  P1 = exp(Qt_matrix %*% coef[1, ]) * P0
  P2 = exp(Qt_matrix %*% coef[2, ]) * P0
  prob = cbind(P0, P1, P2)
  colnames(prob) = c('P0', 'P1', 'P2')
  
  list(prob = prob, coef = coef, pool = pool, nn = nn, m = m)
}


# This function calculates the K-S test statistic for assessing the goodness-of-fit of the DRM.

cal.KS.stat = function(obs, qt) {
  
  # Point estimation under DRM
  out = DRM.est(obs, qt)
  prob = out$prob
  coef = out$coef
  nn = out$nn
  pool = out$pool
  n = sum(nn) # pooled sample size
  M = length(obs) # number of classes
  
  # Construct function strings for CDFs
  F_hat = lapply(1:M, function(i) {
    sprintf("function(x){sum(prob[,%d] * ifelse(pool <= x, 1, 0))}", i)
  })
  F_emp = lapply(1:M, function(i) {
    sprintf("function(x){sum(ifelse(unlist(obs[%d]) <= x, 1, 0) / nn[%d])}", i, i)
  })
  
  # Evaluate CDFs at sorted pooled sample points
  pool.sort = sort(pool)
  FF = lapply(1:M, function(i) {
    sapply(pool.sort, eval(parse(text = F_hat[i])))
  })
  FE = lapply(1:M, function(i) {
    sapply(pool.sort, eval(parse(text = F_emp[i])))
  })
  
  # K-S statistics for each class
  Deltas = sapply(1:M, function(i) {
    sqrt(n) * max(abs(unlist(FF[i]) - unlist(FE[i])))
  })
  
  # The overall K-S statistic
  rhos = nn/nn[1] 
  Delta = mean(rhos * Deltas)
  
  return(Delta)
}


# This function obtains the bootstrap p-value of the K-S goodness-of-fit test.

DRM.gof = function(obs, qt, B = 500) {
  
  # Point estimation under DRM
  out = DRM.est(obs, qt)
  prob = out$prob
  coef = out$coef
  nn = out$nn
  pool = out$pool
  M = length(obs) # number of classes
  
  # Observed test statistic
  stat0 = cal.KS.stat(obs, qt)
  
  # Bootstrap resampling
  boot_stat = c()
  for (b in 1:B) {
    boot_obs = lapply(1:M, function(i) {
      sample(pool, nn[i], replace = TRUE, prob = prob[, i])
    })
    boot_stat = c(boot_stat, cal.KS.stat(boot_obs, qt))
  }
  
  # Bootstrap p-value 
  pval = mean(na.omit(boot_stat) > stat0)
  
  list(stat = stat0, pval = pval)
}


#############################
# # This is the main function.

DRM.multiclass.Youden = function(obs, qt = 1, level = 0.95, warn = TRUE) {

    # =========================================================================
    # Step 1: Point estimation of the generalized Youden index and cutoff points
    # =========================================================================

    # Obtain DRM parameter estimates
    out = DRM.est(obs, qt = qt)
    if (is.null(out)) { return(NULL) }
    prob = out$prob
    P0 = prob[, 1]  # (n, 1)
    P1 = prob[, 2]
    P2 = prob[, 3]
    coef = unname(out$coef)
    theta1 = coef[1, ]  # (1, d+1)
    theta2 = coef[2, ]
    pool = out$pool
    m = out$m

    # The proposed estimates of c1, c2 with analytical solutions for selected q(x)
    if (qt %in% c(1, 2)) {
        alpha1 = coef[1, 1]
        beta1 = coef[1, 2]
        alpha2 = coef[2, 1]
        beta2 = coef[2, 2]
        if (qt == 1) {
            c1 = -alpha1 / beta1
            c2 = -(alpha1 - alpha2) / (beta1 - beta2)
        } else if (qt == 2) {
            c1 = exp(-alpha1 / beta1)
            c2 = exp(-(alpha1 - alpha2) / (beta1 - beta2))
        }
    } else if (qt %in% c(3, 4)) {
        alpha1 = coef[1, 1]
        beta11 = coef[1, 2]
        beta12 = coef[1, 3]
        alpha2 = coef[2, 1]
        beta21 = coef[2, 2]       
        beta22 = coef[2, 3]
        if (qt == 3) {
            c1 = (beta11 - sqrt(beta11^2 - 4 * beta12 * alpha1)) / (-2 * beta12)
            c2 = (-(beta11 - beta21) - sqrt((beta11 - beta21)^2 - 4 * (beta12 - beta22) * (alpha1 - alpha2))) / (2 * (beta12 - beta22))       
        } else if (qt == 4) {
            c1 = exp((beta11 - sqrt(beta11^2 - 4 * beta12 * alpha1)) / (-2 * beta12))
            c2 = exp((-(beta11 - beta21) - sqrt((beta11 - beta21)^2 - 4 * (beta12 - beta22) * (alpha1 - alpha2))) / (2 * (beta12 - beta22)))       
        }
    }

    # Validate the ordering of cutoff points
    if (c1 >= c2) {
        warning('cutoff points are not ordered')
        return(NULL)
    }

    # Calculate the proposed estimate of generalized Youden index 
    F0_hat = function(x) { sum(P0 * ifelse(pool <= x, 1, 0)) }
    F1_hat = function(x) { sum(P1 * ifelse(pool <= x, 1, 0)) }
    F2_hat = function(x) { sum(P2 * ifelse(pool <= x, 1, 0)) }
    J = 0.5 * (F0_hat(c1) - F1_hat(c1) + F1_hat(c2) - F2_hat(c2))
    J = ifelse(is.na(J), 0, J)

    # =========================================================================
    # Step 2: Compute Sigma_Theta (covariance matrix for theta)
    # =========================================================================

    # Sample sizes
    n0 = length(unlist(obs[1]))
    n1 = length(unlist(obs[2]))
    n2 = length(unlist(obs[3]))
    n = n0 + n1 + n2
    rho0 = n0 / n
    rho1 = n1 / n
    rho2 = n2 / n
    rhos = c(rho0, rho1, rho2)
    
    # Compute basis functions based on pooled sample
    if (qt == 1){ 
        qt_matrix = cbind(pool)
        qqt_matrix = cbind(pool * pool)
    }else if (qt == 2){ 
        qt_matrix = cbind(log(pool))
        qqt_matrix = cbind(log(pool) * log(pool))
    }else if (qt == 3){ 
        qt_matrix = cbind(pool, pool^2)
        qqt_matrix = cbind(pool*pool, pool*pool^2, pool^2*pool, pool^2*pool^2)
    }else if (qt == 4){ 
        qt_matrix = cbind(log(pool), log(pool)**2) 
        qqt_matrix = cbind(log(pool)*log(pool), log(pool)*log(pool)^2, log(pool)^2*log(pool), log(pool)^2*log(pool)^2)
    }
    Qt_matrix = cbind(1, qt_matrix) 
    d = ncol(qt_matrix)
    
    # Compute QQt_matrix: (n, (d+1)^2)
    QQt_matrix = Qt_matrix
    for (i in 1:d) {
        QQt_matrix = cbind(QQt_matrix, qt_matrix[, i])
        QQt_matrix = cbind(QQt_matrix, qqt_matrix[, (1 + (i - 1) * d):(i * d)])
    }

    # Compute U matrix: (m*(d+1), m*(d+1))
    h0 = rho0
    h1 = rho1 * exp(Qt_matrix %*% theta1)  # (n,1)
    h2 = rho2 * exp(Qt_matrix %*% theta2)
    s = h0 + h1 + h2
    Ht_matrix = cbind(h1 - h1^2 / s, -h1 * h2 / s, -h1 * h2 / s, h2 - h2^2 / s)  # (n,m^2)
    U11 = (Ht_matrix[, 1] * P0) %*% QQt_matrix  # (1, (d+1)^2)
    U12 = (Ht_matrix[, 2] * P0) %*% QQt_matrix
    U21 = (Ht_matrix[, 3] * P0) %*% QQt_matrix
    U22 = (Ht_matrix[, 4] * P0) %*% QQt_matrix
    U = rbind(
        cbind(matrix(U11, nrow = d + 1), matrix(U12, nrow = d + 1)),
        cbind(matrix(U21, nrow = d + 1), matrix(U22, nrow = d + 1))
    )

    # Compute W matrix: (m*(d+1), m*(d+1))
    T = (1 / rho0) * (rep(1, m) %*% t(rep(1, m))) + diag(1 / rhos[2:(m + 1)])
    W = kronecker(T, diag(c(1, rep(0, d))))

    # Compute Sigma_Theta
    U.solve = solve(U)
    Sigma_Theta = U.solve - W

    # =========================================================================
    # Step 3: Compute Sigma_c (covariance matrix for cutoff points (c1, c2) )
    # =========================================================================
    
    A1 = cbind(diag(d + 1), diag(d + 1) * 0)
    A2 = cbind(diag(d + 1) * 0, diag(d + 1))
    
    # Evaluate basis functions at the cutoff points and compute derivatives
    param1 = -theta1[2:(d + 1)]
    param2 = theta1[2:(d + 1)] - theta2[2:(d + 1)]
    if (qt == 1){ 
        qc_matrix = rbind(c1, c2)
        expres_1 = 'param1 * c1'
        expres_2 = 'param2 * c2'
    }else if (qt == 2){ 
        qc_matrix = rbind(log(c1), log(c2))
        expres_1 = 'param1 * log(c1)'
        expres_2 = 'param2 * log(c2)'
    }else if (qt == 3){ 
        qc_matrix = rbind(c(c1, c1^2), c(c2, c2^2))
        expres_1 = paste(param1, c('c1', 'c1^2'), sep = "*", collapse = "+")
        expres_2 = paste(param2, c('c2', 'c2^2'), sep = "*", collapse = "+")
    }else if (qt == 4){ 
        qc_matrix = rbind(c(log(c1), log(c1)^2), c(log(c2), log(c2)^2))
        expres_1 = paste(param1, c('log(c1)', 'log(c1)^2'), sep = "*", collapse = "+")
        expres_2 = paste(param2, c('log(c2)', 'log(c2)^2'), sep = "*", collapse = "+")
    }
    dqc1 = eval(D(parse(text = expres_1), 'c1'))
    dqc2 = eval(D(parse(text = expres_2), 'c2'))

    # Variance estimates for c1 and c2
    Qc_matrix = cbind(1, qc_matrix)  # (m,d+1)
    Qc1 = Qc_matrix[1, ]
    Qc2 = Qc_matrix[2, ]
    sigma2_c1 = (t(Qc1) %*% A1 %*% Sigma_Theta %*% t(A1) %*% Qc1) / (dqc1^2)
    sigma2_c2 = (t(Qc2) %*% (A1 - A2) %*% Sigma_Theta %*% t(A1 - A2) %*% Qc2) / (dqc2^2)
    sigma_c12 = (-t(Qc1) %*% A1 %*% Sigma_Theta %*% t(A1 - A2) %*% Qc2) / (dqc1 * dqc2)

    sigma2_c1 = as.numeric(sigma2_c1)
    sigma2_c2 = as.numeric(sigma2_c2)
    sigma_c12 = as.numeric(sigma_c12)
    Sigma_c = matrix(c(sigma2_c1, sigma_c12, sigma_c12, sigma2_c2), nrow = 2)

    # Warning for invalid variance estimates
    if (warn) {
        if (is.na(sigma2_c1) | is.na(sigma2_c2)) {
            warning("Warning: sigma2_c1 or sigma2_c2 contains missing values.")
        } else if ((sigma2_c1 < 0) | (sigma2_c2 < 0)) {
            warning("Warning: Negative variance estimate detected. Sigma_Theta may not be positive definite.")
        }
    }
    if (is.na(sigma2_c1) || is.na(sigma2_c2) || sigma2_c1 < 0 || sigma2_c2 < 0) {
        youden_vec = matrix(c(J, NULL, NULL, NULL), nrow = 1, ncol = 4)
        cutoff_matrix = matrix(c(c1, NULL, NULL, NULL, c2, NULL, NULL, NULL), ncol = 4, nrow = 2, byrow = TRUE)
        colnames(youden_vec) = c('est', 'se', 'CI.lower', 'CI.upper')
        colnames(cutoff_matrix) = c('est', 'se', 'CI.lower', 'CI.upper')
        return(list(youden = youden_vec, cutoff = cutoff_matrix))
    }

    se_c1 = sqrt(sigma2_c1 / n)
    se_c2 = sqrt(sigma2_c2 / n)
    Sigma_c = Sigma_c / n

    # =========================================================================
    # Step 4: Compute Sigma_J (variance for generalized Youden index)
    # =========================================================================
    
    F0 = function(x) { sum(P0 * ifelse(pool <= x, 1, 0)) }
    F1 = function(x) { sum(P1 * ifelse(pool <= x, 1, 0)) }
    F2 = function(x) { sum(P2 * ifelse(pool <= x, 1, 0)) }
 
    a00 = function(x) { sum((h0 - h0^2 / s) * P0 * ifelse(pool <= x, 1, 0)) }
    a01 = function(x) { sum((-h0 * h1 / s) * P0 * ifelse(pool <= x, 1, 0)) }
    a02 = function(x) { sum((-h0 * h2 / s) * P0 * ifelse(pool <= x, 1, 0)) }
    a11 = function(x) { sum(Ht_matrix[, 1] * P0 * ifelse(pool <= x, 1, 0)) }
    a12 = function(x) { sum(Ht_matrix[, 2] * P0 * ifelse(pool <= x, 1, 0)) }
    a22 = function(x) { sum(Ht_matrix[, 4] * P0 * ifelse(pool <= x, 1, 0)) }

    B0 = function(x) {
        rbind(
            t(Qt_matrix) %*% (-h1 * (h0 / s) * P0 * ifelse(pool <= x, 1, 0)),
            t(Qt_matrix) %*% (-h2 * (h0 / s) * P0 * ifelse(pool <= x, 1, 0))
        )
    }
    B1 = function(x) {
        rbind(
            t(Qt_matrix) %*% (Ht_matrix[, 1] * P0 * ifelse(pool <= x, 1, 0)),
            t(Qt_matrix) %*% (Ht_matrix[, 2] * P0 * ifelse(pool <= x, 1, 0))
        )
    }
    B2 = function(x) {
        rbind(
            t(Qt_matrix) %*% (Ht_matrix[, 3] * P0 * ifelse(pool <= x, 1, 0)),
            t(Qt_matrix) %*% (Ht_matrix[, 4] * P0 * ifelse(pool <= x, 1, 0))
        )
    }

    # Evaluate at the optimal cutoff c1 and c2
    F0_c1 = F0(c1)
    F0_c2 = F0(c2)
    F1_c1 = F1(c1)
    F1_c2 = F1(c2)
    F2_c1 = F2(c1)
    F2_c2 = F2(c2)
    a00_c1 = a00(c1)
    a01_c1 = a01(c1)
    a02_c1 = a02(c1)
    a11_c1 = a11(c1)
    a11_c2 = a11(c2)
    a12_c1 = a12(c1)
    a12_c2 = a12(c2)
    a22_c2 = a22(c2)
    B0_c1 = B0(c1)
    B1_c1 = B1(c1)
    B1_c2 = B1(c2)
    B2_c2 = B2(c2)

    # Compute Omega matrix (4x4)
    O11 = (F0_c1 - F0_c1^2) / rho0 - (a00_c1 - t(B0_c1) %*% U.solve %*% B0_c1) / (rho0^2)
    O22 = (F1_c1 - F1_c1^2) / rho1 - (a11_c1 - t(B1_c1) %*% U.solve %*% B1_c1) / (rho1^2)
    O33 = (F1_c2 - F1_c2^2) / rho1 - (a11_c2 - t(B1_c2) %*% U.solve %*% B1_c2) / (rho1^2)
    O44 = (F2_c2 - F2_c2^2) / rho2 - (a22_c2 - t(B2_c2) %*% U.solve %*% B2_c2) / (rho2^2)
    O12 = -(a01_c1 - t(B0_c1) %*% U.solve %*% B1_c1) / (rho0 * rho1)
    O13 = -(a01_c1 - t(B0_c1) %*% U.solve %*% B1_c2) / (rho0 * rho1)
    O14 = -(a02_c1 - t(B0_c1) %*% U.solve %*% B2_c2) / (rho0 * rho2)
    O23 = -(a11_c1 - t(B1_c1) %*% U.solve %*% B1_c2) / (rho1^2)
    O24 = -(a12_c1 - t(B1_c1) %*% U.solve %*% B2_c2) / (rho1 * rho2)
    O34 = -(a12_c2 - t(B1_c2) %*% U.solve %*% B2_c2) / (rho2^2)
    Omega = diag(4)
    Omega[upper.tri(Omega)] = c(O12, O13, O14, O23, O24, O34)
    Omega = Omega + t(Omega)
    diag(Omega) = c(O11, O22, O33, O44)

    # Compute sigma_J
    avec = matrix(rep(c(1, -1), m) / m, nrow = 1)
    sigma2_J = as.numeric(avec %*% Omega %*% t(avec))
    se_J = as.numeric(sqrt(sigma2_J / n))

    # =========================================================================
    # Step 5: Construct confidence intervals
    # =========================================================================
    
    z_q = -qnorm((1-level)/2)
    
    CI_c1 = c(c1, c1) + c(-1, 1) * z_q * se_c1
    CI_c2 = c(c2, c2) + c(-1, 1) * z_q * se_c2
    CI_J = c(J, J) + c(-1, 1) * z_q * se_J

    # Generate output
    youden_vec = matrix(c(J, se_J, CI_J), nrow = 1)
    colnames(youden_vec) = c('est', 'se', 'CI.lower', 'CI.upper')

    cutoff_matrix = matrix(c(c1, se_c1, CI_c1, c2, se_c2, CI_c2), ncol = 4, byrow = TRUE)
    colnames(cutoff_matrix) = c('est', 'se', 'CI.lower', 'CI.upper')

    list(youden = youden_vec, cutoff = cutoff_matrix, Sigma_c = Sigma_c)
}

#########################
# # A simulation example:
#
# set.seed(3710)
# obs = list(
#     "1" = rgamma(n=100, shape=1, rate=1),
#     "2" = rgamma(n=100, shape=2, rate=1),
#     "3" = rgamma(n=100, shape=2.838, rate=1))
# out = DRM.multiclass.Youden(obs, qt = 2)
# print(out)
#
# # The output looks like:
#
# > $youden
# >            est         se  CI.lower  CI.upper
# > [1,] 0.3422817 0.03666769 0.2704144 0.4141491

# > $cutoff
# >           est         se  CI.lower CI.upper
# > [1,] 1.020739 0.06741537 0.8886077 1.152871
# > [2,] 2.076907 0.10000526 1.8809005 2.272914

# > $Sigma_c
# >             [,1]       [,2]
# > [1,] 0.004544832 0.00341016
# > [2,] 0.003410160 0.01000105


#########################
# # Application to the Parkinson's disease data
# # To analysis the Parkinson's disease data, you need to first download 
#  the "PDdatabook.xlsx" dataset from https://doi.org/10.13140/RG.2.2.32452.86402

# # load data
# data=read.csv('./parkinson_disease.csv')
# n.cols = ncol(data)
# data[, 1:n.cols] = sapply(data[, 1:n.cols], as.numeric)
# group = 3 - data[, 1]

# obs_moca = split(data[, 'MoCA'], group)
# obs_moca = lapply(1:3, function(i) { na.exclude(obs_moca[[i]]) })

# obs_mmse = split(data[, 'S.MMSE'], group)

# # goodness-of-fit test for the DRM based on MoCA and S-MMSE
# set.seed(3710)
# DRM.gof(obs_moca, qt = 2, B = 500)

# # The output looks like:
# > $stat
# > [1] 0.9310599
# > 
# > $pval
# > [1] 0.326

# set.seed(3710)
# DRM.gof(obs_mmse, qt = 2, B = 500)

# # The output looks like:
# > $stat
# > [1] 1.071284
# > 
# > $pval
# > [1] 0.312

# # the proposed estimates of J, c1, c2 and 95% CIs based on MoCA and S-MMSE
# DRM.multiclass.Youden(obs_moca, qt = 2)

# # The output looks like:
# > $youden
# >            est         se  CI.lower CI.upper
# > [1,] 0.5635407 0.05474046 0.4562514  0.67083
# > 
# > $cutoff
# >           est        se CI.lower CI.upper
# > [1,] 20.83053 0.4377166 19.97262 21.68844
# > [2,] 25.20482 0.2500028 24.71483 25.69482
# > 
# > $Sigma_c
# >            [,1]       [,2]
# > [1,] 0.19159578 0.05045444
# > [2,] 0.05045444 0.06250140

# DRM.multiclass.Youden(obs_mmse, qt = 2)

# # The output looks like:
# > $youden
# >            est        se  CI.lower  CI.upper
# > [1,] 0.4233745 0.0599374 0.3058994 0.5408497
# > 
# > $cutoff
# >           est        se CI.lower CI.upper
# > [1,] 24.38964 0.3674328 23.66948 25.10980
# > [2,] 27.13688 0.2067437 26.73167 27.54209
# > 
# > $Sigma_c
# >            [,1]       [,2]
# > [1,] 0.13500687 0.03926525
# > [2,] 0.03926525 0.04274295


