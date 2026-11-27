# multiclass-Youden
This is the R code to implement the proposed semiparametric procedures and reproducing the numerical results in the associated paper ``Semiparametric Estimation of the Generalized Youden Index and the Optimal Cutoff Points for Multiclass Diagnostic Tests''.

The main function is DRM.multiclass.Youden. 
It computes point and interval estimates for the generalized Youden index J_3 and the optimal cutoff points c_1, c_2 
under the density ratio model for three-class diagnostic test.
    
    - Inputs:
       - obs     : List of 3 numeric vectors. E.g. list(class0_data, class1_data, class2_data)
       - qt      : Integer specifying the basis function q(x):
                   1 = x (default); 2 = log(x); 3 = c(x, x^2); 4 = c(log(x), log(x)^2).
                   Note that when qt = 2 or 4, all observations must be positive.
       - level   : confidence level for the confidence intervals (default = 0.95).
       - warn    : Logical; if TRUE, prints warnings (default = TRUE).
    
    - Outputs:
       - youden  : (1,4) matrix of the generalized Youden index estimates. Columns: est, se, CI.lower, CI.upper.
       - cutoff  : (2,4) matrix of optimal cutoff point estimates. Rows: c1, c2. Columns: est, se, CI.lower, CI.upper.
       - Sigma_c : (2,2) asymptotic covariance matrix of vector c.
