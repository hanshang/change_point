#################################
# Regression-based change points
#################################

# Nonstationary functional dynamic factor model
# data: p by n data matrix
# k_0: large enough number of principal components
# test_alpha: level of significance for the independence test

NFDF_model <- function(data, k_0 = 10, n.ahead = 1)
{
    n_row = nrow(data)
    n_col = ncol(data)

    # computing the mean of functional data

    X_bar = rowMeans(data)

    # de-centering functional data

    center_fun_data = t(scale(t(data), center = TRUE, scale = FALSE))

    # take first-order differencing of the original curve time series

    diff_data = t(diff(t(center_fun_data)))

    # compute the estimated long-run covariance

    est_Lambda = long_run_covariance_estimation(dat = diff_data)$C_0_est

    # eigen-decomposition of the estimated long-run covariance

    eigen_Lambda = eigen(est_Lambda)
    eigen_Lambda_ratio = vector("numeric", k_0)
    for(ik in 1:k_0)
    {
        eigen_Lambda_ratio[ik] = eigen_Lambda$values[ik+1]/eigen_Lambda$values[ik]
        rm(ik)
    }
    tau = which.min(eigen_Lambda_ratio)
    loadings = matrix(eigen_Lambda$vectors[,1:tau], ncol = tau)
    scores = crossprod(center_fun_data, loadings)
    colnames(loadings) = colnames(scores) = 1:tau

    # forecast scores

    forecast_scores = matrix(NA, n.ahead, tau)
    for(ik in 1:tau)
    {
        forecast_scores[,ik] = forecast(auto.arima(scores[,ik]), h = n.ahead)$mean
        rm(ik)
    }
    colnames(forecast_scores) = 1:tau
    rownames(forecast_scores) = 1:n.ahead

    # forecast curves

    forecast_curve_main = loadings %*% t(forecast_scores)

    # obtain model residuals

    residuals = center_fun_data - loadings %*% t(scores)

    # long-run covariance of the residuals

    est_Lambda_residuals = long_run_covariance_estimation(dat = residuals)$C_0_est
    eigen_Lambda_residuals = eigen(est_Lambda_residuals)
    eigen_Lambda_residuals_ratio = vector("numeric", k_0)
    for(ik in 1:k_0)
    {
        eigen_Lambda_residuals_ratio[ik] = eigen_Lambda_residuals$values[ik+1]/eigen_Lambda_residuals$values[ik]
        rm(ik)
    }
    K_Z = which.min(eigen_Lambda_residuals_ratio)

    # perform independence test for functional time series

    colnames(residuals) = 1:n_col
    # GK_p_value = GK_test(data = residuals, K = K_Z)
    rm(diff_data); rm(est_Lambda); rm(eigen_Lambda); rm(eigen_Lambda_ratio)
    rm(eigen_Lambda_residuals_ratio)

    loadings_residuals = matrix(eigen_Lambda_residuals$vectors[,1:K_Z], ncol = K_Z)
    scores_residuals = crossprod(residuals, loadings_residuals)
    colnames(loadings_residuals) = colnames(scores_residuals) = 1:K_Z

    # forecast scores

    forecast_scores_residual = matrix(NA, n.ahead, K_Z)
    for(ik in 1:K_Z)
    {
        forecast_scores_residual[,ik] = forecast(auto.arima(scores_residuals[,ik]), h = n.ahead)$mean
    }
    colnames(forecast_scores_residual) = 1:K_Z
    rownames(forecast_scores_residual) = 1:n.ahead

    # forecast curves

    forecast_curve_residuals = loadings_residuals %*% t(forecast_scores_residual)

    # add all components

    forecast_curve_combo = forecast_curve_main + forecast_curve_residuals + matrix(rep(X_bar, n.ahead), n_row, n.ahead)
    rownames(forecast_curve_combo) = 1:n_row
    return(list(forecast_curve = forecast_curve_combo, ncomp_mean = tau, ncomp_resi = K_Z))
}

## female

AUS_mort_fun <- function(data_set, ini_point = 3)
{
    year = data_set$year
    err_female = matrix(NA, 101, length(year) - ini_point)
    for(iw in ini_point:(length(year) - 1))
    {
        # extract years

        AUS_female_obj = extract.years(data_set, years = year[1:iw])

        # apply NFDF model

        NFDF_female = NFDF_model(data = log(AUS_female_obj$rate$Female, base = 10))

        # compute errors

        err_female[,iw - (ini_point - 1)] = log(extract.years(data_set, years = year[(iw + 1)])$rate$Female, base = 10) - NFDF_female$forecast_curve
        print(iw+1); rm(iw); rm(NFDF_female); rm(AUS_female_obj)
    }
    colnames(err_female) = (ini_point + 1):length(year)
    ISFE_female = ts(colSums(err_female^2), start = (ini_point + 1), end = length(year))
    
    # classical method
    strucchange_detect_female = year[(4:length(year))[breakpoints(ISFE_female ~ 1)$breakpoints[1]]]
    
    # ECP method (ignored)
    ecp_detect_female = year[(4:length(year))[e.divisive(X = matrix(ISFE_female, ncol = 1), k = 1)$estimates[2]]]
    return(c(strucchange_detect_female, ecp_detect_female)) # 1976 1977
}

AUS_mort_fun(data_set = AUS_demo_female_trunc)  # 1976 1977
AUS_mort_fun(data_set = extract.years(AUS_demo_female_trunc, 1921:2010)) # 1974 1977

## male

AUS_mort_male_fun <- function(data_set, ini_point = 3)
{
    year = data_set$year
    err_male = matrix(NA, 101, length(year) - ini_point)
    for(iw in ini_point:(length(year) - 1))
    {
        # extract years

        AUS_male_obj = extract.years(data_set, years = year[1:iw])

        # apply NFDF model
  
        NFDF_male = NFDF_model(data = log(AUS_male_obj$rate$Male, base = 10))

        # compute errors

        err_male[,iw - (ini_point - 1)] = log(extract.years(data_set, years = year[(iw + 1)])$rate$Male, base = 10) - NFDF_male$forecast_curve
        print(iw + 1); rm(iw); rm(NFDF_male); rm(AUS_male_obj)
    }
    colnames(err_male) = (ini_point + 1):length(year)
    ISFE_male = ts(colSums(err_male^2), start = (ini_point + 1), end = length(year))
    strucchange_detect_male = year[(4:length(year))[max(which((breakfactor(breakpoints(ISFE_male ~ 1), breaks = 1)) == "segment1"))]]
    ecp_detect_male = year[(4:length(year))[e.divisive(X = matrix(ISFE_male, ncol = 1), k = 1)$estimates[2]]]
    return(c(strucchange_detect_male, ecp_detect_male)) 
}

AUS_mort_male_fun(data_set = AUS_demo_male_trunc) # 1982 1983
AUS_mort_male_fun(data_set = extract.years(AUS_demo_male_trunc, 1921:2010)) # 1982 1981

