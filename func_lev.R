##########################################
### This supplementary R file provides ###
###    functions for the main files    ###
##########################################
####      The code is developed by     ###
####     @ Juri Trifonov, NRU HSE      ###
####     @ Bogdan Potanin, NRU HSE     ###
##########################################


#####################################################
##### Flexible data generating process function #####
#####################################################
GARCH.Simulate <- function (n,                         # number of observations
                            mu, omega, alpha, beta,    # parameters
                            gamma_2 = 0,
                            lambda,
                            lambda_2 = 0, 
                            family = 'GARCH-M')                
{
  # simulate standardized shocks
  xi <- rnorm(n)
  # initialize main variables
  eps <- rep(NA, n)                               # shocks
  sigma <- rep(NA, n)                             # conditional volatility
  y <- rep(NA, n)
  
  # set initial values
  ####### Unconditional Variance ########
  if (family == "GARCH-M")
  {
    var_uncond <- omega / (1 - alpha - beta)
    
    ksi <- (omega ^ 2 + 2 * omega * alpha * var_uncond + 2 * omega * beta * var_uncond) / 
              (1 - 3 * alpha ^ 2 - beta ^ 2 - 2 * alpha * beta)
    
    sigma_0 <- lambda ^ 2 * ((alpha ^ 2 * (3 * ksi - var_uncond ^ 2) + 2 * alpha * beta * (ksi - var_uncond ^ 2)) / 
                 (1 - beta ^ 2)) + var_uncond
  }
  
  if (family == "GARCH-M-GJR")
  {
    var_uncond <- omega / (1 - alpha - beta - gamma_2 / 2)
    
    exp_sigma_4 <- (omega ^ 2 + omega * var_uncond * (2 * alpha + 2 * beta + gamma_2)) / 
                      (1 - 3 * alpha ^ 2 - beta ^ 2 - 1.5 * gamma_2 ^ 2 - 2 * alpha * beta - 
                         3 * alpha * gamma_2 - beta * gamma_2)
    
    var_sigma_2 <- (alpha ^ 2 * (3 * exp_sigma_4 - var_uncond^2) + 2 * alpha * beta * (exp_sigma_4 - var_uncond^2) + 
                      gamma_2 ^ 2 * (1.5 * exp_sigma_4 - 0.25 * var_uncond^2) + 
                      2 * (alpha * gamma_2 * (1.5 * exp_sigma_4 - var_uncond ^ 2 / 2) + 
                             beta * gamma_2 * 0.5 * (exp_sigma_4 - var_uncond ^ 2))) / (1 - beta ^ 2)
    
    sigma_0 <- lambda ^ 2 * var_sigma_2 + var_uncond
  }
  
  if (family == "GARCH-M-GJR-LEV")
  {
    var_uncond <- omega / (1 - alpha - beta - gamma_2 / 2)
    exp_sigma_4 <- (omega ^ 2 + omega * var_uncond * (2 * alpha + 2 * beta + gamma_2)) / 
                      (1 - 3 * alpha ^ 2 - beta ^ 2 - 1.5 * gamma_2 ^ 2 - 2 * alpha * beta - 
                         3 * alpha * gamma_2 - beta * gamma_2)
    var_sigma_2 <- (alpha ^ 2 * (3 * exp_sigma_4 - var_uncond^2) + 2 * alpha * beta * (exp_sigma_4 - var_uncond^2) + 
                      gamma_2 ^ 2 * (1.5 * exp_sigma_4 - 0.25 * var_uncond^2) + 
                         2 * (alpha * gamma_2 * (1.5 * exp_sigma_4 - var_uncond ^ 2 / 2) + 
                                beta * gamma_2 * 0.5 * (exp_sigma_4 - var_uncond ^ 2))) / (1 - beta ^ 2)

    sigma_0 <- lambda ^ 2 * var_sigma_2 + lambda_2 ^ 2 * (0.5 * exp_sigma_4 - 0.25 * var_uncond^2) + 
                 lambda * lambda_2 * (exp_sigma_4 - var_uncond^2) + var_uncond
  }

  if(family == "GARCH-M-LEV")
  {
    var_uncond <- omega / (1 - alpha - beta)
    exp_sigma_4 <- (omega ^ 2 + 2 * omega * var_uncond * (alpha + beta)) / 
      (1 - 3 * alpha ^ 2 - beta ^ 2 - 2 * alpha * beta)
    
    var_sigma_2 <- (alpha ^ 2 * (3 * exp_sigma_4 - var_uncond^2) + 2 * alpha * beta * (exp_sigma_4 - var_uncond^2)) / (1 - beta ^ 2)
    
    sigma_0 <- lambda ^ 2 * var_sigma_2 + lambda_2 ^ 2 * (0.5 * exp_sigma_4 - 0.25 * var_uncond^2) + 
      lambda * lambda_2 * (exp_sigma_4 - var_uncond^2) + var_uncond
  }
  
  sigma[1] <- sqrt(sigma_0)                         # conditional volatility at 
                                                    # the first period
  eps[1] <- sigma[1] * xi[1]                        # shock in the first period
  if (family != 'GARCH-M')
  {
  is_lower <- rep(NA, n)
  is_lower[1] <- as.numeric(eps[1] < 0)
  }
  y[1] <- eps[1] + mu
  
  # validation conditions
  if (is.nan(sqrt(sigma_0)))
  {
    return(-9999999999999)
  }
  
  if (alpha <= 0 | beta <= 0 |(alpha + beta) > 1 | omega <= 0)
  {
    return(-9999999999999)
  }
  
  # calculate conditional volatility
  for(t in 2:n)
  {
    if (family == "GARCH-M")
    {
    sigma[t] <- sqrt(omega +
                       alpha * eps[t - 1] ^ 2 + 
                       beta * sigma[t - 1] ^ 2)
    eps[t] <- sigma[t] * xi[t]
    y[t] <- eps[t] + mu + lambda * sigma[t-1] ^ 2 
    }
    
    if (family == "GARCH-M-GJR")
    {
      sigma[t] <- sqrt(omega +
                         alpha * eps[t - 1] ^ 2 + 
                         beta * sigma[t - 1] ^ 2 +
                         gamma_2 * is_lower[t - 1] * eps[t - 1] ^ 2)
      eps[t] <- sigma[t] * xi[t]
      is_lower[t] <- as.numeric(eps[t] < 0)
      y[t] <- eps[t] + mu + lambda * sigma[t-1] ^ 2 
    }
    
    if (family == "GARCH-M-GJR-LEV")
    {
      sigma[t] <- sqrt(omega +
                         alpha * eps[t - 1] ^ 2 + 
                         beta * sigma[t - 1] ^ 2 +
                         gamma_2 * is_lower[t - 1] * eps[t-1] ^ 2)
      eps[t] <- sigma[t] * xi[t]
      is_lower[t] <- as.numeric(eps[t] < 0)
      y[t] <- eps[t] + mu + lambda * sigma[t-1] ^ 2 +
                            lambda_2 * is_lower[t-1] * sigma[t-1] ^ 2
    }
    if (family == "GARCH-M-LEV")
    {
      sigma[t] <- sqrt(omega +
                         alpha * eps[t - 1] ^ 2 + 
                         beta * sigma[t - 1] ^ 2)
      eps[t] <- sigma[t] * xi[t]
      is_lower[t] <- as.numeric(eps[t] < 0)
      y[t] <- eps[t] + mu + lambda * sigma[t-1] ^ 2 +
                            lambda_2 * is_lower[t-1] * sigma[t-1] ^ 2
    }
  }
  
  # aggregate the results
  if(family == 'GARCH-M')
  {
  df <- data.frame(sigma = sigma, 
                   y = y,
                   eps = eps,
                   eps_std = eps / sigma,
                   xi = xi)
  }
  if(family != 'GARCH-M')
  {
  df <- data.frame(sigma = sigma, 
                   y = y,
                   eps = eps,
                   eps_std = eps / sigma,
                   is_lower = is_lower, 
                   xi = xi,
                   var_uncond = sigma_0)
  }
  
  return(df)
}

#####################################################
############ Flexible Likelihood Function ###########
#####################################################

lnL <- function (par,
                 data,
                 type = "likelihood",
                 is_aggregate = TRUE,
                 family = 'GARCH-M',
                 is_eps = FALSE,
                 is_var_uncond = FALSE,
                 is_rp = FALSE)
  
  
{ 
  # estimated parameters
  if (family == "GARCH-M")
  {
    mu <- par[1]                               
    omega <- par[2]
    alpha <- par[3]
    beta <- par[4]
    lambda <- par[5]
  }
  if (family == "GARCH-M-GJR")
  {
    mu <- par[1]                               
    omega <- par[2]
    alpha <- par[3]
    beta <- par[4]
    gamma_2 <- par[5]
    lambda <- par[6]
  }
  if (family == "GARCH-M-GJR-LEV")
  {
    mu <- par[1]                               
    omega <- par[2]
    alpha <- par[3]
    beta <- par[4]
    gamma_2 <- par[5]
    lambda <- par[6]
    lambda_2 <- par[7]
  }
  if (family == "GARCH-M-LEV")
  {
    mu <- par[1]                               
    omega <- par[2]
    alpha <- par[3]
    beta <- par[4]
    lambda <- par[5]
    lambda_2 <- par[6]
  }
  
  y <- data$y
  n = length(y)                            # number of observations
  
  sigma <- rep(NA, n)                      # conditional volatility
  
  if(family != 'GARCH-M')
  {
    is_lower <- rep(NA, n)
  }

  eps <- rep(NA, n)
  risk_prem <- rep(NA, n)
  risk_prem[1] <- 0
  ### specification of unconditional variance ###
  
  if (family == "GARCH-M")
  {
    var_uncond <- omega / (1 - alpha - beta)
    
    ksi <- (omega ^ 2 + 2 * omega * alpha * var_uncond + 2 * omega * beta * var_uncond) / 
              (1 - 3 * alpha ^ 2 - beta ^ 2 - 2 * alpha * beta)
    
    sigma_0 <- lambda ^ 2 * ((alpha ^ 2 * (3 * ksi - var_uncond ^ 2) + 2 * alpha * beta * (ksi - var_uncond ^ 2)) / 
                 (1 - beta ^ 2)) + var_uncond
  }
  
  if (family == "GARCH-M-GJR")
  {
    var_uncond <- omega / (1 - alpha - beta - gamma_2 / 2)
    
    exp_sigma_4 <- (omega ^ 2 + omega * var_uncond * (2 * alpha + 2 * beta + gamma_2)) / 
                      (1 - 3 * alpha ^ 2 - beta ^ 2 - 1.5 * gamma_2 ^ 2 - 2 * alpha * beta - 
                         3 * alpha * gamma_2 - beta * gamma_2)
    
    var_sigma_2 <- (alpha ^ 2 * (3 * exp_sigma_4 - var_uncond^2) + 2 * alpha * beta * (exp_sigma_4 - var_uncond^2) + 
                      gamma_2 ^ 2 * (1.5 * exp_sigma_4 - 0.25 * var_uncond^2) + 
                         2 * (alpha * gamma_2 * (1.5 * exp_sigma_4 - var_uncond ^ 2 / 2) + 
                            beta * gamma_2 * 0.5 * (exp_sigma_4 - var_uncond ^ 2))) / (1 - beta ^ 2)
    
    sigma_0 <- lambda ^ 2 * var_sigma_2 + var_uncond
  }
  
  if (family == "GARCH-M-GJR-LEV")
  {
    var_uncond <- omega / (1 - alpha - beta - gamma_2 / 2)
    exp_sigma_4 <- (omega ^ 2 + omega * var_uncond * (2 * alpha + 2 * beta + gamma_2)) / 
                      (1 - 3 * alpha ^ 2 - beta ^ 2 - 1.5 * gamma_2 ^ 2 - 2 * alpha * beta - 
                         3 * alpha * gamma_2 - beta * gamma_2)
    
    var_sigma_2 <- (alpha ^ 2 * (3 * exp_sigma_4 - var_uncond^2) + 2 * alpha * beta * (exp_sigma_4 - var_uncond^2) + 
                      gamma_2 ^ 2 * (1.5 * exp_sigma_4 - 0.25 * var_uncond^2) + 
                        2 * (alpha * gamma_2 * (1.5 * exp_sigma_4 - var_uncond ^ 2 / 2) + 
                         beta * gamma_2 * 0.5 * (exp_sigma_4 - var_uncond ^ 2))) / (1 - beta ^ 2)
    
    sigma_0 <- lambda ^ 2 * var_sigma_2 + lambda_2 ^ 2 * (0.5 * exp_sigma_4 - 0.25 * var_uncond^2) + 
                lambda * lambda_2 * (exp_sigma_4 - var_uncond^2) + var_uncond
  }
  if(family == "GARCH-M-LEV")
  {
    var_uncond <- omega / (1 - alpha - beta)
    exp_sigma_4 <- (omega ^ 2 + 2 * omega * var_uncond * (alpha + beta)) / 
      (1 - 3 * alpha ^ 2 - beta ^ 2 - 2 * alpha * beta)
    
    var_sigma_2 <- (alpha ^ 2 * (3 * exp_sigma_4 - var_uncond^2) + 2 * alpha * beta * (exp_sigma_4 - var_uncond^2)) / (1 - beta ^ 2)
    
    sigma_0 <- lambda ^ 2 * var_sigma_2 + lambda_2 ^ 2 * (0.5 * exp_sigma_4 - 0.25 * var_uncond^2) + 
      lambda * lambda_2 * (exp_sigma_4 - var_uncond^2) + var_uncond
  }
  
  
  sigma[1] <- sqrt(sigma_0)
  eps[1] <- y[1] - mu
  if (family != "GARCH-M")
  {
  is_lower[1] <- eps[1] < 0
  }
  
  # Validation
  if (is.nan(sigma_0))
  {
    return(-9999999999999)
  }
  
  if (sigma_0 < 0)
  {
    return(-9999999999999)
  }
  
  if ((alpha < 0) | (beta < 0) |((alpha + beta) > 1) | (omega <= 0) | (mu <= 0))
  {
    return(-9999999999999)
  }
  
  # Recursive calculations
  for(t in 2:n)
  {
    if (family == "GARCH-M")
    {
      sigma[t] <- sqrt(omega +
                         alpha * eps[t - 1] ^ 2 + 
                         beta * sigma[t - 1] ^ 2)
      eps[t] <- y[t] - mu - lambda * sigma[t-1] ^ 2
    }
    
    if (family == "GARCH-M-GJR")
    {
      sigma[t] <- sqrt(omega +
                         alpha * eps[t - 1] ^ 2 + 
                         beta * sigma[t - 1] ^ 2 +
                         gamma_2 * is_lower[t - 1] * eps[t - 1] ^ 2)
      eps[t] <- y[t] - mu - lambda * sigma[t-1] ^ 2
      is_lower[t] <- as.numeric(eps[t] < 0)
    }
    
    if (family == "GARCH-M-GJR-LEV")
    {
      sigma[t] <- sqrt(omega +
                         alpha * eps[t - 1] ^ 2 + 
                         beta * sigma[t - 1] ^ 2 +
                         gamma_2 * is_lower[t - 1] * eps[t - 1] ^ 2)
      eps[t] <- y[t] - mu - lambda * sigma[t-1] ^ 2 - 
                            lambda_2 * is_lower[t-1] * sigma[t-1] ^ 2
      is_lower[t] <- as.numeric(eps[t] < 0)
      risk_prem[t] <- lambda * sigma[t-1] ^ 2 + lambda_2 * is_lower[t-1] * sigma[t-1] ^ 2
    }
    
    if (family == "GARCH-M-LEV")
    {
      sigma[t] <- sqrt(omega +
                         alpha * eps[t - 1] ^ 2 + 
                         beta * sigma[t - 1] ^ 2)
      eps[t] <- y[t] - mu - lambda * sigma[t-1] ^ 2 - 
                            lambda_2 * is_lower[t-1] * sigma[t-1] ^ 2
      is_lower[t] <- as.numeric(eps[t] < 0)
    }
  }
  
  L_vector <- dnorm(eps, mean = 0, sd = sigma)
  
  
  lnL_value <- log(L_vector)               
  
  if (type == "sigma")
  {
    return(sigma)
  }
  if(type == 'return')
  {
    y_est <- rep(NA, n)
    y_est[1] <- mu
    for(t in 2:n)
    {
      if (family == "GARCH-M")
      {
        y_est[t] <- mu + lambda * sigma[t-1] ^ 2
      }
      if (family == "GARCH-M-GJR")
      {
        y_est[t] <- mu + lambda * sigma[t-1] ^ 2
      }
      if (family == "GARCH-M-GJR-LEV")
      {
        y_est[t] <- mu + lambda * sigma[t-1] ^ 2  + lambda_2 * is_lower[t-1] * sigma[t-1] ^ 2
      }
      if (family == "GARCH-M-LEV")
      {
        y_est[t] <- mu + lambda * sigma[t-1] ^ 2  + lambda_2 * is_lower[t-1] * sigma[t-1] ^ 2
      }
    }
    
    return(y_est)
  }
  if(!is_aggregate)
  {
    return(lnL_value)
  }
  if(is_eps)
  {
    return(eps)
  }
  if(is_var_uncond)
  {
    return(sigma_0)
  }
  if(is_rp)
  {
    return(risk_prem)
  }
  
  return(sum(lnL_value))                   
}


#####################################################
########## Flexible Optimization Algorithm ##########
#####################################################

GARCH <- function(data, x0, family, ga_iter = 100)             # input values are data, initial vector of points, 
                                                         # model family, and number of ga iterations
  
{
  result <- optim(par = x0,                              # initial point
                  method = "Nelder-Mead",                # optimization algorithm
                  fn = lnL,                              # likelihood function
                  control = list(maxit = 100000,         # in order to transform minimized function
                                 fnscale = -1,           # to maximized, we multiply it by -1
                                 reltol = 1e-10),        # set high accuracy           
                  hessian = FALSE,                       # hessian
                  data = data,
                  family = family)                       # arguments of optimized function
  
  ### Additional genetic algorithm optimization ###
  ga_result <- gena(
    fn = lnL,                               # maximized function
    lower = -abs(result$par) * 10,          # vectors of upper and lower bounds,
    upper = abs(result$par) * 10,           # where optimal parameters are being searched
    pop.initial = result$par,
    hybrid.prob = 0.1, 
    mutation.method = "percent",
    data = data,
    family = family,
    maxiter = ga_iter)                             
  
  gamma_est <- ga_result$par               # estimates of coefficients
  
  sigma_est <- lnL(gamma_est, data, type = 'sigma', family = family)
  y_est <- lnL(gamma_est, data, type = 'return', family = family)
  risk_prem <- lnL(gamma_est, data, is_rp = TRUE, family = family)
  eps <- lnL(gamma_est, data, type = 'likelihood', family = family, is_eps = TRUE)
  
  return_list <- list("gamma" = gamma_est,               # coefficient estimates
                      'sigma' = sigma_est,               # volatility etimates
                      'return' = y_est,                  # return estimates
                      'risk premium' = risk_prem,
                      'eps' = eps,
                      "data" = data,                     # dataframe that was used
                      "lnL" = ga_result$value) # log likelihood value
  
  
  class(return_list) <- "GARCH"                          # class
  
  return(return_list)                                    # overall result 
}


# Calculation of accuracy metrics 
rmse_func <- function(real, predicted)
{
  dif <- (real - predicted) ^ 2
  #row_sum <- sum(dif)
  rmse_value <- sqrt(mean(dif))
  return(rmse_value * 100)
}

mae_func <- function(real, predicted)
{
 # n = length(real)
  dif <- abs(real - predicted)
  #row_sum <- sum(dif)
  mae_value <- mean(dif)
  return(mae_value * 100)
}

mape_func <- function(real, predicted)
{
  #n = length(real)
  dif <- abs((real - predicted)/ real)
  #row_sum <- sum(dif)
  mape_value <- (mean(dif)) * 100
  return(mape_value)
}

aic_func <- function(model)
{
  aic <- 2 * length(model$gamma) - 2 * model$lnL
  return(aic)
}

bic_func <- function(model)
{
  bic <- length(model$gamma) * log(length(model$data[, 1])) - 2 * model$lnL
}


########################################################
################ Flexible Summary Method ###############
########################################################
summary.GARCH <- function(model, cov_mle, family)        # arguments are: model, estimate of cov. matrix, 
{                                                        # and family type
  n <- length(as.matrix(model$data))                     # number of observations
  
  cat(family, "Estimation Results\n")                   # title
  
  cat(paste("Observations:", n, "\n"))                   # print number of observations
  cat(paste("Log-likelihood:",                           # print log likelihood value
            round(model$lnL, 3), "\n"))  
  cat(paste("AIC:",                                      # print aic value
            round(aic_func(model), 3), "\n")) 
  cat(paste("BIC:",                                      # print aic value
            round(bic_func(model), 3), "\n")) 
  
  cat("---\n")                                           # visualization
  
  cat("Coefficients:\n")                                 # visualization
  
  ### List estimates of coefficients depending on the family type ###
  if (family == 'GARCH-M')
  {
  
  gamma_est <- c('mu' = model$gamma[1], 'omega' = model$gamma[2], 'alpha' = model$gamma[3], 
                 'beta' = model$gamma[4], 'lambda' = model$gamma[5])
  }
  
  if (family =="GARCH-M-GJR")
  {
    gamma_est <- c('mu' = model$gamma[1], 'omega' = model$gamma[2], 'alpha' = model$gamma[3], 
                   'beta' = model$gamma[4],'gamma' = model$gamma[5], 'lambda' = model$gamma[6])
  }
  
  if(family == "GARCH-M-GJR-LEV")
  {
    gamma_est <- c('mu' = model$gamma[1], 'omega' = model$gamma[2], 'alpha' = model$gamma[3], 
                   'beta' = model$gamma[4],'gamma' = model$gamma[5], 'lambda' = model$gamma[6], 
                   'lambda_2' = model$gamma[7])
  }
  if(family == "GARCH-M-GJR-LEV")
  {
    gamma_est <- c('mu' = model$gamma[1], 'omega' = model$gamma[2], 'alpha' = model$gamma[3], 
                  'beta' = model$gamma[4], 'gamma' = model$gamma[5], 'lambda' = model$gamma[6], 
                  'lambda_2' = model$gamma[7])
  }
  
  as_std_est <- sqrt(diag(cov_mle))                      # calculate as. estimates of
                                                         # standard errors
                                                         # of coefficient estimates
  
  z <- gamma_est / as_std_est                            # calculate z statistics
  p_value <- 2 * pmin(pnorm(z), 1 - pnorm(z))            # рcalculate p-values
  m <- length(gamma_est)                                 # number of estimated parameters
  stars <- rep("", m)                                    # stars for visualization
  stars[p_value <= 0.001] <- "***"
  stars[(p_value > 0.001) & (p_value < 0.01)] <- "**"
  stars[(p_value > 0.01) & (p_value <= 0.05)] <- "*"
  stars[(p_value > 0.05) & (p_value <= 0.1)] <- "."
  
  df <- data.frame("Estimate" = gamma_est,               # construct the final dataframe
                   "As.Std" = as_std_est,
                   "Test.Statistic" = z,
                   "p-value" = p_value,
                   "Significance" = stars)
  is_df_num_col <- sapply(df, is.numeric)                # point out 'numeric' columns in the dataframe
  df[, is_df_num_col] <- round(df[, is_df_num_col], 5)   # and round them up to a 5th sign
  print(df)                                              # print the dataframe
  cat("---\n")                                           # visualization
  cat(paste("Signif. codes:  0 ‘***’ 0.001 ‘**’",        # legend for the output 
            "0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1\n"))                      
}

forecast_garch <- function (par,
                            data,
                            family = 'GARCH-M',
                            steps_ahead = 1)
{
  # estimated parameters
  if (family == "GARCH-M")
  {
    mu <- par[1]                               
    omega <- par[2]
    alpha <- par[3]
    beta <- par[4]
    lambda <- par[5]
  }
  if (family == "GARCH-M-GJR")
  {
    mu <- par[1]                               
    omega <- par[2]
    alpha <- par[3]
    beta <- par[4]
    gamma_2 <- par[5]
    lambda <- par[6]
  }
  if (family == "GARCH-M-GJR-LEV")
  {
    mu <- par[1]                               
    omega <- par[2]
    alpha <- par[3]
    beta <- par[4]
    gamma_2 <- par[5]
    lambda <- par[6]
    lambda_2 <- par[7]
  }
  
  y <- data$y
  n = length(y)                            # number of observations
  sigma <- rep(NA, n)                      # conditional volatility
  eps <- rep(NA, n)                      # conditional volatility
  if(family != 'GARCH-M')
  {
    is_lower <- rep(NA, n)
  }
  ### specification of unconditional variance ###
  
  if (family == "GARCH-M")
  {
    var_uncond <- omega / (1 - alpha - beta)
    
    ksi <- (omega ^ 2 + 2 * omega * alpha * var_uncond + 2 * omega * beta * var_uncond) / 
      (1 - 3 * alpha ^ 2 - beta ^ 2 - 2 * alpha * beta)
    
    sigma_0 <- lambda ^ 2 * ((alpha ^ 2 * (3 * ksi - var_uncond ^ 2) + 2 * alpha * beta * (ksi - var_uncond ^ 2)) / 
                               (1 - beta ^ 2)) + var_uncond
  }
  
  if (family == "GARCH-M-GJR")
  {
    var_uncond <- omega / (1 - alpha - beta - gamma_2 / 2)
    
    exp_sigma_4 <- (omega ^ 2 + omega * var_uncond * (2 * alpha + 2 * beta + gamma_2)) / 
      (1 - 3 * alpha ^ 2 - beta ^ 2 - 1.5 * gamma_2 ^ 2 - 2 * alpha * beta - 
         3 * alpha * gamma_2 - beta * gamma_2)
    
    var_sigma_2 <- (alpha ^ 2 * (3 * exp_sigma_4 - var_uncond^2) + 2 * alpha * beta * (exp_sigma_4 - var_uncond^2) + 
                      gamma_2 ^ 2 * (1.5 * exp_sigma_4 - 0.25 * var_uncond^2) + 
                      2 * (alpha * gamma_2 * (1.5 * exp_sigma_4 - var_uncond ^ 2 / 2) + 
                             beta * gamma_2 * 0.5 * (exp_sigma_4 - var_uncond ^ 2))) / (1 - beta ^ 2)
    
    sigma_0 <- lambda ^ 2 * var_sigma_2 + var_uncond
  }
  
  if (family == "GARCH-M-GJR-LEV")
  {
    var_uncond <- omega / (1 - alpha - beta - gamma_2 / 2)
    exp_sigma_4 <- (omega ^ 2 + omega * var_uncond * (2 * alpha + 2 * beta + gamma_2)) / 
      (1 - 3 * alpha ^ 2 - beta ^ 2 - 1.5 * gamma_2 ^ 2 - 2 * alpha * beta - 
         3 * alpha * gamma_2 - beta * gamma_2)
    
    var_sigma_2 <- (alpha ^ 2 * (3 * exp_sigma_4 - var_uncond^2) + 2 * alpha * beta * (exp_sigma_4 - var_uncond^2) + 
                      gamma_2 ^ 2 * (1.5 * exp_sigma_4 - 0.25 * var_uncond^2) + 
                      2 * (alpha * gamma_2 * (1.5 * exp_sigma_4 - var_uncond ^ 2 / 2) + 
                             beta * gamma_2 * 0.5 * (exp_sigma_4 - var_uncond ^ 2))) / (1 - beta ^ 2)
    
    sigma_0 <- lambda ^ 2 * var_sigma_2 + lambda_2 ^ 2 * (0.5 * exp_sigma_4 - 0.25 * var_uncond^2) + 
      lambda * lambda_2 * (exp_sigma_4 - var_uncond^2) + var_uncond
  }
  if(family == "GARCH-M-LEV")
  {
    var_uncond <- omega / (1 - alpha - beta)
    exp_sigma_4 <- (omega ^ 2 + 2 * omega * var_uncond * (alpha + beta)) / 
      (1 - 3 * alpha ^ 2 - beta ^ 2 - 2 * alpha * beta)
    
    var_sigma_2 <- (alpha ^ 2 * (3 * exp_sigma_4 - var_uncond^2) + 2 * alpha * beta * (exp_sigma_4 - var_uncond^2)) / (1 - beta ^ 2)
    
    sigma_0 <- lambda ^ 2 * var_sigma_2 + lambda_2 ^ 2 * (0.5 * exp_sigma_4 - 0.25 * var_uncond^2) + 
      lambda * lambda_2 * (exp_sigma_4 - var_uncond^2) + var_uncond
  }
  
  
  sigma[1] <- sqrt(sigma_0)
  eps[1] <- y[1] - mu
  
  if (family != "GARCH-M")
  {
    is_lower[1] <- eps[1] < 0
  }
  
  # Validation
  if (is.nan(sigma_0))
  {
    return(-9999999999999)
  }
  
  if (sigma_0 < 0)
  {
    return(-9999999999999)
  }
  
  if ((alpha < 0) | (beta < 0) |((alpha + beta) > 1) | (omega <= 0) | (mu <= 0))
  {
    return(-9999999999999)
  }
  
  # Recursive calculations
  for(t in 2:n)
  {
    if (family == "GARCH-M")
    {
      sigma[t] <- sqrt(omega +
                         alpha * eps[t - 1] ^ 2 + 
                         beta * sigma[t - 1] ^ 2)
      eps[t] <- y[t] - mu - lambda * sigma[t-1] ^ 2
      
    }
    if (family == "GARCH-M-GJR")
    {
      sigma[t] <- sqrt(omega +
                         alpha * eps[t - 1] ^ 2 + 
                         beta * sigma[t - 1] ^ 2 +
                         gamma_2 * is_lower[t - 1] * eps[t - 1] ^ 2)
      eps[t] <- y[t] - mu - lambda * sigma[t-1] ^ 2
      is_lower[t] <- as.numeric(eps[t] < 0)
    }
    
    if (family == "GARCH-M-GJR-LEV")
    {
      sigma[t] <- sqrt(omega +
                         alpha * eps[t - 1] ^ 2 + 
                         beta * sigma[t - 1] ^ 2 +
                         gamma_2 * is_lower[t - 1] * eps[t - 1] ^ 2)
      eps[t] <- y[t] - mu - lambda * sigma[t-1] ^ 2 - 
        lambda_2 * is_lower[t-1] * sigma[t-1] ^ 2
      is_lower[t] <- as.numeric(eps[t] < 0)
    }
  }

    sigma_forc_2 <- rep(NA, steps_ahead)                                                       # volatility forecasts
    y_forc <- rep(NA, steps_ahead)                                                             # return forecasts
    
    if (family == "GARCH-M")
    {
      y_forc[1] <- mu + lambda * sigma[length(sigma)] ^ 2                                        # first step return forecast
      sigma_forc_2[1] <- omega + alpha * eps[length(eps)] ^ 2 + beta * sigma[length(sigma)] ^ 2  # first step variance forecast
      
      if (steps_ahead > 1)
      {
        for(k in 2:steps_ahead)
        {
          sigma_forc_2[k] <- omega + alpha * sigma_forc_2[k-1] + beta * sigma_forc_2[k-1]
          y_forc[k] <- mu + lambda * sigma_forc_2[k-1]
        }
      }
    }
    
    if (family == "GARCH-M-GJR")
    {
      y_forc[1] <- mu + lambda * sigma[length(sigma)] ^ 2                                        # first step return forecast
      sigma_forc_2[1] <- omega + alpha * eps[length(eps)] ^ 2 + beta * sigma[length(sigma)] ^ 2 +
                         gamma_2 * is_lower[length(is_lower)] * eps[length(eps)] ^ 2             # first step variance forecast
      if (steps_ahead > 1)
      {
        for(k in 2:steps_ahead)
        {
        sigma_forc_2[k] <- omega + alpha * sigma_forc_2[k-1] + beta * sigma_forc_2[k-1] + 0.5 * gamma_2 * sigma_forc_2[k-1]
        y_forc[k] <- mu + lambda * sigma_forc_2[k-1]
        }
      }
    }
    
    if (family == "GARCH-M-GJR-LEV")
    {
      y_forc[1] <- mu + lambda * sigma[length(sigma)] ^ 2 +
                        lambda_2 * is_lower[length(is_lower)] * sigma[length(sigma)] ^ 2         # first step return forecast
      sigma_forc_2[1] <- omega + alpha * eps[length(eps)] ^ 2 + beta * sigma[length(sigma)] ^ 2 +
                         gamma_2 * is_lower[length(is_lower)] * eps[length(eps)] ^ 2             # first step variance forecast
      if (steps_ahead > 1)
      {
        for(k in 2:steps_ahead)
        {
          sigma_forc_2[k] <- omega + alpha * sigma_forc_2[k-1] + beta * sigma_forc_2[k-1] + 0.5 * gamma_2 * sigma_forc_2[k-1]
          y_forc[k] <- mu + lambda * sigma_forc_2[k-1] + 0.5 * lambda_2 * sigma_forc_2[k-1]
        }
      }
    }
      
    
    
    
  
    return_list <- list("gamma" = par,                             # coefficient estimates
                        'sigma_forc' = sqrt(sigma_forc_2),         # volatility forecasts
                        'return_forc' = y_forc,                    # return forecasts
                        'sigma' = sigma,
                        'return' = y
                       )
    return(return_list)

}
                                                                                                            
                                                                                                            
                                                                                                            
                                                                  
