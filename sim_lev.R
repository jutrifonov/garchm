############################################
### This R file provides simulated data ####
### analysis for the different GARCH    ####
###           family models             ####
##########################################
####      The code is developed by     ###
####     @ Juri Trifonov, NRU HSE      ###
####                                   ###
##########################################

########################################
#          Libraries required
########################################
library('devtools')
library('gena')
library('rugarch')
library("writexl")
library('readxl')

######################################################################
######################################################################
#                              !!!NB!!!
#         PLEASE, PROVIDE THE PATH TO THE APPROPRIATE SOURCE
#                  FILE WITH FUNCTIONS ON YOUR PC
#                    ↓↓↓↓↓↓↓↓↓↓↓HERE↓↓↓↓↓↓↓↓↓↓↓
source("/Users/yutrifonov/Desktop/Research/RP-GARCH/GARCH-M-GJR-LEV/func_lev.R")
#######################################################################
#######################################################################

# Disable scientific notation
options(scipen = 999)

################################
### 1. Simulation parameters ###
################################

################################
############ SET I #############
################################
n <- 1000                      # sample size
m <- 100                       # number of simulations
ga_maxiter <- 100              # number of GA interations
mu <- 0.01                     # true parameters
omega <- 0.1                   # true parameters
alpha <- 0.1                   # true parameters
beta <- 0.7                    # true parameters
gamma_2 <- 0.15                # true parameters
lambda <- 0.2                  # true parameters
lambda_2 <- 0.5                # true parameters

################################
########### SET II  ############
################################
#% Please, uncomment this part
#% if you would like to replicate
#% simulations for the Set II

#n <- 1000                      # sample size
#m <- 100                       # number of simulations
#ga_maxiter <- 100              # number of GA iterations
#mu <- 0.06                     # true parameters
#omega <- 0.05                  # true parameters
#alpha <- 0.06                  # true parameters
#beta <- 0.7                    # true parameters
#gamma_2 <- 0.23                # true parameters
#lambda <- -0.07                # true parameters
#lambda_2 <- 0.18               # true parameters


#########################################################
### 2. Create matrices/vectors for storing parameters ###
#########################################################
{
gamma_garch <- matrix(NA, nrow = m, ncol = 4)
gamma_garch_m <- matrix(NA, nrow = m, ncol = 5)
gamma_garch_m_gjr <- matrix(NA, nrow = m, ncol = 6)
gamma_garch_m_gjr_lev <- matrix(NA, nrow = m, ncol = 7)

# vectors for storing information criterias
aic_garch <- rep(NA, m)
aic_garch_m <- rep(NA, m)
aic_garch_m_gjr <- rep(NA, m)
aic_garch_m_gjr_lev <- rep(NA, m)

bic_garch <- rep(NA, m)
bic_garch_m <- rep(NA, m)
bic_garch_m_gjr <- rep(NA, m)
bic_garch_m_gjr_lev <- rep(NA, m)

# vectors for storing accuracy metrics
# volatility accuracy metrics
rmse_sigma_garch <- rep(NA, m)
rmse_sigma_garch_m <- rep(NA, m)
rmse_sigma_garch_m_gjr <- rep(NA, m)
rmse_sigma_garch_m_gjr_lev <- rep(NA, m)

mae_sigma_garch <- rep(NA, m)
mae_sigma_garch_m <- rep(NA, m)
mae_sigma_garch_m_gjr <- rep(NA, m)
mae_sigma_garch_m_gjr_lev <- rep(NA, m)

mape_sigma_garch <- rep(NA, m)
mape_sigma_garch_m <- rep(NA, m)
mape_sigma_garch_m_gjr <- rep(NA, m)
mape_sigma_garch_m_gjr_lev <- rep(NA, m)
#returns accuracy metrics
rmse_ret_garch <- rep(NA, m)
rmse_ret_garch_m <- rep(NA, m)
rmse_ret_garch_m_gjr <- rep(NA, m)
rmse_ret_garch_m_gjr_lev <- rep(NA, m)

mae_ret_garch <- rep(NA, m)
mae_ret_garch_m <- rep(NA, m)
mae_ret_garch_m_gjr <- rep(NA, m)
mae_ret_garch_m_gjr_lev <- rep(NA, m)

mape_ret_garch <- rep(NA, m)
mape_ret_garch_m <- rep(NA, m)
mape_ret_garch_m_gjr <- rep(NA, m)
mape_ret_garch_m_gjr_lev <- rep(NA, m)
}

######################
### 3. Simulations ###
######################
for (i in 1:m)
{
  print("##############################")
  print(paste0('Simulation number:', i))
  print("##############################")
set.seed(i * 100)
# true data generating process
df <- GARCH.Simulate(n = n,
                     mu = mu, omega = omega, alpha = alpha, beta = beta,
                     lambda = lambda,
                     gamma_2 = gamma_2,
                     lambda_2 = lambda_2,
                     family = 'GARCH-M-GJR-LEV')

##############################
### simple GARCH (rugarch) ###
##############################
model <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(0, 0), archm = FALSE),
                    start.pars = c(mu, omega, alpha, beta)) 
model_rugarch <- ugarchfit(spec = model,
                           data = df$y)
gamma_garch[i, ] <- coef(model_rugarch)

AIC_rugarch <- 2 * 4 - 2 * likelihood(model_rugarch)
BIC_rugarch <- 4 * log(n) - 2 * likelihood(model_rugarch)
aic_garch[i] <- AIC_rugarch
bic_garch[i] <- BIC_rugarch

rmse_sigma_rugarch <- rmse_func(df$sigma, sigma(model_rugarch))
rmse_ret_rugarch <- rmse_func(df$y, fitted(model_rugarch))
rmse_sigma_garch[i] <- rmse_sigma_rugarch
rmse_ret_garch[i] <- rmse_ret_rugarch

mae_sigma_rugarch <- mae_func(df$sigma, sigma(model_rugarch))
mae_ret_rugarch <- mae_func(df$y, fitted(model_rugarch))
mae_sigma_garch[i] <- mae_sigma_rugarch
mae_ret_garch[i] <- mae_ret_rugarch

mape_sigma_rugarch <- mape_func(df$sigma, sigma(model_rugarch))
mape_ret_rugarch <- mape_func(df$y, fitted(model_rugarch))
mape_sigma_garch[i] <- mape_sigma_rugarch
mape_ret_garch[i] <- mape_ret_rugarch
##############################

##############################
########### GARCH-M ##########
##############################
model_m <- GARCH(df, c(coef(model_rugarch), 0), family = 'GARCH-M', ga_iter = ga_maxiter)
if ((aic_func(model_m) > 10000) | (rmse_func(df$sigma, model_m$sigma) > 1000))
{
  print("OPTIMIZATION TROUBLE")
  attempt_count <- 0
  while((aic_func(model_m) > 10000) | (rmse_func(df$sigma, model_m$sigma) > 1000))
  {
    attempt_count <- attempt_count + 1
    tryCatch(
      {
        model_m <- GARCH(df, c(coef(model_rugarch), 0.1) * runif(1, 0, 1), family = 'GARCH-M', ga_iter = ga_maxiter)
      },
      error = function(cond)
      {
        message("initial point error")
      }
    )
    if (attempt_count > 10)
    {
      message("too much attempts")
      break
    }
  }
}
gamma_garch_m[i, ] <- model_m$gamma

aic_garch_m[i] <- aic_func(model_m)
bic_garch_m[i] <- bic_func(model_m)

rmse_sigma_garch_m[i] <- rmse_func(df$sigma, model_m$sigma)
rmse_ret_garch_m[i] <- rmse_func(df$y, model_m$return)

mae_sigma_garch_m[i] <- mae_func(df$sigma, model_m$sigma)
mae_ret_garch_m[i] <- mae_func(df$y, model_m$return)

mape_sigma_garch_m[i] <- mape_func(df$sigma, model_m$sigma)
mape_ret_garch_m[i] <- mape_func(df$y, model_m$return)
##############################

##############################
######### GARCH-M-GJR ########
##############################
model_m_gjr <- GARCH(df, c(model_m$gamma, 0), family = 'GARCH-M-GJR', ga_iter = ga_maxiter)
if ((aic_func(model_m_gjr) > 10000) | (rmse_func(df$sigma, model_m_gjr$sigma) > 1000))
{
  print("OPTIMIZATION TROUBLE")
  while((aic_func(model_m_gjr) > 10000) | (rmse_func(df$sigma, model_m_gjr$sigma) > 1000))
  {
    model_m_gjr <- GARCH(df, c(model_m$gamma, 0.01) * runif(1, 0, 1), family = 'GARCH-M-GJR', ga_iter = ga_maxiter)
  }
}
gamma_garch_m_gjr[i, ] <- model_m_gjr$gamma

aic_garch_m_gjr[i] <- aic_func(model_m_gjr)
bic_garch_m_gjr[i] <- bic_func(model_m_gjr)

rmse_sigma_garch_m_gjr[i] <- rmse_func(df$sigma, model_m_gjr$sigma)
rmse_ret_garch_m_gjr[i] <- rmse_func(df$y, model_m_gjr$return)

mae_sigma_garch_m_gjr[i] <- mae_func(df$sigma, model_m_gjr$sigma)
mae_ret_garch_m_gjr[i] <- mae_func(df$y, model_m_gjr$return)

mape_sigma_garch_m_gjr[i] <- mape_func(df$sigma, model_m_gjr$sigma)
mape_ret_garch_m_gjr[i] <- mape_func(df$y, model_m_gjr$return)
##############################

##############################
####### GARCH-M-GJR-LEV ######
##############################
model_m_gjr_lev <- GARCH(df, c(model_m_gjr$gamma, 0), family = 'GARCH-M-GJR-LEV', ga_iter = ga_maxiter)
if ((aic_func(model_m_gjr_lev) > 10000) | (rmse_func(df$sigma, model_m_gjr_lev$sigma) > 1000))
{
  print("OPTIMIZATION TROUBLE")
  while((aic_func(model_m_gjr_lev) > 10000) | (rmse_func(df$sigma, model_m_gjr_lev$sigma) > 1000))
  {
    model_m_gjr_lev <- GARCH(df, c(model_m_gjr$gamma, 0.01) * runif(1, 0, 1), family = 'GARCH-M-GJR-LEV', ga_iter = ga_maxiter)
  }
}
gamma_garch_m_gjr_lev[i, ] <- model_m_gjr_lev$gamma

aic_garch_m_gjr_lev[i] <- aic_func(model_m_gjr_lev)
bic_garch_m_gjr_lev[i] <- bic_func(model_m_gjr_lev)

rmse_sigma_garch_m_gjr_lev[i] <- rmse_func(df$sigma, model_m_gjr_lev$sigma)
rmse_ret_garch_m_gjr_lev[i] <- rmse_func(df$y, model_m_gjr_lev$return)

mae_sigma_garch_m_gjr_lev[i] <- mae_func(df$sigma, model_m_gjr_lev$sigma)
mae_ret_garch_m_gjr_lev[i] <- mae_func(df$y, model_m_gjr_lev$return)

mape_sigma_garch_m_gjr_lev[i] <- mape_func(df$sigma, model_m_gjr_lev$sigma)
mape_ret_garch_m_gjr_lev[i] <- mape_func(df$y, model_m_gjr_lev$return)
}

###############################
### 4. Final table creation ###
############################### 
{
# give names to estimated parameters
colnames(gamma_garch) <- c("rugarch.mu", "rugarch.omega", 
                           "rugarch.alpha", "rugarch.beta")
colnames(gamma_garch_m) <- c("garch_m.mu", "garch_m.omega", 
                             "garch_m.alpha", "garch_m.beta", 
                             "garch_m.lambda")
colnames(gamma_garch_m_gjr) <- c("garch_m_gjr.mu", "garch_m_gjr.omega", 
                                 "garch_m_gjr.alpha", "garch_m_gjr.beta", 
                                 "garch_m_gjr.lambda", "garch_m_gjr.gamma")
colnames(gamma_garch_m_gjr_lev) <- c("garch_m_gjr_lev.mu", "garch_m_gjr_lev.omega",
                                     "garch_m_gjr_lev.alpha", "garch_m_gjr_lev.beta",
                                     "garch_m_gjr_lev.lambda", "garch_m_gjr_lev.gamma",
                                     "garch_m_gjr_lev.lambda_2")
# data frames of all estimates and metrics
true_param <- data.frame(mu, omega, alpha, beta, lambda, gamma_2, lambda_2)
colnames(true_param) <- c("mu", "omega", "alpha", 
                          "beta", "lambda","gamma","lambda_2")
gamma_df <- data.frame(gamma_garch, gamma_garch_m, 
                       gamma_garch_m_gjr, gamma_garch_m_gjr_lev)
aic_df <- data.frame(aic_garch, aic_garch_m, 
                     aic_garch_m_gjr, aic_garch_m_gjr_lev)
bic_df <- data.frame(bic_garch, bic_garch_m, 
                     bic_garch_m_gjr, bic_garch_m_gjr_lev)
rmse_sigma_df <- data.frame(rmse_sigma_garch, rmse_sigma_garch_m,
                            rmse_sigma_garch_m_gjr, rmse_sigma_garch_m_gjr_lev)
mae_sigma_df <- data.frame(mae_sigma_garch, mae_sigma_garch_m,
                           mae_sigma_garch_m_gjr, mae_sigma_garch_m_gjr_lev)
mape_sigma_df <- data.frame(mape_sigma_garch, mape_sigma_garch_m,
                            mape_sigma_garch_m_gjr, mape_sigma_garch_m_gjr_lev)
rmse_ret_df <- data.frame(rmse_ret_garch, rmse_ret_garch_m,
                          rmse_ret_garch_m_gjr, rmse_ret_garch_m_gjr_lev)
mae_ret_df <- data.frame(mae_ret_garch, mae_ret_garch_m,
                         mae_ret_garch_m_gjr, mae_ret_garch_m_gjr_lev)
mape_ret_df <- data.frame(mape_ret_garch, mape_ret_garch_m,
                          mape_ret_garch_m_gjr, mape_ret_garch_m_gjr_lev)

# Final Table
final_sim <- data.frame(true_param, gamma_df, aic_df, bic_df, rmse_sigma_df,
                        mae_sigma_df, mape_sigma_df, rmse_ret_df,
                        mae_ret_df, mape_ret_df)
}
print(final_sim)

# Creation of the excel file with results
#% Please provide a path on your PC 
#% where you want to save the simulation results
write_xlsx(final_sim, '/Users/jurytrifonov/Desktop/sim_results_30.xlsx')

###############################
####### 5. Aggregation ########
### and additional metrics ####
###############################

#% Please provide a path on your PC 
#% to the file with simulations that you
#% have saved earlier on line 302

sim_results <- data.frame(read_excel('/Users/yutrifonov/Desktop/Research/RP-GARCH/Sim_setup_1/sim_results_1.xlsx'))

{
true_full <- sim_results[1:7]
garch_full<- sim_results[8:11]
garch_m_full <- sim_results[12:16]
garch_m_gjr_full <- sim_results[17:22]
garch_m_gjr_lev_full <- sim_results[23:29]
full_table <- as.data.frame(cbind(true_full, garch_full, garch_m_full, 
                                  garch_m_gjr_full, garch_m_gjr_lev_full))

# Calculate mean values of estimates
mean_value <- colMeans(sim_results)
true <- mean_value[1:7]
garch <- c(mean_value[8:11], rep(0, 3))
garch_m <- c(mean_value[12:16], rep(0, 2))
garch_m_gjr <- c(mean_value[17:22], rep(0, 1))
garch_m_gjr_lev <- mean_value[23:29]

# table with mean estimates of coefficients
estimates_table <- as.data.frame(cbind(true, garch, garch_m, garch_m_gjr, garch_m_gjr_lev))

# Calculate rmse/mae/mape for coefficient estimates
{
rmse_garch <- rep(NA, 4)
mae_garch = rmse_garch
mape_garch = rmse_garch
rmse_garch_m <- rep(NA, 5)
mae_garch_m = rmse_garch_m
mape_garch_m = rmse_garch_m
rmse_garch_m_gjr <- rep(NA, 6)
mae_garch_m_gjr = rmse_garch_m_gjr
mape_garch_m_gjr = rmse_garch_m_gjr
rmse_garch_m_gjr_lev <- rep(NA, 7)
mae_garch_m_gjr_lev = rmse_garch_m_gjr_lev
mape_garch_m_gjr_lev = rmse_garch_m_gjr_lev


j = 8
for (i in 1:4)
{
  
  rmse_garch[i] <- rmse_func(full_table[, i], full_table[, j])
  mae_garch[i] <- mae_func(full_table[, i], full_table[, j])
  mape_garch[i] <- mape_func(full_table[, i], full_table[, j]) 
  j = j + 1
}
j = 12
for(i in 1:5)
{
  rmse_garch_m[i] <- rmse_func(full_table[, i], full_table[, j])
  mae_garch_m[i] <- mae_func(full_table[, i], full_table[, j])
  mape_garch_m[i] <- mape_func(full_table[, i], full_table[, j]) 
  j = j + 1
}
j = 17
for(i in 1:6)
{
  rmse_garch_m_gjr[i] <- rmse_func(full_table[, i], full_table[, j])
  mae_garch_m_gjr[i] <- mae_func(full_table[, i], full_table[, j])
  mape_garch_m_gjr[i] <- mape_func(full_table[, i], full_table[, j]) 
  j = j + 1
}
j = 23
for(i in 1:7)
{
  rmse_garch_m_gjr_lev[i] <- rmse_func(full_table[, i], full_table[, j])
  mae_garch_m_gjr_lev[i] <- mae_func(full_table[, i], full_table[, j])
  mape_garch_m_gjr_lev[i] <- mape_func(full_table[, i], full_table[, j]) 
  j = j + 1
}

accuracy_table <- t(data.frame(c(rmse_garch, rep(0, 3)), c(rmse_garch_m, rep(0, 2)), 
                               c(rmse_garch_m_gjr, rep(0,1)), rmse_garch_m_gjr_lev,
                               c(mae_garch, rep(0, 3)), c(mae_garch_m, rep(0, 2)),
                               c(mae_garch_m_gjr, rep(0,1)), mae_garch_m_gjr_lev,
                               c(mape_garch, rep(0, 3)), c(mape_garch_m, rep(0, 2)),
                               c(mape_garch_m_gjr, rep(0, 1)), mape_garch_m_gjr_lev
                               ))
colnames(accuracy_table) <- c('mu', 'omega', 'alpha', 'beta', 'lambda', 'gamma', 'lambda_2')
rownames(accuracy_table) <- c('garch_rmse', 'garch_m_rmse', 'garch_m_gjr_rmse', 'garch_m_gjr_lev_rmse',
                              'garch_mae', 'garch_m_mae', 'garch_m_gjr_mae', 'garch_m_gjr_lev_mae',
                              'garch_mape', 'garch_m_mape', 'garch_m_gjr_mape', 'garch_m_gjr_lev_mape')

# table with accuracy metrics for coefficient estimates
accuracy_table
}


aics <- rep(NA, 4)
bics <- rep(NA, 4)
j = 30
k = 34
for (i in 1:4)
{
  
  aics[i] <- mean(sim_results[, j])
  bics[i]<- mean(sim_results[, k])
  j = j + 1
  k = k + 1
}

rmse_sigma <- rep(NA, 4)
mae_sigma <- rep(NA, 4)
mape_sigma <- rep(NA, 4)
rmse_ret <- rep(NA, 4)
mae_ret <- rep(NA, 4)
mape_ret <- rep(NA, 4)
j = 38
k = 42
l = 46
for (i in 1:4)
{
  rmse_sigma[i] <- mean(sim_results[, j])
  mae_sigma[i] <- mean(sim_results[, k])
  mape_sigma[i] <- mean(sim_results[, l])
  j = j + 1
  k = k + 1
  l = l + 1
}
j = 50
k = 54
l = 58
for(i in 1:4)
{
  rmse_ret[i] <- mean(sim_results[, j])
  mae_ret[i] <- mean(sim_results[, k])
  mape_ret[i] <- mean(sim_results[, l])
  j = j + 1
  k = k + 1
  l = l + 1
}
mean(sim_results$mape_ret_garch_m_gjr_lev)
vol_ret_accuracy <- t(data.frame(rmse_sigma, mae_sigma, mape_sigma, rmse_ret, mae_ret, mape_ret))
colnames(vol_ret_accuracy) <- c('garch', 'garch_m', 'garch_m_gjr', 'garch_m_gjr_lev')
vol_ret_accuracy

# table with AICs, BICs ans accuracy metrics on volatility and returns
total_accuracy <- rbind(aics, bics, vol_ret_accuracy)

### FINAL SUMMARY TABLES ###

# table with mean estimates of coefficients
print(estimates_table)
# table with accuracy metrics for coefficient estimates
print(accuracy_table)
# table with AICs, BICs ans accuracy metrics on volatility and returns
print(total_accuracy)
}

{
rmse_sigma_table <- sim_results[, 38:41]
boolean_garch <- rep(NA, length(sim_results[, 1]))
boolean_garch_m = boolean_garch
boolean_garch_m_gjr = boolean_garch
boolean_garch_m_gjr_lev = boolean_garch

for (i in 1: length(sim_results[, 1]))
{
 boolean_garch[i] <- (rmse_sigma_table[i,1]== min(rmse_sigma_table[i,]))
 boolean_garch_m[i] <- (rmse_sigma_table[i,2]== min(rmse_sigma_table[i,]))
 boolean_garch_m_gjr[i] <- (rmse_sigma_table[i,3]== min(rmse_sigma_table[i,]))
 boolean_garch_m_gjr_lev[i] <- (rmse_sigma_table[i,4]== min(rmse_sigma_table[i,]))
}
victory_perc_sigma_garch <- mean(boolean_garch)
victory_perc_sigma_garch_m <- mean(boolean_garch_m)
victory_perc_sigma_garch_m_gjr <- mean(boolean_garch_m_gjr)
victory_perc_sigma_garch_m_gjr_lev <- mean(boolean_garch_m_gjr_lev)
victory_sigma_table <- data.frame(victory_perc_sigma_garch, victory_perc_sigma_garch_m,
                                  victory_perc_sigma_garch_m_gjr, victory_perc_sigma_garch_m_gjr_lev)
colnames(victory_sigma_table) <- c('garch', 'garch_m', 'garch_m_gjr', 'garch_m_gjr_lev')
victory_sigma_table

rmse_ret_table <- sim_results[, 50:53]
boolean_garch <- rep(NA, length(sim_results[, 1]))
boolean_garch_m = boolean_garch
boolean_garch_m_gjr = boolean_garch
boolean_garch_m_gjr_lev = boolean_garch

for (i in 1: length(sim_results[, 1]))
{
  boolean_garch[i] <- (rmse_ret_table[i,1]== min(rmse_ret_table[i,]))
  boolean_garch_m[i] <- (rmse_ret_table[i,2]== min(rmse_ret_table[i,]))
  boolean_garch_m_gjr[i] <- (rmse_ret_table[i,3]== min(rmse_ret_table[i,]))
  boolean_garch_m_gjr_lev[i] <- (rmse_ret_table[i,4]== min(rmse_ret_table[i,]))
}
rmse_ret_table[1, ]
victory_perc_ret_garch <- mean(boolean_garch)
victory_perc_ret_garch_m <- mean(boolean_garch_m)
victory_perc_ret_garch_m_gjr <- mean(boolean_garch_m_gjr)
victory_perc_ret_garch_m_gjr_lev <- mean(boolean_garch_m_gjr_lev)
victory_ret_table <- data.frame(victory_perc_ret_garch, victory_perc_ret_garch_m,
                                  victory_perc_ret_garch_m_gjr, victory_perc_ret_garch_m_gjr_lev)
colnames(victory_ret_table) <- c('garch', 'garch_m', 'garch_m_gjr', 'garch_m_gjr_lev')
print(victory_sigma_table)
print(victory_ret_table)
}

