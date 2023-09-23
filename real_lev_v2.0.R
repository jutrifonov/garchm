################################################
### This R file provides real data analysis ####
###        for the SP 500 market index      ####
###               (2004-2021)               ####
################################################
#####       The code is developed by     #######
#####      @ Juri Trifonov, NRU HSE      #######
#####      @ Bogdan Potanin, NRU HSE     #######
################################################

########################################
#          Libraries required
########################################
library('gena')
library('rugarch')
library("writexl")
library('readxl')
library("numDeriv")

######################################################################
######################################################################
#                              !!!NB!!!
#         PLEASE, PROVIDE THE PATH TO THE APPROPRIATE SOURCE
#                  FILE WITH FUNCTIONS ON YOUR PC
#                    ↓↓↓↓↓↓↓↓↓↓↓HERE↓↓↓↓↓↓↓↓↓↓↓
source("/Users/jurytrifonov/Desktop/Research/RP-GARCH/GARCH-M-GJR-LEV/func_lev.R")
#######################################################################
#######################################################################

# Disable scientific notation
options(scipen = 999)

####### Upload the xlsx file with data #######
#% Please provide an appropriate 
#% path to the file on your PC
SP500 <- read_excel("Desktop/Research/SP500.xlsx")

####### Please, define a ONE specific period on which you  ########
#######   want to estimate and compare all of the models   ########
sp <- data.frame('y' = SP500$y[3776:4531] * 100)  # years 2019-2021 
sp <- data.frame('y' = SP500$y[3525:4280] * 100)  # years 2018-2020
sp <-  data.frame('y' = SP500$y[3274:4027] * 100) # years 2017-2019

sp <-  data.frame('y' = SP500$y[3022:3775] * 100) # years 2016-2018
sp <-  data.frame('y' = SP500$y[2770:3524] * 100) # years 2015-2017 
sp <-  data.frame('y' = SP500$y[2518:3273] * 100) # years 2014-2016
sp <-  data.frame('y' = SP500$y[2266:3021] * 100) # years 2013-2015

sp <- data.frame('y' = SP500$y * 100)             # whole sample analysis

### 1. GARCH-M-GJR-LEV
get1 <- GARCH(sp, c(get3$gamma, 0), family = "GARCH-M-GJR-LEV", ga_iter = 100)

### 2. GARCH-M-GJR
get3 <- GARCH(sp, c(coef(model_rugarch1), 0, 0), family = "GARCH-M-GJR", ga_iter = 100)

### 3. GARCH-M
get2 <- GARCH(train, c(coef(model_rugarch1), 0), family = "GARCH-M", ga_iter = 100)

### 4. GARCH
model1 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(0, 0))) 
model_rugarch1 <- ugarchfit(spec = model1, data = sp$y)
coef(model_rugarch1)

### Store the estimated parameters for every model
gamma_est_garch_m_gjr_lev <- c('mu' = get1$gamma[1], 'omega' = get1$gamma[2], 'alpha' = get1$gamma[3], 
                              'beta' = get1$gamma[4], 'gamma_2' = get1$gamma[5],
                              'lambda' = get1$gamma[6],'lambda_2' = get1$gamma[7])

gamma_est_garch_m_gjr <- c('mu' = get3$gamma[1], 'omega' = get3$gamma[2], 'alpha' = get3$gamma[3], 
                           'beta' = get3$gamma[4], 'gamma_2' = get3$gamma[5],
                           'lambda' = get3$gamma[6], 'lambda_2' = 0)

gamma_est_garch_m <- c('mu' = get2$gamma[1], 'omega' = get2$gamma[2], 'alpha' = get2$gamma[3], 
                       'beta' = get2$gamma[4], 'gamma_2' = 0,
                       'lambda' = get2$gamma[5], 'lambda_2' = 0)

rbind(gamma_est_garch_m,
      gamma_est_garch_m_gjr,
      gamma_est_garch_m_gjr_lev)

### Calculate information criteria
AIC_garch <- 2 * 4 - 2 * likelihood(model_rugarch1)
AIC_garch_m <- 2 * 5 - 2 * get2$lnL
AIC_garch_m_gjr <- 2 * 6 - 2 * get3$lnL
AIC_garch_m_gjr_lev <- 2 * 7 - 2 * get1$lnL
cbind('AIC_garch' = AIC_garch, 'AIC_garch_m' = AIC_garch_m, 
      'AIC_garch_m_gjr' = AIC_garch_m_gjr, 'AIC_garch_gjr_lev' = AIC_garch_m_gjr_lev)
BIC_garch_m <- bic_func(get2)
BIC_garch_m_gjr <- bic_func(get3)
BIC_garch_m_gjr_lev <- bic_func(get1)

#############################################
### Covariance matrix estimation with GOP ###
#############################################

### 1. GARCH-M-GJR-LEV
gamma_test <- get1$gamma
J <- jacobian(func = lnL,                             
              x = gamma_test,                         
              data = sp, is_aggregate = FALSE, family = 'GARCH-M-GJR-LEV')    

cov_MLE_gop <- solve(t(J) %*% J)    

# Print final result
summary(get1, cov_MLE_gop, family = "GARCH-M-GJR-LEV")

### 2. GARCH-M-GJR
gamma_test_3 <- get3$gamma
J <- jacobian(func = lnL,                             
              x = gamma_test_3,                         
              data = sp, is_aggregate = FALSE, family = 'GARCH-M-GJR')        
cov_MLE_gop_m_gjr <- solve(t(J) %*% J) 

# Print final result
summary(get3, cov_MLE_gop_m_gjr, family = "GARCH-M-GJR")

### 3. GARCH-M
gamma_test_2 <- get2$gamma
J <- jacobian(func = lnL,                             
              x = gamma_test_2,                         
              data = sp, is_aggregate = FALSE, family = 'GARCH-M')        
cov_MLE_gop_m <- solve(t(J) %*% J) 

# Print final result
summary(get2, cov_MLE_gop_m, family = "GARCH-M")






