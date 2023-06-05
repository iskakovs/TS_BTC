# Packages
install.packages('ARDL')
install.packages('tidyverse')
install.packages('fpp2')
install.packages('fpp3')
install.packages('urca')
install.packages('rio')
install.packages('patchwork')
install.packages('rvest')
install.packages('lmtest')
install.packages('tseries')
install.packages('ggplot2')
install.packages('vars')
install.packages('dynamac')
install.packages('tseries')
install.packages('TSstudio')
install.packages('tsDyn')

# Libraries
library(tidyverse)
library(fpp2)
library(fpp3)
library(urca)
library(rio)
library(patchwork)
library(rvest)
library(ARDL)
library(lmtest)
library(tseries)
library(ggplot2)
library(vars)
library(dynamac)
library(tseries)
library(TSstudio)
library(tsDyn)
library(changepoint)

####################

# Set working directory
setwd("C:/Users/adimn/Desktop/Data")

data = import('final_dataset.xlsx')
data

# Change dates to weeks in timeseries
data = mutate(data, date = yearweek(ymd('2015-01-19') +
                                  weeks(0:258))) 
# Declare the dataset is timeseries
data = as_tsibble(data, index = date)

# Make 2 new variables - logs of bitcoin and M1
data = mutate(data, lm1 = log(m1), lbit = log(bitcoin))

# Plot the data to check 
p1 = autoplot(data, snp) + labs(x = '', y = 'S&P500 Index')
p2 = autoplot(data, lm1) + labs(x = '', y = 'M1 Aggregate')
p3 = autoplot(data, lbit) + labs(x = '', y = 'BTC Price')
p4 = autoplot(data, gvz) + labs(x = '', y = 'Gold Price')
p1/p2/p3/p4

# ACF plot for check trend
gg_tsdisplay(data, snp, plot_type = "partial")
gg_tsdisplay(data, lm1, plot_type = "partial")
gg_tsdisplay(data, lbit, plot_type = "partial")
gg_tsdisplay(data, gvz, plot_type = "partial")

#Structure break
snp_break = cpt.mean(data$snp, 
                     method = 'AMOC')
snp_break
plot(snp_break)

m1_break = cpt.mean(data$m1, 
                     method = 'AMOC')
m1_break
plot(m1_break)

bit_break = cpt.mean(data$bitcoin, 
                    method = 'AMOC')
bit_break
plot(bit_break)

gvz_break = cpt.mean(data$gvz, 
                     method = 'AMOC')
gvz_break
plot(gvz_break)

# Step 1. Сheck for stationary with ADF test

# H0: Non-stationary
# H1: stationary
# if ADF < DF_cr, then H0 rejected

summary(ur.df(data$snp, type = "none", selectlags = c("AIC")))
summary(ur.df(data$snp, type = "drift", selectlags = c("AIC")))
summary(ur.df(data$snp, type = "trend", selectlags = c("AIC")))
#Non-stationary

summary(ur.df(data$lm1, type = "none", selectlags = c("AIC")))
summary(ur.df(data$lm1, type = "drift", selectlags = c("AIC")))
summary(ur.df(data$lm1, type = "trend", selectlags = c("AIC")))
#Stationary at 5% with drift and non-st with trend

summary(ur.df(data$lbit, type = "none", selectlags = c("AIC")))
summary(ur.df(data$lbit, type = "drift", selectlags = c("AIC")))
summary(ur.df(data$lbit, type = "trend", selectlags = c("AIC")))
#Non-stationary

summary(ur.df(data$gvz, type = "none", selectlags = c("AIC")))
summary(ur.df(data$gvz, type = "drift", selectlags = c("AIC")))
summary(ur.df(data$gvz, type = "trend", selectlags = c("AIC")))
#Stationary with drift and trend

data = mutate(data, dsnp = difference(snp,1), dlbit = difference(lbit,1),
              dgvz = difference(gvz,1), dlm1 = difference(lm1,1))
data = tail(data,-1)
summary(ur.df(data$dsnp, type = "drift"))
#Stationary
summary(ur.df(data$dlbit, type = "drift"))
#Stationary
summary(ur.df(data$dgvz, type = "drift"))
#Stationary
summary(ur.df(data$dlm1, type = "drift"))
#Stationary

# Steps 2-3. ARDL and UECM
ardl = auto_ardl(snp ~ lm1 + lbit + gvz, data = data, max_order = 5)
ardl$best_order
ardl$best_model

ardl2511 = ardl(snp ~ lm1 + lbit + gvz, data = data, order = c(2, 5, 1, 1))
ardl2511

uecm_model = uecm(ardl2511)
res_uecm = uecm_model$residuals
plot(res_uecm)

# Steps 4-5. Check the absence of autocorrelation in the residuals of the UECM model and their stationarity
dwtest(ardl2511) # lm-test H0: absence of autocorrelation
summary(ur.ers(res_uecm)) # res are stationary

# Step 6. Bounds test
bounds_f_test(uecm_model, case = 3)
bft = bounds_f_test(uecm_model, case = 3, alpha = 0.05)
bft
bft$tab 
# no cointegration

# Steps 2.1-3.1. ARDL and UECM for diff
ardld = auto_ardl(dsnp ~ dlm1 + dlbit + dgvz, data = data, max_order = 5)
ardld$best_order
ardld$best_model
ardl3540 = ardl(dsnp ~ dlm1 + dlbit + dgvz, data = data, order = c(3, 5, 4, 0))
ardl3540

uecm_modeld = uecm(ardl3540)
res_uecmd = uecm_modeld$residuals
plot(res_uecmd)

# Steps 4.1-5.1. Check the absence of autocorrelation in the residuals of the UECM model and their stationarity
dwtest(ardl3540) # lm-test H0: absence of autocorrelation
summary(ur.ers(res_uecmd)) # res are stationary

# Step 6.1. Bounds test
bounds_f_test(uecm_modeld, case = 3)
bftd = bounds_f_test(uecm_modeld, case = 3, alpha = 0.05)
bftd
bftd$tab 
# do have cointegration

# Step 7. Estimation
# Estimation of the cointegration ratio 
# OLS
model_coint1 = lm(dsnp ~ dlm1 + dlbit + dgvz, data = data)
summary(model_coint1)
# DOLS
model_coint2 = lm(dsnp ~ dlm1 + dlbit + dgvz + 
                    lag(dlm1) + lead(dlm1, 1) + 
                    lag(dlbit) + lead(dlbit, 1) +
                    lag(dgvz) + lead(dgvz, 1), 
                  data = data) 
summary(model_coint2)

resid = model_coint2$residuals
plot(resid, type = 'l')
summary(ur.ers(resid, type = "DF-GLS", lag.max = 4))

# ECM
data2 = mutate(data, error.ecm = model_coint2$residuals)

model_ecm = lm(dsnp ~ dlm1 + dlbit + dgvz + lag(error.ecm), data = data2)
summary(model_ecm)

dwtest(model_ecm)
plot(model_ecm$residuals, type = "l")
summary(ur.df(model_ecm$residuals, type = "drift", selectlags = "AIC"))


### VAR
# Due to the fact that R predicts 10 steps division was made at (end week minus 10 weeks)
train = filter(data, date < ymd('2019-10-28'))
test = filter(data, date >= ymd('2019-10-28'))

# VAR(p) model estimation
VARselect(train[,8:11], lag.max = 5, type = "const")[["selection"]]
var_model1 = vars::VAR(train[,8:11],p = 3, type = "const")

summary(var_model1)

serial.test(var_model1, lags.pt = 52, type = "PT.asymptotic")
serial.test(var_model1, lags.bg = 52, type = "BG")

# Predicting with VAR Model
fcst1 = predict(var_model1, h = 10)
fcst1
plot(fcst1)

f = as.data.frame(fcst1$fcst)
f
test = mutate(test, pred_dsnp = f$dsnp.fcst, pred_dlbit = f$dlbit.fcst,
              pred_dgvz = f$dgvz.fcst, pred_dlm1 = f$dlm1.fcst)
test


# Impulse Response Function (IRF)

## SP500 Responses
# 1. With orthogonal shocks
# SP500 to SP500 shock
snp_irf1 = irf(var_model1, impulse = "dsnp", response = "dsnp", n.ahead = 20, boot = TRUE, ortho = TRUE)
plot(snp_irf1, ylab = "SP500", main = "Response of the SP500 to the SP500 shocks (with orthogonal)")

# M1 to SP500 shock
lm12snp_irf1 = irf(var_model1, impulse = "dlm1", response = "dsnp", n.ahead = 20, boot = TRUE, ortho = TRUE)
plot(lm12snp_irf1, ylab = "SP500", main = "Response of the SP500 to the M1 shocks (with orthogonal)")

# BTC to SP500 shock
snp_irf1 = irf(var_model1, impulse = "dlbit", response = "dsnp", n.ahead = 20, boot = TRUE, ortho = TRUE)
plot(snp_irf1, ylab = "SP500", main = "Response of the SP500 to the Bitcoin price shocks (with orthogonal)")

# GVZ to M1 shock
lm12snp_irf1 = irf(var_model, impulse = "dgvz", response = "dsnp", n.ahead = 20, boot = TRUE, ortho = TRUE)
plot(lm12snp_irf1, ylab = "SP500", main = "Response of the SP500 to the GVZ shocks (with orthogonal)")

# 2. Without orthogonal shocks
snp_irf1 = irf(var_model1, impulse = "dsnp", response = "dsnp", n.ahead = 20, boot = TRUE, ortho = FALSE)
plot(snp_irf1, ylab = "SP500", main = "Response of the SP500 to the SP500 shocks (without orthogonal)")

# M1 to SP500 shock
lm12snp_irf1 = irf(var_model1, impulse = "dlm1", response = "dsnp", n.ahead = 20, boot = TRUE, ortho = FALSE)
plot(lm12snp_irf1, ylab = "SP500", main = "Response of the SP500 to the M1 shocks (without orthogonal)")

# BTC to SP500 shock
snp_irf1 = irf(var_model1, impulse = "dlbit", response = "dsnp", n.ahead = 20, boot = TRUE, ortho = FALSE)
plot(snp_irf1, ylab = "SP500", main = "Response of the SP500 to the Bitcoin price shocks (without orthogonal)")

# GVZ to M1 shock
lm12snp_irf1 = irf(var_model1, impulse = "dgvz", response = "dsnp", n.ahead = 20, boot = TRUE, ortho = FALSE)
plot(lm12snp_irf1, ylab = "SP500", main = "Response of the SP500 to the GVZ shocks (without orthogonal)")


# FEVD - Forecast Error Variance Decomposition
# Let us analyze the share of variance of variable shocks in the forecast error
# The composition of variances is nothing but the components of the variance of the forecast error of the studied endogenous variable, 
# due to the shock of the remaining endogenous variables, i.e. the contribution of each of these variables to the variance 
# of the forecast of the studied indicator
var_fevd = fevd(var_model1, n.ahead = 10)
var_fevd
plot(var_fevd)
