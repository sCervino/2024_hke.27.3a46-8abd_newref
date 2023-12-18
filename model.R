# model.R - Running rebuilding MPs
# WKREBUILD_toolset/model.R

# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
# modified: Dorleta Garcia 2023-11-13
# Distributed under the terms of the EUPL-1.2


# output of this script is: model.rda

library(icesTAF)
library(mse)
library(msy)
library(FLSRTMB)
library(FLRef)
library(ggplot2)
library(RColorBrewer)
library(Rfast)
library(tidyverse)
library(grid)
library(gridExtra)
library(patchwork)
library(ggpubr)
library(ggExtra)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Arguments and options ----

# CHOOSE number of cores for doFuture
cores <- 2

setwd('C:/use/OneDrive - AZTI/ICES WK/WKNEWREF/2024_hke.27.3a46-8abd_newref')

# LOAD oem and oem
load('data/data.rda')
load("data/Biomass_refpts.rda")


it <- 1000

# important years
iy <- 2023
fy <- 2100


Blim     <- Blim_segreg
  
  
#### Short cut approach: F and SSB deviances ----
sdevs <- shortcut_devs(om, Fcv = 0.212, Fphi= 0.423, SSBcv = 0)

# Grid of F-s to use in the forward projection to fing Fmsy
opts <- list(lim=rep(0,101), target=seq(0,1, 0.01), min=seq(0,1, 0.01))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Create 'model' folder
mkdir("model")

# Sourve utilities scripts with additional functions/utilities.
source("utilities.R")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####    **** Equilibrium Analysis with eqSim **** ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Historical stock object
stk <- iter(window(om@stock, 1978, 2019),1)

# Conduct a bootstrap with 1000 iterations, 
#    Fit the three SRR to each of the iterations
#    Select the model that best fits to the data based on the AIC.
sr_fit <- eqsr_fit(stk,
                nsamp = 1000,
                models = c("Ricker", "Segreg", "Bevholt"))

# deterministic fit to calculate the AR1 parameter
sr_det  <- srrTMB(as.FLSRs(stk, models=c("segreg")), spr0=mean(spr0y(stk)))
resid   <- sr_det[[1]]@residuals[drop=T]
rho_ar1 <- coef(lm(resid[-1] ~ resid[-length(resid)]-1))

# autocorrelation of order(1) in residuals

# Select Blim from the previous analysis.
Blim <- Blim_segreg

# If the standard error of the SSB calculated by the stock assessment model in 
# the last year is bigger that 0.2 
#     Bpa <- Blim*exp(-1.645*sigma)
#  Otherwise:
#     Bpa <- 1.4*Blim 
Bpa <- Blim*1.4

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Equilibrium analysis to calculate Flim ----~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The population is projected using constant F with perfect observation/implementation.
eqPop_Flim <- eqsim_run(sr_fit,
                Fcv = 0, Fphi = 0, SSBcv = 0,
                            rhologRec = rho_ar1,
                            Btrigger = 0, Blim = Blim, Bpa = Bpa,
                            Nrun = 200, Fscan = seq(0,1.0,0.05), verbose = F)


eqsim_plot_range(sim = eqPop_Flim, type="median")

Flim <- round(eqPop_Flim$Refs2['catF','F50'],3)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Equilibrium analysis to calculate Flim ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The same as in the previous step, but assessment error is added, 
# using the default values.
cvF  <- 0.212                                 # Default = 0.212
phiF <- 0.423                                 # Default = 0.423
# SSB
cvSSB <- 0                                    # Default = 0
phiSSB <- 0                                   # Default = 0

eqPop_Fmsy <- eqsim_run(sr_fit,
                        Fcv=cvF, Fphi=phiF, SSBcv=cvSSB,
                        Btrigger = 0, Blim=Blim,Bpa=Bpa,
                        rhologRec = rho_ar1,
                        Nrun=200, Fscan = seq(0,1,0.05),verbose=F)

Fmsy_tmp <- round(eqPop_Fmsy$Refs2["lanF","medianMSY"],3)
eqsim_plot_range(eqPop_Fmsy, type="median")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Theoretical & 'real' Bmsy ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# In the case of northern stock of hake, the stock has been exploited at levels
# below Fmsy for more than 10 years.

dbEq <-eqPop_Fmsy$rbp

dbEq_ssb <- dbEq[dbEq$variable == 'Spawning stock biomass',]

pos <- which(findInterval(unique(dbEq_ssb$Ftarget), Fmsy_tmp) ==1)[1]

BmsyEq <- (Fmsy_tmp - dbEq_ssb[pos-1,'Ftarget'])/(dbEq_ssb[pos,'Ftarget'] - dbEq_ssb[pos-1,'Ftarget'])*dbEq_ssb[pos-1,'p50'] +
        (dbEq_ssb[pos,'Ftarget']-Fmsy_tmp )/(dbEq_ssb[pos,'Ftarget'] - dbEq_ssb[pos-1,'Ftarget'])*dbEq_ssb[pos,'p50'] 

BmsyEq

mean(ssb(stk)[, ac(2010:2019)])

#  SSB in the period with exploitation below 
# Fmsy is higher than the theoretical Bmsy (BmsyEq), which is consistent
# and supports the right definition of Bmsy.


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### MSY Btrigger ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Definition: A lower bound of the expected range of the point at which F is reduced when 
#      SSB when the stock is fished at F applying the ICES MSY advice rule (AR). 
# Basis: MSY Btrigger = maximum (Bpa, the 5th percentile of the distribution of 
#      SSB when fishing at F), modified according to the scheme for determining 
#      MSY BMSYtrigger(described in the section on MSY reference points). 
# WKREF1/2:

# Theoretical MSY Btrigger based on equilibrium values.
#  Extrapolate equilibrium values in the 5% quantile.
F.05   <- dbEq[dbEq$variable == "Spawning stock biomass", ]$Ftarget
ssb.05 <- dbEq[dbEq$variable == "Spawning stock biomass", ]$p05
plot(b.05~x.05, ylab="SSB", xlab="F")
abline(v=Fmsy_tmp)
i <- which(x.05<Flim)
b.lm <- loess(b.05[i] ~ x.05[i])
lines(x.05[i],c(predict(b.lm)),type='l')

MSYBtrigger_theo <- round(predict(b.lm,Fmsy_tmp))
abline(h=MSYBtrigger_theo)

