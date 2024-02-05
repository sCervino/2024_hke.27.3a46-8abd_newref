# model_robustness_test.R - Testing robustness of the reference points.
# /model_robustness_test.R

# Based on the code developed for WKREBUILD2
# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
# Modified: Dorleta Garcia 2024-01-29 <dorleta.garcia@ices.dk> - <dgarcia@azti.es>
# Distributed under the terms of the EUPL-1.2


# input of this script is: 
#      - model/ICES_RefPts.RData 
#      - om_conditioning.rda
# output of this script is: model.rda



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 0. Libraries, Input data, Code
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(icesTAF)
library(FLSRTMB)
library(mse)
# Load the ICES reference points calculated in model.R script and the
# OM conditioning from model_robustness_conditioning.R script.
load("model/ICES_RefPts.RData")
load("model/om_conditioning.rda")
load('bootstrap/data/wgbie2023_nhke_FLStock_csim.RData')
# Rename the object
stk <- hke.stk


# Source utilities scripts with additional functions/utilities.
source("utilities.R")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 1. Arguments and options ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# *** Number of cores for doFuture *** 
#  The scenarios will run in parallel. 
# The amount of cores you can use will depend on number of processors in your 
# PC and amount on memory (RAM). 
cores <- 3

# *** Number of iterations ***
it <- 500

# *** Important years ***
iy <- 2023 # The first intermediate year in the simulation
fy <- 2050 # Last year in the simulation


opts <- list(lim      = rep(0, nrow(MSYBtriggerF_grid)),
             trigger  = MSYBtriggerF_grid[,2],
             target   = MSYBtriggerF_grid[,1],
             min      = rep(0, nrow(MSYBtriggerF_grid)))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 2. SET UP Management procedure runs ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# *** Short cut approach: F and SSB deviances ***
sdevs <- shortcut_devs(om, Fcv = 0.212, Fphi= 0.423, SSBcv = 0)

# Implementation error 
# No implementation error in the calculation of the reference points
mean.iem <- 0
sd.iem   <- 0
iem <- FLiem(method=noise.iem,
             args=list(noise=rlnorm(it, rec(om) %=% mean.iem, sd.iem)))

# *** SET intermediate year + start of runs, lags and frequency ***
# iy in the mseargs = iy-1 because we need advice in 2023
mseargs <- list(iy=iy-1, fy=fy, data_lag=1, management_lag=1, frq=1)


#### 2.a SETUP standard ICES advice rule ----
AR <- mpCtrl(list(

  # (est)imation method: shortcut.sa + SSB deviances
  est = mseCtrl(method=shortcut.sa,
    args=list(SSBdevs=sdevs$SSB)),

  # hcr: hockeystick (fbar ~ ssb | lim, trigger, target, min)
  hcr = mseCtrl(method=hockeystick.hcr,
    args=list(lim=0, trigger=refpts(om)$Btrigger, target=refpts(om)$Fmsy,
    min=0, metric="ssb", output="fbar")),

  # (i)mplementation (sys)tem: tac.is (C ~ F) + F deviances
  isys = mseCtrl(method=tac.is,
    args=list(recyrs=-2, fmin=0, Fdevs=sdevs$F))
  ))

# plot HCRs
plot_hockeystick.hcr(AR$hcr, labels=c(lim="Blim", trigger="MSYBtrigger",
  min="", target="Ftarget")) +
  xlab(expression(hat(SSB))) + ylab(expression(bar(F)))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 3. Run simulations for [MSY Btrigger, Ftarget] combinations ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Be sure that the number of iterations in om@stock and om@sr are exactly the same. 
# In my laptop with 3 cores and 20 [MSYBtrigger, Ftarget] combinations it 
#  Takes 1 hour
om@sr <- iter(om@sr,1:500)
system.time(
      ARSims <- mps(om, ctrl=AR, args=mseargs, hcr= opts)
)

names(ARSims) <- rownames(MSYBtriggerF_grid)

# Plot with all the HCRs tested.
Reduce('+', Map(function(x, y)
  plot_hockeystick.hcr(control(x)$hcr,
                       labels=c(trigger="Btrigger", target="Ftarget")) +
    xlab("SSB") + ylab(expression(bar(F))) +
    ggtitle(y), x=ARSims, y=names(ARSims)))

plot(lapply(ARSims, function(x) x@om@stock))
plot(lapply(ARSims, function(x) qapply(x@om@stock, function(x) apply(x, 1:5, median)))) +
  ggtitle('[MSY Btrigger, Ftarget] GRID') 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 4. Is MSY Btrigger a real buffer for Blim?  ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 1. Reduce the population to a point between Blim and MSY Btrigger
#    For other stocks, higher or lower F could be required.
om_red <- om
om_red@stock <- fwd(om_red@stock, sr=rec(om_red@stock), fbar=fbar(om_red@stock)[, ac(2010:2022)] * 2.5)
ssb(om_red@stock)[, '2022']/MSYBtrigger

# 2. Test the performance of the AR.
system.time(AR_MSYBtrigger <- mp(om_red, iem=iem, ctrl= AR, args=mseargs))

plot(AR_MSYBtrigger, window=FALSE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 4. AR[MSYBtrigger, Blim, Ftarget] Robust to low productivity period? ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# In the case of HKE we will take [1990, 2003] period to condition the SR 
# relationship.
# We have to run two scenarios:
#    4.a) Low productivity in the OM + Reference Points based on historical productivity 
#    4.b) Historical productivity in the OM + Reference Points based on low productivity.

#### *** 4.a) OM[Low productivity] + RefPts[Historical productivity]    *** ----
stk_low <- stk[,ac(1990:2003)]
sr.fits     <- srrTMB(as.FLSRs(stk_low, models=c("segreg", "ricker", "bevholt")), spr0=mean(spr0y(stk)))
plot(sr.fits)

mixedSR_low_boot_all   <- bootstrapSR_list(stk_low, iters=it, models=c("ricker", "bevholt", 'segreg'), method="best")
mixedSR_low_boot       <- mixedSR_boot_all[[1]]
mixedSR_low_boot_params <- mixedSR_boot_all[[2]]

table(sapply(mixedSR_low_boot, function(x) slot(x, 'desc')))

om_lowProd <- om
om_lowProd@sr@params <- om_lowProd@sr@params[c('a', 'b', 'm'),]

system.time(AR_OMlowProd <- mp(om_lowProd, iem=iem, ctrl= AR, args=mseargs))

plot(AR_OMlowProd, window=FALSE)

#### *** 4.b) OM[Historical productivity] + RefPts[Low productivity]    *** ----
# Set up the AR for the regimen with low productivity.
# In this example we set up some hypothetical reference points dividing 
# Blim and MSYBtrigger by 2.
MSYBtrigger_lowProd <- 55000
Fmsy_lowProd <- 0.27
AR_lowProd <- mpCtrl(list(
  
  # (est)imation method: shortcut.sa + SSB deviances
  est = mseCtrl(method=shortcut.sa,
                args=list(SSBdevs=sdevs$SSB)),
  
  # hcr: hockeystick (fbar ~ ssb | lim, trigger, target, min)
  hcr = mseCtrl(method=hockeystick.hcr,
                args=list(lim=0, trigger=MSYBtrigger_lowProd, target=Fmsy_lowProd,
                          min=0, metric="ssb", output="fbar")),
  
  # (i)mplementation (sys)tem: tac.is (C ~ F) + F deviances
  isys = mseCtrl(method=tac.is,
                 args=list(recyrs=-2, fmin=0, Fdevs=sdevs$F))
))

system.time(ARlowProd_OMHistProd <- mp(om, iem=iem, ctrl= AR_lowProd, args=mseargs))


save(ARSims, AR_MSYBtrigger, AR_lowProd, ARlowProd_OMHistProd, file = 'model/robustness_test.RData' )
# CLOSE cluster
plan(sequential)



#### AR[OPRs] Robust to minimum assessment uncertainty?  ----
# Here we use a short-cut MSE with the same implementation error as in
# eqSim by including the 1 year time lag.


#### AR[OPRs] Robust to observed low productivity in the historical period?  ----


#### AR[OPRs] Btrigger well defined to avoid Blim ----
# Applying a constant high F reduce the SSB to a point in (Bpa, MSYBtrigger)
# and analyse the capacity of the AR[OPRs] to rebuild the stock above Bpa.


### AR[OPRs] Btrigger, Blim well defined to recover the stock ----



