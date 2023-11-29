# model.R - Running rebuilding MPs
# WKREBUILD_toolset/model.R

# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
# modified: Dorleta Garcia 2023-11-13
# Distributed under the terms of the EUPL-1.2


# output of this script is: model.rda

library(icesTAF)
library(mse)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Arguments and options ----

# CHOOSE number of cores for doFuture
cores <- 1

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


#### Implementation error ----
# No implementation error in the calculation of the reference points
mean.iem <- 0
sd.iem   <- 0
iem <- FLiem(method=noise.iem,
             args=list(noise=rlnorm(it, rec(om) %=% mean.iem, sd.iem)))


#### SET UP MP runs ----

# SET intermediate year + start of runs, lags and frequency
# iy in the mseargs = iy-1 because we need advice in 2023
mseargs <- list(iy=iy-1, fy=fy, data_lag=1, management_lag=1, frq=1)

# SETUP standard ICES advice rule
AR <- mpCtrl(list(

  # (est)imation method: shortcut.sa + SSB deviances
  est = mseCtrl(method=shortcut.sa,
    args=list(SSBdevs=sdevs$SSB)),

  # hcr: hockeystick (fbar ~ ssb | lim, trigger, target, min)
  hcr = mseCtrl(method=hockeystick.hcr,
    args=list(lim=0, trigger=refpts(om)$Btrigger, target=refpts(om)$Fmsy,
    min=0, drop=0, metric="ssb", output="fbar")),

  # (i)mplementation (sys)tem: tac.is (C ~ F) + F deviances
  isys = mseCtrl(method=tac.is,
    args=list(recyrs=-2, fmin=0, Fdevs=sdevs$F))
  ))


# SETUP a constant F advice rule tunning the ICES AR.
FcteR <- mpCtrl(list(
  
  # (est)imation method: shortcut.sa + SSB deviances
  est = mseCtrl(method=shortcut.sa,
                args=list(SSBdevs=sdevs$SSB)),
  
  # hcr: hockeystick (fbar ~ ssb | lim, trigger, target, min)
  hcr = mseCtrl(method=hockeystick.hcr,
                args=list(lim=0, trigger=refpts(om)$Btrigger, target=refpts(om)$Fmsy,
                          min= refpts(om)$Fmsy, drop = 0, metric="ssb", output="fbar")),
  
  # (i)mplementation (sys)tem: tac.is (C ~ F) + F deviances
  isys = mseCtrl(method=tac.is,
                 args=list(recyrs=-2, fmin=0, Fdevs=sdevs$F))
))



# plot HCRs
plot_hockeystick.hcr(AR$hcr, labels=c(lim="Blim", trigger="MSYBtrigger",
  min="", target="Ftarget")) +
  xlab(expression(hat(SSB))) + ylab(expression(bar(F)))

plot_hockeystick.hcr(FcteR$hcr, labels=c(lim="Blim", trigger="MSYBtrigger",
                                      min="", target="Ftarget")) +
  xlab(expression(hat(SSB))) + ylab(expression(bar(F)))


# - RUN applying ICES advice rule
# system.time(
#   advice <- mp(om, iem=iem, ctrl= AR, args=mseargs)
# )

# PLOT
# plot(runf0, advice, window=FALSE)
 


#### Simulations with constant F in the range [0,2] ---- 
# system.time(
#   FcteSims <- lapply(1:10, function(j){
#           aux <- mps(window(om, start=2020), ctrl=FcteR, args=mseargs, hcr=lapply(opts, function(x) x[j]))
#           res <- metrics(aux)
#           return(res)}
#           ))

system.time(
  FcteSims <- lapply(31:31, function(j){
    cat('.................  ', j, ' ................../n')
    FcteRtemp <- FcteR
    FcteRtemp$hcr@args$target[] <- opts$target[j]
    FcteRtemp$hcr@args$min[] <- opts$target[j]
  system.time(aux <- mp(om, iem=iem, ctrl=FcteRtemp, args=mseargs, verbose = TRUE, parallel = F))
    res <- metrics(aux)
    return(res)}
  ))



# --- SAVE

# save(FcteSims, file="model/FcteSims.rda", compress="xz")
# 
# #### Simulations with constant F in the range [0,2] ---- 
# system.time(
#       ARSims <- mps(window(om, start=2020), ctrl=AR, args=mseargs, hcr=opts)
# )
# 
# save(FcteSims, file="model/ARSims.rda", compress="xz")

# CLOSE cluster
plan(sequential)
