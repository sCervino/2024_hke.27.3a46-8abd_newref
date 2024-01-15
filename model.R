#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# model.R - 
# 2024_2024_hke.27.3a46-8abd_newref/data.R
#
# Dorleta Garcia
# 2024/01/09
#
# Input: Files created in the data.R script
#           - 'data/data.rda'
#           - 'data/Biomass_refpts.rda'
#
# Output:
#      - 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### Load libraries ----
# FLR related packages
library(FLRef)
library(FLSRTMB)
library(mse)
# ICES packages
library(icesTAF)
library(msy)
# tidyverse packages
library(tidyverse)
# ggplot related packages
library(RColorBrewer)
library(ggplot2)
library(gridExtra) # arrange multiple plots in a page
library(ggpubr)    #  provides some easy-to-use functions for creating and customizing 'ggplot2'
library(ggExtra)   # Marginal histograms
library(ggthemes)  # Themes for ggplot


library(devtools)
install_github("ices-tools-prod/msy")

# library(Rfast)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Arguments and options ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load data and reference points
#load('data/data.rda')
load("data/Biomass_refpts.rda")

it <- 100

#### Load the FLStock object ---- 
load('bootstrap/data/wgbie2023_nhke_FLStock_csim.RData')
# Shorten the FLStock object based on previous analysis.
stk <- window(hke.stk, dimnames(hke.stk)[[2]][1], 2021)

# First year in the projection
iy <- 2022

# Final year in the projection
fy <- 2100

bio.years <- c(2017,2021)
sel.years <- c(2017,2021)
#### Proposed Blim based on previous analysis--- 
Blim     <- Blim_segreg
  
#### Short cut approach: F and SSB deviances ----
#sdevs <- shortcut_devs(om, Fcv = 0.212, Fphi= 0.423, SSBcv = 0)

# Grid of F-s to use in the forward projection to fing Fmsy
#opts <- list(lim=rep(0,101), target=seq(0,1, 0.01), min=seq(0,1, 0.01))

# Create 'model' folder
mkdir("model")

# Source utilities scripts with additional functions/utilities.
source("utilities.R")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  1.  Stock recruitment relationships for eqSim ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Conduct a bootstrap with 1000 iterations, 
#    Fit the three SRR to each of the iterations
#    Select the model that best fits to the data based on the AIC.
sr_fit <- eqsr_fit(stk,
                nsamp = it,
                models = c("Ricker", "Segreg", "Bevholt"))

# Conduct a bootstrap with 1000 iterations fitting only the Segmented regression model
sr_segreg <- eqsr_fit(stk,
                   nsamp = it,
                   models =  "Segreg")

# deterministic fit to calculate the AR1 parameter
# autocorrelation of order(1) in residuals
sr_det  <- srrTMB(as.FLSRs(stk, models=c("segreg")), spr0=mean(spr0y(stk)))
resid   <- sr_det[[1]]@residuals[drop=T]
rho_ar1 <- coef(lm(resid[-1] ~ resid[-length(resid)]-1))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  2.  Bpa ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# If the standard error of the SSB calculated by the stock assessment model in 
# the last year is bigger that 0.2 
#     Bpa <- Blim*exp(-1.645*sigma)
#  Otherwise:
#     Bpa <- 1.4*Blim 
Bpa <- Blim*1.4

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 2.  Equilibrium analysis to calculate ** Flim**  ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  *** ICES guidelines ***
# Simulate a stock with a segmented regression S–R relationship, 
# "with the point of inflection at Blim", thus determining the F = Flim which, 
# at equilibrium, yields a 50% probability of SSB > Blim. Note that this simulation 
# should be conducted based on a fixed F (i.e. without inclusion of a Btrigger) and without 
# inclusion of assessment/advice errors (In the MSY R package (also known as EqSim), 
# this means Btrigger, Fcv, and Fphi should all be set to zero). 

eqPop_Flim <- eqsim_run(sr_segreg,
                        bio.years = bio.years ,
                        sel.years = sel.years ,
                        Fcv = 0, Fphi = 0, SSBcv = 0,
                        rhologRec = rho_ar1,
                        Btrigger = 0, Blim = Blim, Bpa = Bpa,
                        Nrun = 200, Fscan = seq(0,1.0,0.05), verbose = F)


taf.png("report/eqsim_Flim.png")
eqsim_plot_range(sim = eqPop_Flim, type="median")
dev.off()

Flim <- round(eqPop_Flim$Refs2['catF','F50'],3)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 3.  Equilibrium analysis to calculate ** Fmsy ** ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# In this case the stock-recruitment relationship is the one that comes from the
# best stock-recruitment model in each of the bootstrap iterations.
# Assessment error is added, the default in this case. 
# The obtained Fmsy is temporal (Fmsy_tmp) because it could be capped by Fpa = Fp.05

# Default assessment error values. 
cvF  <- 0.212                                 # Default = 0.212
phiF <- 0.423                                 # Default = 0.423
# SSB
cvSSB <- 0                                    # Default = 0
phiSSB <- 0                                   # Default = 0

eqPop_Fmsy <- eqsim_run(sr_fit,
                        bio.years = bio.years,
                        sel.years = sel.years,
                        Fcv=cvF, Fphi=phiF, SSBcv=cvSSB,
                        Btrigger = 0, Blim=Blim,Bpa=Bpa,
                        rhologRec = rho_ar1,
                        Nrun=200, Fscan = seq(0,1,0.05),verbose=F)



Fmsy_tmp <- round(eqPop_Fmsy$Refs2["lanF","medianMSY"],3)


taf.png("report/eqsim_Fmsy.png")
eqsim_plot_range(eqPop_Fmsy, type="median")
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 4.  MSY Btrigger ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Definition: A lower bound of the expected range of the point at which SSB is   
#      reduced when the stock is fished at Fmsy applying the ICES MSY advice rule (AR). 
#      In the ICES MSY approach, MSY Btrigger is defined as the 5% percentile on 
#      the distribution of SSB when fishing at FMSY
#     
#      For most stocks that lack data on fishing at FMSY, MSY Btrigger is set at Bpa. 
#      However, as a stock starts to be fished consistently with FMSY, it is possible 
#      to move towards implementation of a value for MSY Btrigger that reflects the 5%
#      percentile definition of MSY Btrigger. 

# Basis: MSY Btrigger = maximum (Bpa, the 5th percentile of the distribution of 
#      SSB when fishing at F), modified according to the scheme for determining 
#      MSY BMSYtrigger(described in the section on MSY reference points). 
# WKREF1/2:

# Theoretical MSY Btrigger based on equilibrium values.
#  Interpolate equilibrium values in the 5% quantile.
dbEq_Fmsy <- eqPop_Fmsy$rbp
Fs     <- dbEq_Fmsy[dbEq_Fmsy$variable == "Spawning stock biomass", ]$Ftarget
ssb.05 <- dbEq_Fmsy[dbEq_Fmsy$variable == "Spawning stock biomass", ]$p05


plot(ssb.05~Fs, ylab="tonnes", xlab="F", main = '5% percentile of SSB versus F')
abline(v=Fmsy_tmp)
i <- which(Fs<Flim)
b.lm <- loess(ssb.05[i] ~ Fs[i])
lines(Fs[i],c(predict(b.lm)),type='l')

MSYBtrigger_temp <- round(predict(b.lm,Fmsy_tmp))
abline(h=MSYBtrigger_temp, col = 2, lwd = 2)
text(0.1,MSYBtrigger_temp*1.15, expression(MSYB[trigger]) )


# Standard deviation, σ, of ln(SSB) in the final assessment year
# we assume ln(SSB) ~ N(est_ln(SSB), sigma) =>
# calculate the 5th percentile.
sd <- 0.147 
SSB_p05 <- exp(qnorm(0.05, log(ssb(stk)[, '2019']), 0.147))


MSYBtrigger <- max(Bpa, MSYBtrigger_temp)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 5.   Fp.05 = Fpa ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Definition: An exploitation rate reference point below which exploitation is 
#     considered to be sustainable, having accounted for estimation uncertainty. 
#  *** ICES guidelines ***
# Calculation: The fishing mortality including the advice rule that, if applied as a
#   target in the ICES MSY advice rule (AR) would lead to SSB ≥ Blim with a
#   95% probability (also known as Fp05). The derivation of Fp.05 should include 
#   the expected stochastic variability in biology and fishery, as well as advice error. 
#  *** WKREF1/WKREF2 ***
#    Calculation without MSY Btrigger.
eqPop_Fp05_AR <- eqsim_run(sr_fit,
                           bio.years = bio.years,
                           sel.years = sel.years,
                           Fcv=cvF, Fphi=phiF, SSBcv=cvSSB,
                           Btrigger = MSYBtrigger, Blim=Blim,Bpa=Bpa,
                           rhologRec = rho_ar1,
                           Nrun=200, Fscan = seq(0,1,0.05),verbose=F)

Fp05_AR <- eqPop_Fp05_AR$Refs2["catF","F05"]

taf.png("report/eqsim_Fp05_AR.png")
eqsim_plot(eqPop_Fp05_AR)
dev.off()

# The proposal of WKREF1/WKREF2 is derived from the same simulations as Fmsy.
taf.png("report/eqsim_Fp05_NAR.png")
eqsim_plot(eqPop_Fmsy)
dev.off()

Fp05_NAR <- eqPop_Fmsy$Refs2["catF","F05"]

#### Final Fmsy ----
Fmsy <- ifelse(Fmsy_tmp < Fp05_AR, Fmsy_tmp, Fp05_AR)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 6.  Equilibrium analysis to calculate ** Bmsy **  ----
# *** ICES guidelines ***^
# Bmsy corresponding to Fmsy in the equilibrium analysis but without
# introducing assessment error in the simulation. (i.e eqPop_Flim)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Fscan <- seq(0,1,0.05)
SSBFscan_p50 <- eqPop_Flim$rbp$p50[eqPop_Flim$rbp$variable=="Spawning stock biomass"]
## Interpolate to get SSB for more F values, percentile 50 of SSB for (SSB,F) pairs. 
SSBF_p50 <- as.data.frame(approx(Fscan,                             # The F-s for which we have F
                                 SSBFscan_p50,                      # The p50 SSB-s corresponding to Fscan
                                 xout = seq(min(Fscan), max(Fscan),length=1000))) # F values over we want to interpolate
names(SSBF_p50) <- c('F', 'SSB')

# BMSY: The SSB corresponding to the F closest to Fmsy.
Bmsy <- SSBF_p50$SSB[which.min(abs(SSBF_p50$F - Fmsy))]

# In the case of Northern stock of Hake, it has been exploited below or or around Fmsy 
# for 10 years so the theoretical Bmsy can be compared with a bomass corrseponding
# to a period of exploitation below Fmsy.
# The SSB is the period is higher than the theoretical Bmsy (BmsyEq), which is consistent
# and supports the right definition of Bmsy.

mean(ssb(stk)[, ac(2010:2019)])

Bmsy0.8 <- Bmsy*0.8
Bmsy0.5 <- Bmsy*0.5

ratio_MSYBtrigger_Bmsy <- MSYBtrigger_temp/Bmsy

taf.png("report/potential_MSYBtrigger.png")
barplot(c(Bpa = Bpa, q0.05_Bmsy = MSYBtrigger_temp,  q0.05_SSB2019 = SSB_p05, 
          `0.80*Bmsy` =  Bmsy0.8, `0.50*Bmsy` = Bmsy0.5), main = 'Potential MSY Btrigger values',
        ylab = 'tonnes')
dev.off()


save(B0, Blim, Blim_segreg, Blim_segreg_boot, Bmsy, Bmsy0.5, Bmsy0.8, Bpa, Flim , Fmsy, Fp05_AR, Fp05_NAR, 
     MSYBtrigger, file = 'model/refpts.RData')
