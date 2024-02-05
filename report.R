# report.R - DESC
# 2024_hke.27.3a46-8abd_newref/report.R

# Dorleta Garcia
# 2024/01/11


#### Load libraries ----
library(icesTAF)
mkdir("report")
library(mse)
library(FLSRTMB)
library(mseviz)
library(scales)

source("utilities.R")


library(FLCore)
library(tidyverse)
library(icesTAF)
library(mse)
library(gridExtra)

source("utilities.R")
mkdir("output")

# LOAD model.R outputs

#### Load the FLStock object ---- 
load('bootstrap/data/wgbie2023_nhke_FLStock_csim.RData')
# Rename the object
stk <- hke.stk

load('model/robustness_test.RData')
load("data/Biomass_refpts.rda")
load("data/SR_analysis.rda")
load('data/brps.rda')
load('model/ICES_RefPts.RData') 

it <- 500

iy <- 2022
fy <- 2050

short  <- iy:(iy+4)
long   <- (fy-10):(fy-1)


stknm <- 'hke.27.3a46-8abd'


#------------------------------------------------------~~~~~~~~~~~~~~~~
#### Reporting TABLE ----
#------------------------------------------------------~~~~~~~~~~~~~~~~

# Stock recruitment models
srmods <- table(sapply(mixedSR_boot, function(x) x@desc))/length(mixedSR_boot)
aux <- as.data.frame(mixedSR_boot_params)
rk_iters <- aux %>% filter(aa$params == 'm', aa$data == 2) %>% select(iter)
sr_iters <- aux %>% filter(aa$params == 'm', aa$data == 3) %>% select(iter)
bh_iters <- aux %>% filter(aa$params == 'm', aa$data == 1) %>% select(iter)

# spawning per recruit in the absence of fishing.
sprf0 <- spr0(stk)

# R0 and B0 from Ricker model
rk_params <- apply(mixedSR_boot_params[c('a', 'b'),rk_iters$iter],1, median)[drop=T]
rk_B0 <- log(rk_params['a'] * sprf0) / rk_params['b'] 
rk_R0 <- rk_params['a']  * rk_B0 * exp(-rk_params['b']  * rk_B0)

# R0 and B0 from Beverton-Holt model
## R = a * S / (b + S)
bh_params <- apply(mixedSR_boot_params[c('a', 'b'),bh_iters$iter],1, median)[drop=T]
bh_B0 <- bh_params['a'] * sprf0 - bh_params['b']
bh_R0 <- bh_params['a'] * S0 / (bh_params['b'] + S0)

# R0 and B0 from Segmented regression model
## R = a * S / (b + S)
sr_params <- apply(mixedSR_boot_params[c('a', 'b'),sr_iters$iter],1, median)[drop=T]
sr_R0 <- prod(sr_params)
sr_B0 <- NA


reference_points_table <- data.frame(
  stock                 = stknm,
  assessment_model      = 'SS',
  software_used         = 1,
  software_explain      = 'wknewref code',
  stock_status          = 3,
  life_span             = 2,
  assessmet_wg          = 'WGBIE',
  contact_person        = 'Dorleta Garcia',
  contact_person_email  = 'dgarcia@azti.es',
  bevholt               = unname(ifelse(is.na(srmods['bevholt']), 0, srmods['bevholt'])),
  ricker                = unname(ifelse(is.na(srmods['ricker']), 0, srmods['ricker'])),
  segreg                = unname(ifelse(is.na(srmods['segreg']), 0, srmods['segreg'])),
  BP_segreg             = median(segreg_boot['b',drop=T]),
  BP_segreg_ok          = TRUE,
  BP_segrreg_explain    = "Breakpoint sensitive to the exclusion of annual data", 
  SR_stock_assessment   = FALSE,
  R0_sr                 = sr_R0,
  R0_bh                 = bh_R0,
  R0_rk                 = rk_B0,
  R0_ok                 = TRUE,
  R0_explain            = "Asymptote well estimated based on high recruitments at high biomass levels",
  B0                    = median(B0),
  B0_ok                 = TRUE, 
  B0_explain            = "60% higher than maximum biomass that was already high",
  Bloss                 = min(ssb(stk)),
  Bmax                  = max(ssb(stk)),
  Blim_current          = 61563,
  Blim_wknewref         = Blim,
  Blim_rationale        = 1,
  Blim_explain          = "segmented regression breakpoint",
  MSYBtrigger_current   = ,
  MSYBtrigger_wknewref  =,
  MSYBtrigger_rationale =,
  MSYBtrigger_explain   =,
  Fmsy_current          =,
  Fmsy_wknewref         =,
  Ftarget               =,
  Ftarget_rationale     =,
  Ftarget_explain       =,
  Fp.05_current         =,
  Fp.05_AR_wknewref     =,
  Fp.05_NAR_wknewref    =,
  Fp.05_explain         =,
  Flim                  =)

#------------------------------------------------------~~~~~~~~~~~~~~~~
#### [Ftarget, MSYBtrigger] results TABLE ----
#------------------------------------------------------~~~~~~~~~~~~~~~~

