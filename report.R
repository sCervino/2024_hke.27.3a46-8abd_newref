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

srmods <- table(sapply(mixedSR_boot, function(x) x@desc))/length(mixedSR_boot)

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
  R0_sr                 = ,
  R0_bh                 = ,
  R0_rk                 =
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

