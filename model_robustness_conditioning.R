


library(mse)
library(FLSRTMB)

# CHOOSE number of cores for doFuture / doParallel
cores <- 3

source("utilities.R")

# Load the FLStock object
load('bootstrap/data/wgbie2023_nhke_FLStock_csim.RData')
stk <- hke.stk

# Load the FLSRTMB generated in data.R
load("data/SR_analysis.rda")
# Reference points calculated in the model.R script based on the ICES framework.
load('model/ICES_RefPts.RData' )

# Default reference points from previous scripts
refpts <- FLPar(Btrigger = MSYBtrigger, Fmsy = Fmsy, Blim = Blim, Bpa = Bpa, 
                Flim = Flim, Fpa = NA, lFmsy = NA, uFmsy = NA, 
                F05 = Fp05_AR, F05noAR = Fp05_NAR)

# Historical year for recruitment deviations.
hy <- 2010

# INTERMEDIATE year
iy <- 2023

# FINAL year
fy <- 2100

# NUMBER of iterations
it <- 500

# DATA year
dy <- iy - 1



set.seed(987)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  1. CONSTRUCT the OM ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#### *** 1.a GENERATE recruitment deviances ***  ----
# - historical period: lognormal error standardized to mean equal 1.
# - projection period: lognormal autocorrelated

srdevs <- FLQuant(dimnames = list(year = hy:fy), iter = it)
#aux <- FLQuant(rlnorm(it*(dy-hy+1), 0,sd=c(srpars$sigmaR)), dimnames = list(year = hy:dy), iter = it)
aux <- FLQuant(rlnorm(it*(dy-hy+1), 0,sd=0.4), dimnames = list(year = hy:dy), iter = it)
srdevs[,ac(hy:dy)] <- aux/apply(aux, 2, mean)
srdevs[,ac((dy+1):fy)] <- rlnormar1(it, 0,sdlog=mixedSR_boot_params$sigmaR, rho=mixedSR_boot_params$rho, years=seq(dy+1, fy))

taf.png("report/recruitment_deviations.png")
plot(srdevs)
dev.off()



# BUILD FLom
om <- FLom(stock=propagate(stk, it), refpts=refpts, model='mixedsrr',
           params=mixedSR_boot_params, deviances=srdevs)


# SETUP om future: bootstrap of the historical values **
om <- fwdWindow(om, end=fy, nsq=3, years=c(wt=5, mat=5, catch.sel=5))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 2. Hindcast: Introduce Uncertainty in the initial population  ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# If the assessment model allows to introduce uncertainty in a different way
# (e.g bayesian models, bootstrap option...) this step can be skipped.

# SET stochastic rec dy, start deviations is year 'hrdy' so a random population is
# generated for starting conditions: 
om@stock <- hindcast_srdevs(om@stock, srdevs, start = hy, end = dy, process_error = FALSE)
  
#initial population with uncertainty in the initial conditions.
plot(window(om@stock,1999, dy)) 

# compare initial population with the stock assessment output
plot(FLStocks(list(est = window(om@stock,1999, dy), 
                   SS = stk, 
                   mean = window(qapply(om@stock, iterMeans), 1999, dy))))


#### SAVE ----
save(om, file= "model/om_conditioning.rda", compress="xz")
# CLOSE cluster
plan(sequential)
