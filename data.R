# data.R - condition OM(s)
# WKREBUILD_toolset/data.R

# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
# modified: Dorleta Garcia 2023-11-13

# Distributed under the terms of the EUPL-1.2



library(icesTAF)
library(mse)
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

# CHOOSE number of cores for doFuture / doParallel
cores <- 3

# Load the FLStock obtained from the assessment results. 
load('bootstrap/data/wgbie2023_nhke_FLStock_csim.RData')

#### Current refpts ----
refpts <- FLPar(Btrigger = 78405, Fmsy = 0.24, Blim = 61563, Bpa = 78405, 
                Flim = 0.73, Fpa = 0.54, lFmsy = 0.147, uFmsy = 0.37, 
                F05 = 0.54, F05noAR = NA)

# INTERMEDIATE year
iy <- 2023
# FINAL year
fy <- 2100
# Years to fit Stock-Recruitment relationship.
# Selected after the first analysis: segreg_finalYr
recy <- 1978:2020

# NUMBER of iterations
it <- 1000

set.seed(527)

# Rename the object
stk <- hke.stk


# Create the data directory to: 
# Create report directory to save plots for the report.
mkdir("data")
mkdir('report')

# Source some utilities
source("utilities.R")


# DATA year
dy <- iy - 1



#### - Stock-recruitment relationship(s) ----


#### Select the time series to be used to calculate reference points ----
# H0: we select the whole time series
# H1: We select a shorter time series to account for retrospective patterns
#     in the most recent years.

# Is the breakpoint sensitive to the length of the time series?
yrs <- dimnames(stk)$year
segreg_finalYr <- matrix(NA,2, 10, dimnames = list(c('a', 'b'), 2013:2022))

for(y in 2013:2022){
  
  temp <- srrTMB(as.FLSRs(window(stk, 1978, y), models=c("segreg")), spr0=mean(spr0y(stk)))
  
  segreg_finalYr[, ac(y)] <- c(temp[[1]]@params)
}

taf.png("report/recruitment_segreg_Finalyr.png")
plot(segreg_finalYr[1,], segreg_finalYr[2,], main = 'SR params, time series up to year "y"',
     xlab = "a", ylab = "b (breakpoint)", pch = 4, col = 3, xlim = c(14.9, 15.10))
text(segreg_finalYr[1,]+0.005, segreg_finalYr[2,]+100, 2013:2022, cex = 0.8)
dev.off()
# The breakpoint is very sensitive to the last three year points, with much more 
# consistency in previous years => Remove last 3 years from the fit.


# Is the breakpoint sensitive to single (rec, ssb) pairs?
segreg_1yrOut <- matrix(NA,2, length(recy), dimnames = list(c('a', 'b'), recy))

for(i in 1:length(recy)){
  if(i == 1)               subs <- 2:length(recy)
  if(i %in% 2:length(recy)) subs <- c(1:(i-1), (i+1):length(recy))
  if(i == length(recy))     subs <- 1:(length(recy)-1)
  
  temp <- srrTMB(as.FLSRs(stk[,subs], models=c("segreg")), spr0=mean(spr0y(stk)))
  
  segreg_1yrOut[, i] <- c(temp[[1]]@params)
}

taf.png("report/recruitment_segreg_1yrOut.png")
plot(segreg_1yrOut[1,], segreg_1yrOut[2,], main = 'SR params, 1 year out',
     xlab = "a", ylab = "b (breakpoint)", pch = 4, col = 3)
text(segreg_1yrOut[1,]+0.01, segreg_1yrOut[2,]+250, yrs, cex = 0.8)
dev.off()


#### FITSR models ----
# Deterministic fit
sr.fits <- srrTMB(as.FLSRs(window(stk, recy[1], recy[length(recy)]), models=c("segreg", "ricker", "bevholt")), spr0=mean(spr0y(stk)))
Blim_segreg <- sr.fits$segreg@params$b[drop=T]
# PLOT
taf.png("report/recruitment_fits.png")
plotsrs(sr.fits) + ggtitle('Stock recruiment relationships')
dev.off()

# Boostrap fit of the segmented regression model.
segreg_boot <- bootstrapSR(window(stk, recy[1], recy[length(recy)]), iters=it, models=c( "segreg"))
Blim_segreg_boot <- c(segreg_boot['b',])

taf.png("report/recruitment_segreg_boot.png")
plot(segreg_boot[c('a', 'b', 'sigmaR', 'R0', 'rho', 'B0'),]) + ggtitle('Bootstrapped segmented regression parameters')
dev.off()


# stats_brkpt: matrix to store the Blim calculated with different methods.
stats_brkpt <- matrix(NA, 2,4, dimnames = list(c('bootstrap', '1yrOut'), c('mean', 'median', 'cv', 'sd')))

stats_brkpt[1,] <- c(mean(segreg_boot[2,]), median(c(segreg_boot[2,])),
                       sd(segreg_boot[2,])/mean(segreg_boot[2,]), sd(segreg_boot[2,]))
stats_brkpt[2,] <- c(mean(segreg_1yrOut[2,]), median(c(segreg_1yrOut[2,])),
                       sd(segreg_1yrOut[2,])/mean(segreg_1yrOut[2,]), sd(segreg_1yrOut[2,]))


# Plot the bootstrap breakpoint
# taf.png("report/recruitment_segreg_breakpoint.png")
# par(mfrow = c(1,2))
# boxplot(c(segreg_boot[2,]), main = "Boostrap Segmented regression  breakpoint")
# mtext(paste('mean = ', round(stats_brkpt[1,1]), ', median = ', round(stats_brkpt[1,2]), 
#             'cv = ', round(stats_brkpt[1,3],2), ', sd = ', round(stats_brkpt[1,4])), 1, line = 2)
# plot(density(c(segreg_boot[2,])), main = "", xlab = "")
# abline(v = round(stats_brkpt[1,2]), col = 2)
# dev.off()
taf.png("report/recruitment_segreg_breakpoint.png")
segreg_boot_df <- as.data.frame(t(segreg_boot[drop=T]))
plot1 <- ggplot(segreg_boot_df, aes(x = a, y = b)) + 
  geom_point(color = "#619CFF") + 
  stat_smooth(method = "lm", fullrange = TRUE) +
  geom_rug() + 
  scale_y_continuous(name = "b (the breakpoint)",  limits = range(segreg_boot_df$b) ) + 
  scale_x_continuous(name = "a (the slope)",  limits = c(10, 40), expand = c(0, 0)) + 
  theme_pubr()  + ggtitle("Boostraped Segmented regression parameters")

ggMarginal(plot1, fill = "#619CFF") 
dev.off()


#### BOOTSTRAP and SELECT the SR model by largest logLik ----
mixedSR_boot <- bootstrapSR_list(stk, iters=it, models=c("ricker", "bevholt", 'segreg'), method="best")

#### SAVE SR FITS ----
save(sr.fits, segreg_boot, mixedSR_boot, stats_brkpt, file="data/SR_analysis.rda", compress="xz")


#### Calculate reference points using FLRef library ----

brps_det <- lapply(sr.fits, 
                    function(x) computeFbrps(stock = stk, sr = x, proxy = 'sprx', f0.1 = TRUE, verbose = FALSE))

brps_det <- as_tibble(Reduce(rbind, lapply(1:3, function(i) cbind(as.data.frame(brps_det[[i]]@refpts[-c(1, 3,7, 12)][,1:5])[,-3], 
                                                        
                                                               model = names(sr.fits)[i]))))


brps_objs <- lapply(mixedSR_boot, 
               function(x) computeFbrps(stock = stk, sr = x, proxy = 'sprx', f0.1 = TRUE, verbose = FALSE))

brps <- as_tibble(Reduce(rbind, lapply(1:it, function(i) cbind(as.data.frame(brps_objs[[i]]@refpts[-c(1, 3,7, 12)][,1:5])[,-3], 
                                              iter = i,
                                              model = mixedSR_boot[[i]]@desc))))
#### SAVE BRPs ----
save(brps, brps_det, file="data/brps.rda", compress="xz")

#### Analyse B0 and potential Blims ---- 

# using the bootstrap of the != SR relationships

# Biomass relative to B0.
SSBrpRel2B0 <- brps %>% filter(quant == 'ssb') %>% group_by(quant, iter) %>% 
  mutate(data.rel = data/data[refpt == 'B0']) %>% group_by(refpt) %>% 
  summarize_at('data.rel', median) %>% select(refpt, data.rel) %>% mutate(data.rel = paste0(round(data.rel*100), "%"))
SSBrpRel2B0 <- t(SSBrpRel2B0)
rownames(SSBrpRel2B0) <- c('Reference point', '% of B0')


ggplot(brps %>% filter(quant == 'ssb'), aes(data, fill = refpt)) + geom_density() + 
  facet_wrap(refpt~., ncol = 2,)  + ggtitle('Boostrapped SRR: SSB') +
  scale_fill_manual(values = c(brewer.pal(8,'Dark2'), 'black'))

p1 <- ggplot(brps %>% filter(quant == 'ssb'), aes(data,fill = refpt, alpha = 0.3)) + 
  geom_density( colour = 'white') +
  scale_fill_manual(values = c(brewer.pal(8,'Dark2'), 'black')) + xlab(NULL)+ 
  theme(plot.margin=unit(c(1,1,0,1), "cm")) 

p1 <- ggplotGrob(p1) 
p2<-tableGrob(SSBrpRel2B0) 

taf.png("report/Biomas_ReferencePoints.png")
grid.arrange(p1, p2, top = "Boostrapped SRR: SSB")
dev.off()


#### B0 versus historical biomasses
B0 <- (brps %>% filter(quant == 'ssb', refpt == 'B0') %>% select(data))$data
Bmsy <- (brps %>% filter(quant == 'ssb', refpt == 'msy') %>% select(data))$data

taf.png("report/Historical_SSB_&_B0_Blim_Bmsy.png")
plot(ssb(stk)) + geom_line(linewidth = 1.5) + 
  ggtitle('Historical SSB, Bmsy, Blim_segreg and fractions of B0') +
  geom_hline(yintercept = median(B0)*c(0.10, 0.15, 0.20, 0.30, 0.5,1), col = 'blue', lty = 2, linewidth = 1) +
  geom_hline(yintercept = median(Bmsy), col = 3, lty = 2, linewidth = 1) +
  geom_hline(yintercept = Blim_segreg, col = 2, lty = 2, linewidth = 1) +
  geom_text(data = data.frame(year = c(2019, 2019, rep(2022,6)) , data = c(median(Bmsy), Blim_segreg, median(B0)*c(0.10, 0.15, 0.20, 0.30, 0.5,1))+1e4, 
                       value = c('Bmsy', 'Blim_segreg', paste0(c(0.10, 0.15, 0.20, 0.30, 0.5,1)*100, "%"))), aes(year, data, label = value))
dev.off()

hist(B0/max(ssb(stk)), main = 'B0/max(SSB)', xlab = "")

save(Blim_segreg_boot, Blim_segreg, B0, Bmsy, file="data/Biomass_refpts.rda", compress="xz")


# Compare with data driven approaches.


#### Fmsy ---- using new tool 
Fref_med <- brps %>% filter(quant == 'harvest') %>% group_by(refpt, quant) %>% 
  summarize_at('data', median)

p1 <- ggplot(brps %>% filter(quant == 'harvest'), aes(refpt, data, fill = refpt)) + geom_boxplot() + 
  ggtitle('Boostrapped SRR: F') + 
  geom_point(data = Fref_med, aes(color = refpt, size = 2)) +
  geom_hline(yintercept =  refpts['Fmsy',drop=T], size = 1.5, lty = 2) +
  scale_color_manual(values = c(brewer.pal(8, 'Dark2'), 'black')) +
  scale_fill_manual(values = c(brewer.pal(8, 'Dark2'), 'black'))


p2 <- plot(fbar(stk)) + geom_line(linewidth = 1.5) + 
  ggtitle('Historical F') +
  geom_hline(data = Fref_med, aes(yintercept = data, color = refpt), linewidth = 1.5, lty = 2)  +
  scale_color_manual(name = NULL, values = c(brewer.pal(8, 'Dark2'), 'black')) 

taf.png("report/Historical_F_&_Deter_Frfpt.png")
grid.arrange(p1, p2)
dev.off()

#### Condition BC Simulations to calculate Fp05 & Fmsy ----

# Recruitment deviations
srdevs.hist <- rlnormar1(it, 0,sdlog=sd(sr.fits$segreg@residuals), rho=ar1(c(sr.fits$segreg@residuals)), 
                    years=2005:2022)
srdevs <- rlnormar1(it, 0,sdlog=sd(sr.fits$segreg@residuals), rho=ar1(c(sr.fits$segreg@residuals)), 
                            years=iy:fy)
plot(srdevs)

# BUILD FLom
om <- FLom(stock=propagate(stk, it), refpts=refpts, model='segreg',
           params=segreg_boot, deviances=srdevs)

# SETUP om future: average of last 3 years **
om <- fwdWindow(om, end=fy)


#### Condition a Hake-like population with Blim < SSB < Btrigger ----
# Apply a large F during several years, and add some uncertainty.
om2 <- fwd(propagate(stk, it), sr=rec(stk)[, ac(2005:2022)], fbar=fbar(stk)[, ac(2005:2022)] * 3.5, deviances = srdevs.hist)
om2@desc <- 'nhke like stock with 3.5*F applied in 2005:2022 so SSB[2022] < Btrigger'


#### Condition a Hake-like population with SSB < Blim  ----
om3 <- fwd(propagate(stk, it), sr=rec(stk)[, ac(2005:2022)], fbar=fbar(stk)[, ac(2005:2022)] * 4, deviances = srdevs.hist)
om3@desc <- 'nhke like stock with 3.5*F applied in 2005:2022 so SSB[2022] < Blim'

save(om, om2, om3, refpts, file="data/data.rda", compress="xz")


plan(sequential)

