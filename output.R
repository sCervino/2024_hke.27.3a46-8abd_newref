# output.R - Performance evaluation and output tables
# WKREBUILD_toolset/output.R

# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
# modified: Dorleta Garcia 2023-11-13

# Distributed under the terms of the EUPL-1.2

library(tidyverse)
library(icesTAF)
library(mse)
library(gridExtra)

source("utilities.R")
mkdir("output")

# LOAD model.R outputs

load("model/FcteSims.rda")
load("data/Biomass_refpts.rda")


it <- 500

#### Identify [Btrigger, Fmsy] pairs ----
# Can we select [Btrigger, Fmsy] pairs in an integrated way?
# We can search a GRID [Btrigger, Fmsy], then we need to select based on
#  - Precautionary
#  - Biggest catches.
# Could we reduce SSB to levels below Btrigger to select a Btrigger with which the AR is able to 
# recover SSB above Btrigger?
# Do we need to test the AR for robustness?


#### Select Blim ----
# In this case the breakpoint of the segmented regression is well defined.
# Segmented regression Blim corresponds ~ 12% Blim.
Blim_segreg
median(Blim_segreg_boot)

Blim <- Blim_segreg


#### Fp.05 without AR ----

# p(ssb < 0.05) over time and F levels.
pssb <- as_tibble(Reduce('rbind', lapply(1:201, function(i){
                             res <- as_tibble(as.data.frame(ssb(FcteSims[[i]]@om@stock))[,c(2,6:7)]) %>% 
                                      mutate(lowerBlim = (data < Blim)) %>% group_by(year) %>% 
                                      summarize_at('lowerBlim', sum) %>% mutate(lowerBlim = lowerBlim/it)
                             res <- cbind(F = seq(0,2, 0.01)[i], res)
                              }))) %>% filter(year > 2039) %>% group_by(F) %>% summarize_at('lowerBlim', mean)

Fp.05_df <- pssb[which(pssb$lowerBlim < 0.05)[length(which(pssb$lowerBlim < 0.05))],]
Fp.05    <- as.numeric(pssb[which(pssb$lowerBlim < 0.05)[length(which(pssb$lowerBlim < 0.05))],'F'])

taf.png(file="report/pssb_bxp.png")
ggplot(pssb, aes(F, lowerBlim)) + geom_line(linewidth = 1) + ggtitle("p(ssb < Blim) in equilibrium") +
  ylab('%') + geom_hline(yintercept = 0.05, col = 2) + 
  geom_point(data = Fp.05_df, aes(F, lowerBlim), size = 3, col = 'blue', pch = 18) +
  geom_text(data = Fp.05_df, aes(F*1.25, lowerBlim*0.9, label= Fp.05)) 
dev.off()

#### Fmsy ----
# yield over time, iterations and F levels.
yields <- as_tibble(Reduce('rbind', 
                           lapply(1:101, function(i) cbind(F = seq(0,2, 0.01)[i],as.data.frame(FcteSims[[i]]@om@stock@catch)[,c(2,6:7)])))) %>% 
  mutate(F = as.factor(F))

# COMPUTE yearly performance statistics
yields.max <- yields %>% filter(year %in% ac(2040:2050)) %>% group_by(F, iter) %>% summarize_at('data', mean) %>% 
                  group_by(F) %>% summarize_at('data', median)
Fmsy_df <- yields.max[which.max(yields.max$data),]
Fmsy <- as.numeric(as.character(Fmsy_df$F))
 msy <- as.numeric(as.character(Fmsy_df$data))

 taf.png(file="report/yield_equilibrium_bxp.png")
 p1 <- ggplot(yields %>% filter(year %in% ac(2045), F %in% as.character(seq(0,0.8, 0.01))), aes(F, data)) + 
   geom_boxplot() + ggtitle("Yield in equilibrium") + ylab('tonnes') + 
   geom_point(data = Fmsy_df, aes(F, data), size = 3, col = 'blue', pch = 18) +
   geom_text(data = Fmsy_df, aes(as.character(F), data*0.9, label= Fmsy)) 
 p1
 dev.off()
 
 # auxiliar plot to look at stability of the output
aux <- yields %>% filter(year %in% ac(2040:2050)) %>% group_by(F, year) %>% summarize_at('data', median) 
 ggplot(aux %>% filter(F %in% seq(0.2,0.5, 0.01)), aes(year, data, group = F, colour = F)) + geom_line()
 
 
# How flat is the curve?
msy99 <- as.numeric(as.character(yields.max[which(yields.max$data > max(yields.max$data)*0.99), ]$F))

taf.png(file="report/yield_equilibrium_line.png")
p2 <- ggplot(yields.max, aes(as.numeric(as.character(F)), data)) + geom_line(linewidth = 1) + xlab('F') + ylab('t') + 
  ggtitle('Long term median yield') + geom_vline(xintercept = msy99[c(1, length(msy99))], col = 2) +
  geom_vline(xintercept = Fmsy, col = 3, linewidth = 2)
p2
dev.off()

# Biomass 
ssbs <- as_tibble(Reduce('rbind', 
                           lapply(1:101, function(i) cbind(F = seq(0,1, 0.01)[i],as.data.frame(ssb(FcteSims[[i]]@om@stock))[,c(2,6:7)])))) %>% 
  mutate(F = as.factor(F)) %>% group_by(F, year) %>%  summarize_at('data', median)
ssb99msy <- ssbs[as.numeric(as.character(ssbs$F)) %in% msy99,]

taf.png(file="report/ssb_equilibrium_bxp.png")
p1 <- ggplot(ssbs %>% filter(year %in% ac(2045), F %in% as.character(seq(0,0.8, 0.01))), aes(F, data)) + 
  geom_boxplot() + ggtitle("SSB in equilibrium") + ylab('tonnes') + 
  geom_point(data = Fmsy_df, aes(F, data), size = 3, col = 'blue', pch = 18) +
  geom_text(data = Fmsy_df, aes(as.character(F), data*0.9, label= Fmsy)) 
p1
dev.off()

ssb99 <- ssbs[as.numeric(as.character(ssbs$F)) %in% msy99,]

taf.png(file="report/ssb_equilibrium_line.png")
p2 <- ggplot(yields.max, aes(as.numeric(as.character(F)), data)) + geom_line(linewidth = 1) + xlab('F') + ylab('t') + 
  ggtitle('Long term median yield') + geom_vline(xintercept = msy99[c(1, length(msy99))], col = 2) +
  geom_vline(xintercept = Fmsy, col = 3, linewidth = 2)
p2
dev.off()


# [Bmsy, Fmsy] area

#### AR[OPRs] Robust to minimum assessment uncertainty?  ----
# Here we use a short-cut MSE with the same implementation error as in
# eqSim by including the 1 year time lag.


#### AR[OPRs] Robust to observed low productivity in the historical period?  ----


#### AR[OPRs] Btrigger well defined to avoid Blim ----
# Applying a constant high F reduce the SSB to a point in (Bpa, MSYBtrigger)
# and analyse the capacity of the AR[OPRs] to rebuild the stock above Bpa.


### AR[OPRs] Btrigger, Blim well defined to recover the stock ----




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


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






save(yields, Fmsy, pssb, Fp.05, msy99, file="output/output.rda", compress="xz")
