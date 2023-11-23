# output.R - Performance evaluation and output tables
# WKREBUILD_toolset/output.R

# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
# modified: Dorleta Garcia 2023-11-13

# Distributed under the terms of the EUPL-1.2

library(tidyverse)
library(icesTAF)
library(mse)

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
                             res <- as_tibble(as.data.frame(ssb(plans[[i]]@om@stock))[,c(2,6:7)]) %>% 
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
                           lapply(1:201, function(i) cbind(F = seq(0,2, 0.01)[i],as.data.frame(plans[[i]]@om@stock@catch)[,c(2,6:7)])))) %>% 
  mutate(F = as.factor(F))

# COMPUTE yearly performance statistics, ADD F0 as reference
yields.max <- yields %>% filter(year %in% ac(2040:2050)) %>% group_by(F, iter) %>% summarize_at('data', mean) %>% 
                  group_by(F) %>% summarize_at('data', median)
Fmsy_df <- yields.max[which.max(yields.max$data),]
Fmsy <- as.numeric(as.character(Fmsy_df$F))
 msy <- as.numeric(as.character(Fmsy_df$data))

 taf.png(file="report/yield_equilibrium_bxp.png")
 ggplot(yields %>% filter(year %in% ac(2045), F %in% as.character(seq(0,0.8, 0.01))), aes(F, data)) + 
   geom_boxplot() + ggtitle("Yield in equilibrium") + ylab('tonnes') + 
   geom_point(data = Fmsy_df, aes(F, data), size = 3, col = 'blue', pch = 18) +
   geom_text(data = Fmsy_df, aes(as.character(F), data*0.9, label= Fmsy)) 
 dev.off()
 
 
 

perf_year <- performance(c(plans, F0=runf0), statistics=annualstats,
  years=2023:2041)

# COMPUTE performance statistics by periods

perf <- performance(c(plans, F0=runf0), statistics=fullstats,
  years=list(short=2024:2028, medium=2028:2034, long=2034:2041, all=2024:2041))


# --- TABLES

tables <- list()

# WHEN does stock recover (P(SB>Blim) >= 95%) by mp?

tables$recovery <- perf_year[statistic == "PBlim" & data > 0.95, .SD[1], 
  by=mp][order(year),]

perf[statistic=='firstyear' & year == 'all', data, by=.(mp)][order(data),]

# WHEN is P(B>Btrigger) > 50% by mp?

tables$status <- perf_year[statistic == "PBtrigger" & data > 0.50, .SD[1],
  by=mp]

# CREATE table of catch by mp and year (mp ~ year | C)

tables$catch_mp <- dcast(perf_year[statistic == 'C', .(data=mean(data)),
  by=.(year, mp, name, desc)], mp ~ year, value.var='data')

# CREATE table of all statistics by mp and year (statistic + mp ~ year)

tables$stats_mp <- dcast(perf_year[, .(data=mean(data)),
  by=.(year, mp, name, desc)], name + mp ~ year, value.var='data')


# --- TRACK decisions (EXAMPLES)

# TRACK decision for a single iter and year
decisions(advice, year=2024, iter=1)

# TRACk decisions for multiple years and all iters
decisions(advice, year=2024:2025)

# SAVE

save(perf_year, perf, tables, file="output/output.rda", compress="xz")
