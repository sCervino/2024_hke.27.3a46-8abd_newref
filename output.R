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
 
save(yields, Fmsy, pssb, Fp.05, msy99, file="output/output.rda", compress="xz")
