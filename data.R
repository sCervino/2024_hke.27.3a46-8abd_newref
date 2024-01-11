#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# data.R - 
# 2024_2024_hke.27.3a46-8abd_newref/data.R
#
# Dorleta Garcia
# 2024/01/05
#
# Input: R file with the FLStock object
# Output:
#      - "data/SR_analysis.rda": SR relationships
#      - "data/brps.rda": Reference points from FLref package.
#      - "data/Biomass_refpts.rda"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### Install FLR packages from r-universe ----
# The rest of the packages are available in CRAN
# install.packages(c("FLCore",  "mse", "FLSRTMB", "FLRef"), repos = c('https://flr.r-universe.dev', 'https://cloud.r-project.org'))

#### Load libraries ----
# FLR related packages
library(FLRef)
library(FLSRTMB)
library(mse)
# ICES packages
library(icesTAF)
# tidyverse packages
library(tidyverse)
# ggplot related packages
library(RColorBrewer)
library(ggplot2)
library(gridExtra) # arrange multiple plots in a page
library(ggpubr)    #  provides some easy-to-use functions for creating and customizing 'ggplot2'
library(ggExtra)   # Marginal histograms
library(ggthemes)  # Themes for ggplot

# Some extra utilities:bootstrapSR_list
source('utilities.R')
#### Load the FLStock object ---- 
load('bootstrap/data/wgbie2023_nhke_FLStock_csim.RData')
# Rename the object
stk <- hke.stk


#### Current refpts ----
refpts <- FLPar(Btrigger = 78405, 
                Fmsy = 0.24, 
                Blim = 61563, 
                Bpa = 78405, 
                Flim = 0.73, 
                Fpa = 0.54, 
                lFmsy = 0.147, 
                uFmsy = 0.37, 
                F05 = 0.54, 
                F05noAR = NA)

#### Some input data ----

# INTERMEDIATE year
iy <- stk@range['maxyear'] + 1
# FIRST assessment year
y0 <- stk@range['minyear'] 
# FINAL year
fy <- 2100

# Years to fit Stock-Recruitment relationship.
# NEED TO BE UPDATED: Selected after the first analysis!!!
recy <- 1978:2020

# DATA year
dy <- stk@range['maxyear']
# year names
yrs <- dimnames(stk)$year

# NUMBER of iterations
it <- 1000

set.seed(527)


# Create the data directory to: 
mkdir("data")
# Create report directory to save plots for the report.
mkdir('report')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### *** 1. Stock-recruitment relationship(s) ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### ** 1.a Select the time series to be used to calculate reference points ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### Is the breakpoint sensitive to the length of the time series? ----
segreg_finalYr <- matrix(NA,2, 10, dimnames = list(c('a', 'b'), (dy-10+1):dy))

for(y in (dy-10+1):dy){
  
  temp <- srrTMB(as.FLSRs(window(stk, y0, y), models=c("segreg")), spr0=mean(spr0y(stk)))
  
  segreg_finalYr[, ac(y)] <- c(temp[[1]]@params)[1:2]
}

taf.png("report/recruitment_segreg_Finalyr.png")
plot(segreg_finalYr[1,], segreg_finalYr[2,], main = 'SR params, time series up to year "y"',
     xlab = "a", ylab = "b (breakpoint)", pch = 16, col = 2)
text(segreg_finalYr[1,]+(range(segreg_finalYr[1,])[2] -range(segreg_finalYr[1,])[1])*1e-5, 
     segreg_finalYr[2,]+(range(segreg_finalYr[2,])[2] -range(segreg_finalYr[2,])[1])*1e-5, 
     (dy-10+1):dy, cex = 0.8)
dev.off()
# The breakpoint is very sensitive to the last three year points, with much more 
# consistency in previous years => Remove last 3 years from the fit.


#### Is the breakpoint sensitive to single (rec, ssb) pairs? ----
segreg_1yrOut <- matrix(NA,2, dim(stk)[2], dimnames = list(c('a', 'b'), dimnames(stk)[[2]]))

for(i in 1:dim(stk)[2]){
  if(i == 1)               subs <- 2:dim(stk)[2]
  if(i %in% 2:dim(stk)[2]) subs <- c(1:(i-1), (i+1):dim(stk)[2])
  if(i == dim(stk)[2])     subs <- 1:(dim(stk)[2]-1)
  
  temp <- srrTMB(as.FLSRs(stk[,subs], models=c("segreg")), spr0=mean(spr0y(stk)))
  
  segreg_1yrOut[, i] <- c(temp[[1]]@params)[1:2]
}

taf.png("report/recruitment_segreg_1yrOut.png")
plot(segreg_1yrOut[1,], segreg_1yrOut[2,], main = 'SR params, 1 year out',
     xlab = "a", ylab = "b (breakpoint)", pch = 15, col = 2)
text(segreg_1yrOut[1,]+0.01, segreg_1yrOut[2,]+250, yrs, cex = 0.8)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### ** 1.b Fit SR models ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### Deterministic fit ----
sr.fits     <- srrTMB(as.FLSRs(window(stk, recy[1], recy[length(recy)]), models=c("segreg", "ricker", "bevholt")), spr0=mean(spr0y(stk)))
Blim_segreg <- sr.fits$segreg@params$b[drop=T]

# PLOT
taf.png("report/recruitment_fits.png")
plotsrs(sr.fits) + ggtitle('Stock recruiment relationships')
dev.off()

#### Boostrap fit of the segmented regression model ----
segreg_boot <- bootstrapSR(window(stk, recy[1], recy[length(recy)]), iters=it, models=c( "segreg"))
Blim_segreg_boot <- c(segreg_boot['b',])

taf.png("report/recruitment_segreg_boot.png")
plot(segreg_boot[c('a', 'b', 'sigmaR', 'R0', 'rho', 'B0'),]) + ggtitle('Bootstrapped segmented regression parameters')
dev.off()


taf.png("report/recruitment_boostrap_segreg_params.png")
segreg_boot_df <- as.data.frame(t(segreg_boot[drop=T]))[,c('a', 'b', 'logLik', 'sigmaR', 'R0', 'rho', 'B0')]
plot1 <- ggplot(segreg_boot_df, aes(x = a, y = b)) + 
  geom_point(color = "#619CFF") + 
  stat_smooth(method = "lm", fullrange = TRUE) +
  geom_rug() + 
  scale_y_continuous(name = "b (the breakpoint)",  limits = range(segreg_boot_df$b) ) + 
  scale_x_continuous(name = "a (the slope)",  limits = c(10, 40), expand = c(0, 0)) + 
 # scale_x_continuous(name = "a (the slope)",  expand = c(0, 0)) + 
  theme_pubr()  + ggtitle("Boostraped Segmented regression parameters")

ggMarginal(plot1, fill = "#619CFF") 
dev.off()


# stats_brkpt: matrix to store the Blim calculated with different methods.
stats_brkpt <- matrix(NA, 2,4, dimnames = list(c('bootstrap', '1yrOut'), c('mean', 'median', 'cv', 'sd')))

stats_brkpt[1,] <- c(mean(segreg_boot[2,]), median(c(segreg_boot[2,])),
                     sd(segreg_boot[2,])/mean(segreg_boot[2,]), sd(segreg_boot[2,]))
stats_brkpt[2,] <- c(mean(segreg_1yrOut[2,]), median(c(segreg_1yrOut[2,])),
                     sd(segreg_1yrOut[2,])/mean(segreg_1yrOut[2,]), sd(segreg_1yrOut[2,]))


#### Bootstrap of 3 SR model and  selection of the model by largest logLik ----
mixedSR_boot <- bootstrapSR_list(stk, iters=it, models=c("ricker", "bevholt", 'segreg'), method="best")
table(sapply(mixedSR_boot, function(x) slot(x, 'desc')))

#### Save SR fits ----
save(sr.fits, segreg_boot, mixedSR_boot, stats_brkpt, file="data/SR_analysis.rda", compress="xz")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 2. Calculate reference points using FLRef library ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### Compute **deterministic reference points** from the 3 deterministic SR relationships ----
aux <- lapply(sr.fits, 
                function(x) computeFbrps(stock = stk, sr = x, proxy = 'sprx', f0.1 = TRUE, verbose = FALSE))

# Reshape the reference points calculated in previous step in a data frame format.
brps_det <- as_tibble(Reduce(rbind, lapply(1:3, function(i) cbind(as.data.frame(aux[[i]]@refpts[-c(1, 3,7, 12)][,1:5])[,-3], 
                                                        
                                                               model = names(sr.fits)[i]))))

#### Compute **stochastic reference points** from the bootstrapped mixed stock-recruitment relationship ----
brps_objs <- lapply(mixedSR_boot, 
               function(x) computeFbrps(stock = stk, sr = x, proxy = 'sprx', f0.1 = TRUE, verbose = FALSE))

brps <- as_tibble(Reduce(rbind, lapply(1:it, function(i) cbind(as.data.frame(brps_objs[[i]]@refpts[-c(1, 3,7, 12)][,1:5])[,-3], 
                                              iter = i,
                                              model = mixedSR_boot[[i]]@desc))))
#### Save BRPs ----
save(brps, brps_det, file="data/brps.rda", compress="xz")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 3.  Analyse B0 and other SSB reference points ---- 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# using the bootstrap of the != SR relationships

#### Distribution of SSB for different reference points ----
# Biomass relative to B0.
SSBrpRel2B0 <- brps %>% filter(quant == 'ssb') %>% group_by(quant, iter) %>% 
  mutate(data.rel = data/data[refpt == 'B0']) %>% group_by(refpt) %>% 
  summarize_at('data.rel', median) %>% select(refpt, data.rel) %>% mutate(data.rel = paste0(round(data.rel*100), "%"))
SSBrpRel2B0 <- t(SSBrpRel2B0)
rownames(SSBrpRel2B0) <- c('Reference point', '% of B0')

p1 <- ggplot(brps %>% filter(quant == 'ssb'), aes(data,fill = refpt, alpha = 0.3)) + 
  geom_density( colour = 'white') +
  scale_fill_manual(values = c(brewer.pal(8,'Dark2'), 'black')) + xlab(NULL)+ 
  theme(plot.margin=unit(c(1,1,0,1), "cm")) + theme_hc()

p1 <- ggplotGrob(p1) 
p2<-tableGrob(SSBrpRel2B0) 

taf.png("report/Biomas_ReferencePoints.png")
grid.arrange(p1, p2, top = "Boostrapped SRR: SSB")
dev.off()


#### B0 versus historical biomasses ----
B0 <- (brps %>% filter(quant == 'ssb', refpt == 'B0') %>% select(data))$data
Bmsy <- (brps %>% filter(quant == 'ssb', refpt == 'msy') %>% select(data))$data

taf.png("report/Historical_SSB_&_B0_Blim_Bmsy.png")
plot(ssb(stk)) + geom_line(linewidth = 1.5) + 
  ggtitle('Historical SSB, Bmsy, Blim_segreg and fractions of B0') +
  geom_hline(yintercept = median(B0)*c(0.10, 0.15, 0.20, 0.30, 0.5,1), col = 'blue', lty = 2, linewidth = 1) +
  geom_hline(yintercept = median(Bmsy), col = 3, lty = 2, linewidth = 1) +
  geom_hline(yintercept = Blim_segreg, col = 2, lty = 2, linewidth = 1) +
  geom_text(data = data.frame(year = c(2019, 2019, rep(2022,6)) , data = c(median(Bmsy), Blim_segreg, median(B0)*c(0.10, 0.15, 0.20, 0.30, 0.5,1))+1e4, 
                       value = c('Bmsy', 'Blim_segreg', paste0(c(0.10, 0.15, 0.20, 0.30, 0.5,1)*100, "%"))), aes(year, data, label = value)) +
  theme_hc()
dev.off()

#### B0 versus SSBmax ----
taf.png("report/Historical_SSB_&_B0_Blim_Bmsy.png")
hist(B0/max(ssb(stk)), main = 'B0/max(SSB)', xlab = "")
dev.off()

save(Blim_segreg_boot, Blim_segreg, B0, Bmsy, file="data/Biomass_refpts.rda", compress="xz")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 4.  Calculate Blim using data driven approaches ---- 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 5. YPR Fishing mortality reference points ---- 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

taf.png("report/Historical_F_&_YPR_DetFrp.png")
grid.arrange(p1, p2)
dev.off()

