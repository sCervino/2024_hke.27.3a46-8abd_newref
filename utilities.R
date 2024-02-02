# utilities.R - Extra functions
# WKREBUILD_toolset/utilities.R

# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


# PARALLEL setup via doFuture

if(exists("cores")) {
  plan(multisession, workers=cores)
  options(doFuture.rng.onMisuse="ignore")
}


# icesmetrics {{{

# NAME = function ~ refpt, e.g. FMSY = fbar(om) / refpts(om)$Fmsy

icesmetrics <- list(FMSY=fbar~Fmsy, SBMSY=ssb~Btrigger,
  SBPA=ssb~Bpa, SBlim=ssb~Blim)

# }}}

# WKREBUILD2 performance statistics {{{

annualstats <- list(

  # P(SB>SBlim)
  PBlim=list(~iterMeans((SB/Blim) > 1), name="P(SB>SB[lim])",
    desc="Probability that spawner biomass is above Blim"),

  # P(SB>SBtrigger)
  PBtrigger=list(~iterMeans((SB/Btrigger) > 1), name="P(SB>B[trigger])",
    desc="Probability that spawner biomass is above Btrigger"),

  # mean(C)
  C=list(~iterMeans(C), name="mean(C)",
    desc="Mean catch over years"),

  # cv(C)
  cvC=list(~sqrt(iterVars(C)) / iterMeans(C), name="cv(C)",
    desc="CV of catch over years")
)

fullstats <- list(

  # mean(C)
  C=list(~yearMeans(C), name="mean(C)",
    desc="Mean catch over years"),

  # AVVC
  # AVVC
  AAVC=list(~yearMeans(abs(C[, -1] - C[, -dim(C)[2]])/C[, -1]),
    name="AAV(C)", desc="Average annual variability in catch"),
  
  # IACC
  IACC=list(~100 * yearSums(abs(C[, -1] - C[, -dim(C)[2]]))/yearSums(C),
    name="IAC(C)",
    desc="Percentage inter-annual change in catch"),

  # P(SB < SBlim) at least once
  risk2=list(~yearMeans(iterMeans(((SB/Blim) < 1) > 0)),
    name="once(P(SB<B[limit]))",
    desc="ICES Risk 2, probability that spawner biomass is above Blim once"),

  # 1st year
  firstyear=list(~firstYear(iterMeans(SB/Blim > 1) >= 0.95),
  name="recovery", desc="First year in which P(SB/SBlim) >= 0.95")
)

# }}}

# firstyear {{{

# firstYear(iterMeans(SB/Blim > 1) >= 0.95)

firstYear <- function(x) {
  year <- as.numeric(dimnames(x)$year[match(TRUE, x)])
  return(FLQuant(year))
}
# }}}

# decisions {{{

decisions <- function(x, year=1, iter=NULL) {

  trac <- tracking(x)
  args <- args(x)

  year <- as.numeric(year)

  if(is.null(iter))
    iter <- seq(dims(x)$iter)

  .table <- function(d) {

    its <- dims(d)$iter
    dmns <- dimnames(d)

    if(its == 1) {
    data.frame(metric=dmns$metric, year=dmns$year, value=prettyNum(d))
    } else {
      data.frame(metric=dmns$metric, year=dmns$year,
        value=sprintf("%s (%s)", 
          prettyNum(apply(d, 1:5, median, na.rm=TRUE)),
          prettyNum(apply(d, 1:5, mad, na.rm=TRUE))))
    }
  }

  res <- lapply(year, function(y) {
  
    ay  <-  y
    dy <- ay - args$data_lag
    my  <- ay + args$management_lag

    dmet <- c("SB.om", "SB.obs", "SB.est")
    amet <- c("met.hcr", "decision.hcr", "fbar.hcr", "hcr", "fbar.isys", "isys",
      "fwd", "C.om")

    dout <- trac[dmet, ac(dy),,,, iter]
    aout <- trac[amet, ac(ay),,,, iter]
    mout <- trac["SB.om", ac(my),,,,iter] / trac["SB.om", ac(ay),,,,iter]
    dimnames(mout)$metric <- "diff(SB.om)"

    rbind(.table(dout), .table(aout), .table(mout))      
  })

  do.call(cbind, res)
}
# }}}


# TODO: CODE showMethod


#### modify bootstrap FLR to return one object per iteration not only pars ----

bootstrapSR_list <- function (x, iters = 100, method = c("best", "logLik", "relative"), 
          models = c("bevholt", "ricker", "segreg"), verbose = TRUE, 
          ...) 
{
  spr0x <- yearMeans(spr0y(x))
  x <- as.FLSR(x)
  id <- matrix(sample(seq(dim(x)[2]), dim(x)[2] * iters, replace = TRUE), 
               nrow = dim(x)[2], ncol = iters)
  models <- match.arg(models, several.ok = TRUE)
  mod <- list(bevholt = bevholtSV, ricker = rickerSV, segreg = segreg)[models]
  method <- match.arg(method)
  p <- progressor(along = seq(iters), offset = 1)
  if (nbrOfWorkers() > 1) 
    message("Running on ", nbrOfWorkers(), " nodes.")
  res <- foreach(i = seq(iters), .options.future = list(globals = structure(TRUE, 
                                                                            add = c("x", "id", "mod", "spr0x", "method", "models", 
                                                                                    "logLik"), packages = c("FLCore"), seed = TRUE))) %dofuture% 
    {
      y <- x
      rec(y) <- rec(y)[, id[, i]]
      ssb(y) <- ssb(y)[, id[, i]]
      fits <- lapply(mod, function(m) {
        model(y) <- m
        srrTMB(y, spr0 = spr0x, verbose = FALSE, ...)
      })
      if (verbose) 
        p(message = sprintf(paste0("[", i, "]")))
      llkhds <- unlist(lapply(fits, "logLik"))
      aics <- unlist(lapply(fits, "AIC"))
      if (method == "best") {
        best <- which.max(llkhds)
      }
      else if (method == "loglik") {
        probs <- llkhds/sum(llkhds)
        u <- runif(1, 0, 1)
        best <- which(u <= cumsum(probs))[1]
      }
      else if (method == "relative") {
        daic <- aics - min(aics)
        relkhd <- exp(-0.5 * daic)
        aicw <- relkhd/sum(relkhd)
        u <- runif(1, 0, 1)
        best <- which(u <= cumsum(aicw))[1]
      }
      m <- match(models[best], c("bevholt", "ricker", "segreg"))
      fit <- fits[[best]]
      fit@desc <- c("bevholt", "ricker", "segreg")[m]
      list(fit = fit,
      rbind(params(fit), FLPar(m = m, spr0 = spr0x, logLik = llkhds[best]), 
             FLPar(attr(fit, "SV"))))    
    }
  # browser()
  res1 <- lapply(1:length(res), function(i) res[[i]][[1]])
  res2 <- lapply(1:length(res), function(i) res[[i]][[2]])
  
  res <- list(fits = res1, params = Reduce(cbind, res2))
  
  return(res)
  
}


# Function to create a table with the model settings.
scenario_table <- function(settings, scenario_name, description){
  nrow <- length(settings)
  res <- data.frame(Component = rep(NA, nrow), Process = rep(NA, nrow),
                    Model = rep(NA, nrow),
                    `Years used` = rep(NA, nrow),
                    `Source`= rep(NA, nrow),
                    `Rationale` = rep(NA, nrow))
  
  res[,1] <- c("Initial population", "Initial population", "Recruitment", "Recruitment",        
               "Recruitment", "Biological parameters", "Biological parameters", "Biological parameters",                 
               "Fishery parameters", "Fishery parameters", "Short cut approach",  "Short cut approach", "Short cut approach", "Short cut approach", "Harvest_Control_Rules", "Implmentation error") 
  res[,2] <- c("Model", "Uncertainty", "Functional form", "Parametric uncertainty",        
               "Process error", "Natural mortality", "Maturity", "Weight at age",                 
               "Selection pattern", "Discards", "Recruitment", "SSB_cv",  "F_cv", "Fphi", "Harvest control rules", "Implementation Error") 
  for(r in 1:dim(res)[1]) res[r,-(1:2)] <- model_settings[[r]]
  
  
  res <- cbind(Space = c(rep('Operating model', 10), rep('Management procedure', nrow-10)), res)
  
  res <- flextable(res) %>% merge_v(j = 1:2) %>%  
    add_header_row(values = c("Description", description),  colwidths = c(1,dim(res)[2]-1)) %>% 
    add_header_row(values = c("Scenario Name", scenario_name),  colwidths = c(1,dim(res)[2]-1)) %>% 
    hline(c(2,5,8,10,11,14,15,16)) %>% 
    width(c(2, 3, 4,6, 7), c(2,5,5,5,10), unit = 'cm')
  
  return(res)
  
}


# hindcast_srdevs: function to introduce uncertainty in the initial population
hindcast_srdevs <- function(stock, srdevs, start, end, process_error = FALSE){ 
  # The process error is the right one for the first year 'hy', but afterwards it needs to be 
  # recalculated as the numbers at ages in years are corrected, so the population in year
  # y included process error.
  # perry <- perr
  
  start <- hy
  end   <- dy
  
  y0 <- as.numeric(dimnames(stock@m)[[2]][1])
  
  A <- dim(stock@m)[1]
  
  #  * Add uncertainty in recruitment
  for(y in hy:dy) rec(stock)[, ac(y)] <- rec(stock)[1, ac(y)] * srdevs[, ac(y)]
  
  # -- SIMPLE hindcast
  dhind <- fwd(stock, sr=rec(stock), catch=catch(stock)[, ac(hy:dy)])
  
  
  if(process_error == FALSE) return(dhind)
  
  # -- IF process_error => correct it.
  
  #  * Correcting for process error
  #  First the process error is calculated deterministically using the SAM (or other) estimate 
  #  and then it is used to correct the stochastic projection 
  
  # COMPUTE process error, e = y/(x exp(-z)) all ages except recruitment.
  perr        <- stock.n(stock)[-1, ac(hy:dy)]/(stock.n(stock)[-A, ac((hy-1):(dy-1))] * 
                                                  exp(-z(stock)[-A, ac((hy-1):(dy-1))]))
  perr[A, ]   <- stock.n(stock)[A, ac(hy:dy)]/(quantSums(stock.n(stock)[(A-1):A, ac((hy-1):(dy-1))] * 
                                                           exp(-z(stock)[(A-1):A, ac((hy-1):(dy-1))])))
  # The process error is the right one for the first year 'hy', but afterwards it needs to be 
  # recalculated as the numbers at ages in years are corrected, so the population in year
  # y included process error.
  # perry <- perr
  
  dhind_perr <- stock
  
  perry <- perr
  for(y in seq(hy,dy)) {
    
    # FWD(catch[y])
    dhind_perr <- fwd(dhind_perr, sr=rec(stk), f=fbar(stk)[, ac(y)])
    
    stock.n(dhind_perr)[-1, ac(y)] <- stock.n(dhind_perr)[-1, ac(y)]*perry[, ac(y)] 
    
    # perr for next year
    if(y < dy){
      # all ages - recruitment
      perry[,ac(y+1)] <- stock.n(stk)[-1, ac(y+1)]/(stock.n(dhind_perr)[-A, ac(y)] * exp(-z(dhind_perr)[-A, ac(y)]))
      # correct plusgroup
      perry[A,ac(y+1)] <- stock.n(stk)[A, ac(y+1)]/(quantSums(stock.n(dhind_perr)[(A-1):A, ac(y)] *
                                                                exp(-z(dhind_perr)[(A-1):A, ac(y)])))
    }
  }
  
  # Stochastic hindcast applying perry to account for the process error.
  for(y in seq(hy,dy)) {
    stock <- fwd(stock, sr=window(rec(stk), 1999,fy), f=fbar(stk)[, ac(y)], deviances = srdevs[, ac(y)])
    
    stock.n(stock)[-1, ac(y)]     <- stock.n(stock)[-1, ac(y)]*perry[, ac(y)] 
  }
  return(stock)
}
