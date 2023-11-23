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
      fit
    #   rbind(params(fit), FLPar(m = m, spr0 = spr0x, logLik = llkhds[best]), 
    #         FLPar(attr(fit, "SV")))
     }
#  params <- Reduce(combine, res)
  return(res)
}