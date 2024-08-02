# utilities.R - Extra functions
# WKREBUILD_toolset/utilities.R

# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

# SET doFuture to ignore warning on RNG
options(doFuture.rng.onMisuse="ignore")


# icesmetrics {{{

# NAME = function ~ refpt, e.g. FMSY = fbar(om) / refpts(om)$Fmsy

icesmetrics <- list(FMSY=fbar~Fmsy, SBMSY=ssb~Btrigger,
  SBPA=ssb~Bpa, SBlim=ssb~Blim)

# }}}

# WKREBUILD2 performance statistics {{{

annualstats <- list(

  # P(SB>SBlim)
  PBlim=list(~iterMeans((SB/Blim) > 0), name="P(SB>SB[lim])",
    desc="Probability that spawner biomass is above Blim"),

  # P(SB>SBtrigger)
  PBtrigger=list(~iterMeans((SB/Btrigger) > 1), name="P(SB>B[trigger])",
    desc="Probability that spawner biomass is above Btrigger"),

  # mean(C)
  Cy=list(~iterMeans(C), name="mean(C)",
    desc="Mean catch per year"),

  # cv(C)
  cvCy=list(~sqrt(iterVars(C)) / iterMeans(C), name="cv(C)",
    desc="CV of catch per year")
)

fullstats <- list(

  # mean(C)
  C=list(~yearMeans(C), name="mean(C)",
    desc="Mean catch over years"),

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
    desc="ICES Risk 2, probability that SSB is below Blim once")
)

# }}}

# fwd function with process error added {{{

pefwd <- function(object, sr, catch=NULL, fbar = NULL,
  deviances = rec(object) %=% 1) {

  # CHECK inputs
  if(!class(object)== "FLom") {
    stop("make sure that class of object is om")
  }
  
  if(is.null(catch) & is.null(fbar)){
    stop("'catch' or 'fbar' missing")
  }

  stock <- object@stock

  # DIMS
  nages <- dims(stock)$age
  dy <- dims(stock)$maxyear

  if(is.null(fbar))
    hiny <- dims(catch)$minyear
  else
    hiny <- dims(fbar)$minyear

  # COMPUTE process error, e = y/(x exp(-z)) all ages except recruitment.
  perr <- stock.n(stock)[-1, ac(hiny:dy)] /
    (stock.n(stock)[-nages, ac(hiny:dy - 1)] *
     exp(-z(stock)[-nages, ac(hiny:dy - 1)]))
  
#  perr <- log(stock.n(stock)[-1, ac(hiny:dy)]) -
#    log(stock.n(stock)[-nages, ac(hiny:dy - 1)]) +
#    z(stock)[-nages, ac(hiny:dy - 1)]
  # nt+1 = nt * exp(-f-m+e)
  #log(nt+1) =  log(nt) -f -m +e
  # e = log(nt+1) - lo
  
  # PLUSGROUP
  perr[ac(dims(stock)$max), ] <- stock.n(stock)[nages, ac(hiny:dy)] /
    (quantSums(stock.n(stock)[(nages - 1):nages, ac(hiny:dy - 1)] *
    exp(-z(stock)[(nages - 1):nages, ac(hiny:dy - 1)])))
      
  dhind_perr <- object
  perry <- perr

  for(Y in hiny:dy) {

    # FWD(catch[y])
    if(is.null(fbar))
      dhind_perr <- fwd(dhind_perr, sr=sr, catch=catch[, ac(Y)],
        deviances=rec(object) %=% 1)
    else
      dhind_perr <- fwd(dhind_perr, sr=sr, fbar=fbar[, ac(Y)],
        deviances=rec(object) %=% 1)

    stock.n(dhind_perr)[-1, ac(Y)] <- stock.n(dhind_perr)[-1, ac(Y)] *
      perry[, ac(Y)] 

    # perr for next year
    if(Y < dy){
      
      # all ages minus recruitment
      perry[,ac(Y+1)] <- stock.n(stock)[-1, ac(Y+1)] / 
        (stock.n(dhind_perr)[-nages, ac(Y)] *
         exp(-z(dhind_perr@stock)[-nages, ac(Y)]))
      
      # correct plusgroup
      perry[ac(dims(stock)$max), ac(Y+1)] <- stock.n(stock)[nages, ac(Y+1)] /
        (quantSums(stock.n(dhind_perr)[(nages-1):nages, ac(Y)] *
        exp(-z(dhind_perr@stock)[(nages-1):nages, ac(Y)])))
    }
  }

  # Stochastic recruitment + Process error with lognormal deviations
  hind_perr_log <- om@stock

  for(YY in seq(hiny,dy)) {
      
    # FWD(catch[y])
    hind_perr_log <- fwd(hind_perr_log, sr=sr[, ac(hiny:dy)],
      f=fbar(stock)[, ac(YY)], deviances = deviances[, ac(YY)])
    
    hind_perr_log@stock.n[-1, ac(YY)] <- stock.n(hind_perr_log)[-1, ac(YY)] *
      perry[, ac(YY)] 
  }

  om@stock <- hind_perr_log

  return(om)
}

# }}}

# fwd function with process error added {{{

pefwd <- function(object, sr, catch=NULL, fbar = NULL,
  deviances = rec(object) %=% 1) {

  # CHECK inputs
  if(!class(object)== "FLom") {
    stop("make sure that class of object is om")
  }
  
  if(is.null(catch) & is.null(fbar)){
    stop("'catch' or 'fbar' missing")
  }

  stock <- stock(object)

  # DIMS
  nages <- dims(stock)$age
  dy <- dims(stock)$maxyear

  if(is.null(fbar))
    hiny <- dims(catch)$minyear
  else
    hiny <- dims(fbar)$minyear

  # COMPUTE process error, e = y/(x exp(-z)) all ages except recruitment.
  perr <- stock.n(stock)[-1, ac(hiny:dy)] /
    (stock.n(stock)[-nages, ac(hiny:dy-1)] *
     exp(-z(stock)[-nages, ac(hiny:dy-1)]))
  
  # PLUSGROUP
  perr[ac(dims(stock)$max), ] <- stock.n(stock)[nages, ac(hiny:dy)] /
    (quantSums(stock.n(stock)[(nages-1):nages, ac(hiny:dy -1)] *
    exp(-z(stock)[(nages-1):nages, ac(hiny:dy -1)])))

  # BUG:
  old <- stock
  browser()
      
  for(Y in hiny:dy) {

    # FWD(catch[y])
    if(is.null(fbar))
      stock <- fwd(stock, sr=sr, catch=catch[, ac(Y)],
        deviances=deviances[, ac(Y)])
    else
      stock <- fwd(stock, sr=sr, fbar=fbar[, ac(Y)],
        deviances=deviances[, ac(Y)])

    stock.n(stock)[-1, ac(Y + 1)] <- stock.n(stock)[-1, ac(Y + 1)] *
      perr[, ac(Y)]
  }

  stock(om) <- stock

  return(om)
}

# }}}
