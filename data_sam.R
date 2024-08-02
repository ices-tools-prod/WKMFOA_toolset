# data.R - condition OM(s)
# WKMFOA_toolset/data.R

# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


library(icesTAF)
mkdir("data")

library(mse)
library(FLSRTMB)

source("utilities.R")

# CHOOSE number of cores for doFuture
cores <- 2

# SET future plan
if(os.unix()) {
  plan(multicore, workers=cores)
} else {
  plan(multisession, workers=cores)
}

# LOAD SAM SA results, 2022 ICES whg.27.7b-ce-k
load('boot/initial/data/whg.27.7b-ce-k.rda')

# DATA year
dy <- dims(run)$maxyear

# INTERMEDIATE year
iy <- dy + 1

# HINDCAST year
hy <- dy - dims(run)$age

# FINAL year
fy <- 2050

# NUMBER of iterations
it <- 500

# RNG seed
set.seed(987)

# - Bootstrap of stock-recruitment relationship(s) **

# FIT models
fits <- srrTMB(as.FLSRs(run, models=c("bevholt", "segreg")), 
  spr0=mean(spr0y(run)))

# PLOT
plotsrs(fits)

# BOOTSTRAP and SELECT model by largest logLik **
srpars <- bootstrapSR(run, iters=it,
  models=c("bevholt", "segreg"), method="best")

# SAVE
save(fits, srpars, file="data/bootstrap.rda", compress="xz")

# - CONSTRUCT OM

# GENERATE lognormal rec deviances, sigma and rho from SS3
srdevs <- rlnormar1(n=500, sdlog=0.5384, rho=0, years=seq(dy - 10, fy))

plot(srdevs) +
  geom_vline(xintercept=2022, linetype=2)

# BUILD FLom w/ SRR
om <- FLom(stock=propagate(run, it), refpts=refpts, model='mixedsrr',
  params=srpars, deviances=srdevs)

# HINDCAST for last 10 years /without process error
som <- fwd(om, catch=catch(om)[, ac(seq(hy, dy))],
  sr=rec(om)[, ac(seq(hy, dy))])

# HINDCAST for last 10 years /w process error and recruitment deviances
fom <- pefwd(om, catch=fbar(om)[, ac(seq(hy, dy))],
  sr=rec(om)[, ac(seq(hy, dy))], deviances = srdevs %=% 1)

com <- pefwd(om, catch=catch(om)[, ac(seq(hy, dy))],
  sr=rec(om)[, ac(seq(hy, dy))], deviances = srdevs %=% 1)

plot(FLStocks(OM=stock(om), F=stock(fom), C=stock(com)))

# SETUP om future: average of last 3 years **
om <- fwdWindow(om, end=fy)

# PROJECT forward for iy assumption (TAC) **
om <- fwd(om, catch=FLQuant(3675, dimnames=list(year=2024)))

# TODO: ADD constant F
# om <- fwd(om, fbar=expand(fbar(run)[,'2023'], year=2024))

# CREATE F and SSB deviances
sdevs <- shortcut_devs(om, Fcv=0.212, Fphi=0.423, SSBcv=0.10)

# - SAVE

save(om, sdevs, file="data/data.rda", compress="xz")

