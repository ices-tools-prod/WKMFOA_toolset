# trove.R - DESC
# WKREBUILD_toolset/trove.R

# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


# --- data.R

# - Stock-recruitment relationship(s)

# SET future recruitments as random walk from last estimated year

sr(om) <- rwalk(window(rec(om), end=2022), end=fy, sd=0.02, delta=0)

# with negative drift

sr(om) <- rwalk(window(rec(om), end=2022), end=fy, sd=0.02, delta=-0.01)


# - SETUP om future

# use mean of last 5 years but resample 'wt' slots

om <- fwdWindow(window(om, end=2023), end=fy, nsq=5, fun=c("mean", wt="sample"))

# resample 10 years for wt and mat, 5 for catch.sel

om <- fwdWindow(om, end=fy, nsq=3,
  fun=c(wt='sample', mat='sample', catch.sel='sample'),
  years=c(wt=10, mat=10, catch.sel=10))


# - PROJECT forward for iy assumption

# F=Fsq
om <- fwd(om, fbar=FLQuant(c(fbar(run)[,'2023']), dimnames=list(year=2024)))


# - CONSTRUCT FLom

# GENERATE future deviances


# --- model.R

# - SET UP MP runs

# (i)mplementation (sys)tem: tac.is (C ~ F) + F deviances

# - EXAMPLE oem

oem <- FLoem(
  observations=list(stk=stock(om)),
  deviances=list(stk=FLQuants(stock.n=rlnorm(500, catch.n(om) %=% 0, 0.2))),
  method=shortcut.oem)

# RUN applying the OEM
error <- mp(om, oem=oem, ctrl=arule, args=mseargs)
