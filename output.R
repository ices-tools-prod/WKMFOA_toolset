# output.R - Performance evaluation and output tables
# WKREBUILD_toolset/output.R

# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


library(icesTAF)
mkdir("output")

library(mse)

source("utilities.R")

# LOAD model.R outputs

load("model/model.rda")


# --- TABLES

tables <- list()
perf <- performance(runs)

# WHAT is P(SB < Blim once) by mp?
perf[statistic == 'risk2',]

# IS P(SB > Blim) by year?
dcast(perf[statistic == 'PBlim', .(data=mean(data)),
  by=.(year, mp, name, desc)], mp ~ year, value.var='data')

# WHAT are the yearly average catch levels?
tables$catch_mp <- dcast(perf[statistic == 'Cy', .(data=mean(data)),
  by=.(year, mp, name, desc)], mp ~ year, value.var='data')

# WHAT are the long-term average annual catch levels?
perf[statistic == 'C', .(data=mean(data)), by=mp]

# WHAT is the average catch variability?
dcast(perf[statistic == 'AAVC', .(data=mean(data)),
  by=.(year, mp, name, desc)], mp ~ year, value.var='data')

# WHAT is the inter-annual catch variability?
dcast(perf[statistic == 'IACC', .(data=mean(data)),
  by=.(year, mp, name, desc)], mp ~ year, value.var='data')

# CREATE table of all statistics by mp and year (statistic + mp ~ year)
tables$stats_mp <- dcast(perf[, .(data=mean(data)),
  by=.(year, mp, name, desc)], name + mp ~ year, value.var='data')

# SAVE
save(tables, file="output/output.rda", compress="xz")
