
# FLR toolset for evaluating multi-annual management plans for [ICES WKMFOA](https://www.ices.dk/advice/Advice-activities/Lists/Posts/Post.aspx?ID=32)

- Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
- Jasper BLEIJENBERG (WMR) <jasper.bleijenberg@wur.nl>

This TAF repository contains a template for conducting an evaluation of multi-annual management plans. A single management plan, assuming a shortcut of the stock assessment, the standard ICES advice rule and short-term forecast, can be run with different *frequencies*, that is, how often advice is provide.

## Installation

The latest version of the required [FLR](https:://flr-project.org) packages, and all their dependencies, can be installed from the [FLR R-universe page](https://flr.r-universe.dev) by calling:

```r
install.packages("icesTAF")

install.packages(icesTAF::deps(), repos=c(
  FLR="https://flr.r-universe.dev",
  CRAN="https://cloud.r-project.org/"))
```

## data.R

There are three version of the data.R script, setup different according to the stock assessment model used to condition the OM.

- data_ss3.R, uses an Stock Synthesis stock assessment (sol.27.4 as example). The OM contains the stock-recruitment relationship estimated by SS3.

- data_sam.R, takes the result of a SAM run (). It builds the SRR uncertainty by bootstrap, and add the computed process error in the hindcast projection.

- data_a4a.R bases the OM on the result of an FLa4a run (). SRR uncertainty is also computed using bootstrap.
