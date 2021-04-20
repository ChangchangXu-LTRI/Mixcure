# Mixcure

R language statistical function package for jointly modelling time-to-event and probability of being event-free. The package consists of three functions for applying mixture cure model on a time-to-event dataset with long term survivors. THe Firth-type penalized likelihood approach was introduced as an alternative to the MLE approach, to handle heavily censored data with imbalanced categorical covariates, where biased or non-converged estimates are likely to be observed as a result of small/sparse sample, known as complete separation or monotone likelihood.

Created by Changchang Xu

Contact:changchang.xu@mail.utoronto.ca

This package can be installed via the following R code:

devtools::install_github("ChangchangXu-LTRI/mixcure", build = TRUE, build_opts = c())

library(mixcure)
