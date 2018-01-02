arg <- commandArgs(trailingOnly=T)
load(arg)
require(unmarked)
summary(all_covs_m)

