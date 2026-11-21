# Teaching examples printed in the manuscript. Run from the package root.
library(dplyr)
library(minpack.lm)
source("Scripts/Study2/simulation_helpers.R")
dat <- read.csv("Data/Study2/Input_Data54.csv")
SCALE_MIN <- 0
SCALE_MAX <- 1
SCALE_MID <- 0.5

dat$pred_visual <- dat$R_V
dat$pred_auditory <- dat$R_A

nls_3p <- nlsLM(
  R_VA ~ wavg_3(S0, W_V, W_A, R_V, R_A),
  data = dat, start = start_params[["3p"]],
  control = nls.lm.control(maxiter = 500)
)
