# Teaching examples printed in the manuscript. Run from the package root.
library(dplyr)
library(minpack.lm)
dat <- read.csv("Data/Study3/Input_Data54.csv")

# Correct the component proportions
p_V <- (dat$R_V * dat$n_V + 0.5) / (dat$n_V + 1)
p_A <- (dat$R_A * dat$n_A + 0.5) / (dat$n_A + 1)
dat$L_V <- qlogis(p_V)
dat$L_A <- qlogis(p_A)
SCALE_MID <- 0

# Load model predictors, starts and fitting helpers.
source("Scripts/Study3/help_function_choice.R")
SCALE_MID <- 0
start_7p <- start_params[["7p"]]
# Evaluate the starting-value prediction.
eta <- choice_eta_functions[["7p"]](start_7p, dat)
predicted_probability <- plogis(eta)

# Stable binomial negative log-likelihood.
choice_nll <- function(par, data, model) {
  eta <- choice_eta_functions[[model]](par, data)
  k <- data$VA_disgust
  n <- data$n_VA
  sum((n - k)*pmax(eta, 0) + k*pmax(-eta, 0) +
        n*log1p(exp(-abs(eta))) - lchoose(n, k))
}

optim_7p <- optim(
  par = start_7p, fn = choice_nll,
  data = dat, model = "7p", method = "L-BFGS-B",
  control = list(maxit = 10000, factr = 1e7, pgtol = 1e-8),
  hessian = TRUE
)
# Store fitted predictions and Hessian-based uncertainty.
fit_7p <- choice_result(optim_7p, "7p", dat)

k_theta <- length(optim_7p$par)
K <- k_theta
N <- nrow(dat)
loglik <- -optim_7p$value
AIC <- -2*loglik + 2*K
c(AICc = AIC + 2*K*(K + 1)/(N - K - 1),
  BIC = -2*loglik + K*log(N))
dat$eta <- predict(fit_7p, type = "link")
dat$predicted_probability <- predict(fit_7p, type = "probability")

choice_wald_parameter_table(fit_7p)

dat$VV <- factor(dat$VV, levels = c(-1, 0, 1))
dat$AV <- factor(dat$AV, levels = c(-1, 0, 1))
dat$VA_anger <- dat$n_VA - dat$VA_disgust

gA <- glm(
  cbind(VA_disgust, VA_anger) ~ AV + offset(eta),
  family = binomial(), data = dat
)
gAdd <- update(gA, . ~ . + VV)
anova(gA, gAdd, test = "LRT")
