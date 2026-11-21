## Study 3: grouped-binomial fitting of anger-disgust choices.
## The working model assumes independent choices with one probability per pair.
## Fit six candidates with one optim start each and no parameter bounds.

## 0. Environment and dependencies -----------------------------------------
rm(list = ls())
# Load shared paths relative to this script.
.script_files <- unlist(lapply(sys.frames(), function(frame) frame$ofile))
.script_file <- if (length(.script_files)) tail(.script_files, 1) else
  sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
if (!length(.script_file)) .script_file <- if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) rstudioapi::getSourceEditorContext()$path else "Step2_ModelFitting.R"
source(file.path(dirname(normalizePath(.script_file)), "..", "Project_Config.R"), encoding = "UTF-8")
setwd(file.path(PROJECT_ROOT, "Scripts", "Study3"))
packages <- c("dplyr", "ggplot2", "tidyr", "tibble")
to_install <- setdiff(packages, rownames(installed.packages()))
if (length(to_install)) stop("Install the missing R packages first: ", paste(to_install, collapse = ", "))

library(dplyr)
dir.create("results", showWarnings = FALSE)
library(ggplot2)
library(tidyr)
library(tibble)

# Load integration predictors, binomial likelihood and fitted-result methods.
source("help_function_choice.R")

# On the log-odds scale, probability .5 corresponds to 0.
SCALE_MID <- 0


# For each candidate, set starting values, fit the model and inspect its summary.
# Use component log odds as inputs and binomial negative log-likelihood for estimation.
## 1. Load empirical data --------------------------------------------------
# Step 1 provides 50 audiovisual responses per pair, matched by stimulus identity.
# Correct component proportions using their observed denominators.
dat <- read.csv(file.path(PROJECT_ROOT, "Data", "Study3", "Input_Data54.csv"))

dat <- choice_prepare_inputs(dat)

dat <- dat %>%
  mutate(
    VV = factor(VV,
      levels = c(-1, 0, 1),
      labels = c("anger", "mixed", "disgust")
    ),
    AV = factor(AV,
      levels = c(-1, 0, 1),
      labels = c("anger", "mixed", "disgust")
    )
  )

# Store corrected component log odds in L_V and L_A.
str(dat)

## 2. Fit each candidate independently with optim --------------------------
# Use L-BFGS-B without parameter bounds; the Hessian supplies approximate Wald uncertainty.

start_3p <- start_params[["3p"]]
optim_3p <- optim(
  par = start_3p,
  fn = choice_nll,
  data = dat,
  model = "3p",
  method = "L-BFGS-B",
  control = list(maxit = 10000, factr = 1e7, pgtol = 1e-8),
  hessian = TRUE
)
diagnostics_3p <- choice_optimizer_diagnostics(optim_3p, "3p", dat)
print(diagnostics_3p)
fit_3p <- choice_result(optim_3p, "3p", dat)
summary(fit_3p)

start_4p <- start_params[["4p"]]
optim_4p <- optim(
  par = start_4p,
  fn = choice_nll,
  data = dat,
  model = "4p",
  method = "L-BFGS-B",
  control = list(maxit = 10000, factr = 1e7, pgtol = 1e-8),
  hessian = TRUE
)
diagnostics_4p <- choice_optimizer_diagnostics(optim_4p, "4p", dat)
print(diagnostics_4p)
fit_4p <- choice_result(optim_4p, "4p", dat)
summary(fit_4p)

start_4plq <- start_params[["4plq"]]
optim_4plq <- optim(
  par = start_4plq,
  fn = choice_nll,
  data = dat,
  model = "4plq",
  method = "L-BFGS-B",
  control = list(maxit = 10000, factr = 1e7, pgtol = 1e-8),
  hessian = TRUE
)
diagnostics_4plq <- choice_optimizer_diagnostics(optim_4plq, "4plq", dat)
print(diagnostics_4plq)
fit_4plq <- choice_result(optim_4plq, "4plq", dat)
summary(fit_4plq)

start_5p <- start_params[["5p"]]
optim_5p <- optim(
  par = start_5p,
  fn = choice_nll,
  data = dat,
  model = "5p",
  method = "L-BFGS-B",
  control = list(maxit = 10000, factr = 1e7, pgtol = 1e-8),
  hessian = TRUE
)
diagnostics_5p <- choice_optimizer_diagnostics(optim_5p, "5p", dat)
print(diagnostics_5p)
fit_5p <- choice_result(optim_5p, "5p", dat)
summary(fit_5p)

start_5plq <- start_params[["5plq"]]
optim_5plq <- optim(
  par = start_5plq,
  fn = choice_nll,
  data = dat,
  model = "5plq",
  method = "L-BFGS-B",
  control = list(maxit = 10000, factr = 1e7, pgtol = 1e-8),
  hessian = TRUE
)
diagnostics_5plq <- choice_optimizer_diagnostics(optim_5plq, "5plq", dat)
print(diagnostics_5plq)
fit_5plq <- choice_result(optim_5plq, "5plq", dat)
summary(fit_5plq)

start_7p <- start_params[["7p"]]
optim_7p <- optim(
  par = start_7p,
  fn = choice_nll,
  data = dat,
  model = "7p",
  method = "L-BFGS-B",
  control = list(maxit = 10000, factr = 1e7, pgtol = 1e-8),
  hessian = TRUE
)
diagnostics_7p <- choice_optimizer_diagnostics(optim_7p, "7p", dat)
print(diagnostics_7p)
fit_7p <- choice_result(optim_7p, "7p", dat)
summary(fit_7p)

## 3. Model comparison and selection ---------------------------------------
fit_lookup <- list(
  `3p` = fit_3p,
  `4p` = fit_4p,
  `4plq` = fit_4plq,
  `5p` = fit_5p,
  `5plq` = fit_5plq,
  `7p` = fit_7p
)

# Single-modality benchmarks predict directly from corrected component log odds.
# An empty parameter vector gives a fixed benchmark with zero fitted coefficients.
single_cue_statistics <- function(model, eta) {
  p_hat <- plogis(eta)
  baseline <- list(
    model = model, par = numeric(0),
    nll = -sum(dbinom(dat$VA_disgust, dat$n_VA, p_hat, log = TRUE)),
    predicted_probability = p_hat,
    diagnostics = data.frame(converged = logical(0))
  )
  model_statistics(baseline, dat) %>%
    mutate(converged_starts = NA_integer_, total_starts = NA_integer_)
}

cmp_tbl <- bind_rows(
  model_statistics(fit_3p, dat),
  model_statistics(fit_4p, dat),
  model_statistics(fit_4plq, dat),
  model_statistics(fit_5p, dat),
  model_statistics(fit_5plq, dat),
  model_statistics(fit_7p, dat),
  single_cue_statistics("Visual-only", dat$L_V),
  single_cue_statistics("Auditory-only", dat$L_A)
) %>%
  mutate(
    delta_AIC = AIC - min(AIC),
    delta_AICc = AICc - min(AICc),
    delta_BIC = BIC - min(BIC),
    AIC_weight = exp(-0.5 * delta_AIC) /
      sum(exp(-0.5 * delta_AIC))
  ) %>%
  arrange(AICc)

print(cmp_tbl)

## 4. Save fitted models and Wald parameter inference ----------------------
# Wald summaries test each coefficient against zero.
parameter_tbl <- bind_rows(
  choice_wald_parameter_table(fit_3p),
  choice_wald_parameter_table(fit_4p),
  choice_wald_parameter_table(fit_4plq),
  choice_wald_parameter_table(fit_5p),
  choice_wald_parameter_table(fit_5plq),
  choice_wald_parameter_table(fit_7p)
)

print(parameter_tbl)

OUTPUT_TAG <- "BinomialLikelihood_54"
# Retain a specification when AICc and BIC agree; stop if they differ.
best_AICc <- cmp_tbl$model[which.min(cmp_tbl$AICc)]
best_BIC <- cmp_tbl$model[which.min(cmp_tbl$BIC)]
if (best_AICc != best_BIC) stop("AICc and BIC prefer different specifications; review the comparison before continuing.")
best_model_code <- best_AICc
best_fit <- fit_lookup[[best_model_code]]

best_predictions <- dat %>%
  mutate(
    # These columns describe fitted weights; they are not extra fitted coefficients.
    Stimcomb = row_number(),
    VA_anger = n_VA - VA_disgust,
    # This auxiliary joint logit is descriptive only; the likelihood uses original counts.
    L_VA = choice_count_to_logodds(VA_disgust, n_VA),
    predicted_eta = best_fit$eta,
    predicted_disgust = best_fit$predicted_probability,
    probability_residual = R_VA - predicted_disgust,
    latent_logit_residual = L_VA - predicted_eta
  )

fit_objects <- list(
  fit_3p = fit_3p,
  fit_4p = fit_4p,
  fit_4plq = fit_4plq,
  fit_5p = fit_5p,
  fit_5plq = fit_5plq,
  fit_7p = fit_7p,
  cmp_tbl = cmp_tbl,
  input_scale = "empirical_logodds",
  response_distribution = "grouped_binomial",
  input_denominators = dat[c("k_V", "n_V", "k_A", "n_A")],
  scale_mid = SCALE_MID,
  best_model = best_model_code,
  best_AICc = best_AICc, best_BIC = best_BIC
)

saveRDS(
  fit_objects,
  file = paste0("results/model_fits_weighted_averaging_", OUTPUT_TAG, ".rds")
)
write.csv(
  cmp_tbl,
  paste0("results/Study3_ItemLevel_", OUTPUT_TAG, "_ModelComparison.csv"),
  row.names = FALSE
)
write.csv(
  parameter_tbl,
  paste0("results/Study3_ItemLevel_", OUTPUT_TAG, "_ParameterEstimates.csv"),
  row.names = FALSE
)
choice_save_diagnostics(fit_lookup, OUTPUT_TAG)
write.csv(
  best_predictions,
  paste0("results/Study3_ItemLevel_", OUTPUT_TAG, "_BestModelPredictions.csv"),
  row.names = FALSE
)

cat("\nBinomial-likelihood fitting and model comparison completed.\n")
cat("AICc/BIC-best model: ", best_model_code, "\n", sep = "")
## 5. Plot model comparison -----------------------------------------------
# Subtract each criterion's own minimum for display.
dir.create("figures", showWarnings = FALSE)
comparison <- cmp_tbl %>% mutate(model = sub(" .*", "", model)) %>%
  select(model, AIC, AICc, BIC) %>%
  pivot_longer(-model, names_to = "Criterion", values_to = "IC") %>%
  group_by(Criterion) %>% mutate(Delta = IC - min(IC)) %>% ungroup() %>%
  mutate(model = factor(model, levels = rev(c("Visual-only", "Auditory-only", "3p", "4p", "4plq", "5p", "5plq", "7p"))))
comparison_display <- comparison %>% filter(Criterion %in% c("AICc", "BIC"))
p_cmp <- ggplot(comparison_display, aes(Delta, model, color = Criterion)) +
  geom_point(position = position_dodge(width = .45), size = 4, alpha = .9) +
  scale_color_manual(name = NULL, values = c(AICc = FIGURE_ORANGE, BIC = FIGURE_GREEN)) +
  scale_y_discrete(limits = rev(sub(" .*", "", cmp_tbl$model))) +
  labs(x = "\u0394 information criterion", y = "Model") + theme_classic(base_family = FIGURE_FONT) + theme(legend.position = "top")
# Display every model while omitting an empty interval of the x axis.
# The two displayed segments have separate plotting scales.
stopifnot(!any(comparison$Delta > 25 & comparison$Delta < 950))
p_cmp <- p_cmp +
  scale_x_continuous(limits = c(-1, 1210),
                     breaks = c(0, 10, 20, 1000, 1200), expand = c(0, 0)) +
  ggbreak::scale_x_break(c(25, 950), scales = 0.55,
                        space = 0.18, symbol = "slash", expand = FALSE)
save_figure("figures/03_ModelComparison.png", p_cmp, width = FIGURE_WIDTH_IN, height = FIGURE_HEIGHT_IN, dpi = FIGURE_DPI)
write.csv(comparison, "results/Plot_ModelComparison.csv", row.names = FALSE)
