# Study 2: fit six continuous integration models and two fixed benchmarks.
# Input_Data54.csv supplies 54 pair-level mean ratings; compare AICc/BIC.

## 0. Environment and dependencies -----------------------------------------
rm(list = ls())
# Load shared paths relative to this script.
.script_files <- unlist(lapply(sys.frames(), function(frame) frame$ofile))
.script_file <- if (length(.script_files)) tail(.script_files, 1) else
  sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
if (!length(.script_file)) .script_file <- if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) rstudioapi::getSourceEditorContext()$path else "Step2_ModelFitting.R"
source(file.path(dirname(normalizePath(.script_file)), "..", "Project_Config.R"), encoding = "UTF-8")
setwd(file.path(PROJECT_ROOT, "Scripts", "Study2"))

source("simulation_helpers.R")

# Check required packages.
packages <- c("dplyr", "ggplot2",
              "minpack.lm", "broom", "tidyr", "tibble")
to_install <- setdiff(packages, rownames(installed.packages()))
if (length(to_install)) stop("请先安装R包：", paste(to_install, collapse = ", "))

library(dplyr)
dir.create("results", showWarnings = FALSE)
library(ggplot2)
library(minpack.lm)
library(broom)
library(tidyr)
library(tibble)

## 1. Load empirical data --------------------------------------------------
# Load the 54 pair-level mean ratings prepared by Step 1.
dat <- read.csv(file.path(PROJECT_ROOT, "Data", "Study2", "Input_Data54.csv"))

dat <- dat %>%
  mutate(
    VV = factor(VV, levels = c(-1, 0, 1), labels = c("anger","mixed","disgust")),
    AV = factor(AV, levels = c(-1, 0, 1), labels = c("anger","mixed","disgust"))
  )

str(dat)

## 2. Model fitting (nonlinear least squares) ------------------------------
## Models: 3p, 4p, 4plq, 5p, 5plq, 7p
## This is Wilson & Collins' "Fit the parameters" step.

## 2.1 3p constant-weight averaging: wavg_3 ---------------------
start_3p <- start_params[["3p"]]

nls_3p <- nlsLM(
  R_VA ~ wavg_3(S0, W_V, W_A, R_V, R_A),
  data    = dat,
  start   = start_3p,
  control = nls.lm.control(maxiter = 500)
)
summary(nls_3p)
rss_3p <- sum(residuals(nls_3p)^2)

## 2.2 4p linear weights (shared slope): wavg_4 ----------------------------
start_4p <- start_params[["4p"]]

nls_4p <- nlsLM(
  R_VA ~ wavg_4(S0, A_V, A_A, B_VA, R_V, R_A),
  data    = dat,
  start   = start_4p,
  control = nls.lm.control(maxiter = 1000)
)
summary(nls_4p)
rss_4p <- sum(residuals(nls_4p)^2)

## 2.3 4plq lin+quad with A_V = A_A = A_VA: wavg_4plq ----------------------
start_4plq <- start_params[["4plq"]]

nls_4plq <- nlsLM(
  R_VA ~ wavg_4plq(S0, A_VA, B_VA, B_VA2, R_V, R_A),
  data    = dat,
  start   = start_4plq,
  control = nls.lm.control(maxiter = 1000)
)
summary(nls_4plq)
rss_4plq <- sum(residuals(nls_4plq)^2)

## 2.4 5p linear weights (separate slopes): wavg_5 -------------------------
start_5p <- start_params[["5p"]]

nls_5p <- nlsLM(
  R_VA ~ wavg_5(S0, A_V, B_V, A_A, B_A, R_V, R_A),
  data    = dat,
  start   = start_5p,
  control = nls.lm.control(maxiter = 1000)
)
summary(nls_5p)
rss_5p <- sum(residuals(nls_5p)^2)

## 2.5 5plq lin+quad with separate intercepts, shared slope & curvature: wavg_5plq ---
start_5plq <- start_params[["5plq"]]

nls_5plq <- nlsLM(
  R_VA ~ wavg_5plq(S0, A_V, A_A, B_VA, B_VA2, R_V, R_A),
  data    = dat,
  start   = start_5plq,
  control = nls.lm.control(maxiter = 1000)
)
summary(nls_5plq)
rss_5plq <- sum(residuals(nls_5plq)^2)

## 2.6 7p lin+quad with separate slopes/curvatures: wavg_7 -----------------
start_7p <- start_params[["7p"]]

nls_7p <- nlsLM(
  R_VA ~ wavg_7(S0, A_V, B_V, B_V2, A_A, B_A, B_A2, R_V, R_A),
  data    = dat,
  start   = start_7p,
  control = nls.lm.control(maxiter = 1000)
)
summary(nls_7p)
rss_7p <- sum(residuals(nls_7p)^2)

## 3. Model comparison and selection (AIC/AICc/BIC) ------------------------

n <- nrow(dat)

p_3p   <- length(coef(nls_3p))
p_4p   <- length(coef(nls_4p))
p_4plq <- length(coef(nls_4plq))
p_5p   <- length(coef(nls_5p))
p_5plq <- length(coef(nls_5plq))
p_7p   <- length(coef(nls_7p))

sst <- sum((dat$R_VA - mean(dat$R_VA))^2)

aic_bic_from_rss <- function(n, rss, k_par) {
  k_eff <- k_par + 1
  aic   <- n * log(rss / n) + 2 * k_eff
  bic   <- n * log(rss / n) + k_eff * log(n)
  r2    <- 1 - rss / sst
  aicc <- if (n > k_eff + 1) aic + 2 * k_eff * (k_eff + 1) / (n - k_eff - 1) else NA_real_
  c(AIC = aic, AICc = aicc, BIC = bic, R2 = r2)
}

# Single-modality benchmarks use observed component ratings without fitted coefficients.
rss_visual <- sum((dat$R_VA - dat$R_V)^2)
rss_auditory <- sum((dat$R_VA - dat$R_A)^2)

cmp_tbl <- tibble::tibble(
  model = c("3p (constant weights)",
            "4p (linear, shared slope)",
            "4plq (lin+quad, A_V=A_A)",
            "5p (linear, separate slopes)",
            "5plq (lin+quad, shared slope/quad)",
            "7p (lin+quad, separate slopes/curvatures)",
            "Visual-only", "Auditory-only"),
  k_par = c(p_3p, p_4p, p_4plq, p_5p, p_5plq, p_7p, 0, 0),
  rss   = c(rss_3p, rss_4p, rss_4plq, rss_5p, rss_5plq, rss_7p,
            rss_visual, rss_auditory)
) %>%
  rowwise() %>%
  mutate(stats = list(aic_bic_from_rss(n, rss, k_par))) %>%
  tidyr::unnest_wider(stats) %>%
  ungroup() %>%
  mutate(n_AICc = n, k_AICc = k_par + 1,
         delta_AICc = AICc - min(AICc),
         delta_BIC = BIC - min(BIC)) %>%
  arrange(AICc)

print(cmp_tbl)
# Retain a specification when AICc and BIC agree; stop if they differ.
best_AICc <- sub(" .*", "", cmp_tbl$model[which.min(cmp_tbl$AICc)])
best_BIC <- sub(" .*", "", cmp_tbl$model[which.min(cmp_tbl$BIC)])
cat("AICc-best:", best_AICc, "; BIC-best:", best_BIC, "\n")
if (best_AICc != best_BIC) stop("AICc与BIC选择不同，请先确定后续分析模型。")

## 4. Save fitted models for later validation ------------------------------

parameter_tbl <- dplyr::bind_rows(
  broom::tidy(nls_3p)   %>% mutate(model = "3p"),
  broom::tidy(nls_4p)   %>% mutate(model = "4p"),
  broom::tidy(nls_4plq) %>% mutate(model = "4plq"),
  broom::tidy(nls_5p)   %>% mutate(model = "5p"),
  broom::tidy(nls_5plq) %>% mutate(model = "5plq"),
  broom::tidy(nls_7p)   %>% mutate(model = "7p")
) %>%
  mutate(
    CI_low_95 = estimate - qt(0.975, df = 54 -
      dplyr::case_when(model == "3p" ~ 3, model %in% c("4p", "4plq") ~ 4,
                       model %in% c("5p", "5plq") ~ 5, model == "7p" ~ 7)) * std.error,
    CI_high_95 = estimate + qt(0.975, df = 54 -
      dplyr::case_when(model == "3p" ~ 3, model %in% c("4p", "4plq") ~ 4,
                       model %in% c("5p", "5plq") ~ 5, model == "7p" ~ 7)) * std.error
  ) %>%
  select(model, everything())

fit_objects <- list(
  nls_3p   = nls_3p,
  nls_4p   = nls_4p,
  nls_4plq = nls_4plq,
  nls_5p   = nls_5p,
  nls_5plq = nls_5plq,
  nls_7p   = nls_7p,
  best_model = best_AICc,
  best_AICc = best_AICc, best_BIC = best_BIC,
  cmp_tbl  = cmp_tbl
)

saveRDS(fit_objects, file = "results/model_fits_weighted_averaging.rds")
write.csv(cmp_tbl, "results/Study2_ItemLevel_ModelComparison.csv", row.names = FALSE)
write.csv(parameter_tbl, "results/Study2_ItemLevel_ParameterEstimates.csv", row.names = FALSE)

cat("\nModel fitting, AIC/AICc/BIC comparison, and saving fitted models completed.\n")

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
stopifnot(!any(comparison$Delta > 10 & comparison$Delta < 120))
p_cmp <- p_cmp +
  scale_x_continuous(limits = c(-0.5, 152),
                     breaks = c(0, 5, 10, 135, 150), expand = c(0, 0)) +
  ggbreak::scale_x_break(c(10, 120), scales = 0.55,
                        space = 0.18, symbol = "slash", expand = FALSE)
save_figure("figures/03_ModelComparison.png", p_cmp, width = FIGURE_WIDTH_IN, height = FIGURE_HEIGHT_IN, dpi = FIGURE_DPI)
write.csv(comparison, "results/Plot_ModelComparison.csv", row.names = FALSE)
