# Gao (2018) Information integration theory.
# Model fitting & model comparison
# Candidate models: 3p, 4p, 4plq, 5p, 5plq, 7p
# Data file: MODELING_Exp1_Valence.xlsx

## 0. Environment and dependencies -----------------------------------------
rm(list = ls())
# Load shared paths relative to this script.
.script_files <- unlist(lapply(sys.frames(), function(frame) frame$ofile))
.script_file <- if (length(.script_files)) tail(.script_files, 1) else
  sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
if (!length(.script_file)) .script_file <- if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) rstudioapi::getSourceEditorContext()$path else "Step3_ModelFitting.R"
source(file.path(dirname(normalizePath(.script_file)), "..", "Project_Config.R"), encoding = "UTF-8")
setwd(file.path(PROJECT_ROOT, "Scripts", "Study1"))
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# One-time installation of required packages
packages <- c("readxl", "dplyr", "ggplot2", "car",
              "minpack.lm", "broom", "tidyr", "tibble")
to_install <- setdiff(packages, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)

library(readxl)
library(dplyr)
library(ggplot2)
library(car)
library(minpack.lm)
library(broom)
library(tidyr)
library(tibble)

# Load helper functions (wavg_3/4/4plq/5/5plq/7p, reference_params, etc.)
source("simulation_helpers.R")

## 1. Load empirical data --------------------------------------------------
dat <- read_excel(file.path(PROJECT_ROOT, "Data", "Study1", "MODELING_Exp1_Valence.xlsx"))

dat <- dat %>% mutate(Stimcomb = dplyr::row_number())

dat <- dat %>%
  mutate(
    VV = factor(VV, levels = c(-1, 0, 1), labels = c("Neg","Neu","Pos")),
    AV = factor(AV, levels = c(-1, 0, 1), labels = c("Neg","Neu","Pos"))
  )

str(dat)

## 2. Model fitting (maximum likelihood / least squares) -------------------
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
  control = nls.lm.control(maxiter = 2000)
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
  control = nls.lm.control(maxiter = 2000)
)
summary(nls_5plq)
coef_5plq <- summary(nls_5plq)$coefficients
ci_margin_5plq <- qt(0.975, df.residual(nls_5plq)) * coef_5plq[, "Std. Error"]
parameter_ci_5plq <- data.frame(
  parameter = rownames(coef_5plq),
  estimate = coef_5plq[, "Estimate"],
  SE = coef_5plq[, "Std. Error"],
  CI_low_95 = coef_5plq[, "Estimate"] - ci_margin_5plq,
  CI_high_95 = coef_5plq[, "Estimate"] + ci_margin_5plq,
  row.names = NULL
)
write.csv(parameter_ci_5plq, "results/Study1_5plq_ParameterCI.csv", row.names = FALSE)
rss_5plq <- sum(residuals(nls_5plq)^2)

## 2.6 7p lin+quad with separate slopes/curvatures: wavg_7 -----------------
start_7p <- start_params[["7p"]]

nls_7p <- nlsLM(
  R_VA ~ wavg_7(S0, A_V, B_V, B_V2, A_A, B_A, B_A2, R_V, R_A),
  data    = dat,
  start   = start_7p,
  control = nls.lm.control(maxiter = 2000)
)
summary(nls_7p)
rss_7p <- sum(residuals(nls_7p)^2)

## 3. Model comparison and selection (AIC/BIC only) ------------------------

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
  c(AIC = aic, BIC = bic, R2 = r2)
}

cmp_tbl <- tibble::tibble(
  model = c("3p (constant weights)",
            "4p (linear, shared slope)",
            "4plq (lin+quad, A_V=A_A)",
            "5p (linear, separate slopes)",
            "5plq (lin+quad, shared slope/quad)",
            "7p (lin+quad, separate slopes/curvatures)"),
  k_par = c(p_3p, p_4p, p_4plq, p_5p, p_5plq, p_7p),
  rss   = c(rss_3p, rss_4p, rss_4plq, rss_5p, rss_5plq, rss_7p)
) %>%
  rowwise() %>%
  mutate(stats = list(aic_bic_from_rss(n, rss, k_par))) %>%
  tidyr::unnest_wider(stats) %>%
  ungroup() %>%
  arrange(AIC)

print(cmp_tbl)

## 4. Save fitted models for later validation ------------------------------

fit_objects <- list(
  nls_3p   = nls_3p,
  nls_4p   = nls_4p,
  nls_4plq = nls_4plq,
  nls_5p   = nls_5p,
  nls_5plq = nls_5plq,
  nls_7p   = nls_7p,
  cmp_tbl  = cmp_tbl
)

saveRDS(fit_objects, file = "results/model_fits_weighted_averaging.rds")

cat("\nModel fitting, AIC/BIC comparison, and saving fitted models completed.\n")

## 5. Plotting -------------------------------------------------------------

cmp_tbl2 <- cmp_tbl %>%
  mutate(
    model_code = sub("\\s*\\(.*$", "", model),
    model_num  = case_when(
      model_code == "3p"   ~ 1L,
      model_code == "4p"   ~ 2L,
      model_code == "4plq" ~ 3L,
      model_code == "5p"   ~ 4L,
      model_code == "5plq" ~ 5L,
      model_code == "7p"   ~ 6L
    ),
    model_label = paste0("Model ", model_num, " (", model_code, ")")
  )

plot_df <- cmp_tbl2 %>%
  select(model_label, AIC, BIC) %>%
  pivot_longer(cols = c(AIC, BIC),
               names_to = "metric",
               values_to = "value")

order_by_aic <- cmp_tbl2 %>%
  arrange(AIC) %>%
  pull(model_label)

plot_df <- plot_df %>%
  mutate(
    model_label = factor(model_label,
                         levels = order_by_aic),   # Order models by AIC.
    metric = factor(metric, levels = c("AIC", "BIC"))
  )

p_ic <- ggplot(plot_df,
               aes(x = value, y = model_label, color = metric)) +
  geom_point(size = 4, alpha = 0.9) +
  scale_color_manual(
    values = c("AIC" = FIGURE_ORANGE,
               "BIC" = FIGURE_GREEN)
  ) +
  scale_y_discrete(limits = rev(order_by_aic)) +   # Lowest AIC at the top.
  labs(
    x = "Information criterion",
    y = "Model",
    color = ""
  ) +
  theme_bw(base_size = FIGURE_TEXT_PT, base_family = FIGURE_FONT) +
  theme(
    plot.title = element_blank(),
    axis.title.x = element_text(size = FIGURE_AXIS_PT, face = "plain", margin = margin(t = 10)),
    axis.title.y = element_text(size = FIGURE_AXIS_PT, face = "plain", margin = margin(r = 10)),
    axis.text.x  = element_text(size = FIGURE_TEXT_PT),
    axis.text.y  = element_text(size = FIGURE_TEXT_PT),
    legend.text  = element_text(size = FIGURE_TEXT_PT),
    legend.position = "top",
    panel.grid = element_blank()
  )

print(p_ic)

save_figure("figures/Fig2_ModelFitting.png", width = FIGURE_WIDTH_IN, height = FIGURE_HEIGHT_IN, dpi = FIGURE_DPI)
