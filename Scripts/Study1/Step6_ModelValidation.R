# Gao (2018) Information integration theory
# Model validation for the retained model (5plq only)
# - Cell means vs 5plq predictions

# Requires:
# - MODELING_Exp1_Valence.xlsx
# - model_fits_weighted_averaging.rds  (from model_fitting script)
# - simulation_helpers.R               (for wavg_5plq etc., used by predict)

## 0. Environment and dependencies -----------------------------------------
rm(list = ls())
# Load shared paths relative to this script.
.script_files <- unlist(lapply(sys.frames(), function(frame) frame$ofile))
.script_file <- if (length(.script_files)) tail(.script_files, 1) else
  sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
if (!length(.script_file)) .script_file <- if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) rstudioapi::getSourceEditorContext()$path else "Step6_ModelValidation.R"
source(file.path(dirname(normalizePath(.script_file)), "..", "Project_Config.R"), encoding = "UTF-8")
setwd(file.path(PROJECT_ROOT, "Scripts", "Study1"))
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

packages <- c("readxl", "dplyr", "ggplot2", "car", "tidyr", "tibble")
to_install <- setdiff(packages, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)

library(readxl)
library(dplyr)
library(ggplot2)
library(car)
library(tidyr)
library(tibble)

# Load helper functions (for wavg_5plq etc., which are needed by predict())
source("simulation_helpers.R")

## 1. Load empirical data --------------------------------------------------

dat <- read_excel(file.path(PROJECT_ROOT, "Data", "Study1", "MODELING_Exp1_Valence.xlsx"))

dat <- dat %>%
  mutate(
    Stimcomb = dplyr::row_number(),
    VV = factor(VV, levels = c(-1, 0, 1), labels = c("Neg","Neu","Pos")),
    AV = factor(AV, levels = c(-1, 0, 1), labels = c("Neg","Neu","Pos"))
  )

str(dat)

## 2. Load fitted 5plq model from disk -------------------------------------

fit_objects <- readRDS("results/model_fits_weighted_averaging.rds")

# Extract the retained model
nls_5plq <- fit_objects$nls_5plq

cat("\nLoaded 5plq fitted model. Coefficients:\n")
print(coef(nls_5plq))

## 3. Add fitted values and residuals for 5plq -----------------------------

dat <- dat %>%
  mutate(
    est_5plq = fitted(nls_5plq),
    res_5plq = residuals(nls_5plq)
  )

write.csv(dat, "results/Validation_PairPredictions.csv", row.names = FALSE)

## 4. Residual ANOVA: VV × AV effects in residuals -------------------------
## If 5plq captures the VV × AV structure, residuals should not show
## strong main effects or interactions.

# Sum coding applies to these Type III diagnostics, not IIT estimation.
contrasts(dat$VV) <- contr.sum(nlevels(dat$VV))
contrasts(dat$AV) <- contr.sum(nlevels(dat$AV))
aov_res_5plq <- aov(res_5plq ~ VV * AV, data = dat)
cat("\nType-III ANOVA on 5plq residuals (VV × AV):\n")
residual_anova <- car::Anova(aov_res_5plq, type = "III")
print(residual_anova)
term_names <- setdiff(rownames(residual_anova), c("(Intercept)", "Residuals"))
df_residual <- residual_anova["Residuals", "Df"]
residual_effects <- data.frame(
  effect = term_names,
  df1 = residual_anova[term_names, "Df"],
  df2 = df_residual,
  F = residual_anova[term_names, "F value"],
  partial_eta2 = residual_anova[term_names, "F value"] *
    residual_anova[term_names, "Df"] /
    (residual_anova[term_names, "F value"] *
       residual_anova[term_names, "Df"] + df_residual),
  row.names = NULL
)
write.csv(residual_effects, "results/Study1_ResidualANOVA_EffectSizes.csv", row.names = FALSE)

## 5. Cell means and 5plq predictions under VV × AV ------------------------
## Compare observed and predicted means across the nine conditions.

# Predict each pair first, then average observed and predicted values by cell.
plot_tbl <- dat %>%
  group_by(VV, AV) %>%
  summarise(
    SE = sd(R_VA) / sqrt(n()),
    R_VA = mean(R_VA),
    Model5plq = mean(est_5plq),
    .groups = "drop"
  )

## 5.4 Observed vs 5plq overlay --------------------------------------------
p_cell_5plq <- ggplot(
  plot_tbl,
  aes(x = VV, group = AV, color = AV)
) +
  geom_line(aes(y = R_VA), linewidth = 1.0) +
  geom_point(aes(y = R_VA), size = 3) +
  geom_errorbar(aes(ymin = R_VA - SE, ymax = R_VA + SE),
                width = 0.08, linewidth = 0.6) +
  geom_point(aes(y = Model5plq),
             shape = 1, size = 3, colour = "black") +
  scale_y_continuous(
    name   = "Valence Rating",
    breaks = SCALE_BREAKS,
    expand = c(0, 0)
  ) +
  coord_cartesian(ylim = c(SCALE_MIN, SCALE_MAX)) +
  scale_x_discrete(
    name   = "Visual Condition",
    labels = c("Neg" = "Negative",
               "Neu" = "Neutral",
               "Pos" = "Positive")
  ) +
  scale_color_manual(
    name   = "Auditory Condition",
    values = c("Neg" = FIGURE_TEAL,
               "Neu" = FIGURE_ORANGE,
               "Pos" = FIGURE_GREEN),
    labels = c("Neg" = "Negative",
               "Neu" = "Neutral",
               "Pos" = "Positive")
  ) +
  theme_bw(base_size = FIGURE_TEXT_PT, base_family = FIGURE_FONT) +
  theme(
    plot.title        = element_blank(),
    axis.title.x      = element_text(size = FIGURE_AXIS_PT, face = "plain",
                                     margin = margin(t = 10)),
    axis.title.y      = element_text(size = FIGURE_AXIS_PT, face = "plain",
                                     margin = margin(r = 10)),
    axis.text.x       = element_text(size = FIGURE_TEXT_PT, face = "plain"),
    axis.text.y       = element_text(size = FIGURE_TEXT_PT),
    legend.title      = element_text(size = FIGURE_TEXT_PT, face = "plain"),
    legend.text       = element_text(size = FIGURE_TEXT_PT),
    legend.position   = "right",
    panel.grid       = element_blank()
  )

p_cell_5plq

save_figure("figures/Fig3_model_validation_cell_means.png",
       p_cell_5plq, width = FIGURE_WIDTH_IN, height = FIGURE_HEIGHT_IN, dpi = FIGURE_DPI)
