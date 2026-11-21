# Gao (2018) Information integration theory
# Data file: MODELING_Exp1_Valence.xlsx

## 0. Environment and dependencies -----------------------------------------
rm(list = ls())
# Load shared paths relative to this script.
.script_files <- unlist(lapply(sys.frames(), function(frame) frame$ofile))
.script_file <- if (length(.script_files)) tail(.script_files, 1) else
  sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
if (!length(.script_file)) .script_file <- if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) rstudioapi::getSourceEditorContext()$path else "Step2_AveragingVSaddition.R"
source(file.path(dirname(normalizePath(.script_file)), "..", "Project_Config.R"), encoding = "UTF-8")
setwd(file.path(PROJECT_ROOT, "Scripts", "Study1"))
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# One-time installation of required packages
packages <- c("readxl", "dplyr", "ggplot2", "car", "broom")
to_install <- setdiff(packages, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)

library(readxl)
library(dplyr)
library(ggplot2)
library(car)      # Type-III ANOVA
library(broom)

# Load the common rating-scale settings
source("simulation_helpers.R")

## 1. Load empirical data --------------------------------------------------
# NOTE: This script only does model-independent analyses (no simulations, no model fitting).

dat <- read_excel(file.path(PROJECT_ROOT, "Data", "Study1", "MODELING_Exp1_Valence.xlsx"))

# Assign a pair index.
dat <- dat %>% mutate(Stimcomb = dplyr::row_number())

# Factorize VV / AV (-1/0/1 -> Neg/Neu/Pos) for ANOVA and plotting
dat <- dat %>%
  mutate(
    VV = factor(VV, levels = c(-1, 0, 1), labels = c("Neg","Neu","Pos")),
    AV = factor(AV, levels = c(-1, 0, 1), labels = c("Neg","Neu","Pos"))
  )

# Numeric coding of VV/AV for the "slope method"
dat_num <- dat %>%
  mutate(
    VVn = dplyr::case_when(
      VV == "Neg" ~ -1,
      VV == "Neu" ~  0,
      VV == "Pos" ~  1
    ),
    AVn = dplyr::case_when(
      AV == "Neg" ~ -1,
      AV == "Neu" ~  0,
      AV == "Pos" ~  1
    )
  )

str(dat)
head(dat)

## 2. Model-independent analysis -------------------------------------------
## This corresponds to Wilson & Collins' recommendation to first analyze behavior
## without models, to check that the process of interest is engaged and that
## qualitative signatures are present in the data.

## 2.1 pair-level ANOVA: R_VA ~ VV * AV ------------------------------------

# Local sum coding for the Type III descriptive ANOVA.
contrasts(dat$VV) <- contr.sum(nlevels(dat$VV))
contrasts(dat$AV) <- contr.sum(nlevels(dat$AV))
aov_item <- aov(R_VA ~ VV * AV, data = dat)
anova_res <- car::Anova(aov_item, type = "III")  # Describe visual, auditory and interaction effects before model fitting.
print(anova_res)

## 2.2 Paired slope differences at the pair level --------------------------
# Each component and joint rating refer to the same audiovisual-pair row.
# Repeated component stimuli can induce dependence across rows; these are
# pair-level OLS tests, not participant-level or cluster-adjusted inference.
lm_R_V_VVn  <- lm(R_V  ~ VVn, data = dat_num)
lm_R_VA_VVn <- lm(R_VA ~ VVn, data = dat_num)
lm_R_A_AVn  <- lm(R_A  ~ AVn, data = dat_num)
lm_R_VA_AVn <- lm(R_VA ~ AVn, data = dat_num)

lm_diff_V <- lm(I(R_V - R_VA) ~ VVn, data = dat_num)
lm_diff_A <- lm(I(R_A - R_VA) ~ AVn, data = dat_num)

visual <- summary(lm_diff_V)$coefficients["VVn", ]
auditory <- summary(lm_diff_A)$coefficients["AVn", ]
slope_comparison <- data.frame(
  channel = c("Visual", "Auditory"),
  slope_unimodal = c(coef(lm_R_V_VVn)["VVn"], coef(lm_R_A_AVn)["AVn"]),
  slope_audiovisual = c(coef(lm_R_VA_VVn)["VVn"], coef(lm_R_VA_AVn)["AVn"]),
  difference = c(visual["Estimate"], auditory["Estimate"]),
  SE = c(visual["Std. Error"], auditory["Std. Error"]),
  t = c(visual["t value"], auditory["t value"]),
  df = c(df.residual(lm_diff_V), df.residual(lm_diff_A)),
  p = c(visual["Pr(>|t|)"], auditory["Pr(>|t|)"])
)
print(slope_comparison, row.names = FALSE)
write.csv(slope_comparison, "results/Study1_SlopeDifference_Regression.csv",
          row.names = FALSE)

## 2.4 quick visualization of slopes ---------------------------------------
# Simple diagnostic plots: mean R_VA / R_V / R_A across VV / AV

# Visual channel: R_V vs R_VA across VV
vv_plot_dat <- dat_num %>%
  group_by(VV) %>%
  summarise(
    R_V_mean  = mean(R_V,  na.rm = TRUE),
    R_VA_mean = mean(R_VA, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  tidyr::pivot_longer(cols = c(R_V_mean, R_VA_mean),
                      names_to = "series",
                      values_to = "mean_rating")

vv_plot_dat2 <- vv_plot_dat %>%
  mutate(series = dplyr::recode(series,
                                "R_V_mean" = "Visual",
                                "R_VA_mean" = "Audiovisual"))

vv_plot_dat2$VV <- factor(vv_plot_dat2$VV,
                          levels = c("Neg", "Neu", "Pos"),
                          labels = c("Negative", "Neutral", "Positive"))

ggplot(vv_plot_dat2,
       aes(x = VV, y = mean_rating,
           group = series, color = series)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3) +
  scale_y_continuous(limits = c(SCALE_MIN, SCALE_MAX),
                     breaks = SCALE_BREAKS, expand = c(0, 0)) +
  scale_color_manual(values = c("Visual" = FIGURE_ORANGE,
                                "Audiovisual" = FIGURE_TEAL)) +
  labs(x = "Visual Condition",
       y = "Valence Rating",
       color = "") +
  theme_bw(base_size = FIGURE_TEXT_PT, base_family = FIGURE_FONT) +
  theme(
    plot.title      = element_blank(),
    axis.title.x    = element_text(size = FIGURE_AXIS_PT, face = "plain",
                                   margin = margin(t = 10)),
    axis.title.y    = element_text(size = FIGURE_AXIS_PT, face = "plain",
                                   margin = margin(r = 10)),
    axis.text.x     = element_text(size = FIGURE_TEXT_PT, face = "plain"),
    axis.text.y     = element_text(size = FIGURE_TEXT_PT),
    legend.text     = element_text(size = FIGURE_TEXT_PT),
    legend.position = "top",
    panel.grid = element_blank()
  )

save_figure("figures/Fig2_Slope_Visual.png", width = FIGURE_WIDTH_IN, height = FIGURE_HEIGHT_IN, dpi = FIGURE_DPI)

# Auditory channel: R_A vs R_VA across AV
av_plot_dat <- dat_num %>%
  dplyr::group_by(AV) %>%
  dplyr::summarise(
    R_A_mean  = mean(R_A,  na.rm = TRUE),
    R_VA_mean = mean(R_VA, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  tidyr::pivot_longer(
    cols      = c(R_A_mean, R_VA_mean),
    names_to  = "series",
    values_to = "mean_rating"
  )

# Rename for publication-friendly legend
av_plot_dat2 <- av_plot_dat %>%
  dplyr::mutate(
    series = dplyr::recode(
      series,
      "R_A_mean"  = "Auditory",
      "R_VA_mean" = "Audiovisual"
    )
  )

ggplot(av_plot_dat2,
       aes(x = AV, y = mean_rating,
           group = series, color = series)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3) +
  scale_y_continuous(limits = c(SCALE_MIN, SCALE_MAX),
                     breaks = SCALE_BREAKS, expand = c(0, 0)) +
  scale_color_manual(
    values = c("Auditory"   = FIGURE_GREEN,
               "Audiovisual" = FIGURE_TEAL),
    breaks = c("Audiovisual", "Auditory")  # Legend order
  ) +
  scale_x_discrete(
    labels = c("Neg" = "Negative",
               "Neu" = "Neutral",
               "Pos" = "Positive")
  ) +
  labs(x = "Auditory Condition",
       y = "Valence Rating",
       color = "") +
  theme_bw(base_size = FIGURE_TEXT_PT, base_family = FIGURE_FONT) +
  theme(
    plot.title       = element_blank(),
    axis.title.x     = element_text(size = FIGURE_AXIS_PT, face = "plain",
                                    margin = margin(t = 10)),
    axis.title.y     = element_text(size = FIGURE_AXIS_PT, face = "plain",
                                    margin = margin(r = 10)),
    axis.text.x      = element_text(size = FIGURE_TEXT_PT, face = "plain"),
    axis.text.y      = element_text(size = FIGURE_TEXT_PT),
    legend.text      = element_text(size = FIGURE_TEXT_PT),
    legend.position  = "top",
    panel.grid = element_blank()
  )

save_figure("figures/Fig2_Slope_Auditory.png", width = FIGURE_WIDTH_IN, height = FIGURE_HEIGHT_IN, dpi = FIGURE_DPI)

