# Gao (2018) Information integration theory. Code written following Wilson & Collins (2019).
# Data file: MODELING_Exp1_Valence.xlsx

## 0. Environment and dependencies -----------------------------------------
rm(list = ls())
# Load shared paths relative to this script.
.script_files <- unlist(lapply(sys.frames(), function(frame) frame$ofile))
.script_file <- if (length(.script_files)) tail(.script_files, 1) else
  sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
if (!length(.script_file)) .script_file <- if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) rstudioapi::getSourceEditorContext()$path else "Step1_ModelSimulation.R"
source(file.path(dirname(normalizePath(.script_file)), "..", "Project_Config.R"), encoding = "UTF-8")
setwd(file.path(PROJECT_ROOT, "Scripts", "Study1"))
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# One-time installation of required packages
packages <- c("readxl", "dplyr", "ggplot2", "car", "minpack.lm", "broom", "tidyr")
to_install <- setdiff(packages, rownames(installed.packages()))
if(length(to_install)) install.packages(to_install)

library(readxl)
library(dplyr)
library(ggplot2)
library(car)          # Type-III ANOVA
library(minpack.lm)   # nlsLM for nonlinear least squares
library(broom)
library(tidyr)

## 1. Experimental design --------------------------------------------------
# Load the empirical design used for simulation.
# 180 audiovisual pairs: nine conditions with 20 pair-level means each.
# R_V and R_A are observed components; R_VA is the observed audiovisual outcome.
# Inspect how linear and quadratic weight terms affect simulated judgments.

## Load empirical inputs
dat <- read_excel(file.path(PROJECT_ROOT, "Data", "Study1", "MODELING_Exp1_Valence.xlsx"))

# Assign a pair index.
dat <- dat %>% mutate(Stimcomb = row_number())

# Factorize VV / AV (-1/0/1 -> Neg/Neu/Pos) for ANOVA and plotting
dat <- dat %>%
  mutate(
    VV = factor(VV, levels = c(-1,0,1), labels = c("Neg","Neu","Pos")),
    AV = factor(AV, levels = c(-1,0,1), labels = c("Neg","Neu","Pos"))
  )

# Numeric coding of VV/AV for "slope method" (model-independent sanity check,
# not relying on contrast coding)
dat_num <- dat %>%
  mutate(
    VVn = case_when(VV == "Neg" ~ -1, VV == "Neu" ~ 0, VV == "Pos" ~ 1),
    AVn = case_when(AV == "Neg" ~ -1, AV == "Neu" ~ 0, AV == "Pos" ~ 1)
  )

str(dat); head(dat)

## 2. Specify the candidate model space (Build models) ---------------------
# NOTE: This section enumerates the model family to be compared.

## 3. Simulation – qualitative signatures for 5plq only --------------------
## We focus on the 5plq averaging model and:
## (A) inspect slope patterns (uni vs multi) under different B_VA / B_VA2 settings;
## (B) inspect how NDOM changes with B_VA / B_VA2;
## (C) inspect how VDOM changes with A_V – A_A.
source("simulation_helpers.R")

RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(20251118)

# Design for simulation: keep full structure including VV/AV numeric codes
design_df <- dat_num %>%
  dplyr::select(R_V, R_A, VV, AV, VVn, AVn)

reference_5plq <- reference_params[["5plq"]]

# Parameter grids

# 1) B_VA: seven values within twice the reference magnitude.
B_VA_seq_all <- abs(reference_5plq$B_VA) * c(-2, -1, -0.5, 0, 0.5, 1, 2)

# 2) B_VA2: zero and four positive multiples of the reference.
B_VA2_seq_all <- reference_5plq$B_VA2 * c(0, 0.5, 1, 2, 3)

# Use explicit representative values for the slope plots.
B_VA_levels_slope  <- abs(reference_5plq$B_VA) * c(-2, 0, 2)
B_VA2_levels_slope <- reference_5plq$B_VA2 * c(0, 1, 3)

# For A_V - A_A sweeps, use 5 levels (two negative, zero, two positive)
A_V0        <- reference_5plq$A_V
A_A0        <- reference_5plq$A_A
mean_A_V_A_A <- (A_V0 + A_A0) / 2
max_delta  <- 2  # can adjust if needed
delta_seq  <- seq(-max_delta, max_delta, length.out = 5)  ## 2 negative, 0, 2 positive -----------------------------------------------

noise_sd_slope   <- 0.39
noise_sd_param   <- 0.39
nsim_per_param   <- 100   # for parameter sweeps
nsim_per_slope   <- 100   # for slope patterns

## 3A. Slope patterns (uni vs multi) under different B_VA/B_VA2 settings ---

## Helper to compute four slopes from one simulated dataset
compute_slopes_4 <- function(sim_df) {
  # Visual: R_V vs VVn (uni), R_VA_sim vs VVn (multi)
  lm_R_V_VVn   <- lm(R_V      ~ VVn, data = sim_df)
  lm_R_VA_VVn  <- lm(R_VA_sim ~ VVn, data = sim_df)
  slope_V_uni   <- coef(lm_R_V_VVn)["VVn"]
  slope_V_multi <- coef(lm_R_VA_VVn)["VVn"]

  # Auditory: R_A vs AVn (uni), R_VA_sim vs AVn (multi)
  lm_R_A_AVn   <- lm(R_A      ~ AVn, data = sim_df)
  lm_R_VA_AVn2 <- lm(R_VA_sim ~ AVn, data = sim_df)
  slope_A_uni   <- coef(lm_R_A_AVn)["AVn"]
  slope_A_multi <- coef(lm_R_VA_AVn2)["AVn"]

  tibble::tibble(
    slope_V_uni   = slope_V_uni,
    slope_V_multi = slope_V_multi,
    slope_A_uni   = slope_A_uni,
    slope_A_multi = slope_A_multi
  )
}

slope_results <- list()
idx <- 1L

# --- 3A.1 B_VA manipulation (three levels: negative / zero / positive)

for (j in seq_along(B_VA_levels_slope)) {
  b_val <- B_VA_levels_slope[j]
  cur_par <- reference_5plq
  cur_par$B_VA <- b_val

  label <- dplyr::case_when(
    j == 1 ~ "B_VA_neg",
    j == 2 ~ "B_VA_zero",
    j == 3 ~ "B_VA_pos"
  )

  for (rep_id in seq_len(nsim_per_slope)) {
    sim_df <- simulate_once(
      true_model    = "5plq",
      design_df     = design_df,
      truth_params  = cur_par,
      noise_sd      = noise_sd_slope,
      clip_to_scale = TRUE,
      lower         = SCALE_MIN, upper = SCALE_MAX
    )

    slopes <- compute_slopes_4(sim_df)
    slope_results[[idx]] <- slopes %>%
      dplyr::mutate(
        manip_type   = "B_VA",
        param_level  = label,
        param_value  = b_val,
        rep          = rep_id
      )
    idx <- idx + 1L
  }
}

# --- 3A.2 B_VA2 manipulation (Zero / Low / High)

for (j in seq_along(B_VA2_levels_slope)) {
  b2_val <- B_VA2_levels_slope[j]
  cur_par <- reference_5plq
  cur_par$B_VA2 <- b2_val

  label <- dplyr::case_when(
    j == 1 ~ "B_VA2_low",
    j == 2 ~ "B_VA2_medium",
    j == 3 ~ "B_VA2_high"
  )

  for (rep_id in seq_len(nsim_per_slope)) {
    sim_df <- simulate_once(
      true_model    = "5plq",
      design_df     = design_df,
      truth_params  = cur_par,
      noise_sd      = noise_sd_slope,
      clip_to_scale = TRUE,
      lower         = SCALE_MIN, upper = SCALE_MAX
    )

    slopes <- compute_slopes_4(sim_df)
    slope_results[[idx]] <- slopes %>%
      dplyr::mutate(
        manip_type   = "B_VA2",
        param_level  = label,
        param_value  = b2_val,
        rep          = rep_id
      )
    idx <- idx + 1L
  }
}

slope_results <- dplyr::bind_rows(slope_results)

# Reshape to long format: modality × (uni/multi)
slope_long <- slope_results %>%
  tidyr::pivot_longer(
    cols = dplyr::starts_with("slope_"),
    names_to = "which_slope",
    values_to = "slope"
  ) %>%
  dplyr::mutate(
    modality = dplyr::case_when(
      which_slope %in% c("slope_V_uni","slope_V_multi") ~ "Visual",
      which_slope %in% c("slope_A_uni","slope_A_multi") ~ "Auditory"
    ),
    type = dplyr::case_when(
      grepl("_uni$",   which_slope) ~ "uni",
      grepl("_multi$", which_slope) ~ "multi"
    )
  ) %>%
  dplyr::select(manip_type, param_level, param_value, rep, modality, type, slope)

# Summary across simulation replicates
slope_summary <- slope_long %>%
  dplyr::group_by(manip_type, param_level, modality, type) %>%
  dplyr::summarise(
    mean_slope = mean(slope, na.rm = TRUE),
    sd_slope = sd(slope, na.rm = TRUE),
    .groups    = "drop"
  )

# Plot 1: B_VA manipulation – 3 versions of 5plq, slopes for uni vs multi
library(dplyr)
library(ggplot2)

## --------- Plot 1: b1 manipulation

plot_B_VA <- slope_summary %>%
  filter(manip_type == "B_VA") %>%
  mutate(
    # 4) multi / uni → Multimodal / Unimodal
    type = dplyr::recode(type,
                         "multi" = "Multimodal",
                         "uni"   = "Unimodal"),
    type = factor(type, levels = c("Unimodal", "Multimodal")),
    param_level = factor(
      param_level,
      levels = c("B_VA_neg", "B_VA_zero", "B_VA_pos"),
      labels = c("Negative", "Zero", "Positive")
    )
  )

ggplot(
  plot_B_VA,
  aes(x = type, y = mean_slope,
      color = param_level, linetype = param_level,
      group = param_level)
) +
  geom_line(
    position = position_dodge(width = 0.25),
    linewidth = 1.1
  ) +
  geom_point(
    position = position_dodge(width = 0.25),
    size = 3
  ) +
  facet_wrap(~ modality) +
  labs(
    x        = NULL,
    y        = "Slope",
    color    = expression(b[1] ~ "level"),
    linetype = expression(b[1] ~ "level")
  ) +
  scale_color_manual(values = c(FIGURE_TEAL, FIGURE_ORANGE, FIGURE_GREEN)) +
  theme_bw(base_size = FIGURE_TEXT_PT, base_family = FIGURE_FONT) +
  theme(
    plot.title      = element_blank(),
    axis.text.x     = element_text(size = FIGURE_TEXT_PT),
    axis.text.y     = element_text(size = FIGURE_TEXT_PT),
    axis.title.y    = element_text(size = FIGURE_AXIS_PT, face = "plain",
                                   margin = margin(r = 10)),
    strip.text      = element_text(size = FIGURE_TEXT_PT, face = "plain"),
    strip.background = element_blank(),
    legend.title    = element_text(size = FIGURE_TEXT_PT),
    legend.text     = element_text(size = FIGURE_TEXT_PT),
    panel.grid      = element_blank(),
    legend.position = "right"
  )

file_name <- "figures/Fig1_slope_b1.png"
save_figure(file_name, width = 10, height = 5, dpi = FIGURE_DPI)

## --------- Plot 2: b2 manipulation
plot_B_VA2 <- slope_summary %>%
  filter(manip_type == "B_VA2") %>%
  mutate(
    type = dplyr::recode(type,
                         "multi" = "Multimodal",
                         "uni"   = "Unimodal"),
    type = factor(type, levels = c("Unimodal", "Multimodal")),
    param_level = factor(
      param_level,
      levels = c("B_VA2_low", "B_VA2_medium", "B_VA2_high"),
      labels = c("Zero", "Low", "High")
    )
  )

ggplot(
  plot_B_VA2,
  aes(x = type, y = mean_slope,
      color = param_level, linetype = param_level,
      group = param_level)
) +
  geom_line(
    position = position_dodge(width = 0.25),
    linewidth = 1.1
  ) +
  geom_point(
    position = position_dodge(width = 0.25),
    size = 3
  ) +
  facet_wrap(~ modality) +
  labs(
    x        = NULL,
    y        = "Slope",
    color    = expression(b[2] ~ "level"),
    linetype = expression(b[2] ~ "level")
  ) +
  scale_color_manual(values = c(FIGURE_TEAL, FIGURE_ORANGE, FIGURE_GREEN)) +
  theme_bw(base_size = FIGURE_TEXT_PT, base_family = FIGURE_FONT) +
  theme(
    plot.title      = element_blank(),
    axis.text.x     = element_text(size = FIGURE_TEXT_PT),
    axis.text.y     = element_text(size = FIGURE_TEXT_PT),
    axis.title.y    = element_text(size = FIGURE_AXIS_PT, face = "plain",
                                   margin = margin(r = 10)),
    strip.text      = element_text(size = FIGURE_TEXT_PT, face = "plain"),
    strip.background = element_blank(),
    legend.title    = element_text(size = FIGURE_TEXT_PT),
    legend.text     = element_text(size = FIGURE_TEXT_PT),
    panel.grid      = element_blank(),
    legend.position = "right"
  )

file_name <- "figures/Fig1_slope_b2.png"
save_figure(file_name, width = 10, height = 5, dpi = FIGURE_DPI)

## 3B. Negativity dominance under seven B_VA and five B_VA2 values ---

compute_NDOM <- function(sim_df) {
  V_pos_A_neg <- mean(sim_df$R_VA_sim[sim_df$VV == "Pos" & sim_df$AV == "Neg"], na.rm = TRUE)
  V_neg_A_pos <- mean(sim_df$R_VA_sim[sim_df$VV == "Neg" & sim_df$AV == "Pos"], na.rm = TRUE)
  V_pos_A_pos <- mean(sim_df$R_VA_sim[sim_df$VV == "Pos" & sim_df$AV == "Pos"], na.rm = TRUE)
  V_neg_A_neg <- mean(sim_df$R_VA_sim[sim_df$VV == "Neg" & sim_df$AV == "Neg"], na.rm = TRUE)
  num_ndom <- 0.5 * (V_pos_A_neg + V_neg_A_pos) - V_pos_A_pos
  den_ndom <- V_neg_A_neg - V_pos_A_pos
  if (!is.na(den_ndom) && abs(den_ndom) > .Machine$double.eps) {
    num_ndom / den_ndom
  } else {
    NA_real_
  }
}

ndom_B_VA_res  <- list()
ndom_B_VA2_res <- list()
idx <- 1L

# Sweep B_VA (7 levels: 3 negative, 0, 3 positive)
for (b in B_VA_seq_all) {
  nd_vals <- numeric(nsim_per_param)
  cur_par <- reference_5plq
  cur_par$B_VA <- b

  for (rep_id in seq_len(nsim_per_param)) {
    sim_df <- simulate_once(
      true_model    = "5plq",
      design_df     = design_df,
      truth_params  = cur_par,
      noise_sd      = noise_sd_param,
      clip_to_scale = TRUE,
      lower         = SCALE_MIN, upper = SCALE_MAX
    )
    nd_vals[rep_id] <- compute_NDOM(sim_df)
  }

  ndom_B_VA_res[[idx]] <- tibble::tibble(
    B_VA  = b,
    NDOM  = nd_vals
  )
  idx <- idx + 1L
}
ndom_B_VA_res <- dplyr::bind_rows(ndom_B_VA_res)

# Sweep B_VA2 (five levels: 0, 0.5, 1, 2, and 3 times the reference).
idx <- 1L
for (b2 in B_VA2_seq_all) {
  nd_vals <- numeric(nsim_per_param)
  cur_par <- reference_5plq
  cur_par$B_VA2 <- b2

  for (rep_id in seq_len(nsim_per_param)) {
    sim_df <- simulate_once(
      true_model    = "5plq",
      design_df     = design_df,
      truth_params  = cur_par,
      noise_sd      = noise_sd_param,
      clip_to_scale = TRUE,
      lower         = SCALE_MIN, upper = SCALE_MAX
    )
    nd_vals[rep_id] <- compute_NDOM(sim_df)
  }

  ndom_B_VA2_res[[idx]] <- tibble::tibble(
    B_VA2 = b2,
    NDOM  = nd_vals
  )
  idx <- idx + 1L
}
ndom_B_VA2_res <- dplyr::bind_rows(ndom_B_VA2_res)

# Summaries for plotting; error bars show ±1 SD across simulation replicates.
ndom_B_VA_summary <- ndom_B_VA_res %>%
  dplyr::group_by(B_VA) %>%
  dplyr::summarise(
    mean_NDOM = mean(NDOM, na.rm = TRUE),
    sd_NDOM = sd(NDOM, na.rm = TRUE),
    .groups   = "drop"
  )

ndom_B_VA2_summary <- ndom_B_VA2_res %>%
  dplyr::group_by(B_VA2) %>%
  dplyr::summarise(
    mean_NDOM = mean(NDOM, na.rm = TRUE),
    sd_NDOM = sd(NDOM, na.rm = TRUE),
    .groups   = "drop"
  )

## 1. Negativity dominance as a function of b1 -----------------------------

p_ndom_b1 <- ggplot(ndom_B_VA_summary,
                    aes(x = B_VA, y = mean_NDOM)) +
  geom_hline(yintercept = 0.5,
             linetype = "dashed",
             linewidth = 0.6) +
  geom_line(linewidth = 1.1, color = FIGURE_TEAL) +
  geom_point(size = 3, color = FIGURE_TEAL) +
  geom_errorbar(
    aes(ymin = mean_NDOM - sd_NDOM,
        ymax = mean_NDOM + sd_NDOM),
    width = 0.06,
    linewidth = 0.7,
    color = FIGURE_TEAL
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = expression(b[1]),
    y = "NDOM"
  ) +
  theme_bw(base_size = FIGURE_TEXT_PT, base_family = FIGURE_FONT) +
  theme(
    plot.title = element_blank(),
    axis.title.x = element_text(size = FIGURE_AXIS_PT, face = "plain", margin = margin(t = 10)),
    axis.title.y = element_text(size = FIGURE_AXIS_PT, face = "plain", margin = margin(r = 10)),
    axis.text.x  = element_text(size = FIGURE_TEXT_PT),
    axis.text.y  = element_text(size = FIGURE_TEXT_PT),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 0.7),
    legend.position = "none"
  )

print(p_ndom_b1)


## 2. Negativity dominance as a function of b2 (with error bar) ------------

p_ndom_b2 <- ggplot(ndom_B_VA2_summary,
                    aes(x = B_VA2, y = mean_NDOM)) +
  geom_hline(yintercept = 0.5,
             linetype = "dashed",
             linewidth = 0.6) +
  geom_line(linewidth = 1.1, color = FIGURE_TEAL) +
  geom_point(size = 3, color = FIGURE_TEAL) +
  geom_errorbar(
    aes(ymin = mean_NDOM - sd_NDOM,
        ymax = mean_NDOM + sd_NDOM),
    width = 0.06 * diff(range(B_VA2_seq_all)) / diff(range(B_VA_seq_all)),
    linewidth = 0.7,
    color = FIGURE_TEAL
  ) +
  scale_y_continuous(breaks = seq(0.5, 1, by = 0.1)) +
  coord_cartesian(ylim = c(0.5, 1), expand = FALSE) +
  labs(
    x = expression(b[2]),
    y = "NDOM"
  ) +
  theme_bw(base_size = FIGURE_TEXT_PT, base_family = FIGURE_FONT) +
  theme(
    plot.title = element_blank(),
    axis.title.x = element_text(size = FIGURE_AXIS_PT, face = "plain", margin = margin(t = 10)),
    axis.title.y = element_text(size = FIGURE_AXIS_PT, face = "plain", margin = margin(r = 10)),
    axis.text.x  = element_text(size = FIGURE_TEXT_PT),
    axis.text.y  = element_text(size = FIGURE_TEXT_PT),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 0.7),
    legend.position = "none"
  )

print(p_ndom_b2)


## 3C. Visual dominance as a function of A_V – A_A (5 levels) --------------

compute_VDOM <- function(sim_df) {
  V_pos_A_all  <- mean(sim_df$R_VA_sim[sim_df$VV == "Pos"], na.rm = TRUE)
  V_neg_A_all  <- mean(sim_df$R_VA_sim[sim_df$VV == "Neg"], na.rm = TRUE)
  V_all_A_pos <- mean(sim_df$R_VA_sim[sim_df$AV == "Pos"], na.rm = TRUE)
  V_all_A_neg <- mean(sim_df$R_VA_sim[sim_df$AV == "Neg"], na.rm = TRUE)
  num_vdom <- (V_pos_A_all - V_neg_A_all)
  den_vdom <- (V_pos_A_all - V_neg_A_all) + (V_all_A_pos - V_all_A_neg)
  if (!is.na(den_vdom) && abs(den_vdom) > .Machine$double.eps) {
    num_vdom / den_vdom
  } else {
    NA_real_
  }
}

vdom_res <- list()
idx <- 1L

for (d in delta_seq) {
  vdom_vals <- numeric(nsim_per_param)
  cur_par    <- reference_5plq

  # Keep average(A_V, A_A) constant, vary A_V - A_A = d
  cur_par$A_V <- mean_A_V_A_A + d/2
  cur_par$A_A <- mean_A_V_A_A - d/2

  for (rep_id in seq_len(nsim_per_param)) {
    sim_df <- simulate_once(
      true_model    = "5plq",
      design_df     = design_df,
      truth_params  = cur_par,
      noise_sd      = noise_sd_param,
      clip_to_scale = TRUE,
      lower         = SCALE_MIN, upper = SCALE_MAX
    )
    vdom_vals[rep_id] <- compute_VDOM(sim_df)
  }

  vdom_res[[idx]] <- tibble::tibble(
    delta_A_V_A_A = d,
    VDOM        = vdom_vals
  )
  idx <- idx + 1L
}

vdom_res <- dplyr::bind_rows(vdom_res)

# Error bars show ±1 SD across simulation replicates.
vdom_summary <- vdom_res %>%
  dplyr::group_by(delta_A_V_A_A) %>%
  dplyr::summarise(
    mean_VDOM = mean(VDOM, na.rm = TRUE),
    sd_VDOM = sd(VDOM, na.rm = TRUE),
    .groups   = "drop"
  )

p_vdom <- ggplot(vdom_summary, aes(x = delta_A_V_A_A, y = mean_VDOM)) +
  geom_hline(yintercept = 0.5, linetype = 2, color = "gray50", linewidth = 0.7) +
  geom_line(linewidth = 1.1, color = FIGURE_TEAL) +
  geom_point(size = 3, color = FIGURE_TEAL) +
  geom_errorbar(
    aes(ymin = mean_VDOM - sd_VDOM,
        ymax = mean_VDOM + sd_VDOM),
    width = 0.08,
    linewidth = 0.8,
    color = FIGURE_TEAL
  ) +
  labs(
    x = expression(a[V] - a[A]),
    y = "VDOM"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw(base_size = FIGURE_TEXT_PT, base_family = FIGURE_FONT) +
  theme(
    plot.title   = element_blank(),
    axis.title.x = element_text(size = FIGURE_AXIS_PT, face = "plain", margin = margin(t = 10)),
    axis.title.y = element_text(size = FIGURE_AXIS_PT, face = "plain", margin = margin(r = 10)),
    axis.text    = element_text(size = FIGURE_TEXT_PT),
    panel.grid = element_blank()
  )


## 3D. Weight functions for the representative parameter settings ---------

# Use the same raw visual-weight function as Step7_PlotWeightFunctions.R.
W_V_5plq <- function(x, A_V, B_VA, B_VA2) {
  x_c <- x - SCALE_MID
  A_V + B_VA * x_c + B_VA2 * x_c^2
}

newx <- tibble::tibble(x = seq(SCALE_MIN, SCALE_MAX, length.out = 161))

# Solid lines cover observed inputs; dashed lines show extrapolation.
make_weight_segments <- function(modality, weight_function, observed_range) {
  bounds <- c(SCALE_MIN, observed_range, SCALE_MAX)
  bind_rows(lapply(seq_len(3), function(segment) {
    lower <- bounds[segment]
    upper <- bounds[segment + 1]
    if (upper <= lower) return(NULL)
    x <- sort(unique(c(lower, newx$x[newx$x > lower & newx$x < upper], upper)))
    tibble::tibble(x = x, weight = weight_function(x), modality = modality,
                  input_range = if (segment == 2) "Observed range" else "Extrapolation",
                  segment = segment)
  }))
}

plot_simulated_weights <- function(parameter, values, limits, legend_labels) {
  observed <- range(design_df$R_V)
  curves <- list()
  minima <- list()
  for (i in seq_along(values)) {
    par <- reference_5plq
    par[[parameter]] <- values[i]
    weight_function <- function(x) W_V_5plq(x, par$A_V, par$B_VA, par$B_VA2)
    curve <- make_weight_segments("Visual", weight_function, observed)
    curve$param_value <- values[i]
    curves[[i]] <- curve
    turning <- if (par$B_VA2 == 0) NA_real_ else SCALE_MID - par$B_VA / (2 * par$B_VA2)
    minima[[i]] <- tibble::tibble(
      param_value = values[i], turning_rating = turning,
      turning_weight = if (is.finite(turning)) weight_function(turning) else NA_real_,
      within_scale = is.finite(turning) && turning >= SCALE_MIN && turning <= SCALE_MAX,
      within_observed = is.finite(turning) && turning >= observed[1] && turning <= observed[2]
    )
  }
  curves <- bind_rows(curves)
  minima <- bind_rows(minima)
  curves$param_level <- factor(curves$param_value, levels = values)
  curves$input_range <- factor(curves$input_range, levels = c("Observed range", "Extrapolation"))
  marks <- minima[minima$within_scale, ]
  marks$param_level <- factor(marks$param_value, levels = values)
  marks$input_range <- ifelse(marks$within_observed, "Observed range", "Extrapolation")
  plot <- ggplot(curves,
    aes(x = x, y = weight, colour = param_level, linetype = input_range,
        group = interaction(param_level, segment))) +
    geom_line(linewidth = 1.2) +
    geom_point(data = marks,
      aes(x = turning_rating, y = turning_weight, colour = param_level, shape = input_range),
      inherit.aes = FALSE, size = 3, stroke = 1, show.legend = FALSE) +
    scale_colour_manual(values = c(FIGURE_TEAL, FIGURE_ORANGE, FIGURE_GREEN),
                        labels = legend_labels) +
    scale_linetype_manual(values = c("Observed range" = "solid", Extrapolation = "dashed")) +
    scale_shape_manual(values = c("Observed range" = 16, Extrapolation = 1)) +
    scale_x_continuous(breaks = seq(SCALE_MIN, SCALE_MAX, by = 2), expand = c(0, 0)) +
    coord_cartesian(xlim = c(SCALE_MIN, SCALE_MAX), ylim = limits) +
    labs(x = "Unimodal Valence Rating", y = "Weight") +
    theme_bw(base_size = FIGURE_TEXT_PT, base_family = FIGURE_FONT) +
    theme(legend.position = "top", legend.direction = "horizontal",
          panel.grid = element_blank())
  plot <- style_figure(plot)
  plot <- plot + guides(
    colour = guide_legend(title = NULL, nrow = 1,
                         override.aes = list(linetype = "solid", shape = NA)),
    linetype = "none", shape = "none")
  list(plot = plot, curves = curves, minima = minima)
}

weights_b1 <- plot_simulated_weights("B_VA", B_VA_levels_slope, c(0, 7.5),
                                     c("Negative", "Zero", "Positive"))
weights_b2 <- plot_simulated_weights("B_VA2", B_VA2_levels_slope, c(-1, 10),
                                     c("Zero", "Low", "High"))

## 3E. Baseline weights across the modality-difference settings -----------

baseline_weights <- bind_rows(
  tibble(delta = delta_seq, modality = "Visual", weight = mean_A_V_A_A + delta_seq / 2),
  tibble(delta = delta_seq, modality = "Auditory", weight = mean_A_V_A_A - delta_seq / 2)
)
baseline_weights$modality <- factor(baseline_weights$modality,
                                    levels = c("Visual", "Auditory"))

# At the scale midpoint, the raw weights equal the baseline coefficients.
p_baseline_weights <- ggplot(baseline_weights,
                            aes(x = delta, y = weight, colour = modality)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_colour_manual(values = FIGURE_MODALITY_COLOURS, name = NULL) +
  scale_x_continuous(breaks = delta_seq) +
  scale_y_continuous(breaks = seq(0, 2.5, by = 0.5)) +
  coord_cartesian(ylim = c(0, 2.5)) +
  labs(x = expression(a[V] - a[A]), y = "Weight") +
  theme_bw(base_size = FIGURE_TEXT_PT, base_family = FIGURE_FONT) +
  theme(legend.position = "top", legend.direction = "horizontal",
          panel.grid = element_blank()) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(linetype = "solid", shape = NA)))


# Align plotting regions before assembling the panels in Visio.
save_simulation_panels <- function(plots) {
  plots <- lapply(plots, style_figure)
  grobs <- lapply(plots, ggplotGrob)
  shared_widths <- do.call(grid::unit.pmax, lapply(grobs, function(grob) grob$widths))
  shared_heights <- do.call(grid::unit.pmax, lapply(grobs, function(grob) grob$heights))
  for (name in names(plots)) {
    grob <- grobs[[name]]
    grob$widths <- shared_widths
    grob$heights <- shared_heights
    filename <- file.path("figures", paste0(name, ".png"))
    export_figure(filename, grob, width = FIGURE_WIDTH_IN,
                   height = FIGURE_HEIGHT_IN, dpi = FIGURE_DPI)
    saveRDS(plots[[name]], sub("[.]png$", ".rds", filename))
  }
  invisible(plots)
}

save_simulation_panels(list(
  Fig1_weights_b1 = weights_b1$plot,
  Fig1_NDOM_b1 = p_ndom_b1,
  Fig1_weights_b2 = weights_b2$plot,
  Fig1_NDOM_b2 = p_ndom_b2,
  Fig1_baseline_weights = p_baseline_weights,
  Fig1_VDOM = p_vdom
))

## 4. Save the simulation outputs ----------------------------------------
write.csv(slope_results, "results/simulation_slopes_replicates.csv", row.names = FALSE)
write.csv(slope_summary, "results/simulation_slopes_summary.csv", row.names = FALSE)
write.csv(ndom_B_VA_summary, "results/simulation_NDOM_b1_summary.csv", row.names = FALSE)
write.csv(ndom_B_VA2_summary, "results/simulation_NDOM_b2_summary.csv", row.names = FALSE)
write.csv(vdom_summary, "results/simulation_VDOM_summary.csv", row.names = FALSE)
write.csv(baseline_weights, "results/simulation_baseline_weights.csv", row.names = FALSE)
write.csv(weights_b1$curves, "results/simulation_weights_b1.csv", row.names = FALSE)
write.csv(weights_b2$curves, "results/simulation_weights_b2.csv", row.names = FALSE)
write.csv(weights_b1$minima, "results/simulation_weight_minima_b1.csv", row.names = FALSE)
write.csv(weights_b2$minima, "results/simulation_weight_minima_b2.csv", row.names = FALSE)
saveRDS(list(reference = reference_5plq,
             grids = list(b1 = B_VA_seq_all, b2 = B_VA2_seq_all, delta = delta_seq),
             slope_grids = list(b1 = B_VA_levels_slope, b2 = B_VA2_levels_slope),
             slope_results = slope_results, slope_summary = slope_summary,
             ndom_b1 = ndom_B_VA_res, ndom_b2 = ndom_B_VA2_res, vdom = vdom_res,
             ndom_b1_summary = ndom_B_VA_summary, ndom_b2_summary = ndom_B_VA2_summary,
             vdom_summary = vdom_summary), "results/simulation_outputs.rds")
