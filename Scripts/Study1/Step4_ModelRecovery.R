# Gao (2018) Information integration theory
# Model recovery
# Candidate models: 3p, 4p, 4plq, 5p, 5plq, 7p
# Model recovery procedure:
# - For each model, simulate many datasets using random
# parameter draws sampled uniformly from broad ranges.
# - Fit all candidate models to each simulated dataset.
# - Build confusion matrices (AIC/BIC) showing how often
# each model is recovered when it was the true generative model.

## 0. Environment and dependencies -----------------------------------------
rm(list = ls())
# Load shared paths relative to this script.
.script_files <- unlist(lapply(sys.frames(), function(frame) frame$ofile))
.script_file <- if (length(.script_files)) tail(.script_files, 1) else
  sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
if (!length(.script_file)) .script_file <- if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) rstudioapi::getSourceEditorContext()$path else "Step4_ModelRecovery.R"
source(file.path(dirname(normalizePath(.script_file)), "..", "Project_Config.R"), encoding = "UTF-8")
setwd(file.path(PROJECT_ROOT, "Scripts", "Study1"))
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# Required packages
packages <- c("readxl", "dplyr", "tidyr", "purrr", "tibble", "ggplot2")
to_install <- setdiff(packages, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)

library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ggplot2)

# Load helper functions:
# - wavg_* models
# - reference_params
# - generator_dict
# - simulate_once()
# - fit_all_models()
# - sample_uniform_params()
source("simulation_helpers.R")

RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(20251105)

## 1. Load empirical data and construct design matrix ----------------------

dat <- read_excel(file.path(PROJECT_ROOT, "Data", "Study1", "MODELING_Exp1_Valence.xlsx"))

# Preserve original stimulus structure (used as design for simulation)
design_df <- dat %>%
  dplyr::select(R_V, R_A, VV, AV)

str(design_df)

## 2. Configuration for model recovery -------------------------------------

# List of generative models to be tested
true_models <- c("3p", "4p", "4plq", "5p", "5plq", "7p")

# Number of synthetic experiments per generative model
nsim_per_model <- 200

# Measurement noise (kept constant; not treated as an experimental factor)
noise_sd <- 0.39  # Main-text SD; use 0.20 for the supplementary recovery analysis.

# Jitter scale for model-fitting initial values
jitter_scale <- 0.2

# Width of uniform parameter sampling range (relative to empirical parameters)
width_factor <- 0.2   # e.g., ±20% around default values
abs_min      <- 0.05   # coefficient-scale floor; minimum half-width is 0.01

## 3. Simulation loop: generate datasets + fit all models ------------------

sim_rows <- list()
idx <- 1L

for (m in true_models) {
  message("Simulating for true model: ", m)

  for (r in seq_len(nsim_per_model)) {

    ## 3.1 Draw a random parameter vector for the generative model -------------
    truth_par <- sample_uniform_params(
      model                = m,
      reference_params = reference_params,
      width_factor         = width_factor,
      abs_min              = abs_min
    )

    ## 3.2 Simulate one dataset using the generative model + random parameters ---
    dat_sim <- simulate_once(
      true_model    = m,
      design_df     = design_df,
      truth_params  = truth_par,
      noise_sd      = noise_sd,
      clip_to_scale = TRUE,
      lower         = SCALE_MIN,
      upper         = SCALE_MAX
    )

    ## 3.3 Fit all candidate models to this synthetic dataset ------------------
    fits_tbl <- fit_all_models(
      df               = dat_sim,
      start_params = start_params,
      jitter_scale     = jitter_scale
    ) %>%
      dplyr::mutate(
        true_model = m,
        sim_rep    = r
      )

    sim_rows[[idx]] <- fits_tbl
    idx <- idx + 1L
  }
}

# Combine results across all simulations
sim_summary <- dplyr::bind_rows(sim_rows)
saveRDS(sim_summary, "results/Study1_ModelRecovery_AllFits.rds")
convergence <- sim_summary %>% group_by(true_model, model) %>%
  summarise(attempted = n(), converged = sum(converged),
            convergence_rate = mean(converged), .groups = "drop")
write.csv(convergence, "results/Study1_ModelRecovery_Convergence.csv", row.names = FALSE)

# Standardize factor levels for consistent ordering in confusion matrices
sim_summary <- sim_summary %>%
  dplyr::mutate(
    true_model = factor(true_model, levels = true_models),
    model      = factor(model,      levels = true_models)
  )

## 4. Model recovery: confusion matrices -----------------------------------

# Confusion matrix: rows = recovered (best-fitting) model,
# columns = true (generative) model

## 4.1 AIC-based recovery --------------------------------------------------

model_recovery_AIC <- sim_summary %>%
  dplyr::filter(is.finite(AIC)) %>%
  dplyr::group_by(true_model, sim_rep) %>%
  dplyr::slice_min(AIC, with_ties = FALSE) %>%   # Best-fitting model per simulation
  dplyr::ungroup() %>%
  dplyr::count(true_model, selected_model = model) %>%
  # Fill in all true × selected combinations (missing cells = 0)
  tidyr::complete(
    true_model     = factor(true_models, levels = true_models),
    selected_model = factor(true_models, levels = true_models),
    fill = list(n = 0L)
  ) %>%
  dplyr::group_by(true_model) %>%
  dplyr::mutate(prop = n / nsim_per_model) %>%           # Normalize within each true model
  dplyr::ungroup()

model_recovery_AIC

## 4.2 BIC-based recovery --------------------------------------------------

model_recovery_BIC <- sim_summary %>%
  dplyr::filter(is.finite(BIC)) %>%
  dplyr::group_by(true_model, sim_rep) %>%
  dplyr::slice_min(BIC, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::count(true_model, selected_model = model) %>%
  tidyr::complete(
    true_model     = factor(true_models, levels = true_models),
    selected_model = factor(true_models, levels = true_models),
    fill = list(n = 0L)
  ) %>%
  dplyr::group_by(true_model) %>%
  dplyr::mutate(prop = n / nsim_per_model) %>%
  dplyr::ungroup()

model_recovery_BIC

## 5. Plot confusion matrices ----------------------------------------------

# Save model-recovery proportions for plotting.
write.csv(model_recovery_AIC, "results/Plot_ModelRecovery_AIC.csv", row.names = FALSE)
write.csv(model_recovery_BIC, "results/Plot_ModelRecovery_BIC.csv", row.names = FALSE)

# AIC model-recovery matrix

p_conf_AIC <- ggplot(model_recovery_AIC,
                     aes(x = true_model,
                         y = selected_model,
                         fill = prop)) +
  geom_tile(color = "grey40", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.2f", prop),
                colour = ifelse(prop >= .65, "white", "black")),
            size = 4.0, fontface = "plain") +
  scale_color_identity(guide = "none") +
  scale_fill_gradientn(
    colours = c("white", FIGURE_TEAL),
    limits  = c(0, 1),
    breaks  = seq(0, 1, by = 0.25),
    name    = "p(Rec|True)",
    guide   = guide_colourbar(display = "rectangles")
  ) +
  labs(x = "True Model", y = "Recovered Model") +
  coord_equal(ratio = 0.85) +
  theme_bw(base_size = FIGURE_TEXT_PT, base_family = FIGURE_FONT) +
  theme(
    plot.title      = element_blank(),
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 13, face = "plain"),
    axis.text.y     = element_text(size = FIGURE_TEXT_PT, face = "plain"),
    axis.title      = element_text(size = FIGURE_AXIS_PT, face = "plain"),
    panel.grid      = element_blank(),
    legend.title    = element_text(size = FIGURE_TEXT_PT, face = "plain"),
    legend.text     = element_text(size = FIGURE_TEXT_PT),
    legend.position = "right"
  )

p_conf_AIC
save_figure("figures/Fig2_confusion_matrix_AIC.png", p_conf_AIC,
       width = FIGURE_WIDTH_IN, height = FIGURE_HEIGHT_IN, dpi = FIGURE_DPI)

# BIC model-recovery matrix

p_conf_BIC <- ggplot(model_recovery_BIC,
                     aes(x = true_model,
                         y = selected_model,
                         fill = prop)) +
  geom_tile(color = "grey40", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.2f", prop),
                colour = ifelse(prop >= .65, "white", "black")),
            size = 4.0, fontface = "plain") +
  scale_color_identity(guide = "none") +
  scale_fill_gradientn(
    colours = c("white", FIGURE_TEAL),
    limits  = c(0, 1),
    breaks  = seq(0, 1, by = 0.25),
    name    = "p(Rec|True)",
    guide   = guide_colourbar(display = "rectangles")
  ) +
  labs(x = "True Model", y = "Recovered Model") +
  coord_equal(ratio = 0.85) +
  theme_bw(base_size = FIGURE_TEXT_PT, base_family = FIGURE_FONT) +
  theme(
    plot.title      = element_blank(),
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 13, face = "plain"),
    axis.text.y     = element_text(size = FIGURE_TEXT_PT, face = "plain"),
    axis.title      = element_text(size = FIGURE_AXIS_PT, face = "plain"),
    panel.grid      = element_blank(),
    legend.title    = element_text(size = FIGURE_TEXT_PT, face = "plain"),
    legend.text     = element_text(size = FIGURE_TEXT_PT),
    legend.position = "right"
  )

p_conf_BIC
save_figure("figures/Fig2_confusion_matrix_BIC.png", p_conf_BIC,
       width = FIGURE_WIDTH_IN, height = FIGURE_HEIGHT_IN, dpi = FIGURE_DPI)
