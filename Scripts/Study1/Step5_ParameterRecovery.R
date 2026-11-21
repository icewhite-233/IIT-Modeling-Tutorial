# Study 1: parameter recovery for the retained 5plq specification.
# Sample reference-centered parameters, simulate, refit and summarize accuracy.

## 0. Environment and dependencies -----------------------------------------
rm(list = ls())
# Load shared paths relative to this script.
.script_files <- unlist(lapply(sys.frames(), function(frame) frame$ofile))
.script_file <- if (length(.script_files)) tail(.script_files, 1) else
  sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
if (!length(.script_file)) .script_file <- if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) rstudioapi::getSourceEditorContext()$path else "Step5_ParameterRecovery.R"
source(file.path(dirname(normalizePath(.script_file)), "..", "Project_Config.R"), encoding = "UTF-8")
setwd(file.path(PROJECT_ROOT, "Scripts", "Study1"))
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# Packages
packages <- c("readxl", "dplyr", "tibble", "ggplot2", "tidyr", "purrr", "minpack.lm")
to_install <- setdiff(packages, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)

library(readxl)
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyr)
library(purrr)
library(minpack.lm)

# Load helper functions and empirical reference parameters
# (wavg_5plq, simulate_once(), reference_params, etc.)
source("simulation_helpers.R")

RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(20251120)

## 1. Load empirical data and define design --------------------------------

dat <- read_excel(file.path(PROJECT_ROOT, "Data", "Study1", "MODELING_Exp1_Valence.xlsx"))

# Design for simulation: keep R_V, R_A (and VV/AV if you want to carry them through)
design_df <- dat %>%
  dplyr::select(R_V, R_A, VV, AV)

str(design_df)

## 2. Define reference-centered parameter ranges ---------------------------
## Center the sampling ranges on the 5plq reference values.

## For each parameter θ*, we define:
## θ_min = θ* - range_factor * max(|θ*|, abs_min)
## θ_max = θ* + range_factor * max(|θ*|, abs_min)

reference_5plq <- reference_params[["5plq"]]  # list(S0, A_V, A_A, B_VA, B_VA2)

range_factor <- 0.2  # 20% of max(abs(theta), abs_min) on each side
abs_min    <- 0.05  # magnitude floor before multiplying by range_factor

param_ranges_5plq <- lapply(reference_5plq, function(theta) {
  width <- range_factor * max(abs(theta), abs_min)
  c(theta - width, theta + width)
})

# Inspect the parameter ranges.
cat("Parameter ranges for 5plq (derived from empirical best-fit):\n")
print(param_ranges_5plq)

## 3. Configuration for parameter recovery ---------------------------------

nsim_param_rec <- 200   # number of simulated pair-level datasets
noise_sd        <- 0.39  # Main-text SD; use 0.20 for the supplementary recovery analysis.
jitter_scale    <- 0.2 # additive jitter around generating parameters

## 4. Main loop: simulate with random parameters, then refit 5plq ----------

param_rec_list <- vector("list", nsim_param_rec)
recovery_attempts <- vector("list", nsim_param_rec)

for (i in seq_len(nsim_param_rec)) {

  ## 4.1 Sample generative parameters uniformly within the empirically-derived ranges ---
  S0_true   <- runif(1, param_ranges_5plq$S0[1],   param_ranges_5plq$S0[2])
  A_V_true   <- runif(1, param_ranges_5plq$A_V[1],   param_ranges_5plq$A_V[2])
  A_A_true   <- runif(1, param_ranges_5plq$A_A[1],   param_ranges_5plq$A_A[2])
  B_VA_true  <- runif(1, param_ranges_5plq$B_VA[1],  param_ranges_5plq$B_VA[2])
  B_VA2_true <- runif(1, param_ranges_5plq$B_VA2[1], param_ranges_5plq$B_VA2[2])

  truth_par <- list(
    S0   = S0_true,
    A_V   = A_V_true,
    A_A   = A_A_true,
    B_VA  = B_VA_true,
    B_VA2 = B_VA2_true
  )

  ## 4.2 Simulate one dataset under 5plq with these generative parameters ----
  sim_df <- simulate_once(
    true_model    = "5plq",
    design_df     = design_df,
    truth_params  = truth_par,
    noise_sd      = noise_sd,
    clip_to_scale = TRUE,
    lower         = SCALE_MIN,
    upper         = SCALE_MAX
  )

  ## 4.3 Refit the same model (5plq) to the simulated data -------------------
  ## Starts are jittered generating parameters for this recovery simulation.
  start_list <- jitter_list(truth_par, jitter_scale = jitter_scale, abs_min = abs_min)

  fit_5plq <- try(
    nlsLM(
      R_VA_sim ~ wavg_5plq(S0, A_V, A_A, B_VA, B_VA2, R_V, R_A),
      data    = sim_df,
      start   = start_list,
      control = nls.lm.control(maxiter = 2000)
    ),
    silent = TRUE
  )

  converged <- !inherits(fit_5plq, "try-error") && isTRUE(fit_5plq$convInfo$isConv)
  recovery_attempts[[i]] <- data.frame(sim_rep = i, converged = converged)
  if (!converged) {
    warning(sprintf("Parameter-recovery fit %d did not converge; retained in the convergence denominator.", i))
    next
  }

  hat_par <- as.list(coef(fit_5plq))

  ## 4.4 Store generative and recovered parameter values ---------------------
  param_rec_list[[i]] <- tibble(
    sim_rep   = i,
    S0_true   = S0_true,
    A_V_true   = A_V_true,
    A_A_true   = A_A_true,
    B_VA_true  = B_VA_true,
    B_VA2_true = B_VA2_true,
    S0_hat    = hat_par$S0,
    A_V_hat    = hat_par$A_V,
    A_A_hat    = hat_par$A_A,
    B_VA_hat   = hat_par$B_VA,
    B_VA2_hat  = hat_par$B_VA2
  )
}

attempts <- bind_rows(recovery_attempts)
write.csv(attempts, "results/Study1_ParameterRecovery_Attempts.csv", row.names = FALSE)

# Combine all valid rows
param_rec_df <- dplyr::bind_rows(param_rec_list)
write.csv(param_rec_df, "results/Study1_ParameterRecovery_Pairs.csv", row.names = FALSE)

cat("\nNumber of successful fits in parameter recovery:", nrow(param_rec_df), "out of", nsim_param_rec, "\n")

## 5. Correlation analysis: generative vs recovered ------------------------

true_mat <- as.matrix(param_rec_df[, c("S0_true", "A_V_true", "A_A_true", "B_VA_true", "B_VA2_true")])
hat_mat  <- as.matrix(param_rec_df[, c("S0_hat",  "A_V_hat",  "A_A_hat",  "B_VA_hat",  "B_VA2_hat")])

# 5x5 matrix: rows = generative parameters, columns = recovered parameters
cor_mat <- cor(true_mat, hat_mat, use = "pairwise.complete.obs")

cat("\nCross-correlation matrix between generative (rows) and recovered (columns) parameters:\n")
print(round(cor_mat, 3))

# Diagonal = same-parameter recovery (S0, A_V, A_A, B_VA, B_VA2)
diag_cor    <- diag(cor_mat)
param_names <- c("S0", "A_V", "A_A", "B_VA", "B_VA2")

cat("\nParameter-wise recovery (Pearson r for true vs recovered):\n")
for (j in seq_along(param_names)) {
  cat(sprintf("  %s: r = %.3f\n", param_names[j], diag_cor[j]))
}

# Summarize accuracy from the same 200 generating/recovered parameter pairs.
recovery_summary <- do.call(rbind, lapply(param_names, function(parameter) {
  generating <- param_rec_df[[paste0(parameter, "_true")]]
  recovered <- param_rec_df[[paste0(parameter, "_hat")]]
  paired <- complete.cases(generating, recovered)
  generating <- generating[paired]
  recovered <- recovered[paired]
  correlation <- cor.test(generating, recovered, method = "pearson")
  data.frame(
    parameter = parameter,
    n = length(generating),
    r = unname(correlation$estimate),
    CI_low_95 = correlation$conf.int[1],
    CI_high_95 = correlation$conf.int[2],
    mean_generating = mean(generating),
    mean_recovered = mean(recovered),
    bias = mean(recovered - generating),
    RMSE = sqrt(mean((recovered - generating)^2))
  )
}))
print(recovery_summary)
write.csv(recovery_summary, "results/Study1_5plq_ParameterRecovery_Summary.csv",
          row.names = FALSE)

# Parameter recovery heatmap (plotmath version)

library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)

cor_mat

cor_long <- as.data.frame(cor_mat) %>%
  mutate(TrueParam = rownames(.)) %>%
  pivot_longer(
    cols      = -TrueParam,
    names_to  = "RecoveredParam",
    values_to = "Correlation"
  ) %>%
  mutate(
    TrueParam      = gsub("_true", "", TrueParam),
    RecoveredParam = gsub("_hat",  "", RecoveredParam)
  )

param_levels <- c("S0", "A_V", "A_A", "B_VA", "B_VA2")

cor_long <- cor_long %>%
  mutate(
    TrueParam      = factor(TrueParam,      levels = param_levels),
    RecoveredParam = factor(RecoveredParam, levels = param_levels)
  )

head(cor_long)

pretty_labs <- c(
  S0   = expression(S[0]),
  A_V   = expression(a[V]),
  A_A   = expression(a[A]),
  B_VA  = expression(b[1]),
  B_VA2 = expression(b[2])
)

# Preserve signed correlations in the figure source table.
write.csv(cor_long, "results/Plot_ParameterRecovery_Correlations.csv", row.names = FALSE)

p_param_cor <- ggplot(cor_long,
                      aes(x = TrueParam,
                          y = RecoveredParam,
                          fill = Correlation)) +
  geom_tile(color = "grey70", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", Correlation),
                colour = ifelse(Correlation >= .65, "white", "black")),
            size = 4, fontface = "plain") +
  scale_color_identity(guide = "none") +
  scale_fill_gradient2(
    low       = FIGURE_ORANGE,
    mid       = "white",
    high      = FIGURE_GREEN,
    midpoint  = 0,
    limits    = c(-1, 1),
    name      = "Correlation",
    guide     = guide_colourbar(display = "rectangles")
  ) +
  scale_x_discrete(labels = pretty_labs) +
  scale_y_discrete(labels = pretty_labs) +
  labs(
    x = "True Parameter",
    y = "Recovered Parameter"
  ) +
  coord_equal() +
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

p_param_cor

save_figure("figures/Fig2_param_recovery.png", p_param_cor,
       width = FIGURE_WIDTH_IN, height = FIGURE_HEIGHT_IN, dpi = FIGURE_DPI)

## 6. Parameter recovery scatterplots (true vs recovered) ------------------
## These are the typical "parameter recovery" plots for each free parameter.

plot_df <- purrr::map_dfr(
  param_names,
  function(p) {
    tibble(
      param = p,
      true  = param_rec_df[[paste0(p, "_true")]],
      hat   = param_rec_df[[paste0(p, "_hat")]]
    )
  }
)

# Per-parameter Pearson r for annotation (optional)
cor_annot <- plot_df %>%
  dplyr::group_by(param) %>%
  dplyr::summarise(
    r = cor(true, hat, use = "pairwise.complete.obs"),
    .groups = "drop"
  )

cat("\nPer-parameter correlations from plot_df (for checking):\n")
print(cor_annot)

# Scatterplots
p_param_rec <- ggplot(plot_df, aes(x = true, y = hat)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  facet_wrap(~ param, scales = "free") +
  labs(
    title = "Parameter recovery for 5plq (true vs recovered)",
    x     = "Generative parameter value",
    y     = "Recovered parameter estimate"
  ) +
  theme_bw(base_size = FIGURE_TEXT_PT, base_family = FIGURE_FONT)

print(p_param_rec)

cat("\nParameter recovery simulation for 5plq completed.\n")
