# Study 1: fitted 5plq weights and observed/predicted ratings.

rm(list = ls())
# Load shared paths relative to this script.
.script_files <- unlist(lapply(sys.frames(), function(frame) frame$ofile))
.script_file <- if (length(.script_files)) tail(.script_files, 1) else
  sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
if (!length(.script_file)) .script_file <- if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) rstudioapi::getSourceEditorContext()$path else "Step7_PlotWeightFunctions.R"
source(file.path(dirname(normalizePath(.script_file)), "..", "Project_Config.R"), encoding = "UTF-8")
setwd(file.path(PROJECT_ROOT, "Scripts", "Study1"))
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

library(readxl)
library(dplyr)
library(ggplot2)
library(tibble)
source("simulation_helpers.R")

## 1. Load data and fitted model -------------------------------------------

dat <- read_excel(file.path(PROJECT_ROOT, "Data", "Study1", "MODELING_Exp1_Valence.xlsx")) %>%
  mutate(
    Stimcomb = row_number(),
    VV = factor(VV, levels = c(-1,0,1), labels=c("Neg","Neu","Pos")),
    AV = factor(AV, levels = c(-1,0,1), labels=c("Neg","Neu","Pos"))
  )

fit_objects <- readRDS("results/model_fits_weighted_averaging.rds")
nls_5plq    <- fit_objects$nls_5plq
coef_5plq   <- as.list(coef(nls_5plq))

## 2. Define W_V / W_A for 5plq (raw, unnormalized) ------------------------

W_V_5plq <- function(x, A_V, B_VA, B_VA2) {
  x_c <- x - SCALE_MID
  A_V + B_VA * x_c + B_VA2 * x_c^2
}

W_A_5plq <- function(x, A_A, B_VA, B_VA2) {
  x_c <- x - SCALE_MID
  A_A + B_VA * x_c + B_VA2 * x_c^2
}

## 3. Weighting functions (global) -----------------------------------------

newx <- tibble(x = seq(SCALE_MIN, SCALE_MAX,
                       length.out = 161))

# Solid segments span the observed component-input range for each modality.
# Dashed segments show the same fitted function outside that range.
make_weight_segments <- function(modality, weight_function, observed_range) {
  bounds <- c(SCALE_MIN, observed_range, SCALE_MAX)
  bind_rows(lapply(seq_len(3), function(segment) {
    lower <- bounds[segment]
    upper <- bounds[segment + 1]
    if (upper <= lower) return(NULL)
    x <- sort(unique(c(lower, newx$x[newx$x > lower & newx$x < upper], upper)))
    tibble(x = x, weight = weight_function(x), modality = modality,
           input_range = if (segment == 2) "Observed range" else "Extrapolation",
           segment = segment)
  }))
}

weight_segments <- bind_rows(
  make_weight_segments(
    "Visual",
    function(x) W_V_5plq(x, coef_5plq$A_V, coef_5plq$B_VA, coef_5plq$B_VA2),
    range(dat$R_V, na.rm = TRUE)
  ),
  make_weight_segments(
    "Auditory",
    function(x) W_A_5plq(x, coef_5plq$A_A, coef_5plq$B_VA, coef_5plq$B_VA2),
    range(dat$R_A, na.rm = TRUE)
  )
)
weight_segments$input_range <- factor(
  weight_segments$input_range, levels = c("Observed range", "Extrapolation")
)

plot_weightfunc <- ggplot(
  weight_segments,
  aes(x = x, y = weight, color = modality, linetype = input_range,
      group = interaction(modality, segment))
) +
  geom_line(linewidth = 1.2) +
  labs(x = "Unimodal Valence Rating", y = "Weight") +
  scale_color_manual(name = NULL,
                     values = c(Visual = FIGURE_ORANGE, Auditory = FIGURE_GREEN)) +
  scale_linetype_manual(name = NULL,
                        values = c("Observed range" = "solid",
                                   Extrapolation = "dashed")) +
  guides(color = guide_legend(order = 1, override.aes = list(linetype = "solid")),
         linetype = "none") +
  coord_cartesian(xlim = c(SCALE_MIN, SCALE_MAX)) +
  scale_x_continuous(
    breaks = SCALE_BREAKS,
    expand = c(0, 0)
  ) +
  theme_bw(base_size = FIGURE_TEXT_PT, base_family = FIGURE_FONT) +
  theme(
    axis.title.x    = element_text(size = FIGURE_AXIS_PT, face = "plain",
                                   margin = margin(t = 10)),
    axis.title.y    = element_text(size = FIGURE_AXIS_PT, face = "plain",
                                   margin = margin(r = 10)),
    axis.text.x     = element_text(size = FIGURE_TEXT_PT, face = "plain"),
    axis.text.y     = element_text(size = FIGURE_TEXT_PT),
    legend.title    = element_text(size = FIGURE_TEXT_PT, face = "plain"),
    legend.text     = element_text(size = FIGURE_TEXT_PT),
    legend.position = "right",
    panel.grid      = element_blank()
  )

# Display
plot_weightfunc

# Save high-resolution figure
save_figure("figures/Fig3_weighting_functions_5plq.png",
       plot_weightfunc, width = FIGURE_WIDTH_IN, height = FIGURE_HEIGHT_IN, dpi = FIGURE_DPI)

## 4. Pair-level fitted weights and predictions ----------------------------

dat_latent <- dat %>%
  mutate(
    W_V = W_V_5plq(R_V, coef_5plq$A_V, coef_5plq$B_VA, coef_5plq$B_VA2),
    W_A = W_A_5plq(R_A, coef_5plq$A_A, coef_5plq$B_VA, coef_5plq$B_VA2),
    R_VA_hat = fitted(nls_5plq),
    diff_W = W_V - W_A  # Raw weight difference at each pair's component inputs.
  )

## 5. Observed vs predicted R_VA -------------------------------------------

plot_fit <- ggplot(dat_latent, aes(x = R_VA, y = R_VA_hat)) +

  geom_point(color = FIGURE_TEAL, alpha = 0.55, size = 3) +

  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", linewidth = 1, color = "grey40") +

  geom_smooth(method = "lm", se = FALSE,
              color = "black", linewidth = 1.2) +

  scale_x_continuous(limits = c(SCALE_MIN, SCALE_MAX),
                     breaks = SCALE_BREAKS) +
  scale_y_continuous(limits = c(SCALE_MIN, SCALE_MAX),
                     breaks = SCALE_BREAKS) +

  labs(
    x = "Observed rating",
    y = "Predicted rating"
  ) +

  theme_bw(base_size = FIGURE_TEXT_PT, base_family = FIGURE_FONT) +
  theme(
    plot.title      = element_blank(),
    axis.title.x    = element_text(size = FIGURE_AXIS_PT, face = "plain",
                                   margin = margin(t = 10)),
    axis.title.y    = element_text(size = FIGURE_AXIS_PT, face = "plain",
                                   margin = margin(r = 10)),
    axis.text.x     = element_text(size = FIGURE_TEXT_PT, face = "plain"),
    axis.text.y     = element_text(size = FIGURE_TEXT_PT, face = "plain"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )

plot_fit

save_figure("figures/Fig3_fit_observed_vs_predicted.png",
       plot_fit, width = FIGURE_WIDTH_IN, height = FIGURE_HEIGHT_IN, dpi = FIGURE_DPI)

cat("\nAnalyze winning model completed.\n")
