# Plot visual and auditory weighting functions, as in Study 1.
# Read the existing fit; weights are raw, not normalized probabilities.

rm(list = ls())
# Load shared paths relative to this script.
.script_files <- unlist(lapply(sys.frames(), function(frame) frame$ofile))
.script_file <- if (length(.script_files)) tail(.script_files, 1) else
  sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
if (!length(.script_file)) .script_file <- if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) rstudioapi::getSourceEditorContext()$path else "Step5_PlotWeightFunctions.R"
source(file.path(dirname(normalizePath(.script_file)), "..", "Project_Config.R"), encoding = "UTF-8")
setwd(file.path(PROJECT_ROOT, "Scripts", "Study3"))
library(dplyr)
library(ggplot2)
library(tibble)
source("help_function_choice.R")
dir.create("figures", showWarnings = FALSE)

## 1. Load fitted model ----------------------------------------------------
fit_objects <- readRDS("results/model_fits_weighted_averaging_BinomialLikelihood_54.rds")
# Display the focal specification from Step 2.
fit_7p <- fit_objects$fit_7p
coef_7p <- as.list(coef(fit_7p))
dat <- choice_prepare_inputs(read.csv(file.path(PROJECT_ROOT, "Data", "Study3", "Input_Data54.csv")))

# Weighting functions receive corrected component log odds.
# Use the corrected log odds for proportions 0 and 1 with n = 50 as endpoints.
n_display <- max(c(dat$n_V, dat$n_A))
SCALE_MIN <- choice_count_to_logodds(0, n_display)
SCALE_MAX <- choice_count_to_logodds(n_display, n_display)
SCALE_MID <- 0
SCALE_BREAKS <- c(-4, -2, 0, 2, 4)

## 2. Define W_V / W_A (raw, unnormalized) ---------------------------------
W_V_7p <- function(x, A_V, B_V, B_V2) {
  x_c <- x - SCALE_MID
  A_V + B_V * x_c + B_V2 * x_c^2
}
W_A_7p <- function(x, A_A, B_A, B_A2) {
  x_c <- x - SCALE_MID
  A_A + B_A * x_c + B_A2 * x_c^2
}
## 3. Weighting functions (global) -----------------------------------------
newx <- tibble(x = seq(SCALE_MIN, SCALE_MAX, length.out = 161))
# Solid lines cover observed corrected component inputs; dashed lines show extrapolation.
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
    function(x) W_V_7p(x, coef_7p$A_V, coef_7p$B_V, coef_7p$B_V2),
    range(dat$L_V, na.rm = TRUE)
  ),
  make_weight_segments(
    "Auditory",
    function(x) W_A_7p(x, coef_7p$A_A, coef_7p$B_A, coef_7p$B_A2),
    range(dat$L_A, na.rm = TRUE)
  )
)
weight_segments$input_range <- factor(
  weight_segments$input_range, levels = c("Observed range", "Extrapolation")
)
weight_segments$probability <- ((n_display + 1) * plogis(weight_segments$x) - .5) / n_display
weight_plot <- ggplot(
  weight_segments,
  aes(y = weight, color = modality, linetype = input_range,
      group = interaction(modality, segment))
) +
  geom_line(linewidth = 1.2) +
  labs(y = "Weight") +
  scale_color_manual(name = NULL,
                     values = c(Visual = FIGURE_ORANGE, Auditory = FIGURE_GREEN)) +
  scale_linetype_manual(name = NULL,
                        values = c("Observed range" = "solid",
                                   Extrapolation = "dashed")) +
  guides(color = guide_legend(order = 1, override.aes = list(linetype = "solid")),
         linetype = "none") +
  theme_bw(base_size = FIGURE_TEXT_PT, base_family = FIGURE_FONT) +
  theme(
    plot.title = element_blank(),
    axis.title.x = element_text(size = FIGURE_AXIS_PT, face = "plain", margin = margin(t = 10)),
    axis.title.y = element_text(size = FIGURE_AXIS_PT, face = "plain", margin = margin(r = 10)),
    axis.text.x = element_text(size = FIGURE_TEXT_PT, face = "plain"),
    axis.text.y = element_text(size = FIGURE_TEXT_PT),
    legend.title = element_text(size = FIGURE_TEXT_PT, face = "plain"),
    legend.text = element_text(size = FIGURE_TEXT_PT),
    legend.position = "right",
    panel.grid = element_blank()
  )
plot_weightfunc <- weight_plot +
  aes(x = x) +
  labs(x = "Unimodal Log Odds") +
  coord_cartesian(xlim = c(SCALE_MIN, SCALE_MAX)) +
  scale_x_continuous(breaks = SCALE_BREAKS, expand = c(0, 0))
plot_weightfunc
save_figure("figures/Fig3_weighting_functions_7p.png",
       plot_weightfunc, width = FIGURE_WIDTH_IN, height = FIGURE_HEIGHT_IN, dpi = FIGURE_DPI)

# Save the values used to draw the fitted weight curves.
newx$W_V <- W_V_7p(newx$x, coef_7p$A_V, coef_7p$B_V, coef_7p$B_V2)
newx$W_A <- W_A_7p(newx$x, coef_7p$A_A, coef_7p$B_A, coef_7p$B_A2)
write.csv(newx, "results/Plot_WeightFunctions_7p.csv", row.names = FALSE)

## 4. Display the same weights against observed component proportions ------
# The weights still use log-odds inputs; only the displayed x axis changes.
# Invert the half-count correction: p = (51 * plogis(L) - 0.5) / 50.
# Weights are not transformed to probabilities.
# One uncorrected-proportion axis is available only with a common denominator.
if (length(unique(c(dat$n_V, dat$n_A))) == 1L) {
probability_grid <- newx %>%
  mutate(probability = ((n_display + 1) * plogis(x) - .5) / n_display)
plot_weightfunc_probability <- weight_plot +
  aes(x = probability) +
  scale_x_continuous(breaks = seq(0, 1, .25), expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = "Unimodal disgust-choice probability")
save_figure("figures/Fig3_weighting_functions_7p_probability.png",
       plot_weightfunc_probability, width = FIGURE_WIDTH_IN, height = FIGURE_HEIGHT_IN, dpi = FIGURE_DPI)
write.csv(probability_grid, "results/Plot_WeightFunctions_7p_probability.csv", row.names = FALSE)

} else {
  message("Unequal component denominators: use the log-odds plot; no unique inverse correction exists for one probability axis.")
}
