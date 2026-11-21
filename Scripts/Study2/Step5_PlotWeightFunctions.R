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
setwd(file.path(PROJECT_ROOT, "Scripts", "Study2"))
library(dplyr)
library(ggplot2)
library(tibble)
source("simulation_helpers.R")
# Tick marks for the 0-1 component-rating scale.
SCALE_BREAKS <- seq(0, 1, .25)
dir.create("figures", showWarnings = FALSE)

## 1. Load fitted model ----------------------------------------------------
fit_objects <- readRDS("results/model_fits_weighted_averaging.rds")
# Display the focal specification from Step 2.
nls_3p <- fit_objects$nls_3p
coef_3p <- as.list(coef(nls_3p))
dat <- read.csv(file.path(PROJECT_ROOT, "Data", "Study2", "Input_Data54.csv"))

## 2. Define W_V / W_A (raw, unnormalized) ---------------------------------
# The 3p weights are constant, so both curves are horizontal.
W_V_3p <- function(x, W_V) {
  rep(W_V, length(x))
}
W_A_3p <- function(x, W_A) {
  rep(W_A, length(x))
}
## 3. Weighting functions (global) -----------------------------------------
newx <- tibble(x = seq(SCALE_MIN, SCALE_MAX, length.out = 161))
# Solid lines cover observed component inputs; dashed lines show extrapolation.
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
  make_weight_segments("Visual", function(x) W_V_3p(x, coef_3p$W_V),
                       range(dat$R_V, na.rm = TRUE)),
  make_weight_segments("Auditory", function(x) W_A_3p(x, coef_3p$W_A),
                       range(dat$R_A, na.rm = TRUE))
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
  labs(x = "Unimodal Disgust Rating",
       y = "Weight") +
  scale_color_manual(name = NULL,
                     values = c(Visual = FIGURE_ORANGE, Auditory = FIGURE_GREEN)) +
  scale_linetype_manual(name = NULL,
                        values = c("Observed range" = "solid",
                                   Extrapolation = "dashed")) +
  guides(color = guide_legend(order = 1, override.aes = list(linetype = "solid")),
         linetype = "none") +
  coord_cartesian(xlim = c(SCALE_MIN, SCALE_MAX),
                  ylim = c(0, 1.1 * max(coef_3p$W_V, coef_3p$W_A))) +
  scale_x_continuous(breaks = SCALE_BREAKS, expand = c(0, 0)) +
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
plot_weightfunc
save_figure("figures/Fig3_weighting_functions_3p.png",
       plot_weightfunc, width = FIGURE_WIDTH_IN, height = FIGURE_HEIGHT_IN, dpi = FIGURE_DPI)

# Save the values used to draw the fitted weight curves.
newx$W_V <- W_V_3p(newx$x, coef_3p$W_V)
newx$W_A <- W_A_3p(newx$x, coef_3p$W_A)
write.csv(newx, "results/Plot_WeightFunctions_3p.csv", row.names = FALSE)
