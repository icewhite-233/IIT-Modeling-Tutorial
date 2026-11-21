# Candidate-model teaching figure: one illustrative parameter set per model.
# Use the shared prediction functions; no fitting, noise or clipping.
.script_arg <- grep("^--file=", commandArgs(), value = TRUE)
.script_file <- if (length(.script_arg)) sub("^--file=", "", .script_arg[1]) else
  tail(unlist(lapply(sys.frames(), function(x) x$ofile)), 1)
source(file.path(dirname(normalizePath(.script_file)), "..", "Project_Config.R"), encoding = "UTF-8")
suppressPackageStartupMessages({library(ggplot2); library(cowplot)})
cowplot::set_null_device("cairo")
models <- new.env(parent = globalenv())
sys.source(file.path(PROJECT_ROOT, "Scripts", "Study1", "simulation_helpers.R"), envir = models)
study_style <- list(style_figure = style_figure, FIGURE_TEXT_PT = FIGURE_TEXT_PT,
                    FIGURE_AXIS_PT = FIGURE_AXIS_PT)
font <- FIGURE_FONT
modality_colours <- c(Visual = FIGURE_ORANGE, Auditory = FIGURE_GREEN)
response_colours <- c(Low = FIGURE_TEAL, Midpoint = FIGURE_ORANGE, High = FIGURE_GREEN)
CURVE_LINEWIDTH <- 1.4
PANEL_ASPECT_RATIO <- 1
MAIN_WIDTH_IN <- 14.5
MAIN_HEIGHT_IN <- 18
weight_limits <- c(0, 6)
response_limits <- c(1, 9)
model_order <- c("3p", "4p", "4plq", "5p", "5plq", "7p")

# These values illustrate model hypotheses; they are not empirical estimates.
parameters <- data.frame(
  model = model_order, S0 = 5,
  aV = c(2.4, 2.4, 2 - .25*16/3, 2.4, 2.4 - .18*16/3, 2.4 - .2*16/3),
  aA = c(1.6, 1.6, 2 - .25*16/3, 1.6, 1.6 - .18*16/3, 2),
  b1V = c(0, -.3, 0, -.3, -.15, -.25),
  b1A = c(0, -.3, 0, -.05, -.15, .1),
  b2V = c(0, 0, .25, 0, .18, .2),
  b2A = c(0, 0, .25, 0, .18, -.075))
as_helper_parameters <- function(row) {
  with(row, switch(model,
    `3p` = list(S0 = S0, W_V = aV, W_A = aA),
    `4p` = list(S0 = S0, A_V = aV, A_A = aA, B_VA = b1V),
    `4plq` = list(S0 = S0, A_VA = aV, B_VA = b1V, B_VA2 = b2V),
    `5p` = list(S0 = S0, A_V = aV, B_V = b1V, A_A = aA, B_A = b1A),
    `5plq` = list(S0 = S0, A_V = aV, A_A = aA, B_VA = b1V, B_VA2 = b2V),
    `7p` = list(S0 = S0, A_V = aV, B_V = b1V, B_V2 = b2V,
                A_A = aA, B_A = b1A, B_A2 = b2A)))
}

# This polynomial representation is for weight plots and independent checks.
# All plotted response predictions come from Study 1's actual generator.
get_weights <- function(row, R_V, R_A) {
  x_v <- R_V - 5
  x_a <- R_A - 5
  list(Visual = row$aV + row$b1V*x_v + row$b2V*x_v^2,
       Auditory = row$aA + row$b1A*x_a + row$b2A*x_a^2)
}


# Compute continuous weights and nine joint predictions for each specification.
rating <- seq(1, 9, length.out = 401)
design <- expand.grid(R_V = c(1, 5, 9), R_A = c(3, 5, 7))
weight_list <- response_list <- list()
for (i in seq_len(nrow(parameters))) {
  row <- as.list(parameters[i, ])
  w <- get_weights(row, rating, rating)
  weight_list[[i]] <- data.frame(set = "Contrast", model = row$model,
    rating = rep(rating, 2), centered = rep(rating - 5, 2),
    modality = rep(c("Visual", "Auditory"), each = length(rating)),
    weight = c(w$Visual, w$Auditory))
  d <- design
  d$auditory <- factor(c("Low", "Midpoint", "High")[match(d$R_A, c(3, 5, 7))],
                       levels = c("Low", "Midpoint", "High"))
  d$prediction <- models$generator_dict[[row$model]](
    as_helper_parameters(row), d$R_V, d$R_A)
  d$model <- row$model
  d$set <- "Contrast"
  response_list[[i]] <- d
}
weights <- do.call(rbind, weight_list)
weights$modality <- factor(weights$modality, levels = c("Visual", "Auditory"))
responses <- do.call(rbind, response_list)
theme_small <- function() {
  # Extract only the shared visual theme, so the teaching plot keeps its inputs
  # and semantics. One legend per column avoids repeating it in all six rows.
  reference <- study_style$style_figure(ggplot() + theme_bw(
    base_size = study_style$FIGURE_TEXT_PT, base_family = font))
  reference$theme + theme(
    legend.position = "none",
    aspect.ratio = PANEL_ASPECT_RATIO,
    axis.text = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(),
    axis.ticks = element_blank())
}
plot_weights <- function(set_name, model_name, centered = FALSE, bottom = FALSE, local_limits = weight_limits) {
  d <- subset(weights, set == set_name & model == model_name)
  d$x <- if (centered) d$centered else d$rating
  p <- ggplot(d, aes(x, weight, colour = modality)) +
    geom_line(linewidth = CURVE_LINEWIDTH) +
    scale_colour_manual(values = modality_colours) +
    scale_x_continuous(limits = if (centered) c(-4, 4) else c(1, 9),
                       breaks = if (centered) c(-4, 0, 4) else c(1, 5, 9),
                       expand = expansion(mult = c(.08, .04))) +
    scale_y_continuous(limits = local_limits, breaks = seq(0, local_limits[2], by = 2),
                       expand = expansion(mult = .015)) +
    labs(x = if (bottom) if (centered) "Centered Component\nRating" else "Component Rating" else NULL,
         y = "Weight") + theme_small()
  if (model_name == "4plq") p <- p + annotate(
    "text", x = if (centered) 0 else 5, y = local_limits[2] * .92,
    label = "Shared function", size = 24 / ggplot2::.pt,
    family = font, colour = "#333333")
  p
}
plot_response <- function(set_name, model_name, bottom = FALSE, local_limits = response_limits) {
  d <- subset(responses, set == set_name & model == model_name)
  # A 3 x 3 factorial display; segments are guides between computed points,
  # not interpolated claims about the continuous prediction function.
  d <- subset(d, R_V %in% c(1, 5, 9))
  ggplot(d, aes(R_V, prediction, colour = auditory, group = auditory)) +
    geom_line(linewidth = CURVE_LINEWIDTH) + geom_point(size = 3, shape = 16) +
    scale_colour_manual(values = response_colours,
                        labels = c(Low = "Negative", Midpoint = "Neutral", High = "Positive")) +
    scale_x_continuous(limits = c(1, 9), breaks = c(1, 5, 9),
                       expand = expansion(mult = c(.08, .04))) +
    scale_y_continuous(limits = local_limits, breaks = seq(1, 9, 2), expand = expansion(mult = .02)) +
    labs(x = if (bottom) "Visual Component\nRating" else NULL,
         y = "Predicted Joint\nResponse") + theme_small()
}

parameter_counts <- c(`3p` = 3, `4p` = 4, `4plq` = 4, `5p` = 5, `5plq` = 5, `7p` = 7)
descriptions <- c(`3p` = "Constant weights",
                  `4p` = "Separate baselines\nShared linear slope",
                  `4plq` = "Shared quadratic\nweighting function",
                  `5p` = "Separate baselines\nand linear slopes",
                  `5plq` = "Separate baselines\nShared slope and\ncurvature",
                  `7p` = "Separate baselines,\nslopes, and curvatures")
plot_model_hypothesis <- function(model_name) {
  ggdraw() +
    draw_label(model_name, x = .035, y = .81, hjust = 0, size = 31, fontface = "bold", fontfamily = font) +
    draw_label(paste(parameter_counts[[model_name]], "free parameters"),
               x = .035, y = .62, hjust = 0, size = study_style$FIGURE_TEXT_PT,
               colour = "#555555", fontfamily = font) +
    draw_label(descriptions[[model_name]],
               x = .035, y = .30, hjust = 0, size = study_style$FIGURE_TEXT_PT,
               colour = "#333333", fontfamily = font)
}
header <- function(label, items = NULL, colours = NULL, shapes = NULL, legend_title = NULL) {
  p <- ggdraw() + draw_label(label, x = .03, y = .85, hjust = 0, size = 26,
                            fontface = "bold", fontfamily = font)
  if (!is.null(legend_title)) p <- p + draw_label(
    legend_title, x = .03, y = .52, hjust = 0, size = study_style$FIGURE_TEXT_PT,
    colour = "#333333", fontfamily = font)
  if (!is.null(items)) {
    positions <- if (length(items) == 2) c(.07, .51) else c(.02, .36, .69)
    key_length <- if (length(items) == 2) .09 else .055
    for (i in seq_along(items)) {
      p <- p + draw_line(c(positions[i], positions[i] + key_length), c(.16, .16),
                         colour = colours[i], linewidth = CURVE_LINEWIDTH) +
        draw_label(items[i], x = positions[i] + key_length + .02, y = .16, hjust = 0,
                   size = study_style$FIGURE_TEXT_PT, colour = "#333333", fontfamily = font)
      if (!is.null(shapes)) p <- p + geom_point(
        data = data.frame(x = positions[i] + key_length/2, y = .16), aes(x = x, y = y),
        inherit.aes = FALSE, shape = shapes[i], colour = colours[i], size = 3)
    }
  }
  p
}
make_full_figure <- function(set_name, centered = TRUE, y_range = response_limits,
                             row_separators = TRUE) {
  # Common row and column dimensions allow direct comparison.
  col_widths <- c(.72, 1, 1.28)
  headings <- plot_grid(
    header("Model Hypothesis"),
    header("Weight Functions", c("Visual", "Auditory"), modality_colours),
    header("Illustrative Joint Predictions", c("Negative", "Neutral", "Positive"),
           response_colours, shapes = rep(16, 3), legend_title = "Auditory Component Rating"),
    nrow = 1, rel_widths = col_widths)
  rows <- lapply(seq_along(model_order), function(i) {
    m <- model_order[i]
    # Align the two charts, not the hypothesis text. Shared x labels below the
    # grid keep the sixth row's plotting area the same size as the other rows.
    charts <- plot_grid(plot_weights(set_name, m, centered, FALSE),
                        plot_response(set_name, m, FALSE, y_range), nrow = 1,
                        rel_widths = col_widths[2:3], align = "h", axis = "tb")
    plot_grid(plot_model_hypothesis(m), charts, nrow = 1,
              rel_widths = c(col_widths[1], sum(col_widths[2:3])))
  })
  valence_axis_title <- function(label) {
    ggdraw() +
      draw_label(label, x = .57, y = .73,
                 size = study_style$FIGURE_AXIS_PT, fontfamily = font) +
      draw_label("Negative \u2192 Positive", x = .57, y = .27,
                 size = study_style$FIGURE_TEXT_PT, colour = "#333333", fontfamily = font)
  }
  footer <- plot_grid(ggdraw(),
    valence_axis_title("Component Valence"),
    valence_axis_title("Visual Component Valence"),
    nrow = 1, rel_widths = col_widths)
  # As in the study figures, the figure title and setup belong in the caption.
  row_heights <- c(.9, rep(1.85, 6), .6)
  figure <- plot_grid(plotlist = c(list(headings), rows, list(footer)), ncol = 1,
                      rel_heights = row_heights)
  if (row_separators) {
    # Group each model across the three columns without adding data gridlines.
    boundary_y <- 1 - cumsum(row_heights)[2:6]/sum(row_heights)
    figure <- ggdraw(figure)
    for (y in boundary_y) figure <- figure + draw_line(
      x = c(.01, .99), y = c(y, y), colour = "#D2D2D2", linewidth = .35)
  }
  figure
}

# The plotting functions share limits within columns and omit numeric ticks.
figure <- make_full_figure("Contrast", y_range = c(1, 9))
export_figure(file.path(PROJECT_ROOT, "Figures", "Model_Space.png"), figure,
              width = MAIN_WIDTH_IN, height = MAIN_HEIGHT_IN, dpi = FIGURE_DPI)
result_dir <- file.path(PROJECT_ROOT, "Scripts", "Study1", "results")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(parameters, file.path(result_dir, "ModelSpace_Parameters.csv"), row.names = FALSE)
write.csv(responses[c("model", "R_V", "R_A", "auditory", "prediction")],
          file.path(result_dir, "ModelSpace_Predictions.csv"), row.names = FALSE)
