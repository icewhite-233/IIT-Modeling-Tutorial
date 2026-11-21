# Portable project root; do not enter a machine-specific path.
PROJECT_ROOT <- ""
if (!nzchar(PROJECT_ROOT)) {
  .config_files <- unlist(lapply(sys.frames(), function(frame) frame$ofile))
  PROJECT_ROOT <- dirname(dirname(normalizePath(tail(.config_files, 1), winslash = "/")))
}
PROJECT_ROOT <- normalizePath(PROJECT_ROOT, winslash = "/", mustWork = TRUE)

# Shared plotting/export settings; these do not change model calculations.
FIGURE_WIDTH_IN <- 7.5
FIGURE_HEIGHT_IN <- 6.5
FIGURE_DPI <- 300
FIGURE_FONT <- "Arial"
FIGURE_TEXT_PT <- 28
FIGURE_AXIS_PT <- 32
FIGURE_CELL_PT <- 26  # matrix annotations in points
# Approved reference palette. These settings affect drawing only.
FIGURE_TEAL <- "#319C9B"
FIGURE_ORANGE <- "#E65E23"
FIGURE_GREEN <- "#56973D"
FIGURE_CONDITION_COLOURS <- c(FIGURE_TEAL, FIGURE_ORANGE, FIGURE_GREEN)
FIGURE_MODALITY_COLOURS <- c(Visual = FIGURE_ORANGE, Auditory = FIGURE_GREEN)
FIGURE_SERIES_COLOURS <- c(Audiovisual = FIGURE_TEAL,
                          Visual = FIGURE_ORANGE, Auditory = FIGURE_GREEN)
if (.Platform$OS.type == "windows") {
  do.call(grDevices::windowsFonts,
          setNames(list(grDevices::windowsFont(FIGURE_FONT)), FIGURE_FONT))
}
# Stand-alone Rscript calls also print plots before exporting them. A Cairo
# preview device supports the figure font without leaving Rplots.pdf in the
# package; it does not affect fitted models or numerical outputs.
if (!interactive() && capabilities("cairo")) {
  options(device = function(...) grDevices::cairo_pdf(
    file = tempfile("iit-plot-preview-", fileext = ".pdf"), ...))
}

reference_recovery_layout <- function(plot) {
  # Match the reference matrices: corresponding entries run from bottom-left
  # to top-right. Axis order changes only the display, never matrix values.
  data_names <- names(plot$data)
  if (all(c("true_model", "selected_model", "prop") %in% data_names)) {
    model_order <- levels(plot$data$true_model)
    if (is.null(model_order)) model_order <- unique(as.character(plot$data$true_model))
    plot <- plot + ggplot2::scale_y_discrete(limits = model_order) +
      ggplot2::labs(x = "True Model", y = "Recovered Model")
    fill_scale <- plot$scales$get_scales("fill")
    fill_scale$name <- "p(Rec|True)"
  }
  if (all(c("TrueParam", "RecoveredParam", "Correlation") %in% data_names)) {
    parameter_order <- levels(plot$data$TrueParam)
    if (is.null(parameter_order)) parameter_order <- unique(as.character(plot$data$TrueParam))
    labels <- plot$scales$get_scales("x")$labels
    plot <- plot + ggplot2::scale_y_discrete(limits = parameter_order, labels = labels) +
      ggplot2::labs(x = "True Parameter", y = "Recovered Parameter")
    fill_scale <- plot$scales$get_scales("fill")
    fill_scale$name <- "Correlation"
  }
  plot
}

reference_comparison_layout <- function(plot) {
  if (!all(c("model", "Criterion", "Delta") %in% names(plot$data))) return(plot)
  plot <- unserialize(serialize(plot, NULL))
  # As in Study 1, put the lowest primary information criterion at the top.
  primary <- plot$data[as.character(plot$data$Criterion) == "AICc", , drop = FALSE]
  ranked_models <- as.character(primary$model[order(primary$Delta)])
  short_model_labels <- function(x) {
    x <- as.character(x)
    x[x == "Visual-only"] <- "V-only"
    x[x == "Auditory-only"] <- "A-only"
    x
  }
  plot <- plot + ggplot2::scale_y_discrete(
    limits = rev(ranked_models), labels = short_model_labels)
  # Remove the former retained-model outline from older cached plots as well.
  plot$layers <- Filter(function(layer) {
    !(inherits(layer$geom, "GeomPoint") && isTRUE(layer$aes_params$shape == 21) &&
        isFALSE(layer$inherit.aes) && identical(layer$aes_params$colour, "black"))
  }, plot$layers)
  plot$mapping$shape <- NULL
  plot$scales$scales <- Filter(function(scale) !"shape" %in% scale$aesthetics,
                              plot$scales$scales)
  for (j in seq_along(plot$layers)) {
    if (inherits(plot$layers[[j]]$geom, "GeomPoint")) {
      plot$layers[[j]]$mapping$shape <- NULL
      plot$layers[[j]]$aes_params$shape <- 16
      plot$layers[[j]]$aes_params$size <- 4
      plot$layers[[j]]$aes_params$alpha <- .9
      # Small vertical offsets keep coincident zero-delta points visible.
      plot$layers[[j]]$position <- ggplot2::position_dodge(width = .45)
    }
  }
  plot
}

review_figure <- function(plot) {
  if (identical(plot$meta$figure_review, "20261006"))
    return(reference_recovery_layout(plot))
  # Presentation only: retain all original plot data, fitted values and SEs.
  plot <- unserialize(serialize(plot, NULL))
  data_names <- names(plot$data)

  # Ascending row and column orders reproduce the reference matrix orientation.
  if (all(c("true_model", "selected_model", "prop") %in% data_names)) {
    model_order <- levels(plot$data$true_model)
    if (is.null(model_order)) model_order <- unique(as.character(plot$data$true_model))
    plot <- plot + ggplot2::scale_y_discrete(limits = model_order) +
      ggplot2::coord_equal(ratio = 0.85)
  }
  if (all(c("TrueParam", "RecoveredParam", "Correlation") %in% data_names)) {
    parameter_order <- levels(plot$data$TrueParam)
    if (is.null(parameter_order)) parameter_order <- unique(as.character(plot$data$TrueParam))
    old_y <- plot$scales$get_scales("y")
    parameter_labels <- if (is.null(old_y)) ggplot2::waiver() else old_y$labels
    plot <- plot +
      ggplot2::scale_y_discrete(limits = parameter_order, labels = parameter_labels) +
      ggplot2::scale_fill_gradient2(
        low = FIGURE_ORANGE, mid = "#FFFFFF", high = FIGURE_GREEN, midpoint = 0,
        limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1),
        name = "Correlation",
        labels = scales::label_number(accuracy = .5, style_negative = "minus"),
        guide = ggplot2::guide_colourbar(display = "rectangles")) +
      ggplot2::coord_equal(ratio = 0.85)
  }

  # Keep solid/dashed curves; explain their meaning in the caption.
  if (all(c("modality", "input_range", "weight") %in% data_names)) {
    plot <- plot + ggplot2::guides(
      colour = ggplot2::guide_legend(title = NULL, order = 1,
                                    override.aes = list(linetype = "solid")),
      linetype = "none")
    for (scale in plot$scales$scales) {
      if (any(scale$aesthetics %in% c("colour", "color"))) scale$name <- NULL
      if ("linetype" %in% scale$aesthetics) scale$guide <- "none"
    }
  }

  # Preserve the previously agreed broken-axis spacing.
  if (all(c("model", "Criterion", "Delta") %in% data_names)) {
    break_settings <- attr(plot, "axis_break_x", exact = TRUE)
    if (!is.null(break_settings) && inherits(break_settings, "ggbreak_params")) {
      break_settings$space <- 0.25
      attr(plot, "axis_break_x") <- break_settings
    }
  }
  if (all(c("model_label", "metric", "value") %in% data_names)) {
    model_scale <- plot$scales$get_scales("y")
    model_scale$labels <- function(x) sub("^Model [0-9]+ \\((.*)\\)$", "\\1", as.character(x))
  }

  # Give long condition labels room without reducing their font size.
  if (all(c("VV", "AV", "SE", "R_VA", "Model5plq") %in% data_names)) {
    condition_scale <- plot$scales$get_scales("x")
    condition_scale$expand <- ggplot2::expansion(add = .35)
  }

  # Existing observed-versus-predicted scatters use a true identity line.
  if (all(c("R_VA", "predicted") %in% data_names) &&
      identical(rlang::as_label(plot$mapping$x), "R_VA") &&
      identical(rlang::as_label(plot$mapping$y), "predicted")) {
    binary <- all(c("VA_disgust", "n_VA") %in% data_names)
    plot <- plot + ggplot2::labs(
      x = if (binary) "Observed Proportion" else "Observed Rating",
      y = if (binary) "Predicted Probability" else "Predicted Rating") +
      ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1))
  }

  for (j in seq_along(plot$layers)) {
    layer <- plot$layers[[j]]
    if (inherits(layer$geom, "GeomPoint") && identical(layer$aes_params$shape, 1)) {
      plot$layers[[j]]$aes_params$size <- 4
      plot$layers[[j]]$aes_params$stroke <- .6
      plot$layers[[j]]$show.legend <- FALSE
    }
    if (inherits(layer$geom, "GeomErrorbar"))
      plot$layers[[j]]$aes_params$linewidth <- .45
    if (inherits(layer$geom, "GeomText") && "Correlation" %in% data_names) {
      plot$layers[[j]]$mapping$label <- ggplot2::aes(
        label = sub("-", "\u2212", sprintf("%.2f", Correlation), fixed = TRUE))$label
    }
  }
  # Map colors to the displayed quantity without changing any source values.
  colour_mapping <- if (is.null(plot$mapping$colour)) "" else
    rlang::as_label(plot$mapping$colour)
  if (colour_mapping %in% c("AV", "emotion")) {
    condition_keys <- unique(as.character(plot$data[[colour_mapping]]))
    condition_order <- if (any(tolower(condition_keys) %in% c("anger", "mixed", "disgust")))
      c("anger", "mixed", "disgust") else c("neg", "neu", "pos")
    normalized_keys <- tolower(condition_keys)
    normalized_keys <- c(negative = "neg", neutral = "neu", positive = "pos")[normalized_keys] |>
      as.character()
    # Abbreviations and full emotion labels are both present in the saved plots.
    missing_keys <- is.na(normalized_keys)
    normalized_keys[missing_keys] <- tolower(condition_keys[missing_keys])
    values <- setNames(FIGURE_CONDITION_COLOURS[match(normalized_keys, condition_order)],
                       condition_keys)
    stopifnot(!anyNA(values))
    labels <- if (identical(condition_order[1], "anger")) c("Anger", "Mixed", "Disgust") else
      c("Negative", "Neutral", "Positive")
    label_values <- setNames(labels[match(normalized_keys, condition_order)], condition_keys)
    plot <- plot + ggplot2::scale_colour_manual(
      values = values, labels = label_values,
      name = if (identical(colour_mapping, "AV")) "Auditory\nCondition" else NULL)
  } else if (identical(colour_mapping, "modality")) {
    plot <- plot + ggplot2::scale_colour_manual(values = FIGURE_MODALITY_COLOURS, name = NULL)
  } else if (colour_mapping %in% c("Criterion", "metric")) {
    plot <- plot + ggplot2::scale_colour_manual(
      values = c(AIC = FIGURE_ORANGE, AICc = FIGURE_ORANGE, BIC = FIGURE_GREEN), name = NULL)
  } else if (identical(colour_mapping, "series")) {
    # Keep the visual/auditory assignment consistent with the weight panels.
    old_scale <- plot$scales$get_scales("colour")
    stopifnot(all(unique(as.character(plot$data$series)) %in% names(FIGURE_SERIES_COLOURS)))
    plot <- plot + ggplot2::scale_colour_manual(
      values = FIGURE_SERIES_COLOURS,
      labels = old_scale$labels, name = old_scale$name)
  } else if (identical(colour_mapping, "param_level")) {
    old_scale <- plot$scales$get_scales("colour")
    plot <- plot + ggplot2::scale_colour_manual(
      values = FIGURE_CONDITION_COLOURS, labels = old_scale$labels, name = old_scale$name)
  }
  if (all(c("true_model", "selected_model", "prop") %in% data_names)) {
    plot <- plot + ggplot2::scale_fill_gradient(
      low = "#FFFFFF", high = FIGURE_TEAL, limits = c(0, 1),
      name = "p(Rec|True)",
      guide = ggplot2::guide_colourbar(display = "rectangles"))
  }
  old_colours <- c("#6D28D9" = FIGURE_TEAL, "#8B5CF6" = FIGURE_ORANGE,
                   "#C4B5FD" = FIGURE_GREEN)
  for (j in seq_along(plot$layers)) {
    colour <- plot$layers[[j]]$aes_params$colour
    if (is.character(colour) && length(colour) == 1L && toupper(colour) %in% names(old_colours))
      plot$layers[[j]]$aes_params$colour <- unname(old_colours[toupper(colour)])
    # Dark text remains readable on the lighter reference heatmap colors.
    if (inherits(plot$layers[[j]]$geom, "GeomText") &&
        any(c("prop", "Correlation") %in% data_names)) {
      plot$layers[[j]]$mapping$colour <- NULL
      plot$layers[[j]]$aes_params$colour <- "#111111"
    }
  }
  plot$meta$figure_review <- "20261006"
  reference_recovery_layout(plot)
}

style_figure <- function(plot) {
  plot <- review_figure(plot)
  plot <- reference_comparison_layout(plot)
  # Keep the condition title on one line above the legend entries.
  for (scale in plot$scales$scales) {
    if (is.character(scale$name) && length(scale$name) == 1L) {
      if (!identical(scale$name, "p(Rec|True)")) scale$name <- tools::toTitleCase(scale$name)
      if (gsub("\n", " ", scale$name, fixed = TRUE) %in%
          c("Auditory Condition", "Auditory Emotion", "Auditory Valence"))
        scale$name <- "Auditory Condition"
    }
  }
  for (key in c("x", "y")) {
    if (is.character(plot$labels[[key]]) && length(plot$labels[[key]]) == 1L)
      plot$labels[[key]] <- tools::toTitleCase(plot$labels[[key]])
  }
  if (isTRUE(plot$labels$x %in% c("Unimodal Response", "Unimodal Disgust Rating",
                                "Unimodal Disgust-Choice Probability"))) {
    plot <- plot + ggplot2::scale_x_continuous(
      breaks = c(0, 0.5, 1), labels = c("0.0", "0.5", "1.0"), expand = c(0, 0))
  }
  unit_interval_data <- "R_VA" %in% names(plot$data) &&
    is.numeric(plot$data$R_VA) && all(is.finite(plot$data$R_VA)) &&
    all(plot$data$R_VA >= 0 & plot$data$R_VA <= 1)
  if (unit_interval_data &&
      isTRUE(plot$labels$x %in% c("Observed Response", "Observed Rating", "Observed Proportion")) &&
      isTRUE(plot$labels$y %in% c("Predicted Response", "Predicted Rating", "Predicted Probability"))) {
    plot <- plot +
      ggplot2::scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0.0", "0.5", "1.0")) +
      ggplot2::scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0.0", "0.5", "1.0"))
  }
  # One naming framework across studies, with the measured quantity retained.
  axis_names <- c("Visual Valence" = "Visual Condition",
                  "Auditory Valence" = "Auditory Condition",
                  "Visual Emotion" = "Visual Condition",
                  "Rating" = "Valence Rating",
                  "Raw Weight" = "Weight", "Modality Weight" = "Weight",
                  "Disgust-Choice Probability" = "Proportion of\nDisgust Choices",
                  "Proportion of Disgust Choices" = "Proportion of\nDisgust Choices")
  for (key in c("x", "y")) {
    label <- plot$labels[[key]]
    if (is.character(label) && length(label) == 1L) {
      label <- trimws(gsub("\n", " ", label, fixed = TRUE))
      if (label %in% names(axis_names)) plot$labels[[key]] <- unname(axis_names[label])
    }
    axis_scale <- plot$scales$get_scales(key)
    if (!is.null(axis_scale) && is.character(axis_scale$name) && length(axis_scale$name) == 1L) {
      label <- trimws(gsub("\n", " ", axis_scale$name, fixed = TRUE))
      if (label %in% names(axis_names)) axis_scale$name <- unname(axis_names[label])
    }
  }
  if (all(c("modality", "input_range", "weight") %in% names(plot$data))) {
    label <- tolower(trimws(gsub("\n", " ", plot$labels$x, fixed = TRUE)))
    plot$labels$y <- "Weight"
    if (identical(label, "valence rating")) plot$labels$x <- "Unimodal Valence\nRating"
    if (identical(label, "unimodal response")) plot$labels$x <- "Unimodal Disgust Rating"
    if (label %in% c("unimodal log-odds", "unimodal log odds")) plot$labels$x <- "Unimodal Log Odds"
  }
  plot <- plot + ggplot2::theme(
    text = ggplot2::element_text(family = FIGURE_FONT, face = "plain"),
    axis.title = ggplot2::element_text(size = FIGURE_AXIS_PT, face = "plain"),
    axis.title.x = ggplot2::element_text(size = FIGURE_AXIS_PT, face = "plain",
                                       margin = ggplot2::margin(t = 18), colour = "#111111"),
    axis.title.y = ggplot2::element_text(size = FIGURE_AXIS_PT, face = "plain",
                                       margin = ggplot2::margin(r = 10), colour = "#111111"),
    axis.text = ggplot2::element_text(size = FIGURE_TEXT_PT, face = "plain", colour = "#555555"),
    axis.text.x = ggplot2::element_text(size = FIGURE_TEXT_PT, face = "plain", colour = "#555555"),
    axis.text.y = ggplot2::element_text(size = FIGURE_TEXT_PT, face = "plain", colour = "#555555"),
    axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.6),
    legend.title = ggplot2::element_text(size = 26, face = "plain", colour = "#111111"),
    legend.text = ggplot2::element_text(size = 26, face = "plain", colour = "#333333",
                                        margin = ggplot2::margin(l = 5, r = 8)),
    legend.position = "top", legend.direction = "horizontal",
    legend.title.position = "top",
    legend.key.width = grid::unit(28, "pt"),
    legend.key.spacing.x = grid::unit(14, "pt"),
    legend.key.height = grid::unit(24, "pt"),
    legend.spacing.y = grid::unit(6, "pt"),
    legend.margin = ggplot2::margin(0, 0, 10, 0),
    legend.box.spacing = grid::unit(12, "pt"),
    legend.box.margin = ggplot2::margin(l = 8),
    strip.text = ggplot2::element_text(size = FIGURE_TEXT_PT, face = "plain"),
    plot.title = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(t = 12, r = 24, b = 12, l = 12),
    panel.border = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(colour = "black", linewidth = 0.4)
  )
  # Separate the two presentation labels within each simulation facet.
  if (all(c("type", "mean_slope", "param_level") %in% names(plot$data))) {
    plot <- plot + ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = .45)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = 24),
                     axis.text.x.bottom = ggplot2::element_text(margin = ggplot2::margin(t = 8)))
  }
  dense <- any(vapply(plot$layers, function(layer) inherits(layer$geom, "GeomTile"), logical(1)))
  # Text annotations use millimetres in ggplot2, unlike theme text in points.
  for (j in seq_along(plot$layers)) {
    if (inherits(plot$layers[[j]]$geom, "GeomText")) {
      plot$layers[[j]]$aes_params$family <- FIGURE_FONT
      plot$layers[[j]]$aes_params$fontface <- "plain"
      plot$layers[[j]]$aes_params$size <- (if (dense) FIGURE_CELL_PT else FIGURE_TEXT_PT) * 25.4 / 72.27
    }
  }
  plot
}

export_figure <- function(filename, plot, width, height, dpi = FIGURE_DPI, ...) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  if (!capabilities("cairo")) stop("Vector export requires R with Cairo support.")
  ggplot2::ggsave(filename, plot, width = width, height = height, dpi = dpi, bg = "white", ...)
  stem <- sub("\\.png$", "", filename)
  svg_file <- paste0(stem, ".svg")
  ggplot2::ggsave(svg_file, plot, width = width, height = height,
                  device = svglite::svglite, bg = "white", ...)
  svg_lines <- readLines(svg_file, warn = FALSE)
  header <- grep("^<svg ", svg_lines)
  svg_lines[header] <- sub('width="[^"]+"', paste0('width="', width, 'in"'), svg_lines[header])
  svg_lines[header] <- sub('height="[^"]+"', paste0('height="', height, 'in"'), svg_lines[header])
  writeLines(svg_lines, svg_file, useBytes = TRUE)
  ggplot2::ggsave(paste0(stem, ".pdf"), plot, width = width, height = height,
                  device = grDevices::cairo_pdf, bg = "white", ...)
  invisible(plot)
}

save_figure <- function(filename, plot = ggplot2::last_plot(),
                        width = FIGURE_WIDTH_IN, height = FIGURE_HEIGHT_IN,
                        dpi = FIGURE_DPI, ...) {
  plot <- style_figure(plot)
  # Give recovery matrices room for readable cell values and axis labels.
  dense <- any(vapply(plot$layers, function(layer) inherits(layer$geom, "GeomTile"), logical(1)))
  if (dense) {
    width <- max(width, 7.5)
    height <- max(height, 6.5)
    correlation <- "Correlation" %in% names(plot$data)
    values <- if (correlation) plot$data$Correlation else plot$data$prop
    plot$data$.matrix_label <- sub("^(-?)0\\.", "\\1.", sprintf("%.2f", values))
    fill_scale <- plot$scales$get_scales("fill")
    fill_scale$breaks <- if (correlation) c(-1, 0, 1) else c(0, .5, 1)
    for (j in seq_along(plot$layers)) {
      if (inherits(plot$layers[[j]]$geom, "GeomText")) {
        plot$layers[[j]]$mapping$label <- ggplot2::aes(label = .matrix_label)$label
      }
    }
    plot <- plot + ggplot2::labs(
      y = if (correlation) "Recovered\nParameter" else "Recovered\nModel") +
      ggplot2::theme(
        legend.position = "bottom", legend.direction = "horizontal",
        legend.title = ggplot2::element_text(size = 26, margin = ggplot2::margin(b = 8)),
        legend.text = ggplot2::element_text(size = 26, margin = ggplot2::margin(t = 9)),
        legend.ticks = ggplot2::element_blank(),
        legend.margin = ggplot2::margin(t = 8, b = 6)) +
      ggplot2::guides(fill = ggplot2::guide_colourbar(
        title.position = "top", label.position = "bottom", display = "rectangles", nbin = 50,
        barwidth = grid::unit(145, "pt"), barheight = grid::unit(13, "pt")))
  }
  export_figure(filename, plot, width, height, dpi, ...)
  # The assembled results figures reuse these styled plots; no fits are rerun.
  saveRDS(plot, sub("\\.png$", ".rds", filename))
  invisible(plot)
}

# Plot model selections on a common count scale.
plot_selection_frequency <- function(summary) {
  if (!"Study" %in% names(summary)) summary$Study <- ""
  selection_counts <- summary %>% filter(omitted_participant != "None") %>%
    select(Study, best_AICc, best_BIC) %>%
    tidyr::pivot_longer(c(best_AICc, best_BIC), names_prefix = "best_",
                        names_to = "Criterion", values_to = "model") %>%
    count(Study, Criterion, model, name = "n_selected") %>%
    mutate(Criterion = factor(Criterion, levels = c("AICc", "BIC")),
           model = factor(model, levels = c("Visual-only", "Auditory-only", "3p", "4p", "4plq", "5p", "5plq", "7p"))) %>%
    arrange(Study, Criterion, model) %>% group_by(Study, Criterion) %>%
    mutate(proportion = n_selected / sum(n_selected),
           label_position = cumsum(n_selected) - n_selected / 2,
           label = paste0(model, "\n", n_selected, "/", sum(n_selected)),
           row = trimws(paste(Study, Criterion))) %>% ungroup()
  selection_counts$row <- factor(selection_counts$row,
                                levels = rev(unique(selection_counts$row)))
  model_colours <- c("Visual-only" = "#E0E0E0", "Auditory-only" = "#CCCCCC",
                     "4p" = "#8D9DAB", "5p" = "#B9A0BD",
                     "3p" = FIGURE_TEAL, "4plq" = FIGURE_GREEN,
                     "5plq" = "#B0B0B0", "7p" = FIGURE_ORANGE)
  p_selection <- ggplot(selection_counts, aes(x = n_selected, y = row, fill = model)) +
    geom_col(width = .64, colour = "white", linewidth = .6,
             position = position_stack(reverse = TRUE)) +
    geom_text(aes(x = label_position, label = label),
              family = FIGURE_FONT, colour = "#111111", lineheight = .95) +
    scale_fill_manual(values = model_colours, guide = "none") +
    scale_x_continuous(breaks = 0:10, expand = expansion(mult = c(0, .015))) +
    labs(x = "Selection Count", y = NULL) +
    theme_classic(base_size = FIGURE_TEXT_PT, base_family = FIGURE_FONT)
  style_figure(p_selection) +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.y = element_text(colour = "#111111", hjust = 0),
          plot.margin = margin(12, 18, 12, 12))
}

# Combine the saved study summaries after both analyses have been run.
export_combined_selection_frequency <- function() {
  summary_files <- file.path(PROJECT_ROOT, "Results", c("Study2", "Study3"),
                              "ParticipantSensitivity", "Summary.csv")
  if (all(file.exists(summary_files))) {
    combined <- bind_rows(lapply(seq_along(summary_files), function(i) {
      read.csv(summary_files[i]) %>% mutate(Study = paste("Study", i + 1))
    }))
    plot <- plot_selection_frequency(combined)
    filename <- file.path(PROJECT_ROOT, "Figures", "ParticipantSensitivity",
                          "SelectionFrequency.png")
    export_figure(filename, plot, width = 11, height = 6.5)
    saveRDS(plot, sub("\\.png$", ".rds", filename))
  }
}
