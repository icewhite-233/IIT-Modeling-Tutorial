# Step 4: descriptive summaries, model comparison and fit plots.
# Read existing fitted objects; no refitting.
# Colored lines show observations; open circles show predictions.

rm(list = ls())
# Load shared paths relative to this script.
.script_files <- unlist(lapply(sys.frames(), function(frame) frame$ofile))
.script_file <- if (length(.script_files)) tail(.script_files, 1) else
  sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
if (!length(.script_file)) .script_file <- if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) rstudioapi::getSourceEditorContext()$path else "Step4_Plotting.R"
source(file.path(dirname(normalizePath(.script_file)), "..", "Project_Config.R"), encoding = "UTF-8")
setwd(file.path(PROJECT_ROOT, "Scripts", "Study3"))
library(dplyr)
library(tidyr)
library(ggplot2)
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

## 1. Load data and fitted models ------------------------------------------
FIT_FILE <- "results/model_fits_weighted_averaging_BinomialLikelihood_54.rds"
Y_LABEL <- "Proportion of Disgust Choices"
FIT_PREFIX <- "fit_"
dat <- read.csv(file.path(PROJECT_ROOT, "Data", "Study3", "Input_Data54.csv"))
objects <- readRDS(FIT_FILE)
source("help_function_choice.R")
SCALE_MID <- 0
models <- c("3p", "4p", "4plq", "5p", "5plq", "7p")
category <- function(x) factor(x, levels = c(-1, 0, 1),
                              labels = c("Anger", "Mixed", "Disgust"))
palette <- c(Anger = FIGURE_TEAL, Mixed = FIGURE_ORANGE, Disgust = FIGURE_GREEN)
plot_theme <- theme_bw(base_size = FIGURE_TEXT_PT, base_family = FIGURE_FONT) +
  theme(legend.position = "right",
        legend.title = element_text(size = FIGURE_TEXT_PT, face = "plain"),
        legend.text = element_text(size = FIGURE_TEXT_PT),
        axis.title.x = element_text(size = FIGURE_AXIS_PT, face = "plain", margin = margin(t = 10)),
        axis.title.y = element_text(size = FIGURE_AXIS_PT, face = "plain", margin = margin(r = 10)),
        axis.text.x = element_text(size = FIGURE_TEXT_PT, face = "plain"),
        axis.text.y = element_text(size = FIGURE_TEXT_PT),
        panel.grid = element_blank())

## 2. Descriptive summaries; count each component stimulus once ------------
visual <- dat %>% distinct(visual_label, VV, R_V) %>%
  transmute(modality = "Visual", emotion = VV, rating = R_V)
auditory <- dat %>% distinct(auditory_label, AV, R_A) %>%
  transmute(modality = "Auditory", emotion = AV, rating = R_A)
single <- bind_rows(visual, auditory) %>% group_by(modality, emotion) %>%
  summarise(observed = mean(rating), SE = sd(rating)/sqrt(n()),
            n_items = n(), .groups = "drop") %>%
  mutate(emotion = category(emotion),
         modality = factor(modality, levels = c("Visual", "Auditory")))
observed <- dat %>% group_by(VV, AV) %>%
  summarise(observed = mean(R_VA), SE = sd(R_VA)/sqrt(n()),
            n_items = n(), .groups = "drop")
observed <- observed %>% mutate(VV = category(VV), AV = category(AV))
write.csv(single, "results/Plot_Unimodal.csv", row.names = FALSE)
write.csv(observed, "results/Plot_AV_Observed.csv", row.names = FALSE)

p_single <- ggplot(single, aes(emotion, observed, color = emotion)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = observed-SE, ymax = observed+SE), width = .08) +
  facet_wrap(~modality) + scale_color_manual(values = palette) +
  scale_y_continuous(breaks = seq(0,1,.2)) +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Emotion category", y = Y_LABEL) +
  plot_theme + theme(legend.position = "none")
p_av <- ggplot(observed, aes(VV, observed, group = AV, color = AV)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  geom_errorbar(aes(ymin = observed-SE, ymax = observed+SE), width = .08) +
  scale_color_manual(name = "Auditory Condition", values = palette) +
  scale_y_continuous(breaks = seq(0,1,.2)) + coord_cartesian(ylim = c(0,1)) +
  labs(x = "Visual Condition", y = Y_LABEL) + plot_theme
save_figure("figures/01_Unimodal.png", p_single, width = FIGURE_WIDTH_IN, height = FIGURE_HEIGHT_IN, dpi = FIGURE_DPI)
save_figure("figures/02_AV_Observed.png", p_av, width = FIGURE_WIDTH_IN, height = FIGURE_HEIGHT_IN, dpi = FIGURE_DPI)

## 3. Compare criterion-specific differences from the minimum --------------
cmp <- objects$cmp_tbl %>% mutate(model = sub(" .*", "", model))
comparison <- cmp %>% select(model, AIC, AICc, BIC) %>%
  pivot_longer(-model, names_to = "Criterion", values_to = "IC") %>%
  group_by(Criterion) %>% mutate(Delta = IC-min(IC)) %>% ungroup() %>%
  mutate(model = factor(model, levels = rev(c("Visual-only", "Auditory-only", models))))
comparison_display <- comparison %>% filter(Criterion %in% c("AICc", "BIC"))
p_cmp <- ggplot(comparison_display, aes(Delta, model, color = Criterion, shape = Criterion)) +
  geom_point(position = position_dodge(width = .45), size = 3) +
  scale_color_manual(name = NULL, values = c(AICc = FIGURE_ORANGE, BIC = FIGURE_GREEN)) +
  scale_shape_manual(name = NULL, values = c(AICc = 16, BIC = 17)) +
  labs(x = "\u0394 information criterion", y = "Model") + plot_theme + theme(legend.position = "top")
# Display every model while omitting an empty interval of the x axis.
# The two displayed segments have separate plotting scales.
stopifnot(!any(comparison$Delta > 25 & comparison$Delta < 950))
p_cmp <- p_cmp +
  scale_x_continuous(limits = c(-1, 1210),
                     breaks = c(0, 10, 20, 1000, 1200), expand = c(0, 0)) +
  ggbreak::scale_x_break(c(25, 950), scales = 0.55,
                        space = 0.18, symbol = "slash", expand = FALSE)
save_figure("figures/03_ModelComparison.png", p_cmp, width = FIGURE_WIDTH_IN, height = FIGURE_HEIGHT_IN, dpi = FIGURE_DPI)
write.csv(comparison, "results/Plot_ModelComparison.csv", row.names = FALSE)

## 4. Overlay each model's predictions on the observed summaries -----------
dat$L_V <- choice_prob_to_logodds(dat$R_V, dat$n_V)
dat$L_A <- choice_prob_to_logodds(dat$R_A, dat$n_A)
predictions <- list()
cells <- list()
for (model in models) {
  fit <- objects[[paste0(FIT_PREFIX, model)]]
  pred <- as.numeric(predict(fit, newdata = dat, type = "probability"))
  one <- dat %>% mutate(model = model, predicted = pred,
                         residual = R_VA-pred)
  predictions[[model]] <- one
  means <- one %>% group_by(VV, AV) %>%
    summarise(predicted = mean(predicted), .groups = "drop") %>%
    mutate(VV = category(VV), AV = category(AV))
  cell <- left_join(observed, means, by = c("VV", "AV"))
  cells[[model]] <- cell %>% mutate(model = model)
  # Reuse the observed condition means and SEs for each model.
  p_fit <- p_av +
    geom_point(data = cell, aes(y = predicted), shape = 1, size = 3, color = "black", show.legend = FALSE)
  p_scatter <- ggplot(one, aes(R_VA, predicted)) +
    geom_abline(slope = 1, intercept = 0, color = "grey55", linetype = "dashed") +
    geom_point(color = FIGURE_TEAL, size = 2.5) + coord_equal(xlim = c(0,1), ylim = c(0,1)) +
    labs(x = "Observed Proportion", y = "Predicted Probability") + plot_theme
  save_figure(paste0("figures/04_Fit_", model, ".png"), p_fit, width = FIGURE_WIDTH_IN, height = FIGURE_HEIGHT_IN, dpi = FIGURE_DPI)
  save_figure(paste0("figures/05_Scatter_", model, ".png"), p_scatter, width = FIGURE_WIDTH_IN, height = FIGURE_HEIGHT_IN, dpi = FIGURE_DPI)
}
write.csv(bind_rows(cells), "results/Plot_AV_AllModels.csv", row.names = FALSE)
write.csv(bind_rows(predictions), "results/Plot_Predictions_AllModels.csv", row.names = FALSE)
cat("绘图完成：figures内含单通道、联合描述、模型比较和六模型拟合图。\n")
