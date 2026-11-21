# Study 2: predictions, residual ANOVA and fit plots for the retained specification.
rm(list = ls())
# Load shared paths relative to this script.
.script_files <- unlist(lapply(sys.frames(), function(frame) frame$ofile))
.script_file <- if (length(.script_files)) tail(.script_files, 1) else
  sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
if (!length(.script_file)) .script_file <- if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) rstudioapi::getSourceEditorContext()$path else "Step3_ModelValidation.R"
source(file.path(dirname(normalizePath(.script_file)), "..", "Project_Config.R"), encoding = "UTF-8")
setwd(file.path(PROJECT_ROOT, "Scripts", "Study2"))
library(dplyr)
library(ggplot2)
library(car)
source("simulation_helpers.R")
dir.create("figures", showWarnings = FALSE)

## 1. Load the retained fit and compute predictions and residuals ----------
dat <- read.csv(file.path(PROJECT_ROOT, "Data", "Study2", "Input_Data54.csv")) %>%
  mutate(VV = factor(VV, levels = c(-1, 0, 1), labels = c("anger", "mixed", "disgust")),
         AV = factor(AV, levels = c(-1, 0, 1), labels = c("anger", "mixed", "disgust")))
objects <- readRDS("results/model_fits_weighted_averaging.rds")
# Select the stored fit with lowest BIC (3p for these data).
best_model <- objects$best_model
fit <- objects[[paste0("nls_", best_model)]]
dat$predicted <- as.numeric(predict(fit, newdata = dat))
dat$residual <- dat$R_VA - dat$predicted
write.csv(dat, "results/Validation_BestModel_Predictions.csv", row.names = FALSE)

## 2. Type-III residual ANOVA ----------------------------------------------
# Sum coding applies to these Type III diagnostics, not IIT estimation.
contrasts(dat$VV) <- contr.sum(nlevels(dat$VV))
contrasts(dat$AV) <- contr.sum(nlevels(dat$AV))
aov_residual <- aov(residual ~ VV * AV, data = dat)
validation <- car::Anova(aov_residual, type = "III")
print(validation)
validation_tbl <- data.frame(effect = rownames(validation), unclass(validation),
                             row.names = NULL, check.names = FALSE)
write.csv(validation_tbl, "results/Validation_BestModel_ResidualANOVA.csv", row.names = FALSE)

## 3. Observed and predicted condition means with pair-level SEs -----------
# Predict each pair first, then average predictions within conditions.
observed <- dat %>% group_by(VV, AV) %>%
  summarise(observed = mean(R_VA), se = sd(R_VA)/sqrt(n()),
            n_items = n(), .groups = "drop")
means <- dat %>% group_by(VV, AV) %>%
  summarise(predicted = mean(predicted), .groups = "drop")
cells <- left_join(observed, means, by = c("VV", "AV"))
write.csv(cells, "results/Validation_BestModel_CellMeans.csv", row.names = FALSE)
p_fit <- ggplot(cells, aes(VV, observed, group = AV, color = AV)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  geom_errorbar(aes(ymin = observed-se, ymax = observed+se), width = .08) +
  geom_point(aes(y = predicted), shape = 1, size = 3, color = "black", show.legend = FALSE) +
  scale_color_manual(name = "Auditory Condition", values = c(anger = FIGURE_TEAL, mixed = FIGURE_ORANGE, disgust = FIGURE_GREEN),
                     labels = c("Anger", "Mixed", "Disgust")) +
  scale_x_discrete(labels = c("Anger", "Mixed", "Disgust")) +
  scale_y_continuous(breaks = seq(0,1,.2)) + coord_cartesian(ylim = c(0,1)) +
  labs(x = "Visual Condition", y = "Disgust Rating") +
  theme_bw(base_size = FIGURE_TEXT_PT, base_family = FIGURE_FONT) +
  theme(legend.position = "right",
        legend.title = element_text(size = FIGURE_TEXT_PT, face = "plain"),
        legend.text = element_text(size = FIGURE_TEXT_PT),
        axis.title.x = element_text(size = FIGURE_AXIS_PT, face = "plain", margin = margin(t = 10)),
        axis.title.y = element_text(size = FIGURE_AXIS_PT, face = "plain", margin = margin(r = 10)),
        axis.text.x = element_text(size = FIGURE_TEXT_PT, face = "plain"),
        axis.text.y = element_text(size = FIGURE_TEXT_PT),
        panel.grid = element_blank())
save_figure("figures/Validation_BestModel.png", p_fit, width = FIGURE_WIDTH_IN, height = FIGURE_HEIGHT_IN, dpi = FIGURE_DPI)
