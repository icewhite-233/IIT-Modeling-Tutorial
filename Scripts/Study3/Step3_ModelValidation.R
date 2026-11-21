# Study 3: predictions and remaining condition effects under the retained specification.
rm(list = ls())
# Load shared paths relative to this script.
.script_files <- unlist(lapply(sys.frames(), function(frame) frame$ofile))
.script_file <- if (length(.script_files)) tail(.script_files, 1) else
  sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
if (!length(.script_file)) .script_file <- if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) rstudioapi::getSourceEditorContext()$path else "Step3_ModelValidation.R"
source(file.path(dirname(normalizePath(.script_file)), "..", "Project_Config.R"), encoding = "UTF-8")
setwd(file.path(PROJECT_ROOT, "Scripts", "Study3"))
library(dplyr)
library(ggplot2)
source("help_function_choice.R")
SCALE_MID <- 0
dir.create("figures", showWarnings = FALSE)

## 1. Load the retained fit and compute predictions and residuals ----------
dat <- choice_prepare_inputs(read.csv(file.path(PROJECT_ROOT, "Data", "Study3", "Input_Data54.csv"))) %>% mutate(
 VA_anger = n_VA - VA_disgust,
 VV = factor(VV, levels = c(-1, 0, 1), labels = c("anger", "mixed", "disgust")),
 AV = factor(AV, levels = c(-1, 0, 1), labels = c("anger", "mixed", "disgust")))
objects <- readRDS("results/model_fits_weighted_averaging_BinomialLikelihood_54.rds")
best_model <- objects$best_model
fit <- objects[[paste0("fit_", best_model)]]
dat$eta <- predict(fit, dat, type = "link")
dat$predicted <- plogis(dat$eta)
dat$residual <- dat$R_VA - dat$predicted
write.csv(dat, "results/Validation_BestModel_Predictions.csv", row.names = FALSE)

## 2. Test condition corrections beyond the fixed IIT prediction -----------
# Keep original counts as responses; offset(eta) fixes the IIT coefficient at one.
# This is an approximate in-sample diagnostic with the IIT parameters fixed.
compare_glm <- function(reduced, full, effect) {
 lr <- 2 * as.numeric(logLik(full) - logLik(reduced))
 df <- attr(logLik(full), "df") - attr(logLik(reduced), "df")
 data.frame(effect, LR = lr, df, p = pchisq(max(lr, 0), df, lower.tail = FALSE))
}
base_ll <- sum(dbinom(dat$VA_disgust, dat$n_VA, dat$predicted, log = TRUE))
gV <- glm(cbind(VA_disgust, VA_anger) ~ VV + offset(eta), family = binomial(), data = dat)
gA <- glm(cbind(VA_disgust, VA_anger) ~ AV + offset(eta), family = binomial(), data = dat)
gAdd <- glm(cbind(VA_disgust, VA_anger) ~ VV + AV + offset(eta), family = binomial(), data = dat)
gFull <- glm(cbind(VA_disgust, VA_anger) ~ VV * AV + offset(eta), family = binomial(), data = dat)
lr <- 2 * (as.numeric(logLik(gFull)) - base_ll)
validation <- bind_rows(
  compare_glm(gA, gAdd, "visual_given_auditory"),
  compare_glm(gV, gAdd, "auditory_given_visual"),
  compare_glm(gAdd, gFull, "interaction"),
  data.frame(effect = "global9_including_intercept", LR = lr, df = 9,
             p = pchisq(max(lr, 0), 9, lower.tail = FALSE))
) %>% mutate(model = best_model)
print(validation)
write.csv(validation, "results/Validation_BestModel_OffsetGLM.csv", row.names = FALSE)

## 3. Observed proportions, predicted probabilities and pair-level SEs -----
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
  labs(x = "Visual Condition", y = "Proportion of Disgust Choices") +
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
