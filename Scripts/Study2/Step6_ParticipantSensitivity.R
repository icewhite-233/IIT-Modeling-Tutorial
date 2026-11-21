# Study 2: omit each participant and repeat preparation and least-squares fitting.
# Compare within-sample information criteria and changes in fitted predictions.
# Run Steps 1-2 first, then omit participants and compare with the full fit.
rm(list = ls())

# Select the omissions to run. The full sample is always fitted first.
omit_participants <- sprintf("P%02d", 1:10)
# Use "P01" for one omission, or character(0) for the full sample only.
.script_files <- unlist(lapply(sys.frames(), function(frame) frame$ofile))
.script_file <- if (length(.script_files)) tail(.script_files, 1) else
  sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
if (!length(.script_file)) .script_file <- if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) rstudioapi::getSourceEditorContext()$path else "Step6_ParticipantSensitivity.R"
source(file.path(dirname(normalizePath(.script_file)), "..", "Project_Config.R"), encoding = "UTF-8")
setwd(file.path(PROJECT_ROOT, "Scripts", "Study2"))
source("simulation_helpers.R")
library(dplyr)
library(minpack.lm)
library(ggplot2)

out_dir <- file.path(PROJECT_ROOT, "Results", "Study2", "ParticipantSensitivity")
figure_dir <- file.path(PROJECT_ROOT, "Figures", "Study2", "ParticipantSensitivity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## 1. Read the same formal trials as Step 1 --------------------------------
files <- list.files(file.path(PROJECT_ROOT, "Data", "Study2", "rawdata"),
                    pattern = "_Continuum_.*\\.csv$", full.names = TRUE)
raw <- bind_rows(lapply(files, function(file) {
  x <- read.csv(file, colClasses = "character", na.strings = c("", "None"))
  x %>% transmute(
    participant_id = sub("_.*", "", basename(file)),
    block = as.integer(block), trial = as.integer(trial), modality,
    face_file, voice_file, face_model_id,
    face_morph = as.numeric(face_morph),
    voice_model_id = sub("\\.0$", "", voice_model_id),
    voice_morph = as.numeric(voice_morph), response = as.numeric(resp)
  )
})) %>%
  filter(!is.na(block), !is.na(trial), !is.na(response),
         modality %in% c("vision", "auditory", "av")) %>%
  mutate(
    visual_label = paste(face_model_id, face_morph, sep = "_"),
    auditory_label = paste(voice_model_id, voice_morph, sep = "_"),
    VV = match(face_morph, c(20, 50, 80)) - 2L,
    AV = match(voice_morph, c(20, 50, 80)) - 2L
  )

prepare_pairs <- function(trials) {
  # Repeat Step 1 after omission, including the unimodal predictors.
  visual <- trials %>% filter(modality == "vision") %>%
    group_by(face_file) %>% summarise(
      R_V = mean(response), n_V = n(), n_participants_V = n_distinct(participant_id),
      .groups = "drop")
  auditory <- trials %>% filter(modality == "auditory") %>%
    group_by(voice_file) %>% summarise(
      R_A = mean(response), n_A = n(), n_participants_A = n_distinct(participant_id),
      .groups = "drop")
  trials %>% filter(modality == "av") %>%
    group_by(face_file, voice_file, visual_label, auditory_label, VV, AV) %>%
    summarise(R_VA = mean(response), n_VA = n(),
              n_participants_VA = n_distinct(participant_id), .groups = "drop") %>%
    left_join(visual, by = "face_file") %>%
    left_join(auditory, by = "voice_file") %>%
    arrange(visual_label, auditory_label) %>%
    mutate(av_pair_id = sprintf("P%02d", row_number())) %>%
    select(av_pair_id, visual_label, auditory_label, R_V, R_A, R_VA, VV, AV,
           n_V, n_A, n_VA, n_participants_V, n_participants_A, n_participants_VA)
}

## 2. Evaluate fitted predictions and weights ------------------------------
model_predictions <- function(model, fit, dat) {
  if (model %in% c("Visual-only", "Auditory-only")) return(tibble(
    predicted = if (model == "Visual-only") dat$R_V else dat$R_A,
    W_V = NA_real_, W_A = NA_real_, denominator = NA_real_))
  p <- coef(fit)
  x_V <- dat$R_V - SCALE_MID
  x_A <- dat$R_A - SCALE_MID
  # Predictions use the model helpers; weights are also recorded for diagnostics.
  weights <- switch(model,
    `3p` = list(rep(p["W_V"], nrow(dat)), rep(p["W_A"], nrow(dat))),
    `4p` = list(p["A_V"] + p["B_VA"] * x_V, p["A_A"] + p["B_VA"] * x_A),
    `4plq` = list(p["A_VA"] + p["B_VA"] * x_V + p["B_VA2"] * x_V^2,
                  p["A_VA"] + p["B_VA"] * x_A + p["B_VA2"] * x_A^2),
    `5p` = list(p["A_V"] + p["B_V"] * x_V, p["A_A"] + p["B_A"] * x_A),
    `5plq` = list(p["A_V"] + p["B_VA"] * x_V + p["B_VA2"] * x_V^2,
                  p["A_A"] + p["B_VA"] * x_A + p["B_VA2"] * x_A^2),
    `7p` = list(p["A_V"] + p["B_V"] * x_V + p["B_V2"] * x_V^2,
                p["A_A"] + p["B_A"] * x_A + p["B_A2"] * x_A^2)
  )
  tibble(predicted = as.numeric(predict(fit, newdata = dat)),
         W_V = as.numeric(weights[[1]]), W_A = as.numeric(weights[[2]]),
         denominator = 1 + W_V + W_A)
}

## 3. Fit the full-sample reference ----------------------------------------
full_recomputed <- prepare_pairs(raw)
saved_dat <- read.csv(file.path(PROJECT_ROOT, "Data", "Study2", "Input_Data54.csv"))
# Match pair identities, then retain the prepared CSV's exact input values.
stopifnot(identical(full_recomputed$av_pair_id, saved_dat$av_pair_id),
          identical(full_recomputed$visual_label, saved_dat$visual_label),
          identical(full_recomputed$auditory_label, saved_dat$auditory_label))
full_dat <- bind_cols(saved_dat, full_recomputed %>%
  select(n_V, n_A, n_VA, n_participants_V, n_participants_A, n_participants_VA))

## 4. Select data, fit the models, and compare -------------------
inputs <- comparisons <- coefficients <- predictions <- fold_counts <- list()
for (omitted in c("None", omit_participants)) {
  trials <- if (omitted == "None") raw else filter(raw, participant_id != omitted)
  dat <- if (omitted == "None") full_dat else prepare_pairs(trials)
  # Fit the six models with the same explicit calls as Step 2.
  nls_3p <- nlsLM(
    R_VA ~ wavg_3(S0, W_V, W_A, R_V, R_A),
    data = dat, start = start_params[["3p"]],
    control = nls.lm.control(maxiter = 500)
  )
  nls_4p <- nlsLM(
    R_VA ~ wavg_4(S0, A_V, A_A, B_VA, R_V, R_A),
    data = dat, start = start_params[["4p"]],
    control = nls.lm.control(maxiter = 1000)
  )
  nls_4plq <- nlsLM(
    R_VA ~ wavg_4plq(S0, A_VA, B_VA, B_VA2, R_V, R_A),
    data = dat, start = start_params[["4plq"]],
    control = nls.lm.control(maxiter = 1000)
  )
  nls_5p <- nlsLM(
    R_VA ~ wavg_5(S0, A_V, B_V, A_A, B_A, R_V, R_A),
    data = dat, start = start_params[["5p"]],
    control = nls.lm.control(maxiter = 1000)
  )
  nls_5plq <- nlsLM(
    R_VA ~ wavg_5plq(S0, A_V, A_A, B_VA, B_VA2, R_V, R_A),
    data = dat, start = start_params[["5plq"]],
    control = nls.lm.control(maxiter = 1000)
  )
  nls_7p <- nlsLM(
    R_VA ~ wavg_7(S0, A_V, B_V, B_V2, A_A, B_A, B_A2, R_V, R_A),
    data = dat, start = start_params[["7p"]],
    control = nls.lm.control(maxiter = 1000)
  )
  models <- list(
    `3p` = nls_3p, `4p` = nls_4p,
    `4plq` = nls_4plq, `5p` = nls_5p,
    `5plq` = nls_5plq, `7p` = nls_7p
  )
  comparison <- bind_rows(lapply(names(models), function(model) {
    fitted_model <- models[[model]]
    tibble(model = model, k_par = length(coef(fitted_model)),
           rss = sum(residuals(fitted_model)^2), converged = fitted_model$convInfo$isConv)
  }), tibble(
    model = c("Visual-only", "Auditory-only"), k_par = 0L, converged = TRUE,
    rss = c(sum((dat$R_VA - dat$R_V)^2), sum((dat$R_VA - dat$R_A)^2))
  )) %>% mutate(
    n_pairs = nrow(dat), k_IC = k_par + 1L,
    AIC = n_pairs * log(rss / n_pairs) + 2 * k_IC,
    AICc = AIC + 2 * k_IC * (k_IC + 1) / (n_pairs - k_IC - 1),
    BIC = n_pairs * log(rss / n_pairs) + k_IC * log(n_pairs),
    R2 = 1 - rss / sum((dat$R_VA - mean(dat$R_VA))^2),
    delta_AICc = AICc - min(AICc), delta_BIC = BIC - min(BIC)
  )
  fit <- list(models = models, comparison = comparison)
  if (omitted == "None") full_fit <- fit
  stopifnot(identical(dat$av_pair_id, full_dat$av_pair_id),
            identical(dat$visual_label, full_dat$visual_label),
            identical(dat$auditory_label, full_dat$auditory_label))
  best_AICc <- fit$comparison$model[which.min(fit$comparison$delta_AICc)]
  best_BIC <- fit$comparison$model[which.min(fit$comparison$delta_BIC)]
  inputs[[omitted]] <- dat %>% mutate(omitted_participant = omitted, .before = 1)
  fold_counts[[omitted]] <- tibble(
    omitted_participant = omitted, n_participants = n_distinct(trials$participant_id),
    n_formal_trials = nrow(trials), n_visual_trials = sum(trials$modality == "vision"),
    n_auditory_trials = sum(trials$modality == "auditory"),
    n_av_trials = sum(trials$modality == "av"), n_pairs = nrow(dat),
    best_AICc = best_AICc, best_BIC = best_BIC, criteria_agree = best_AICc == best_BIC,
    n_failed_candidates = sum(!fit$comparison$converged)
  )
  comparisons[[omitted]] <- fit$comparison %>% mutate(
    omitted_participant = omitted, n_participants = n_distinct(trials$participant_id),
    focal_3p = model == "3p", preferred_AICc = model == best_AICc,
    preferred_BIC = model == best_BIC, .before = 1)
  coefficients[[omitted]] <- bind_rows(lapply(names(start_params), function(model) {
    fit_model <- fit$models[[model]]
    terms <- names(start_params[[model]])
    tibble(omitted_participant = omitted, model = model, term = terms,
           estimate = unname(coef(fit_model)[terms]))
  }))
  predictions[[omitted]] <- bind_rows(lapply(fit$comparison$model, function(model) {
    fit_model <- fit$models[[model]]
    actual <- model_predictions(model, fit_model, dat)
    reference <- model_predictions(model, fit_model, full_dat)
    full_reference <- model_predictions(model, full_fit$models[[model]], full_dat)
    bind_cols(dat %>% select(av_pair_id, visual_label, auditory_label, VV, AV, R_V, R_A, R_VA),
              actual) %>% mutate(
      omitted_participant = omitted, model = model, focal_3p = model == "3p",
      preferred_AICc = model == best_AICc, preferred_BIC = model == best_BIC,
      reference_R_V = full_dat$R_V, reference_R_A = full_dat$R_A,
      predicted_reference_inputs = reference$predicted,
      W_V_reference_inputs = reference$W_V, W_A_reference_inputs = reference$W_A,
      denominator_reference_inputs = reference$denominator,
      full_fit_prediction = full_reference$predicted,
      prediction_change = predicted - full_fit_prediction,
      fixed_input_prediction_change = predicted_reference_inputs - full_fit_prediction,
      .before = 1)
  }))
  cat("Study 2", omitted, " AICc:", best_AICc, " BIC:", best_BIC, "\n")
}
inputs <- bind_rows(inputs)
comparisons <- bind_rows(comparisons)
coefficients <- bind_rows(coefficients)
predictions <- bind_rows(predictions)
fold_counts <- bind_rows(fold_counts)

## 5. Separate changes in inputs from changes in fitted parameters ---------
# Compare recomputed-input predictions and predictions at fixed full-sample inputs.
prediction_changes <- predictions %>% group_by(omitted_participant, model) %>%
  summarise(
    prediction_rmse_change = sqrt(mean(prediction_change^2)),
    prediction_max_abs_change = max(abs(prediction_change)),
    fixed_input_rmse_change = sqrt(mean(fixed_input_prediction_change^2)),
    fixed_input_max_abs_change = max(abs(fixed_input_prediction_change)),
    min_W_V = min(W_V), max_W_V = max(W_V), min_W_A = min(W_A), max_W_A = max(W_A),
    min_abs_denominator = min(abs(denominator)),
    .groups = "drop")
comparisons <- left_join(comparisons, prediction_changes, by = c("omitted_participant", "model"))

# Summarize changes in the focal 3p predictions reported in the supplement.
focal_changes <- prediction_changes %>% filter(model == "3p") %>%
  select(omitted_participant, prediction_rmse_change, prediction_max_abs_change,
         fixed_input_rmse_change, fixed_input_max_abs_change) %>%
  rename_with(~ paste0("focal_3p_", .x), -omitted_participant)
summary <- left_join(fold_counts, focal_changes, by = "omitted_participant")
weight_diagnostics <- predictions %>% filter(model %in% names(start_params)) %>%
  select(omitted_participant, model, av_pair_id, R_V, R_A, W_V, W_A, denominator,
         reference_R_V, reference_R_A, W_V_reference_inputs, W_A_reference_inputs,
         denominator_reference_inputs) %>%
  mutate(negative_W_V = W_V < 0, negative_W_A = W_A < 0,
         singular_denominator = abs(denominator) <= 1e-10 * (1 + abs(W_V) + abs(W_A)))
write.csv(inputs, file.path(out_dir, "RecomputedInputs.csv"), row.names = FALSE)
write.csv(summary, file.path(out_dir, "Summary.csv"), row.names = FALSE)
write.csv(comparisons, file.path(out_dir, "ModelComparison.csv"), row.names = FALSE)
write.csv(coefficients, file.path(out_dir, "Parameters.csv"), row.names = FALSE)
write.csv(predictions %>% select(-W_V, -W_A, -denominator, -W_V_reference_inputs,
                                -W_A_reference_inputs, -denominator_reference_inputs),
          file.path(out_dir, "Predictions.csv"), row.names = FALSE)
write.csv(weight_diagnostics, file.path(out_dir, "WeightDiagnostics.csv"), row.names = FALSE)
cat("Participant sensitivity outputs saved in", out_dir, "\n")

## 6. Summarize model selection across participant omissions --------------
# Positive margins favor 3p over the closest alternative.
selection_margins <- bind_rows(lapply(fold_counts$omitted_participant, function(omitted) {
  x <- comparisons %>% filter(omitted_participant == omitted)
  focal <- x %>% filter(model == "3p")
  other <- x %>% filter(model != "3p")
  tibble(
    omitted_participant = omitted,
    AICc_margin = min(other$AICc) - focal$AICc,
    BIC_margin = min(other$BIC) - focal$BIC,
    other_AICc = other$model[which.min(other$AICc)],
    other_BIC = other$model[which.min(other$BIC)]
  )
}))
write.csv(selection_margins, file.path(out_dir, "SelectionMargins.csv"), row.names = FALSE)
# The shared plot helper summarizes selections across the ten omissions.
if (length(omit_participants) > 0) {
  p_selection <- plot_selection_frequency(summary)
  write.csv(p_selection$data, file.path(out_dir, "SelectionFrequency.csv"), row.names = FALSE)
  figure_file <- file.path(figure_dir, "SelectionFrequency.png")
  export_figure(figure_file, p_selection, FIGURE_WIDTH_IN, FIGURE_HEIGHT_IN)
  saveRDS(p_selection, sub("\\.png$", ".rds", figure_file))
}

export_combined_selection_frequency()
