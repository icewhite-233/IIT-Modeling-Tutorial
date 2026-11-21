# Study 3: omit each participant and repeat preparation and binomial fitting.
# Compare within-sample information criteria and changes in fitted predictions.
# Run Steps 1-2 first, then omit participants and compare with the full fit.

## 0. Environment and dependencies ----------------------------------------
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
setwd(file.path(PROJECT_ROOT, "Scripts", "Study3"))
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyr)
source("help_function_choice.R")
SCALE_MID <- 0
output_dir <- file.path(PROJECT_ROOT, "Results", "Study3", "ParticipantSensitivity")
figure_dir <- file.path(PROJECT_ROOT, "Figures", "Study3", "ParticipantSensitivity")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## 1. Read the same formal trials as Step 1 -------------------------------
files <- list.files(file.path(PROJECT_ROOT, "Data", "Study3", "rawdata"),
                    pattern = "_Category_.*\\.csv$", full.names = TRUE)
raw <- bind_rows(lapply(files, function(file) {
  x <- read.csv(file, colClasses = "character", na.strings = c("", "None"))
  x %>% transmute(
    participant_id = sub("_.*", "", basename(file)),
    block = as.integer(block), trial = as.integer(trial), modality,
    face_file, voice_file, face_model_id,
    face_morph = as.numeric(face_morph),
    voice_model_id = sub("\\.0$", "", voice_model_id),
    voice_morph = as.numeric(voice_morph),
    response = case_when(resp_label == "disgust" ~ 1, resp_label == "angry" ~ 0)
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

# Remove a participant before recomputing both components and joint counts.
prepare_inputs <- function(trials) {
  visual <- trials %>% filter(modality == "vision") %>%
    group_by(face_file) %>%
    summarise(k_V = sum(response), n_V = n(), R_V = k_V/n_V, .groups = "drop")
  auditory <- trials %>% filter(modality == "auditory") %>%
    group_by(voice_file) %>%
    summarise(k_A = sum(response), n_A = n(), R_A = k_A/n_A, .groups = "drop")
  trials %>% filter(modality == "av") %>%
    group_by(face_file, voice_file, visual_label, auditory_label, VV, AV) %>%
    summarise(VA_disgust = sum(response), n_VA = n(), .groups = "drop") %>%
    left_join(visual, by = "face_file") %>%
    left_join(auditory, by = "voice_file") %>%
    arrange(visual_label, auditory_label) %>%
    mutate(av_pair_id = sprintf("P%02d", row_number()), R_VA = VA_disgust/n_VA) %>%
    select(av_pair_id, visual_label, auditory_label, R_V, R_A, R_VA,
           VV, AV, VA_disgust, n_VA, k_V, n_V, k_A, n_A)
}

## 3. Fit the full-sample reference ----------------------------------------
reference_inputs <- read.csv(file.path(PROJECT_ROOT, "Data", "Study3", "Input_Data54.csv"))
prepared_full <- prepare_inputs(raw)
stopifnot(identical(prepared_full$av_pair_id, reference_inputs$av_pair_id),
          identical(prepared_full$visual_label, reference_inputs$visual_label),
          identical(prepared_full$auditory_label, reference_inputs$auditory_label))
reference_inputs <- choice_prepare_inputs(reference_inputs)

## 4. Refit after each participant omission -------------------------------
comparison_rows <- parameter_rows <- prediction_rows <- weight_rows <- input_rows <- fold_rows <- list()
for (omitted in c("None", omit_participants)) {
  trials <- raw[raw$participant_id != omitted, ]
  dat <- if (omitted == "None") reference_inputs else
    choice_prepare_inputs(prepare_inputs(trials))
  stopifnot(identical(dat$av_pair_id, reference_inputs$av_pair_id),
            identical(dat$visual_label, reference_inputs$visual_label),
            identical(dat$auditory_label, reference_inputs$auditory_label))
  # Fit each binomial model as in Step 2, using the existing likelihood helper.
  optim_3p <- optim(
    par = start_params[["3p"]], fn = choice_nll, data = dat, model = "3p",
    method = "L-BFGS-B",
    control = list(maxit = 10000, factr = 1e7, pgtol = 1e-8), hessian = TRUE
  )
  fit_3p <- choice_result(optim_3p, "3p", dat)
  optim_4p <- optim(
    par = start_params[["4p"]], fn = choice_nll, data = dat, model = "4p",
    method = "L-BFGS-B",
    control = list(maxit = 10000, factr = 1e7, pgtol = 1e-8), hessian = TRUE
  )
  fit_4p <- choice_result(optim_4p, "4p", dat)
  optim_4plq <- optim(
    par = start_params[["4plq"]], fn = choice_nll, data = dat, model = "4plq",
    method = "L-BFGS-B",
    control = list(maxit = 10000, factr = 1e7, pgtol = 1e-8), hessian = TRUE
  )
  fit_4plq <- choice_result(optim_4plq, "4plq", dat)
  optim_5p <- optim(
    par = start_params[["5p"]], fn = choice_nll, data = dat, model = "5p",
    method = "L-BFGS-B",
    control = list(maxit = 10000, factr = 1e7, pgtol = 1e-8), hessian = TRUE
  )
  fit_5p <- choice_result(optim_5p, "5p", dat)
  optim_5plq <- optim(
    par = start_params[["5plq"]], fn = choice_nll, data = dat, model = "5plq",
    method = "L-BFGS-B",
    control = list(maxit = 10000, factr = 1e7, pgtol = 1e-8), hessian = TRUE
  )
  fit_5plq <- choice_result(optim_5plq, "5plq", dat)
  optim_7p <- optim(
    par = start_params[["7p"]], fn = choice_nll, data = dat, model = "7p",
    method = "L-BFGS-B",
    control = list(maxit = 10000, factr = 1e7, pgtol = 1e-8), hessian = TRUE
  )
  fit_7p <- choice_result(optim_7p, "7p", dat)
  fits <- list(
    `3p` = fit_3p, `4p` = fit_4p,
    `4plq` = fit_4plq, `5p` = fit_5p,
    `5plq` = fit_5plq, `7p` = fit_7p
  )
  comparison <- bind_rows(lapply(fits, function(fit) {
    model_statistics(fit, dat) %>%
      mutate(converged = fit$diagnostics$converged,
             hessian_positive_definite = fit$hessian_positive_definite)
  }))
  for (model in c("Visual-only", "Auditory-only")) {
    probability <- plogis(if (model == "Visual-only") dat$L_V else dat$L_A)
    benchmark <- list(model = model, par = numeric(0),
      nll = -sum(dbinom(dat$VA_disgust, dat$n_VA, probability, log = TRUE)),
      predicted_probability = probability,
      diagnostics = data.frame(converged = logical(0)))
    comparison <- bind_rows(comparison, model_statistics(benchmark, dat) %>%
      mutate(converged = NA))
  }
  comparison <- comparison %>%
    mutate(delta_AICc = AICc - min(AICc), delta_BIC = BIC - min(BIC)) %>%
    arrange(AICc)
  fitted <- list(fits = fits, comparison = comparison)
  if (omitted == "None") full <- fitted
  best_AICc <- comparison$model[which.min(comparison$AICc)]
  best_BIC <- comparison$model[which.min(comparison$BIC)]
  comparison_rows[[omitted]] <- comparison %>% mutate(omitted_participant = omitted, .before = 1)
  input_rows[[omitted]] <- dat %>% mutate(omitted_participant = omitted, .before = 1)
  fold_rows[[omitted]] <- tibble(
    omitted_participant = omitted, n_participants = n_distinct(trials$participant_id),
    n_formal_trials = nrow(trials), n_visual_trials = sum(trials$modality == "vision"),
    n_auditory_trials = sum(trials$modality == "auditory"),
    n_av_trials = sum(trials$modality == "av"), n_pairs = nrow(dat),
    component_n_min = min(c(dat$n_V, dat$n_A)), component_n_max = max(c(dat$n_V, dat$n_A)),
    joint_n_min = min(dat$n_VA), joint_n_max = max(dat$n_VA),
    best_AICc = best_AICc, best_BIC = best_BIC, criteria_agree = best_AICc == best_BIC,
    n_failed_candidates = sum(comparison$converged == FALSE, na.rm = TRUE))

  for (model in comparison$model) {
    fit <- fitted$fits[[model]]
    if (model %in% names(start_params)) {
      predicted <- fit$predicted_probability
      predicted_reference <- predict(fit, reference_inputs)
      full_prediction <- full$fits[[model]]$predicted_probability
      weights <- choice_weight_functions[[model]](fit$par, dat)
      reference_weights <- choice_weight_functions[[model]](fit$par, reference_inputs)
      parameter_rows[[paste(omitted, model)]] <- tibble(
        omitted_participant = omitted, model = model,
        term = names(fit$par), estimate = unname(fit$par))
      weight_rows[[paste(omitted, model)]] <- choice_weight_diagnostics(fit) %>%
        mutate(omitted_participant = omitted, .before = 1)
    } else {
      component <- if (model == "Visual-only") "L_V" else "L_A"
      predicted <- plogis(dat[[component]])
      predicted_reference <- full_prediction <- plogis(reference_inputs[[component]])
      weights <- list(W_V = rep(NA_real_, nrow(dat)), W_A = rep(NA_real_, nrow(dat)))
      reference_weights <- weights
    }
    prediction_rows[[paste(omitted, model)]] <- tibble(
      omitted_participant = omitted, model = model, av_pair_id = dat$av_pair_id,
      VV = dat$VV, AV = dat$AV, L_V = dat$L_V, L_A = dat$L_A,
      observed = dat$R_VA, predicted = predicted, full_prediction = full_prediction,
      prediction_change = predicted - full_prediction,
      predicted_reference = predicted_reference,
      reference_prediction_change = predicted_reference - full_prediction,
      W_V = weights$W_V, W_A = weights$W_A,
      denominator = 1 + weights$W_V + weights$W_A,
      reference_L_V = reference_inputs$L_V, reference_L_A = reference_inputs$L_A,
      W_V_reference_inputs = reference_weights$W_V,
      W_A_reference_inputs = reference_weights$W_A,
      denominator_reference_inputs = 1 + reference_weights$W_V + reference_weights$W_A)
  }
  cat("Omitted participant:", omitted, " AICc:", best_AICc, " BIC:", best_BIC, "\n")
}

## 5. Save the full reference and all ten omissions ------------------------
# Summarize changes in the focal 7p predictions reported in the supplement.
prediction_changes <- bind_rows(prediction_rows) %>% filter(model == "7p") %>%
  group_by(omitted_participant) %>% summarise(
    prediction_rmse_change = sqrt(mean(prediction_change^2)),
    prediction_max_abs_change = max(abs(prediction_change)),
    fixed_input_rmse_change = sqrt(mean(reference_prediction_change^2)),
    fixed_input_max_abs_change = max(abs(reference_prediction_change)),
    .groups = "drop") %>%
  rename_with(~ paste0("focal_7p_", .x), -omitted_participant)
summary <- left_join(bind_rows(fold_rows), prediction_changes, by = "omitted_participant")
write.csv(bind_rows(comparison_rows), file.path(output_dir, "ModelComparison.csv"), row.names = FALSE)
write.csv(bind_rows(parameter_rows), file.path(output_dir, "Parameters.csv"), row.names = FALSE)
write.csv(bind_rows(prediction_rows), file.path(output_dir, "Predictions.csv"), row.names = FALSE)
write.csv(bind_rows(weight_rows), file.path(output_dir, "WeightDiagnostics.csv"), row.names = FALSE)
write.csv(bind_rows(input_rows), file.path(output_dir, "RecomputedInputs.csv"), row.names = FALSE)
write.csv(summary, file.path(output_dir, "Summary.csv"), row.names = FALSE)
cat("Participant sensitivity results saved to", output_dir, "\n")

## 6. Summarize model selection across participant omissions --------------
# Positive margins favor 7p; negative margins favor another candidate.
margins <- bind_rows(comparison_rows) %>%
  group_by(omitted_participant) %>%
  summarise(
    AICc_margin = min(AICc[model != "7p"]) - AICc[match("7p", model)],
    BIC_margin = min(BIC[model != "7p"]) - BIC[match("7p", model)],
    other_AICc = model[model != "7p"][which.min(AICc[model != "7p"])],
    other_BIC = model[model != "7p"][which.min(BIC[model != "7p"])],
    .groups = "drop")
write.csv(margins, file.path(output_dir, "SelectionMargins.csv"), row.names = FALSE)
# The shared plot helper summarizes selections across the ten omissions.
if (length(omit_participants) > 0) {
  p_selection <- plot_selection_frequency(summary)
  write.csv(p_selection$data, file.path(output_dir, "SelectionFrequency.csv"), row.names = FALSE)
  figure_file <- file.path(figure_dir, "SelectionFrequency.png")
  export_figure(figure_file, p_selection, FIGURE_WIDTH_IN, FIGURE_HEIGHT_IN)
  saveRDS(p_selection, sub("\\.png$", ".rds", figure_file))
}

export_combined_selection_frequency()
