# Study 3: prepare 54 pair-level count summaries from formal trials.
rm(list = ls())
# Load shared paths relative to this script.
.script_files <- unlist(lapply(sys.frames(), function(frame) frame$ofile))
.script_file <- if (length(.script_files)) tail(.script_files, 1) else
  sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
if (!length(.script_file)) .script_file <- if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) rstudioapi::getSourceEditorContext()$path else "Step1_PrepareData.R"
source(file.path(dirname(normalizePath(.script_file)), "..", "Project_Config.R"), encoding = "UTF-8")
setwd(file.path(PROJECT_ROOT, "Scripts", "Study3"))
library(dplyr)
dir.create("results", showWarnings = FALSE)

## 1. Read participant trial files -----------------------------------------
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
}))

## 2. Retain formal responses and identify stimuli -------------------------
# Block and trial fields distinguish formal trials from practice/instructions.
# Code disgust as 1 and anger as 0, independent of response position.
raw <- raw %>%
  filter(!is.na(block), !is.na(trial), !is.na(response),
         modality %in% c("vision", "auditory", "av")) %>%
  mutate(
    visual_label = paste(face_model_id, face_morph, sep = "_"),
    auditory_label = paste(voice_model_id, voice_morph, sep = "_"),
    VV = match(face_morph, c(20, 50, 80)) - 2L,
    AV = match(voice_morph, c(20, 50, 80)) - 2L
  )

## 3. Summarize each component stimulus across 50 responses ----------------
# Summarize components before joining them to audiovisual pairs.
visual <- raw %>% filter(modality == "vision") %>%
  group_by(face_file) %>% summarise(k_V = sum(response), n_V = n(), R_V = k_V/n_V, .groups = "drop")
auditory <- raw %>% filter(modality == "auditory") %>%
  group_by(voice_file) %>% summarise(k_A = sum(response), n_A = n(), R_A = k_A/n_A, .groups = "drop")

## 4. Summarize audiovisual pairs and match component inputs ---------------
dat <- raw %>% filter(modality == "av") %>%
  group_by(face_file, voice_file, visual_label, auditory_label, VV, AV) %>%
  summarise(VA_disgust = sum(response), n_VA = n(), .groups = "drop") %>%
  left_join(visual, by = "face_file") %>%
  left_join(auditory, by = "voice_file") %>%
  arrange(visual_label, auditory_label) %>%
  mutate(av_pair_id = sprintf("P%02d", row_number()), R_VA = VA_disgust / n_VA) %>%
  select(av_pair_id, visual_label, auditory_label, R_V, R_A, R_VA, VV, AV, VA_disgust, n_VA, k_V, n_V, k_A, n_A)
write.csv(dat, file.path(PROJECT_ROOT, "Data", "Study3", "Input_Data54.csv"), row.names = FALSE)

## 5. Condition means and SEs across pair-level summaries ------------------
cells <- dat %>% group_by(VV, AV) %>%
  summarise(observed = mean(R_VA), se = sd(R_VA)/sqrt(n()),
            n_items = n(), .groups = "drop")
write.csv(cells, "results/Observed_CellMeans.csv", row.names = FALSE)
cat("Formal responses:", nrow(raw), "; audiovisual pairs:", nrow(dat), "\n")
