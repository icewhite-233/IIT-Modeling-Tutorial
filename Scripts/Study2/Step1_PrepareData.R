# Study 2: prepare 54 pair-level inputs from recorded formal trials.
rm(list = ls())
# Load shared paths relative to this script.
.script_files <- unlist(lapply(sys.frames(), function(frame) frame$ofile))
.script_file <- if (length(.script_files)) tail(.script_files, 1) else
  sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
if (!length(.script_file)) .script_file <- if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) rstudioapi::getSourceEditorContext()$path else "Step1_PrepareData.R"
source(file.path(dirname(normalizePath(.script_file)), "..", "Project_Config.R"), encoding = "UTF-8")
setwd(file.path(PROJECT_ROOT, "Scripts", "Study2"))
library(dplyr)
dir.create("results", showWarnings = FALSE)

## 1. Read participant trial files -----------------------------------------
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
    voice_morph = as.numeric(voice_morph),
    response = as.numeric(resp)
  )
}))

## 2. Retain formal responses and identify stimuli -------------------------
# Block and trial fields distinguish formal trials from practice/instructions.
# Responses are already on the 0-1 disgust-rating scale.
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
  group_by(face_file) %>% summarise(R_V = mean(response), .groups = "drop")
auditory <- raw %>% filter(modality == "auditory") %>%
  group_by(voice_file) %>% summarise(R_A = mean(response), .groups = "drop")

## 4. Summarize audiovisual pairs and match component inputs ---------------
dat <- raw %>% filter(modality == "av") %>%
  group_by(face_file, voice_file, visual_label, auditory_label, VV, AV) %>%
  summarise(R_VA = mean(response), .groups = "drop") %>%
  left_join(visual, by = "face_file") %>%
  left_join(auditory, by = "voice_file") %>%
  arrange(visual_label, auditory_label) %>%
  mutate(av_pair_id = sprintf("P%02d", row_number())) %>%
  select(av_pair_id, visual_label, auditory_label, R_V, R_A, R_VA, VV, AV)
write.csv(dat, file.path(PROJECT_ROOT, "Data", "Study2", "Input_Data54.csv"), row.names = FALSE)

## 5. Condition means and SEs across pair-level summaries ------------------
cells <- dat %>% group_by(VV, AV) %>%
  summarise(observed = mean(R_VA), se = sd(R_VA)/sqrt(n()),
            n_items = n(), .groups = "drop")
write.csv(cells, "results/Observed_CellMeans.csv", row.names = FALSE)
cat("原始正式试次：", nrow(raw), "；视听配对：", nrow(dat), "\n")
