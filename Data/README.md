# Data dictionary

Study 1 contains 180 pair-level records. Studies 2–3 each include 54 pair-level records and 10 de-identified participant files, with 450 formal trials per participant (270 audiovisual, 90 visual, and 90 auditory).

## Trial-level files: Studies 2–3

Files in `Study2/rawdata` and `Study3/rawdata` retain the names expected by the scripts. Here, `rawdata` refers to de-identified trial-level tables, rather than complete experiment logs. Participant codes P01–P10 are read from the filenames; codes are scoped to each study.

| Variable | Meaning |
| --- | --- |
| `block`, `trial` | Block and trial indices |
| `modality` | `vision`, `auditory`, or `av` |
| `face_file`, `voice_file` | Relative stimulus filenames; stimulus media are not included |
| `face_model_id`, `voice_model_id` | Stimulus actor identifiers, not participant identifiers |
| `face_morph`, `voice_morph` | Disgust morph level: 20, 50, or 80 (%) |
| `resp` (Study 2) | Continuous disgust rating, from 0 to 1 |
| `resp_label` (Study 3) | `angry` or `disgust`; analysis codes these as 0 and 1 |

`None` indicates an unavailable or inapplicable field, such as a face identifier on an auditory-only trial. The shared tables contain no participant names, contact details, demographic fields, or collection timestamps.

## Pair-level files

| Variable | Meaning |
| --- | --- |
| `R_V`, `R_A`, `R_VA` | Mean visual, auditory, and audiovisual responses; Study 1 uses 1–9 valence ratings, Study 2 uses 0–1 disgust ratings, and Study 3 uses disgust-choice proportions |
| `VV`, `AV` | Visual and auditory conditions: −1/0/1 for negative/neutral/positive valence in Study 1, and anger/mixed/disgust in Studies 2–3 |
| `Visual$`, `Auditory$` (Study 1) | Component stimulus identifiers |
| `av_pair_id` (Studies 2–3) | Pair identifier; P01–P54 denotes pairs, not participants |
| `visual_label`, `auditory_label` (Studies 2–3) | Component identifiers combining actor and morph level |
| `VA_disgust`, `n_VA` (Study 3) | Audiovisual disgust-choice count and total response count |
| `k_V`, `n_V`, `k_A`, `n_A` (Study 3) | Component disgust-choice counts and total response counts |

Study 2–3 Step 1 scripts regenerate the pair-level files from the trial-level tables. Step 6 scripts recompute these summaries after omitting each participant.

## License

These data are shared under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). See the repository [LICENSE](../LICENSE).
