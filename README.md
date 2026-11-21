# IIT tutorial: data and R code

Data and R code for the Information Integration Theory tutorial.

## Setup

The scripts were checked with R 4.5.0 on 64-bit Windows. RStudio is optional. Install the required packages once:

```r
install.packages(c(
  "readxl", "dplyr", "ggplot2", "car", "minpack.lm", "broom",
  "tidyr", "purrr", "tibble", "RColorBrewer", "rlang", "scales",
  "ggbreak", "svglite", "cowplot"
))
```

For RStudio's script-location support, also install `rstudioapi`. Figure exports use Arial and require R with Cairo support; check this with `capabilities("cairo")`. Tested package versions are listed in [Environment.txt](Environment.txt).

## Data

| Study | Analysis input | Contents |
| --- | --- | --- |
| 1 | `Data/Study1/MODELING_Exp1_Valence.xlsx` | Ratings for 180 video–music pairs |
| 2 | `Data/Study2/Input_Data54.csv` | Prepared continuous-rating data for 54 face–voice pairs |
| 3 | `Data/Study3/Input_Data54.csv` | Prepared binary-choice counts for 54 face–voice pairs |

De-identified trial-level files for Studies 2–3 are in each study's `Data/.../rawdata` folder. These contain stimulus identifiers, trial information, and responses. See [Data/README.md](Data/README.md) for variable definitions.

## Run the examples

Keep the folder structure intact. In RStudio, open and source each complete script. Alternatively, open a terminal in the repository root and run, for example:

```sh
Rscript Scripts/Study2/Step2_ModelFitting.R
```

Run the numbered scripts in order within each study:

- **Study 1:** Steps 1–7 cover simulation, slope attenuation, fitting, model recovery, parameter recovery, validation, and weight plots.
- **Studies 2–3:** Steps 2–6 cover fitting, validation, result plots, weight plots, and participant-deletion analyses. Step 1 regenerates the supplied analysis input from participant-level files.

Scripts load their helpers and `Scripts/Project_Config.R` automatically. This configuration locates the project and defines shared figure sizes, fonts, colours, and export settings. No local path needs to be entered.

Numerical results and individual figures are saved in `Scripts/Study*/results` and `Scripts/Study*/figures`. Participant-deletion outputs are saved under `Results/Study*` and `Figures/Study*`, with a combined comparison under `Figures/ParticipantSensitivity`.

## Figures and simulation settings

After the numbered scripts finish, run `Scripts/Export_ManuscriptPanels.R` for the panels in Figures 3–7. `Scripts/Study1/Plot_ModelSpace.R` generates Figure 2. R produces the individual plots; final multi-panel layouts were assembled in Visio. Visio is not required for the analyses.

Random seeds and simulation settings are specified in the scripts. Main analyses use noise SD = 0.39. The SD = 0.20 alternative for supplementary recovery results is marked in Study 1 Steps 4–5.

The `manuscript_examples_study2.R` and `manuscript_examples_study3.R` files provide condensed examples; run them with the repository root as the working directory.

## Access

Experimental procedures and stimulus sources are described in the accompanying article and supplementary material. Experiment-collection scripts and stimulus media are not included. The relative stimulus paths in the trial files identify stimuli; the analysis scripts read the recorded responses without loading these media.

## License

Except where otherwise noted, the code, data, and documentation in this repository are licensed under the [Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/). Please credit the original work and indicate any changes. See [LICENSE](LICENSE) for the full terms.

Third-party software and materials retain their original licenses.
