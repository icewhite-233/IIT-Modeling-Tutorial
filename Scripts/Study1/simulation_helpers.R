# Continuous integration models and simulation/fitting helpers.
# Simulation adds noise to model predictions; recovery refits the simulated data.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(minpack.lm)
  library(broom)
  library(purrr)
  library(tibble)
})

# Rating-scale settings. Change these two values when using another scale.
SCALE_MIN <- 1
SCALE_MAX <- 9
SCALE_MID <- (SCALE_MIN + SCALE_MAX) / 2
SCALE_BREAKS <- seq(SCALE_MIN, SCALE_MAX, length.out = 9)

## A. Weighted-averaging family (wavg_*): ----------------------------------
# Map component ratings R_V/R_A through weighting functions W_V/W_A and integrate.

# Apply the integration rule with the initial-information weight fixed at W0 = 1.

# 3p: constant-weight averaging with separate visual and auditory weights.

# Algebraically reduced prediction; W0 remains fixed at 1.
# A signed weight is allowed by the existing examples; a singular sum is not.
integrate_components <- function(S0, W_V, W_A, R_V, R_A,
                                 denominator_tol = 1e-10) {
  denominator <- 1 + W_V + W_A
  threshold <- denominator_tol * (1 + abs(W_V) + abs(W_A))
  invalid <- !is.finite(denominator) | abs(denominator) <= threshold
  R_VA_hat <- (R_V * (1 + W_V) + R_A * (1 + W_A) - S0) /
    denominator
  R_VA_hat[rep_len(invalid, length(R_VA_hat))] <- NA_real_
  R_VA_hat
}

wavg_3 <- function(S0, W_V, W_A, R_V, R_A) {
  integrate_components(S0, W_V, W_A, R_V, R_A)
}

# 4p: separate intercepts and a shared linear slope.
# W_V = A_V + B_VA*(R_V - SCALE_MID); W_A = A_A + B_VA*(R_A - SCALE_MID)
wavg_4 <- function(S0, A_V, A_A, B_VA, R_V, R_A) {
  W_V <- A_V + B_VA*(R_V - SCALE_MID)
  W_A <- A_A + B_VA*(R_A - SCALE_MID)
  integrate_components(S0, W_V, W_A, R_V, R_A)
}

# 4plq: shared intercept, slope, and curvature parameter.
wavg_4plq <- function(S0, A_VA, B_VA, B_VA2, R_V, R_A) {
  W_V <- A_VA + B_VA*(R_V - SCALE_MID) + B_VA2*(R_V - SCALE_MID)^2
  W_A <- A_VA + B_VA*(R_A - SCALE_MID) + B_VA2*(R_A - SCALE_MID)^2
  integrate_components(S0, W_V, W_A, R_V, R_A)
}

# 5p: separate intercepts and linear slopes.
# W_V = A_V + B_V*(R_V - SCALE_MID); W_A = A_A + B_A*(R_A - SCALE_MID)
wavg_5 <- function(S0, A_V, B_V, A_A, B_A, R_V, R_A) {
  W_V <- A_V + B_V*(R_V - SCALE_MID)
  W_A <- A_A + B_A*(R_A - SCALE_MID)
  integrate_components(S0, W_V, W_A, R_V, R_A)
}

# 5plq: separate intercepts, shared slope and curvature parameter.
wavg_5plq <- function(S0, A_V, A_A, B_VA, B_VA2, R_V, R_A) {
  W_V <- A_V + B_VA*(R_V - SCALE_MID) + B_VA2*(R_V - SCALE_MID)^2
  W_A <- A_A + B_VA*(R_A - SCALE_MID) + B_VA2*(R_A - SCALE_MID)^2
  integrate_components(S0, W_V, W_A, R_V, R_A)
}

# 7p: separate intercepts, slopes, and curvature parameters.
wavg_7 <- function(S0, A_V, B_V, B_V2, A_A, B_A, B_A2, R_V, R_A) {
  W_V <- A_V + B_V*(R_V - SCALE_MID) + B_V2*(R_V - SCALE_MID)^2
  W_A <- A_A + B_A*(R_A - SCALE_MID) + B_A2*(R_A - SCALE_MID)^2
  integrate_components(S0, W_V, W_A, R_V, R_A)
}

## B. Reference locations for simulation and initial fitting ------------------
# Order: 3p, 4p, 4plq, 5p, 5plq, 7p

reference_params <- list(
  `3p`   = list(S0 = 6.6673, W_V = 2.109582, W_A = 1.552700),
  `4p`   = list(S0 = 5.770255, A_V = 1.73551, A_A = 1.32501, B_VA = -0.255595),
  `4plq` = list(S0 = 5.782606, A_VA = 1.47692, B_VA = -0.37464, B_VA2 = 0.129709),
  `5p`   = list(S0 = 5.77261, A_V = 1.74375, B_V = -0.26429, A_A = 1.32156, B_A = -0.24859),
  `5plq` = list(S0 = 5.780349, A_V = 1.64939, A_A = 1.17742,
                B_VA = -0.35569, B_VA2 = 0.127104),
  `7p`   = list(S0 = 5.83261, A_V = 1.68313, B_V = -0.69688, B_V2 = 0.21354,
                A_A = 1.25840, B_A = -0.12651, B_A2 = 0.082705)
)

start_params <- reference_params

## C. Dictionary of generative functions: true parameters + (R_V, R_A) → μ_R_VA ---

generator_dict <- list(
  `3p`   = function(par, R_V, R_A) wavg_3(par$S0, par$W_V, par$W_A, R_V, R_A),
  `4p`   = function(par, R_V, R_A) wavg_4(par$S0, par$A_V, par$A_A, par$B_VA, R_V, R_A),
  `4plq` = function(par, R_V, R_A) wavg_4plq(par$S0, par$A_VA, par$B_VA, par$B_VA2, R_V, R_A),
  `5p`   = function(par, R_V, R_A) wavg_5(par$S0, par$A_V, par$B_V, par$A_A, par$B_A, R_V, R_A),
  `5plq` = function(par, R_V, R_A) wavg_5plq(par$S0, par$A_V, par$A_A,
                                            par$B_VA, par$B_VA2, R_V, R_A),
  `7p`   = function(par, R_V, R_A) wavg_7(par$S0, par$A_V, par$B_V, par$B_V2,
                                          par$A_A, par$B_A, par$B_A2, R_V, R_A)
)

## D. Jitter starting values -----------------------------------------------

# abs_min is a minimum coefficient scale, before multiplying by jitter_scale.
jitter_list <- function(par_list, jitter_scale = 0.10, abs_min = 0.05) {
  stopifnot(jitter_scale >= 0, abs_min > 0)
  lapply(par_list, function(v) {
    if (!is.numeric(v)) return(v)
    jitter_sd <- jitter_scale * max(abs(v), abs_min)
    as.numeric(v) + rnorm(1, mean = 0, sd = jitter_sd)
  })
}

## E. Wrapper for fitting a single nlsLM model -----------------------------

add_nls <- function(name, formula, df, start_list,
                    jitter_scale = 0.10, abs_min = 0.05) {
  st <- jitter_list(start_list, jitter_scale = jitter_scale, abs_min = abs_min)
  fit <- tryCatch(nlsLM(formula, data = df, start = st,
                        control = nls.lm.control(maxiter = 1024)),
                  error = function(e) e)
  ok <- !inherits(fit, "error") && isTRUE(fit$convInfo$isConv)
  if (!ok) return(tibble(
    model = name, converged = FALSE, k = length(start_list),
    rss = Inf, AIC = Inf, BIC = Inf, params = list(NULL)))
  tibble(
    model = name, converged = TRUE, k = length(coef(fit)),
    rss = sum(residuals(fit)^2), AIC = AIC(fit), BIC = BIC(fit),
    params = list(as.list(coef(fit))))
}

## F. Fit all candidate models to one simulated dataset --------------------
# Order: 3p, 4p, 4plq, 5p, 5plq, 7p

fit_all_models <- function(df,
                           start_params = reference_params,
                           jitter_scale = 0.10, abs_min = 0.05) {
  stopifnot(all(c("R_V","R_A","R_VA_sim") %in% names(df)))

  out_3p   <- add_nls("3p",
                      R_VA_sim ~ wavg_3(S0, W_V, W_A, R_V, R_A),
                      df, start_params[["3p"]], jitter_scale, abs_min)
  out_4p   <- add_nls("4p",
                      R_VA_sim ~ wavg_4(S0, A_V, A_A, B_VA, R_V, R_A),
                      df, start_params[["4p"]], jitter_scale, abs_min)
  out_4plq <- add_nls("4plq",
                      R_VA_sim ~ wavg_4plq(S0, A_VA, B_VA, B_VA2, R_V, R_A),
                      df, start_params[["4plq"]], jitter_scale, abs_min)
  out_5p   <- add_nls("5p",
                      R_VA_sim ~ wavg_5(S0, A_V, B_V, A_A, B_A, R_V, R_A),
                      df, start_params[["5p"]], jitter_scale, abs_min)
  out_5plq <- add_nls("5plq",
                      R_VA_sim ~ wavg_5plq(S0, A_V, A_A, B_VA, B_VA2, R_V, R_A),
                      df, start_params[["5plq"]], jitter_scale, abs_min)
  out_7p   <- add_nls("7p",
                      R_VA_sim ~ wavg_7(S0, A_V, B_V, B_V2, A_A, B_A, B_A2, R_V, R_A),
                      df, start_params[["7p"]], jitter_scale, abs_min)

  dplyr::bind_rows(out_3p, out_4p, out_4plq, out_5p, out_5plq, out_7p)
}

## G. Generate one simulated dataset: true model + design + noise → R_VA_sim ---

simulate_once <- function(true_model,
                          design_df,
                          truth_params = reference_params[[true_model]],
                          noise_sd = 0.39,
                          clip_to_scale = TRUE,
                          lower = SCALE_MIN, upper = SCALE_MAX) {
  stopifnot(all(c("R_V","R_A") %in% names(design_df)))
  mu <- generator_dict[[true_model]](truth_params, design_df$R_V, design_df$R_A)
  y  <- mu + rnorm(length(mu), 0, noise_sd)
  if (clip_to_scale) y <- pmin(pmax(y, lower), upper)
  dplyr::mutate(design_df, R_VA_sim = y)
}

## H. Sample random parameter sets for model recovery ----------------------
# Draw each coefficient uniformly around its model's reference value.

sample_uniform_params <- function(model,
                                  reference_params,
                                  width_factor = 0.5,
                                  abs_min = 0.2) {
  # width_factor: how wide the uniform range is (relative to |par|)
  # abs_min: coefficient-scale floor, before multiplying by width_factor.

  par0 <- reference_params[[model]]

  if (is.null(par0)) {
    stop("sample_uniform_params: model '", model, "' not found in reference_params.")
  }

  par_rand <- lapply(par0, function(v) {
    if (!is.numeric(v)) return(v)
    spread <- width_factor * max(abs(v), abs_min)
    lower  <- v - spread
    upper  <- v + spread
    runif(1, min = lower, max = upper)
  })

  names(par_rand) <- names(par0)
  par_rand
}
