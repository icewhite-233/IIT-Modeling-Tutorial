# Study 2: six integration functions and default starting values.

# Rating-scale settings. Change these two values when using another scale.
SCALE_MIN <- 0
SCALE_MAX <- 1
SCALE_MID <- (SCALE_MIN + SCALE_MAX) / 2
SCALE_BREAKS <- seq(SCALE_MIN, SCALE_MAX, length.out = 9)

## A. Weighted-averaging family (wavg_*): ----------------------------------
# Map component ratings R_V/R_A through weighting functions W_V/W_A and integrate.

# Infer component scale values from unimodal ratings; integrate with W0 = 1.

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

## B. Default starting values for the six candidates -----------------------
# Order: 3p, 4p, 4plq, 5p, 5plq, 7p

start_params <- list(
  `3p`   = list(S0 = 0.5, W_V = 1, W_A = 1),
  `4p`   = list(S0 = 0.5, A_V = 1, A_A = 1, B_VA = 0),
  `4plq` = list(S0 = 0.5, A_VA = 1, B_VA = 0, B_VA2 = 0),
  `5p`   = list(S0 = 0.5, A_V = 1, B_V = 0, A_A = 1, B_A = 0),
  `5plq` = list(S0 = 0.5, A_V = 1, A_A = 1,
                B_VA = 0, B_VA2 = 0),
  `7p`   = list(S0 = 0.5, A_V = 1, B_V = 0, B_V2 = 0,
                A_A = 1, B_A = 0, B_A2 = 0)
)
