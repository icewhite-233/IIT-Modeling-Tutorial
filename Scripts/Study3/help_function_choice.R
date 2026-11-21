# Study 3: integration predictors, binomial likelihood and fitted-result methods.
# The calling script sets SCALE_MID and loads dplyr and tibble.

# Half-count correction of component proportions using their denominators.
# Audiovisual counts remain unchanged in the binomial likelihood.
choice_prob_to_logodds <- function(p, n) {
  count <- p * n
  if (any(abs(count - round(count)) > 1e-8)) {
    warning("Some proportions are not integer multiples of 1/n; check the denominator.")
  }
  qlogis((count + 0.5) / (n + 1))
}

# Count-to-log-odds mapping for descriptive residuals and plotting limits.
choice_count_to_logodds <- function(k, n) {
  if (any(!is.finite(k) | !is.finite(n) | n <= 0 | k < 0 | k > n |
          abs(k-round(k)) > 1e-8 | abs(n-round(n)) > 1e-8))
    stop("Component counts must be integers between zero and their positive denominators.")
  qlogis((k + 0.5) / (n + 1))
}

choice_prepare_inputs <- function(data) {
  data$L_V <- choice_prob_to_logodds(data$R_V, data$n_V)
  data$L_A <- choice_prob_to_logodds(data$R_A, data$n_A)
  data
}

## A. Six candidate weighted-averaging models ------------------------------

# Recover component scale values, then apply weighted averaging with W0 = 1.
integrate_components <- function(S0, W_V, W_A, R_V, R_A,
                                 denominator_tol = 1e-10) {
  denominator <- 1 + W_V + W_A
  invalid <- !is.finite(denominator) |
    abs(denominator) <= denominator_tol*(1 + abs(W_V) + abs(W_A))
  prediction <- (1 * S0 +
    W_V * ((R_V * (1 + W_V) - 1 * S0) / W_V) +
    W_A * ((R_A * (1 + W_A) - 1 * S0) / W_A)) /
    (1 + W_V + W_A)
  prediction[rep_len(invalid, length(prediction))] <- NA_real_
  prediction
}

# 3p: constant visual and auditory weights.
wavg_3 <- function(S0, W_V, W_A, R_V, R_A) {
  integrate_components(S0, W_V, W_A, R_V, R_A)
}

# 4p: separate intercepts and a shared linear slope.
wavg_4 <- function(S0, A_V, A_A, B_VA, R_V, R_A) {
  W_V <- A_V + B_VA * (R_V - SCALE_MID)
  W_A <- A_A + B_VA * (R_A - SCALE_MID)
  integrate_components(S0, W_V, W_A, R_V, R_A)
}

# 4plq: shared intercept, slope, and curvature parameter.
wavg_4plq <- function(S0, A_VA, B_VA, B_VA2, R_V, R_A) {
  W_V <- A_VA + B_VA * (R_V - SCALE_MID) +
    B_VA2 * (R_V - SCALE_MID)^2
  W_A <- A_VA + B_VA * (R_A - SCALE_MID) +
    B_VA2 * (R_A - SCALE_MID)^2
  integrate_components(S0, W_V, W_A, R_V, R_A)
}

# 5p: separate intercepts and linear slopes.
wavg_5 <- function(S0, A_V, B_V, A_A, B_A, R_V, R_A) {
  W_V <- A_V + B_V * (R_V - SCALE_MID)
  W_A <- A_A + B_A * (R_A - SCALE_MID)
  integrate_components(S0, W_V, W_A, R_V, R_A)
}

# 5plq: separate intercepts, shared slope and curvature parameter.
wavg_5plq <- function(S0, A_V, A_A, B_VA, B_VA2, R_V, R_A) {
  W_V <- A_V + B_VA * (R_V - SCALE_MID) +
    B_VA2 * (R_V - SCALE_MID)^2
  W_A <- A_A + B_VA * (R_A - SCALE_MID) +
    B_VA2 * (R_A - SCALE_MID)^2
  integrate_components(S0, W_V, W_A, R_V, R_A)
}

# 7p: separate intercepts, slopes, and curvature parameters.
wavg_7 <- function(S0, A_V, B_V, B_V2, A_A, B_A, B_A2, R_V, R_A) {
  W_V <- A_V + B_V * (R_V - SCALE_MID) +
    B_V2 * (R_V - SCALE_MID)^2
  W_A <- A_A + B_A * (R_A - SCALE_MID) +
    B_A2 * (R_A - SCALE_MID)^2
  integrate_components(S0, W_V, W_A, R_V, R_A)
}

## B. Model-specific latent predictions and explicit weighting functions ------

choice_eta_functions <- list(
  `3p` = function(par, data) {
    with(data, wavg_3(
      par[["S0"]], par[["W_V"]], par[["W_A"]], L_V, L_A
    ))
  },
  `4p` = function(par, data) {
    with(data, wavg_4(
      par[["S0"]], par[["A_V"]], par[["A_A"]], par[["B_VA"]],
      L_V, L_A
    ))
  },
  `4plq` = function(par, data) {
    with(data, wavg_4plq(
      par[["S0"]], par[["A_VA"]], par[["B_VA"]], par[["B_VA2"]],
      L_V, L_A
    ))
  },
  `5p` = function(par, data) {
    with(data, wavg_5(
      par[["S0"]], par[["A_V"]], par[["B_V"]],
      par[["A_A"]], par[["B_A"]], L_V, L_A
    ))
  },
  `5plq` = function(par, data) {
    with(data, wavg_5plq(
      par[["S0"]], par[["A_V"]], par[["A_A"]],
      par[["B_VA"]], par[["B_VA2"]], L_V, L_A
    ))
  },
  `7p` = function(par, data) {
    with(data, wavg_7(
      par[["S0"]], par[["A_V"]], par[["B_V"]], par[["B_V2"]],
      par[["A_A"]], par[["B_A"]], par[["B_A2"]], L_V, L_A
    ))
  }
)

choice_weight_functions <- list(
  `3p` = function(par, data) {
    list(
      W_V = rep(par[["W_V"]], nrow(data)),
      W_A = rep(par[["W_A"]], nrow(data))
    )
  },
  `4p` = function(par, data) {
    list(
      W_V = par[["A_V"]] + par[["B_VA"]] * (data$L_V - SCALE_MID),
      W_A = par[["A_A"]] + par[["B_VA"]] * (data$L_A - SCALE_MID)
    )
  },
  `4plq` = function(par, data) {
    list(
      W_V = par[["A_VA"]] + par[["B_VA"]] * (data$L_V - SCALE_MID) +
        par[["B_VA2"]] * (data$L_V - SCALE_MID)^2,
      W_A = par[["A_VA"]] + par[["B_VA"]] * (data$L_A - SCALE_MID) +
        par[["B_VA2"]] * (data$L_A - SCALE_MID)^2
    )
  },
  `5p` = function(par, data) {
    list(
      W_V = par[["A_V"]] + par[["B_V"]] * (data$L_V - SCALE_MID),
      W_A = par[["A_A"]] + par[["B_A"]] * (data$L_A - SCALE_MID)
    )
  },
  `5plq` = function(par, data) {
    list(
      W_V = par[["A_V"]] + par[["B_VA"]] * (data$L_V - SCALE_MID) +
        par[["B_VA2"]] * (data$L_V - SCALE_MID)^2,
      W_A = par[["A_A"]] + par[["B_VA"]] * (data$L_A - SCALE_MID) +
        par[["B_VA2"]] * (data$L_A - SCALE_MID)^2
    )
  },
  `7p` = function(par, data) {
    list(
      W_V = par[["A_V"]] + par[["B_V"]] * (data$L_V - SCALE_MID) +
        par[["B_V2"]] * (data$L_V - SCALE_MID)^2,
      W_A = par[["A_A"]] + par[["B_A"]] * (data$L_A - SCALE_MID) +
        par[["B_A2"]] * (data$L_A - SCALE_MID)^2
    )
  }
)

## C. Grouped-binomial maximum-likelihood fitting with optim() -------------

# Binomial negative log-likelihood passed to optim for the specified model.
choice_nll <- function(par, data, model) {
  eta <- choice_eta_functions[[model]](par, data)
  if (any(!is.finite(eta))) return(1e6)

  k <- data$VA_disgust
  n <- data$n_VA

  sum(
    (n - k) * pmax(eta, 0) +
      k * pmax(-eta, 0) +
      n * log1p(exp(-abs(eta))) -
      lchoose(n, k)
  )
}

# Store predictions and Hessian-based covariance without rerunning optimization.
# Diagnose the completed optimization before computing Wald inference.
choice_optimizer_diagnostics <- function(best, model, data,
                                         lower = rep(-Inf, length(best$par)),
                                         upper = rep(Inf, length(best$par)),
                                         boundary_tol = 1e-4) {
  eta <- choice_eta_functions[[model]](best$par, data)
  predicted_probability <- plogis(eta)
  hessian <- best$hessian
  hessian_finite <- is.matrix(hessian) &&
    all(dim(hessian) == length(best$par)) && all(is.finite(hessian))
  eigenvalues <- if (hessian_finite)
    eigen((hessian + t(hessian))/2, symmetric = TRUE, only.values = TRUE)$values else numeric()
  positive_definite <- length(eigenvalues) == length(best$par) && all(eigenvalues > 0)
  at_bound <- abs(best$par-lower) < boundary_tol | abs(best$par-upper) < boundary_tol
  tibble::tibble(
    model = model, start_id = "fixed", valid = is.finite(best$value) && best$value < 1e6,
    converged = best$convergence == 0L, nll = best$value,
    convergence_code = best$convergence, at_optimizer_bound = any(at_bound),
    boundary_parameters = paste(names(best$par)[at_bound], collapse = ";"),
    hessian_finite = hessian_finite,
    hessian_positive_definite = positive_definite,
    probabilities_valid = all(is.finite(eta)) &&
      all(is.finite(predicted_probability) & predicted_probability >= 0 & predicted_probability <= 1),
    selected = TRUE)
}

choice_result <- function(best, model, data,
                          lower = rep(-Inf, length(best$par)),
                          upper = rep(Inf, length(best$par))) {
  diagnostics <- choice_optimizer_diagnostics(best, model, data, lower, upper)
  if (!diagnostics$valid || !diagnostics$probabilities_valid) stop("Invalid fit: ", model)
  if (!diagnostics$converged) warning("Not converged: ", model)
  if (diagnostics$at_optimizer_bound) warning("Boundary estimate: ", model)
  eta <- choice_eta_functions[[model]](best$par, data)
  weights <- choice_weight_functions[[model]](best$par, data)
  hessian <- if (diagnostics$hessian_finite) (best$hessian+t(best$hessian))/2 else NULL
  # Wald inference is withheld at a boundary or after an unsuccessful fit.
  inference_ok <- diagnostics$converged && !diagnostics$at_optimizer_bound &&
    diagnostics$hessian_positive_definite
  covariance_matrix <- if (inference_ok)
    tryCatch(solve(hessian), error = function(e) NULL) else NULL
  if (!is.null(covariance_matrix))
    dimnames(covariance_matrix) <- list(names(best$par), names(best$par))
  eigenvalues <- if (diagnostics$hessian_positive_definite)
    eigen(hessian, symmetric = TRUE, only.values = TRUE)$values else numeric()
  output <- list(
    model = model, par = best$par, nll = best$value, eta = eta,
    predicted_probability = plogis(eta), W_V = weights$W_V, W_A = weights$W_A,
    denominator = 1 + weights$W_V + weights$W_A, diagnostics = diagnostics,
    hessian = hessian, covariance_matrix = covariance_matrix,
    hessian_positive_definite = diagnostics$hessian_positive_definite,
    hessian_condition_number = if (length(eigenvalues)) max(eigenvalues)/min(eigenvalues) else NA_real_,
    lower = lower, upper = upper, observed_probability = data$R_VA)
  class(output) <- "choice_iit_fit"
  output
}

## D. Wald parameter table and weight diagnostics --------------------------

choice_wald_parameter_table <- function(fit) {
  estimate <- as.numeric(fit$par)
  term <- names(fit$par)

  if (!is.null(fit$covariance_matrix)) {
    standard_error <- sqrt(diag(fit$covariance_matrix))
    z_value <- estimate / standard_error
    p_value <- 2 * pnorm(abs(z_value), lower.tail = FALSE)
    confidence_low <- estimate - qnorm(.975) * standard_error
    confidence_high <- estimate + qnorm(.975) * standard_error
  } else {
    standard_error <- rep(NA_real_, length(estimate))
    z_value <- rep(NA_real_, length(estimate))
    p_value <- rep(NA_real_, length(estimate))
    confidence_low <- rep(NA_real_, length(estimate))
    confidence_high <- rep(NA_real_, length(estimate))
  }

  tibble::tibble(
    model = fit$model,
    term = term,
    estimate = estimate,
    standard_error = standard_error,
    z_value = z_value,
    p_value = p_value,
    confidence_low_95 = confidence_low,
    confidence_high_95 = confidence_high,
    at_optimizer_bound = abs(estimate - fit$lower) < 1e-4 | abs(estimate - fit$upper) < 1e-4,
    hessian_positive_definite = fit$hessian_positive_definite,
    hessian_condition_number = fit$hessian_condition_number
  )
}

choice_weight_diagnostics <- function(fit) {
  tibble::tibble(
    model = fit$model,
    W_V_min = min(fit$W_V),
    W_V_max = max(fit$W_V),
    W_A_min = min(fit$W_A),
    W_A_max = max(fit$W_A),
    n_W_V_nonpositive = sum(fit$W_V <= 0),
    n_W_A_nonpositive = sum(fit$W_A <= 0),
    denominator_min = min(fit$denominator),
    denominator_max = max(fit$denominator),
    predicted_probability_min = min(fit$predicted_probability),
    predicted_probability_max = max(fit$predicted_probability)
  )
}

## E. Familiar model methods used by Step 2 and Step 3 ---------------------

# These methods let the forced-choice fit objects be used in the same style as
# the original nlsLM objects: coef(), summary(), broom::tidy(), predict(), and
# residuals(). The underlying estimation remains optim() binomial likelihood.

coef.choice_iit_fit <- function(object, ...) {
  object$par
}

predict.choice_iit_fit <- function(object, newdata = NULL,
                                   type = c("probability", "link"), ...) {
  type <- match.arg(type)

  if (is.null(newdata)) {
    eta <- object$eta
  } else {
    eta_function <- choice_eta_functions[[object$model]]
    eta <- eta_function(object$par, newdata)
  }

  if (type == "link") eta else plogis(eta)
}

residuals.choice_iit_fit <- function(object, ...) {
  object$observed_probability - object$predicted_probability
}

summary.choice_iit_fit <- function(object, ...) {
  out <- choice_wald_parameter_table(object)
  cat("\nModel:", object$model, "\n")
  cat("Grouped-binomial negative log-likelihood:", object$nll, "\n")
  print(out)
  invisible(out)
}

tidy.choice_iit_fit <- function(x, conf.int = TRUE, conf.level = 0.95, ...) {
  out <- choice_wald_parameter_table(x)

  ans <- tibble::tibble(
    term = out$term,
    estimate = out$estimate,
    std.error = out$standard_error,
    statistic = out$z_value,
    p.value = out$p_value
  )

  if (isTRUE(conf.int)) {
    alpha <- 1 - conf.level
    critical_value <- qnorm(1 - alpha / 2)
    ans$conf.low <- ans$estimate - critical_value * ans$std.error
    ans$conf.high <- ans$estimate + critical_value * ans$std.error
  }

  ans
}

## F. Starting values, comparison statistics and result export -------------
start_params <- list(
  `3p` = c(S0 = 0, W_V = 1, W_A = 1),
  `4p` = c(S0 = 0, A_V = 1, A_A = 1, B_VA = 0),
  `4plq` = c(S0 = 0, A_VA = 1, B_VA = 0, B_VA2 = 0),
  `5p` = c(S0 = 0, A_V = 1, B_V = 0, A_A = 1, B_A = 0),
  `5plq` = c(S0 = 0, A_V = 1, A_A = 1, B_VA = 0, B_VA2 = 0),
  `7p` = c(
    S0 = 0, A_V = 1, B_V = 0, B_V2 = 0,
    A_A = 1, B_A = 0, B_A2 = 0
  )
)

model_statistics <- function(fit, dat) {
n_rows <- nrow(dat)
n_trials <- sum(dat$n_VA)
null_probability <- sum(dat$VA_disgust) / n_trials
null_loglik <- sum(dbinom(
  dat$VA_disgust, dat$n_VA, null_probability,
  log = TRUE
))

  k <- length(fit$par)
  loglik <- -fit$nll
  p_hat <- fit$predicted_probability
  aic <- -2 * loglik + 2 * k
  pearson <- sum(
    (dat$VA_disgust - dat$n_VA * p_hat)^2 /
      (dat$n_VA * p_hat * (1 - p_hat))
  )

  tibble::tibble(
    model = fit$model,
    k_par = k,
    logLik = loglik,
    AIC = aic,
    # AICc and BIC use the number of grouped pairs.
    AICc = aic + 2 * k * (k + 1) / (n_rows - k - 1),
    BIC = -2 * loglik + log(n_rows) * k,
    R2_probability = 1 -
      sum((dat$R_VA - p_hat)^2) /
        sum((dat$R_VA - mean(dat$R_VA))^2),
    RMSE_probability = sqrt(mean((dat$R_VA - p_hat)^2)),
    McFadden_R2 = 1 - loglik / null_loglik,
    Pearson_dispersion = pearson / (n_rows - k),
    converged_starts = sum(fit$diagnostics$converged),
    total_starts = nrow(fit$diagnostics)
  )
}

choice_save_diagnostics <- function(fit_lookup, OUTPUT_TAG) {
  start_diagnostics <- dplyr::bind_rows(lapply(fit_lookup, function(fit) fit$diagnostics))
  weight_diagnostics <- dplyr::bind_rows(lapply(fit_lookup, choice_weight_diagnostics))
  write.csv(start_diagnostics,
            paste0("results/Study3_ItemLevel_", OUTPUT_TAG, "_StartDiagnostics.csv"), row.names = FALSE)
  write.csv(weight_diagnostics,
            paste0("results/Study3_ItemLevel_", OUTPUT_TAG, "_WeightDiagnostics.csv"), row.names = FALSE)
}
