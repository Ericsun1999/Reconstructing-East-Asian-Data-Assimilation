here::i_am("Code/SupportingAnalyses/test_observation_variance.R")

# Lightweight numerical checks; no source workbooks or LME data needed.
source(here::here("Code", "Analysis", "extra_observation_variance.R"))

V <- annual_unconditional_variance(
  M = c(0.5, 0.5),
  r2 = c(4, 1, 1)
)
stopifnot(isTRUE(all.equal(V, c(4, 2, 1.5), tolerance = 1e-12)))

# When the baseline predictive variance is 1 and squared residuals
# average to 4, the predictive-normal log score selects u = 3.
residual <- c(-2, 2, -2, 2)
base <- rep(1, length(residual))
objective <- function(u) predictive_log_score(residual, base, u)
fit <- stats::optimize(objective, c(0, 20), tol = 1e-10)
stopifnot(
  abs(fit$minimum - 3) < 1e-5,
  objective(fit$minimum) < objective(0),
  objective(0) == predictive_log_score(residual, base, 0),
  is.finite(objective(5))
)

# Adding u to the measurement error should reduce the Kalman gain and
# never cause a negative filtered variance.
P <- 2
beta <- 1.2
vdelta <- 0.5
u <- 1.75
gain0 <- beta * P / (beta^2 * P + vdelta)
gain1 <- beta * P / (beta^2 * P + vdelta + u)
P0 <- (1 - beta * gain0) * P
P1 <- (1 - beta * gain1) * P
stopifnot(gain1 < gain0, P1 > P0, P0 >= 0, P1 >= 0)

# Test the real measurement-projection function in both city and grid
# pipelines. When kappa=0, kriging reverts to its zero-mean prior,
# so the new projection must have alpha=g(0), beta=0, and no excess
# Monte Carlo residual variance. Preserve the historical alternative.
variance_tolerance <- 1e-10
load_mc <- function(path) {
  expressions <- parse(file = path)
  candidates <- Filter(function(expr) {
    is.call(expr) &&
      identical(expr[[1L]], as.name("<-")) &&
      is.symbol(expr[[2L]]) &&
      as.character(expr[[2L]]) == "compute_measurement_parameters_mc"
  }, as.list(expressions))
  if (length(candidates) != 1L) {
    stop("Unable to locate measurement projection in ", path)
  }
  eval(candidates[[1L]][[3L]], envir = .GlobalEnv)
}
paths <- c(
  here::here("Code", "Figure9d.R"),
  here::here("Code", "Assimilation_grid.R")
)
for (path in paths) {
  mc <- load_mc(path)
  args <- list(
    sigmaY2 = 1, vhat = c(0, 1),
    reference_yhat = c(-0.4, -0.2),
    reference_nu = c(0.8, 0.8),
    target_nu = c(1, 0),
    g = identity, n_mc = 1000L, seed = 7L
  )
  fixed <- do.call(mc, c(args, list(center_mode = "zero")))
  legacy <- do.call(mc, c(args, list(center_mode = "reference_mean")))
  stopifnot(
    abs(fixed$alpha[1L]) < 1e-10,
    abs(fixed$beta[1L]) < 1e-10,
    fixed$vdelta[1L] == 1e-10,
    abs(legacy$alpha[1L] - mean(args$reference_yhat)) < 1e-10,
    abs(fixed$beta[2L] - 1) < 1e-10,
    abs(fixed$alpha[2L]) < 1e-10,
    fixed$vdelta[2L] == 1e-10
  )
}

message("Observation-discrepancy variance numerical checks passed.")
