here::i_am("Code/Analysis/test_observation_variance.R")

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
fit <- stats::optimize(objective, c(0, 20))
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

message("Observation-discrepancy variance numerical checks passed.")
