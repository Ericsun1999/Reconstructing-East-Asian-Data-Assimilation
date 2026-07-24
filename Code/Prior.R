here::i_am("Code/Prior.R")

### prior estimation

Data1 <- read.csv("./d1.csv", row.names=1)
Data2 <- read.csv("./d2.csv", row.names=1)
Data3 <- read.csv("./d3.csv", row.names=1)

#Data1 Hong Kong
#Data2 Shanghai
#Data3 Beijing

#Beijing
Data<-Data3 #Plot Shanghai or hong kong just change it to Data2 or Data1


haa1 <- t(Data[,c(21:564)])
haa1 <- haa1 -273.15


# D: size (n-1) x n
build_D <- function(n) {
  D <- matrix(0, n - 1, n)
  for (i in 1:(n - 1)) {
    D[i, i]   <- -1
    D[i, i+1] <-  1
  }
  D
}

ell_one_series <- function(xj, mu, M, r2) {
  N <- length(mu)
  if (any(!is.finite(r2)) || any(r2 <= 0)) {
    return(Inf)
  }
  val <- 0.5 * sum(log(r2))
  # t = 1
  val <- val + (xj[1] - mu[1])^2 / (2 * r2[1])
  # t >= 2
  if (N >= 2) {
    resid <- xj[2:N] - mu[2:N] - M[1:(N - 1)] * (xj[1:(N - 1)] - mu[1:(N - 1)])
    val <- val + sum(resid^2 / (2 * r2[2:N]))
  }
  val
}

compute_St <- function(X, mu, M) {
  N <- nrow(X)
  J <- ncol(X) # 應為 13
  S <- numeric(N)
  # S1
  S[1] <- mean( (X[1, ] - mu[1])^2 )
  if (N >= 2) {
    for (t in 2:N) {
      resid_t <- X[t, ] - mu[t] - M[t-1] * (X[t-1, ] - mu[t-1])
      S[t] <- mean( resid_t^2 )
    }
  }
  S
}

solve_r2_cubic <- function(S, neighbor_values, lambda3, J, x_init = NULL,
  eps = 1e-10
) {
  d <- length(neighbor_values)
  q <- sum(neighbor_values)

  if (lambda3 == 0 || d == 0) {
    return(max(S, eps))
  }

  # 4 lambda3 d x^3 - 4 lambda3 q x^2 + J x - J S = 0
  roots <- polyroot(c(
    -J * S,
     J,
    -4 * lambda3 * q,
     4 * lambda3 * d
  ))
  is_real <- abs(Im(roots)) < 1e-8
  cand <- Re(roots[is_real])
  cand <- cand[is.finite(cand) & cand > eps]

  coordinate_objective <- function(x) {
    likelihood_part <- (J / 2) * (log(x) + S / x)
    penalty_part <- lambda3 * sum((x - neighbor_values)^2)

    likelihood_part + penalty_part
  }

  if (length(cand) > 0) {
    values <- vapply(cand, coordinate_objective, numeric(1))
    return(cand[which.min(values)])
  }


  center <- max(c(S, neighbor_values, x_init, eps), na.rm = TRUE)
  z_center <- log(center)
  opt <- optimize(
    function(z) coordinate_objective(exp(z)),
    interval = c(z_center - 20, z_center + 20)
  )

  max(exp(opt$minimum), eps)
}


update_M <- function(X, mu, r2, lambda1) {
  N <- nrow(X)
  J <- ncol(X)
  D_M <- build_D(N - 1)          
  # A: diag_t sum_j (X_t^{(j)} - mu_t)^2 / r_{t+1}^2
  A_diag <- numeric(N - 1)
  for (t in 1:(N - 1)) {
    A_diag[t] <- sum( (X[t, ] - mu[t])^2 ) / r2[t+1]
  }
  A <- diag(A_diag, nrow = N - 1)
  # a_t: sum_j (X_t^{(j)} - mu_t)*(X_{t+1}^{(j)} - mu_{t+1}) / r_{t+1}^2
  a <- numeric(N - 1)
  for (t in 1:(N - 1)) {
    a[t] <- sum( (X[t, ] - mu[t]) * (X[t+1, ] - mu[t+1]) ) / r2[t+1]
  }
  # M = (A + 2 λ1 D'D)^{-1} a
  K <- A + 2 * lambda1 * t(D_M) %*% D_M
  as.vector( solve(K, a) )
}


update_mu <- function(X, M, r2, lambda2) {
  N <- nrow(X); J <- ncol(X)
  D_mu <- build_D(N)
  B <- matrix(0, N, N)
  b <- numeric(N)
  B[1, 1] <- B[1, 1] + J / r2[1]
  b[1]    <- b[1]    + sum(X[1, ]) / r2[1]  
  
  for (t in 2:N) {
    w <- J / r2[t]               
    m <- M[t - 1]
    B[t,   t]   <- B[t,   t]   + w
    B[t-1, t-1] <- B[t-1, t-1] + w * m^2
    B[t,   t-1] <- B[t,   t-1] - w * m
    B[t-1, t]   <- B[t-1, t]   - w * m

    S <- sum(X[t, ]) - m * sum(X[t-1, ])     # = ∑_j (x_t^{(j)} - m x_{t-1}^{(j)})
    b[t]   <- b[t]   + S / r2[t]             # + (1/r_t^2) ∑ S
    b[t-1] <- b[t-1] - m * (S / r2[t])       # - (m/r_t^2) ∑ S  
  }

  K  <- B + 2 * lambda2 * crossprod(D_mu)     # 2*λ2 * D'D
  as.vector(solve(K, b))
}


update_r2 <- function(X, mu, M, r2, lambda3) {
  N <- nrow(X)
  J <- ncol(X)  
  S <- compute_St(X, mu, M)
  r2_new <- r2
  for (t in seq_len(N)) {
    # Sequential coordinate update：
    neighbor_values <- numeric(0)
    if (t > 1) {
      neighbor_values <- c(neighbor_values, r2_new[t - 1])
    }
    if (t < N) {
      neighbor_values <- c(neighbor_values, r2[t + 1])
    }
    r2_new[t] <- solve_r2_cubic(
      S = S[t],
      neighbor_values = neighbor_values,
      lambda3 = lambda3,
      J = J,
      x_init = r2[t]
    )
    if (!is.finite(r2_new[t]) || r2_new[t] <= 0) {
      r2_new[t] <- max(S[t], 1e-10)
    }
  }
  r2_new
}


fit_theta_once <- function(X, lambda1, lambda2, lambda3, max_iter = 100, tol = 1e-6, verbose = FALSE) {
  N <- nrow(X); J <- ncol(X)
  mu  <- as.vector(rowMeans(X))
  r2  <- as.vector(apply(X, 1, function(z) mean( (z - mean(z))^2 )))
  M   <- rep(0, N - 1)
  obj <- function() {
    ell_sum <- 0
  for (j in seq_len(J)) {
    ell_sum <- ell_sum +
      ell_one_series(X[, j], mu, M, r2)
  }
  D_M  <- build_D(N - 1)
  D_mu <- build_D(N)
  D_r2 <- build_D(N)

  penalty <- lambda1 * sum((D_M %*% M)^2) +
    lambda2 * sum((D_mu %*% mu)^2) +
    lambda3 * sum((D_r2 %*% r2)^2)

  ell_sum + penalty
  }
  prev <- obj()
  for (it in 1:max_iter) {
    # 2) Update M
    M  <- update_M(X, mu, r2, lambda1)
    # 3) Update mu
    mu <- update_mu(X, M, r2, lambda2)
    # 4) Update r2
    r2 <- update_r2(X, mu, M, r2, lambda3)
    cur <- obj()
    if (verbose) cat(sprintf("Iter %d: obj = %.6f\n", it, cur))
    if (abs(prev - cur) < tol * (1 + abs(prev))) break
    prev <- cur
  }
  list(M = M, mu = mu, r2 = r2, obj = prev, iters = it, converged = (it < max_iter))
}

cv_select_lambdas <- function(
  X,
  lambda1_grid = 10^seq(-3, 1, length.out = 5),
  lambda2_grid = 10^seq(-3, 1, length.out = 5),
  lambda3_grid = 10^seq(-3, 1, length.out = 5),
  max_iter = 50, tol = 1e-3, verbose = FALSE
) {
  J <- ncol(X)
  combos <- expand.grid(lambda1 = lambda1_grid, lambda2 = lambda2_grid, lambda3 = lambda3_grid)
  best <- list(score = Inf)
  for (i in 1:nrow(combos)) {
    l1 <- combos$lambda1[i]; l2 <- combos$lambda2[i]; l3 <- combos$lambda3[i]
    # 13 folds
    scores <- numeric(J)
    for (l in 1:J) {
      X_train <- X[, setdiff(1:J, l), drop = FALSE]
      X_test  <- X[, l]
      fit <- fit_theta_once(X_train, l1, l2, l3, max_iter = max_iter, tol = tol, verbose = FALSE)
      scores[l] <- ell_one_series(X_test, fit$mu, fit$M, fit$r2)
    }
    cv <- mean(scores)
    if (verbose) cat(sprintf("λ1=%.4g λ2=%.4g λ3=%.4g  CV=%.6f\n", l1, l2, l3, cv))
    if (cv < best$score) best <- list(score = cv, lambda1 = l1, lambda2 = l2, lambda3 = l3)
  }
  best
}

fit_LME_fusedridge <- function(
  X,                            # N x 13 matrix (rows=time, cols=series)
  lambda1_grid = 10^seq(-3, 1, length.out = 5),
  lambda2_grid = 10^seq(-3, 1, length.out = 5),
  lambda3_grid = 10^seq(-3, 1, length.out = 5),
  max_iter = 3, tol = 1e-2, verbose = FALSE
) {
  stopifnot(is.matrix(X), ncol(X) == 13)
  best <- cv_select_lambdas(X, lambda1_grid, lambda2_grid, lambda3_grid, max_iter, tol, verbose)
  if (verbose) {
    cat("Best lambdas: ",
        sprintf("lambda1=%.6g, lambda2=%.6g, lambda3=%.6g; CV=%.6f\n",
                best$lambda1, best$lambda2, best$lambda3, best$score))
  }
  fit <- fit_theta_once(X, best$lambda1, best$lambda2, best$lambda3, max_iter, tol, verbose)
  list(
    lambdas = unlist(best[c("lambda1","lambda2","lambda3")]),
    CV = best$score,
    theta = list(M = fit$M, mu = fit$mu, r2 = fit$r2),
    optimization = list(obj = fit$obj, iters = fit$iters, converged = fit$converged)
  )
}

res <- fit_LME_fusedridge(
  haa1,
  lambda1_grid = 10^seq(-3, 3, length.out = 25),  
  lambda2_grid = 10^seq(-3, 3, length.out = 25),
  lambda3_grid = 10^seq(-3, 3, length.out = 25),
  max_iter = 100,                                 
  tol = 1e-2,
  verbose = TRUE
)

res1<-fit_theta_once(haa1,0,0,0,100,1e-2)

library(tidyr)
                         
mt<-res$theta$M
mt0<-res1$M
mu<-res$theta$mu
mu0<-res1$mu
rt<-res$theta$r2
rt0<-res1$r2
dfc<- data.frame(year = c(1368:1910), penalized = mt[1:543], unpenalized = mt0[1:543])
dfc<-tidyr::gather(dfc,key = "coefficient", value = "value", -year)
dfc[1:543,2]<-"Penalized ML"
dfc[544:1086,2]<-"ML"
df2c<- data.frame(year = c(1368:1911), penalized = mu, unpenalized = mu0)
df2c<-tidyr::gather(df2c,key = "coefficient", value = "value", -year)
df2c[1:544,2]<-"Penalized ML"
df2c[545:1088,2]<-"ML"
df3c<- data.frame(year = c(1368:1911), penalized = rt, unpenalized = rt0)
df3c<-tidyr::gather(df3c,key = "coefficient", value = "value", -year)
df3c[1:544,2]<-"Penalized ML"
df3c[545:1088,2]<-"ML"


# write csv
write_excel_csv(dfc,"mtB.csv")
write_excel_csv(df2c,"muB.csv")
write_excel_csv(df3c,"rtB.csv")

write_excel_csv(dfc,"mtS.csv")
write_excel_csv(df2c,"muS.csv")
write_excel_csv(df3c,"rtS.csv")

write_excel_csv(dfc,"mtH.csv")
write_excel_csv(df2c,"muH.csv")
write_excel_csv(df3c,"rtH.csv")

