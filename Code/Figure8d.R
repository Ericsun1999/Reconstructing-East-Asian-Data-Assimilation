#Kriged REACHES

tempe_all_data <- read.csv("~/Downloads/DA3/tempe_all_v3.csv", header=FALSE)
year3<- as.integer(tempe_all_data[1,-c(1,2)]) 
tempe_all<-tempe_all_data[c(2:4),]

Data1 <- read.csv("~/Downloads/DA3/d1.csv", row.names=1)
Data2 <- read.csv("~/Downloads/DA3/d2.csv", row.names=1)
Data3 <- read.csv("~/Downloads/DA3/d3.csv", row.names=1)

#Data1 Hong Kong
#Data2 Shanghai
#Data3 Beijing

#Beijing
DData<-Data3 #Plot Shanghai or hong kong just change it to Data2 or Data1
loc <- 3  #Shanghai or hong kong just change 3 to 2 or 1

tempe_use<-tempe_all[loc,] 
nu1 <- (read.csv("~/Downloads/DA3/tempe_all_std.csv"))
nu <- nu1[loc, -c(1,2)] 

haa<- (DData[,-c(1,2)]) - 273.15
haave <- colMeans(haa)
haave2 <- as.numeric(t(haa))

#QM 
library(ggplot2)
library(np)
  
# ---- Input you need to prepare ----
# x     : vector of observed LME temperatures (historical sample)
# yhat  : vector of kriging predictions Ŷ_t for each year
# nu    : vector of kriging RMSE (ν_t), one per year
  
yhat <- as.numeric(as.matrix(tempe_use[,-c(1,2)]))
nu <- as.numeric(as.matrix(nu))

# ---- Estimate F_X using npudist with bandwidth chosen by LS-CV ----
fx_bw <- npudistbw(dat = c(haave2))       # select bindwidth h via LS-CV

Fx_hat <- function(q) {
  fx_ob <- npudist(bws = fx_bw, edat = data.frame(x = q))
  # Return estimated CDF values for given q using npudist
  fitted(fx_ob) 
}

# ---- Estimate F_Y following equation (9) ----
FY_hat <- function(y, yhat, nu) {
  # Ensure numeric vectors
  yhat <- as.numeric(yhat)
  nu   <- as.numeric(nu)
  stopifnot(length(yhat) == length(nu))
  
  # Avoid division by zero
  nu_safe <- pmax(nu, 1e-8)
  
  # Vectorized computation
  sapply(y, function(yy) {
    mean(pnorm((yy - yhat) / nu_safe))
  })
}

# ---- Inverse of F_X (quantile function) using uniroot ----
qx_min <- min(haave2) - 5 * sd(haave2)
qx_max <- max(haave2) + 5 * sd(haave2)

Fx_inv <- function(u) {
  # u must be between 0 and 1
  u <- pmin(pmax(u, 0 + 1e-8), 1 - 1e-8)
  sapply(u, function(uu) {
    f <- function(q) Fx_hat(q) - uu
    uniroot(f, interval = c(qx_min, qx_max), tol = 1e-6)$root
  })
}

ycorrected <- Fx_inv(FY_hat(as.numeric(tempe_use[,-c(1,2)]), yhat, nu))

## ---- Inputs you already have ----
## x    : sample used to estimate F_X (e.g., haave2 or stacked LME temps)
## Fx_hat(y) and FY_hat(y, yhat, nu) are already defined above

x <- as.numeric(haave2)                 # or whatever vector you used for Fx
h <- fx_bw$bw                           # bandwidth used in npudist for F_X
stopifnot(length(h) == 1 && h > 0)

## ---- Density of Y (derivative of your FY_hat) ----
fY_hat <- function(y, yhat, nu) {
  nu_safe <- pmax(as.numeric(nu), 1e-8)
  yhat    <- as.numeric(yhat)
  stopifnot(length(yhat) == length(nu_safe))
  mean(dnorm((y - yhat) / nu_safe) / nu_safe)
}

## ---- Density of X from the same kernel-CDF estimator for F_X ----
##      derivative: mean[ phi((q - x_i)/h) / h ]
fX_hat <- function(q) {
  mean(dnorm((q - x) / h) / h)
}

## ---- Inverse CDF g(y) = F_X^{-1}(F_Y(y)) via uniroot ----
g_of <- function(y) {
  u <- FY_hat(y, yhat = yhat, nu = nu)      # target probability in (0,1)
  # search bounds: expand the data range a bit (±6h is ~tail-safe for Gaussian)
  lo <- min(x) - 5 * h
  hi <- max(x) + 5 * h
  uniroot(function(q) Fx_hat(q) - u, lower = lo, upper = hi, tol = 1e-8)$root
}

g_fun <- Vectorize(g_of, SIMPLIFY = TRUE)

y_grid <- seq(-2.5, 2.5, length.out = 1000)
g_grid <- sapply(y_grid, g_of)   
g_fun_fast <- approxfun(y_grid, g_grid)

## =========================================================
##  Kalman filter + RTS smoother with measurement equation
##  X*_t = alpha_t + beta_t X_t + delta_t
##  where (alpha_t, beta_t, var(delta_t)) computed numerically
##  from bivariate normal of (Y_t, Yhat_t) in Eq (14).
## =========================================================

library(mvtnorm)

round_and_clip <- function(x, lower = -2, upper = 1, digits = 2) {
  x <- round(x, digits = digits)   
  pmin(pmax(x, lower), upper)      
}

compute_meas_params_mc <- function(
  sigmaY2,
  v,            # length N, v_t = cov(Y_t, Yhat_t) = var(Yhat_t)
  g,
  n_mc = 50000,
  seed = 1
) {
  set.seed(seed)
  N <- length(v)

  alpha <- numeric(N)
  beta  <- numeric(N)
  vdelta <- numeric(N)

  for (t in seq_len(N)) {
    vt <- v[t]

    ## Cov matrix in Eq (14):
    ## [ var(Y)        cov(Y,Yhat) ]
    ## [ cov(Y,Yhat)   var(Yhat)   ]
    ## = [ sigmaY2  vt ]
    ##   [ vt       vt ]
    Sigma <- matrix(c(sigmaY2, vt,
                      vt,      vt), nrow = 2, byrow = TRUE)

    eig <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
    if (min(eig) <= 0) {
      Sigma <- Sigma + diag(abs(min(eig)) + 1e-10, 2)
    }

    samp <- rmvnorm(n = n_mc, mean = c(0, 0), sigma = Sigma)

    Y    <- round_and_clip(samp[, 1], lower = -2, upper = 1, digits = 3)
    Yhat <- round_and_clip(samp[, 2], lower = -2, upper = 1, digits = 3)

    X    <- g(Y)
    Xst  <- g(Yhat)

    EX  <- mean(X)
    EXs <- mean(Xst)
    VX  <- var(X)
    COV <- cov(Xst, X)

    bt <- COV / VX
    at <- EXs - bt * EX

    ## delta_t = X*_t - alpha_t - beta_t X_t
    del <- Xst - at - bt * X

    alpha[t]  <- at
    beta[t]   <- bt
    vdelta[t] <- var(del)
  }

  list(alpha = alpha, beta = beta, vdelta = vdelta)
}


## ---------------------------
## 3) Kalman filter + RTS smoother (Sec 4.3)
##    Eqs (16)-(23) :contentReference[oaicite:3]{index=3}
## ---------------------------
kalman_filter_smoother <- function(mu, M, r2, Xstar, alpha, beta, vdelta) {
  N <- length(mu)
  stopifnot(length(r2) == N,
            length(Xstar) == N,
            length(alpha) == N,
            length(beta) == N,
            length(vdelta) == N,
            length(M) == N - 1)

  ## Storage
  X_pred <- numeric(N)  # X_{t|t-1}
  P_pred <- numeric(N)  # P_{t|t-1}
  X_filt <- numeric(N)  # X_{t|t}
  P_filt <- numeric(N)  # P_{t|t}

  ## ---- t=1 prior (from Eq (12) / AR prior idea)
  X_pred[1] <- mu[1]
  P_pred[1] <- r2[1]

  ## Update at t=1
  K1 <- beta[1] * P_pred[1] / (beta[1]^2 * P_pred[1] + vdelta[1])  # Eq (18)
  X_filt[1] <- X_pred[1] + K1 * (Xstar[1] - alpha[1] - beta[1] * X_pred[1]) # Eq (19)
  P_filt[1] <- (1 - beta[1] * K1) * P_pred[1] # Eq (20)

  ## ---- t=2..N
  for (t in 2:N) {
    Mt1 <- M[t - 1]

    ## Prediction: Eq (16)(17)
    X_pred[t] <- mu[t] + Mt1 * (X_filt[t - 1] - mu[t - 1])
    P_pred[t] <- Mt1^2 * P_filt[t - 1] + r2[t]

    ## Update: Eq (18)(19)(20)
    Kt <- beta[t] * P_pred[t] / (beta[t]^2 * P_pred[t] + vdelta[t])
    X_filt[t] <- X_pred[t] + Kt * (Xstar[t] - alpha[t] - beta[t] * X_pred[t])
    P_filt[t] <- (1 - beta[t] * Kt) * P_pred[t]
  }

  ## ---- RTS smoother: Eq (21)(22)(23)
  X_smooth <- numeric(N)
  P_smooth <- numeric(N)
  J <- numeric(N - 1)

  X_smooth[N] <- X_filt[N]
  P_smooth[N] <- P_filt[N]

  for (t in (N - 1):1) {
    Mt <- M[t]
    ## Eq (21): J_t = P_{t|t} M_t P^{-1}_{t+1|t}
    J[t] <- P_filt[t] * Mt / P_pred[t + 1]

    ## Eq (22)(23)
    X_smooth[t] <- X_filt[t] + J[t] * (X_smooth[t + 1] - X_pred[t + 1])
    P_smooth[t] <- P_filt[t] + J[t]^2 * (P_smooth[t + 1] - P_pred[t + 1])
  }

  list(
    X_pred = X_pred, P_pred = P_pred,
    X_filt = X_filt, P_filt = P_filt,
    X_smooth = X_smooth, P_smooth = P_smooth,
    J = J
  )
}


  dfc <- read.csv("./mtB.csv")
  df2c <- read.csv("./muB.csv")
  df3c <- read.csv("./rtB.csv")
  
  
  sigma<- read.csv("~/Downloads/DA3/tempe_all_std.csv")
  stdd1<- as.numeric(sigma[loc, -c(1,2)]) 
  
  sigmaY22 <- as.numeric(array(0.887, 524)) 
  varr2 <- sigmaY22 - stdd1
  
  year4<- year3[1:523]
  
  mumu <- df2c$value[year3 - 1367]
  mtmt <- dfc$value[year4 - 1367]
  rtrt <- df3c$value[year3 - 1367]


meas <- compute_meas_params_mc(sigmaY2 = sigmaY22[1], v = varr2, g = g_fun_fast, n_mc = 10000, seed = 1)

out <- kalman_filter_smoother(mu = mumu, M = mtmt, r2 = rtrt,
                               Xstar = ycorrected,
                               alpha = meas$alpha, beta = meas$beta, vdelta = meas$vdelta)

 XtN = out$X_smooth
 PtN = out$P_smooth
 
    

df87<-data.frame(year = c(year4), predicted = out$X_smooth[1:523], 
                 sigmasmooth2 = out$P_smooth[1:523], REACHES = ycorrected[1:523], 
                 alpha = meas$alpha[1:523], beta = meas$beta[1:523], 
                 vdelata = meas$vdelta[1:523],
                 LME = haave[year4-1350])

#Figure8d
#jpeg("Figure8d.png",width=6,height=3, res=300, units = "in")
      print(
      ggplot(data = df87, aes(x = year, y = predicted)) +
      geom_line(color = "red") +
      theme(text=element_text(size=11),legend.position="right",
            legend.key.height=unit(1.5,"cm"),
            plot.title = element_text(hjust = 0.5)) +
      geom_ribbon(aes(ymin = predicted - sqrt(sigmasmooth2),
                      ymax = predicted + sqrt(sigmasmooth2)),
                  alpha = 0.5, fill = "grey3") +
      xlab("year") +
      #ylim(10.4,12.4) +
      #ylim(8.85,12.8) +
      ylab("temperature") 
      )
#   dev.off()

   
