here::i_am("Code/Figure5.R")

###It should be executed after running Figure4.R

# Kriging on 0.5x0.5 grid points

year.start <- range(temp2$year)[1]
year.end <- range(temp2$year)[2]
year.all <- year.start:year.end

year2 <- year.all[c(98,484)]
arr.pred <- arr.std <- array(NA, c(53, 49, length(year2)))
loc <- expand.grid(long = seq(98.25, 124.25, by = 0.5), lat = seq(18.25, 42.25, by = 0.5))
coordinates(loc) <- ~ long + lat
proj4string(loc) <- CRS('+proj=longlat +datum=WGS84')
y2 <- var.fit2$y2

library(RColorBrewer)
library(ggplot2)

#Figure5(a1), 5(b1)

  for(i in 1:length(year2)) {
    temp.dat <- as.data.frame(y2[y2$year==year2[i],])
    #jpeg(paste(getwd(),"/data-temp-",i,".jpg",sep=""),width=5,height=4, res=300, units = "in")
     print(
     ggplot(temp.dat,aes(long,lat)) +
      geom_point(aes(colour=level),cex=2) +
      ggtitle(paste0("Year ",year2[i])) + 
      coord_map(xlim=c(98,124.5),ylim=c(18,42.5)) + 
      scale_colour_gradientn(colours=rev(brewer.pal(n=9,name='RdBu')),limits=c(-2,2),
      na.value="transparent",guide="colourbar") +
      borders(database="world",xlim=c(76,132),ylim=c(19,52),fill=NA,colour="grey50") +
      theme(text=element_text(size=15),legend.title=element_blank(),
        legend.position=c(1.15, 0.7),legend.key.height=unit(2.4,"cm"))
     )
     #dev.off()
  }


library(MASS)
library(mvtnorm)

sigmay = sqrt(vario.fit2$psill2)
sigmae = sqrt(vario.fit2$psill1)

sigma = sqrt(vario.fit2$psill1 + vario.fit2$psill2)
alpha = vario.fit2$range/100
cuts = c(-Inf, -1.5, -0.5, 0.5, Inf)
vals = c(-2, -1, 0, 1)

mu_z <- -pnorm(-1.5, mean = 0, sd = sqrt(sigmay + sigmae))

cov_Y <- function(s1, s2) {
  (sigmay)^2 * exp(-sqrt(sum((s1-s2)^2))/alpha)
}

calc_c_Y <- function(s0, locations) {
  n <- nrow(locations)
  c_Y <- numeric(n)
  for (i in 1:n) {
    si <- locations[i, ]
    c_Y[i] <- cov_Y(s0, si)
  }
  return(c_Y)
}

# E[h(Z*)]
Ez_discrete <- function(sigma, cuts, vals) {
  stopifnot(length(vals) == length(cuts) - 1)
  a <- head(cuts, -1) / sigma
  b <- tail(cuts, -1) / sigma
  probs <- pnorm(b) - pnorm(a)
  sum(vals * probs)
}


# E[Z* h(Z*)]
EZstar_h <- function(sigma, cuts, vals) {
  stopifnot(length(vals) == length(cuts) - 1)
  a <- head(cuts, -1) / sigma
  b <- tail(cuts, -1) / sigma
  part <- sigma * (dnorm(a) - dnorm(b))
  sum(vals * part)
}

# Cov(Z_i, Z_j)
cov_Z_pair <- function(rho, sigma, cuts, vals, Eh = NULL) {
  if (is.null(Eh)) Eh <- Ez_discrete(sigma, cuts, vals)
  K <- length(vals)
  a_std <- head(cuts, -1) / sigma
  b_std <- tail(cuts, -1) / sigma

  Sigma2 <- matrix(c(1, rho, rho, 1), 2, 2)
  Eh2 <- 0
  for (k in seq_len(K)) {
    for (l in seq_len(K)) {
      lower <- c(a_std[k], a_std[l])
      upper <- c(b_std[k], b_std[l])
      pij <- as.numeric(pmvnorm(lower = lower, upper = upper, mean = c(0,0), sigma = Sigma2))
      Eh2 <- Eh2 + vals[k] * vals[l] * pij
    }
  }
  Eh2 - Eh^2
}

cZY_vector <- function(s_coords, s0, sigma_Y2, sigma_E2, alpha,
                       cuts = c(-Inf, -1.5, -0.5, 0.5, Inf),
                       vals = c(-2, -1, 0, 1)) {
  sigma2 <- sigma_Y2 + sigma_E2
  sigma  <- sqrt(sigma2)
  Ezstar_h <- EZstar_h(sigma, cuts, vals)
  dists <- sqrt(rowSums((s_coords - matrix(s0, nrow(s_coords), ncol(s_coords), byrow=TRUE))^2))
  covYY  <- sigma_Y2 * exp(-dists / alpha)
  (covYY / sigma2) * Ezstar_h
}

cZY_matrix <- function(s_coords, s0_mat, sigma_Y2, sigma_E2, alpha,
                       cuts, vals) {
  if (is.null(dim(s0_mat))) s0_mat <- matrix(s0_mat, nrow = 1)
  n <- nrow(s_coords); m <- nrow(s0_mat)
  sigma2 <- sigma_Y2 + sigma_E2
  sigma  <- sqrt(sigma2)
  # E[Z* h(Z*)]
  EZh <- EZstar_h(sigma, cuts, vals)
  dists <- as.matrix(dist(rbind(s_coords, s0_mat)))[1:n, (n+1):(n+m)]
  covYY <- sigma_Y2 * exp(-dists / alpha)          # n × m
  (covYY / sigma2) * EZh                           # n × m
}

SigmaZ_matrix <- function(s_coords, sigma_Y2, sigma_E2, alpha,
                          cuts = c(-Inf, -1.5, -0.5, 0.5, Inf),
                          vals = c(-2, -1, 0, 1)) {
  n <- nrow(s_coords)
  sigma2 <- sigma_Y2 + sigma_E2
  sigma  <- sqrt(sigma2)
  Eh <- Ez_discrete(sigma, cuts, vals)
  
  # d_ij = ||s_i - s_j||
  dmat <- as.matrix(dist(s_coords, method = "euclidean", upper = TRUE, diag = TRUE))
  rho  <- (sigma_Y2 * exp(-dmat / alpha)) / sigma2
  diag(rho) <- 1  

  Sig <- matrix(NA_real_, n, n)
  for (i in 1:n) {
    for (j in i:n) {
      cij <- cov_Z_pair(rho[i,j], sigma, cuts, vals, Eh = Eh)
      Sig[i,j] <- cij
      if (j != i) Sig[j,i] <- cij
    }
  }
  Sig
}

y_pred <- function(s_coords, s0,
                   z_obs = NULL,
                   sigma_Y2, sigma_E2, alpha,
                   cuts = c(-Inf, -1.5, -0.5, 0.5, Inf),
                   vals = c(-2, -1, 0, 1),
                   tol = 1e-2,
                   return = c("mean", "var")) {

  #stabilize <- match.arg(stabilize)
  return <- unique(return)

  if (!is.null(z_obs)) {
    stopifnot(length(z_obs) == nrow(s_coords))
  }

  c_zy   <- cZY_vector(s_coords, s0, sigma_Y2, sigma_E2, alpha, cuts, vals)
  SigmaZ <- SigmaZ_matrix(s_coords, sigma_Y2, sigma_E2, alpha, cuts, vals)
  SigmaZ <- (SigmaZ + t(SigmaZ)) / 2  

  out <- list()
  
  # --- mean: c' Σ^{-1} (z - Ez) ---
  if ("mean" %in% return) {
    sigma2 <- sigma_Y2 + sigma_E2
    sigma  <- sqrt(sigma2)
    Ez <- Ez_discrete(sigma, cuts, vals)

    rhs <- z_obs - Ez
    w_mean <- solve(SigmaZ, rhs, tol = tol)
    out$mean <- drop(crossprod(c_zy, w_mean))
  }

  # --- var: sigma_Y2 - c' Σ^{-1} c ---
  if ("var" %in% return) {
    w_var <- solve(SigmaZ, c_zy, tol = tol)
    var_pred <- sigma_Y2 - drop(crossprod(c_zy, w_var))
    out$var <- max(var_pred, 0)
  }
  if (length(out) == 1) return(out[[1]])
  out
}

#Figure5(a2), 5(a3), 5(b2), 5(b3)

  for(i in 1:length(year2)) {
    temp14 <- y2[y2$year==year2[i],]
    temp15 <- rep(0, length(loc@coords[,1]))
    temp16 <- rep(0, length(loc@coords[,1]))
    
    locations<- temp14@coords
    
    
    c_Y <- matrix(0, nrow = length(loc@coords[,1]), ncol = length(temp14$level))
    Mat1<- matrix(mu_z, nrow = 1, ncol = length(temp14$level))
    
    for (j in 1:length(loc@coords[,1])) {
      
      s0 <- as.numeric(loc@coords[j,])
      c_Y[j,] <- calc_c_Y(s0, locations)
      res <- y_pred(locations, s0,  as.matrix(temp14$level), sigmay, sigmae, alpha, 
                    return = c("mean","var"))
      temp15[j]<- res$mean
      temp16[j] <- res$var
    }
    arr.pred[,,i] <- matrix(temp15,53,49)
    arr.std[,,i] <- matrix(temp16,53,49)
    temp.dat1 <- cbind(as.data.frame(loc@coords)[,1:2],mu=c(arr.pred[,,i]), std = c(arr.std[,,i]))
     #jpeg(paste(getwd(),"/krige4-temp-",i,".jpg",sep=""),width=5,height=4, res=300, units = "in")
    print(
      ggplot(temp.dat1,aes(long,lat)) +
      geom_point(aes(colour=mu),cex=9) +
      ggtitle(paste0("Year ",year2[i])) + 
      coord_map(xlim=c(98,124.5),ylim=c(18,42.5)) + 
      scale_colour_gradientn(colours=rev(brewer.pal(n=9,name='RdBu')),limits=c(-2,2),
      na.value="transparent",guide="colourbar") +
      borders(database="world",xlim=c(76,132),ylim=c(19,52),fill=NA,colour="grey50") +
      theme(text=element_text(size=15),legend.title=element_blank(),
        legend.position=c(1.15, 0.7),legend.key.height=unit(2.4,"cm"))
    )
    #dev.off()
    #jpeg(paste(getwd(),"/krigestd4-temp-",i,".jpg",sep=""),width=5,height=4, res=300, units = "in")
    print(
      ggplot(temp.dat1,aes(long,lat)) +
      geom_point(aes(colour=std),cex=9) +
      ggtitle(paste0("Year ",year2[i])) + 
      coord_map(xlim=c(98,124.5),ylim=c(18,42.5)) + 
      scale_colour_gradient(
        low  = "skyblue3",
        high = "white",
        limits = c(0, 1),
        na.value = "transparent",
        guide = "colourbar"
      ) +
      borders(database="world",xlim=c(76,132),ylim=c(19,52),fill=NA,colour="grey50") +
      theme(text=element_text(size=15),legend.title=element_blank(),
        legend.position=c(1.15, 0.5),legend.key.height=unit(1.745,"cm"))
    )
    #dev.off()
  }
