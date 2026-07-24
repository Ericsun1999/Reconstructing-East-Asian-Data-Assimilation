here::i_am("Code/Figure7e.R")

###It should be executed after running Figure4.R

# Kriging on 2.5x1.89 grid points

year.start <- range(temp2$year)[1]
year.end <- range(temp2$year)[2]
year.all <- year.start:year.end
nyear <- year.end - year.start + 1
ncase1 <- rep(NA, nyear)
for (i in 1:nyear) {
  ncase1[i] <- sum(temp2$year == year.start + i - 1)
  if (ncase1[i] > 0) {
    temp21 <- temp2[temp2$year == (year.start + i - 1), ]
  }
}
m <- 1
year2 <- year.all[which(ncase1 >= m)]

arr.pred <- arr.std <- array(NA, c(12, 14, length(year2)))

loc <- expand.grid(long = seq(97.5, 125, by = 2.5), lat = seq(18, 42.63158, by = 1.89473692308))
coordinates(loc) <- ~ long + lat
proj4string(loc) <- CRS('+proj=longlat +datum=WGS84')

y2 <- var.fit2$y2

library(RColorBrewer)
library(ggplot2)
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

  }

  for(i in 1:length(year2)) {
    temp14 <- y2[y2$year==year2[i],]
    temp15 <- rep(0, length(loc@coords[,1]))
    temp16 <- rep(0, length(loc@coords[,1]))
    
    locations<- temp14@coords
    
    Sigma_zz <- calc_Sigma_zz(locations)
    invsigma_zz<- solve(Sigma_zz)
    c_Y <- matrix(0, nrow = length(loc@coords[,1]), ncol = length(temp14$level))
    Mat1<- matrix(mu_z, nrow = 1, ncol = length(temp14$level))
    
    for (j in 1:length(loc@coords[,1])) {
      
      s0 <- as.numeric(loc@coords[j,])
      c_Y[j,] <- calc_c_Y(s0, locations)
      temp15[j]<- t(as.matrix(c_Y[j,])) %*% as.matrix(invsigma_zz) %*% (as.matrix(temp14$level) - t(Mat1) )
      temp16[j] <- cov_Y(s0,s0) - t(c_Y[j,]) %*% as.matrix(invsigma_zz) %*% (c_Y[j,])
    }
    arr.pred[,,i] <- matrix(temp15,12,14)
    arr.std[,,i] <- matrix(temp16,12,14)
    
    if(i==1) sk.all <- sk.all1 <- as.data.frame(loc@coords)
    sk.all <- cbind(sk.all,c(arr.pred[,,i]))
    sk.all1 <- cbind(sk.all1,c(arr.std[,,i]))
    
  }
  
  names(sk.all)[-c(1,2)] <- year2
  names(sk.all1)[-c(1,2)] <- year2
  
tempe_all <- sk.all

library(sf)
library(rnaturalearth)
library(dplyr)
library(stringr)


world <- ne_countries(scale = "medium", returnclass = "sf")
targets <- world %>%
  filter(str_detect(admin, regex("China|Taiwan|Hong Kong|Macao|Macau", ignore_case = TRUE)))

china_tw_hk_mo <- st_union(targets) 
china_tw_hk_mo <- st_as_sf(data.frame(name = "China+Taiwan+HK+MO"), geometry = china_tw_hk_mo)
st_crs(china_tw_hk_mo) <- 4326

coast_km <- 80
china_tw_hk_mo_valid <- st_make_valid(china_tw_hk_mo)
poly_m <- st_transform(china_tw_hk_mo_valid, 3857)
poly_m_buf <- st_buffer(poly_m, dist = coast_km * 1000)
china_tw_hk_mo_buf <- st_transform(poly_m_buf, 4326)
pts_sf <- st_as_sf(tempe_all, coords = c("long", "lat"), crs = 4326)

tempe_all$in_china_tw_hk_mo_coast <- st_within(
  pts_sf,
  china_tw_hk_mo_buf,
  sparse = FALSE
)[,1]
tempe_all

tempe_all1 <- tempe_all[(tempe_all$in_china_tw_hk_mo_coast == "TRUE"), -c(527)]

nu_reaches <- sk.all1
nu_reaches1 <- nu_reaches[(tempe_all$in_china_tw_hk_mo_coast == "TRUE"),]

#Load LME data
lme_df <- read.csv("~/Downloads/DA3/b0.csv", row.names=1)

#QM 
library(np)
library(dplyr)
library(ggplot2)
library(maps)

n_loc <- 121
n_rep <- 13

get_block <- function(df, r){
  df[((r-1)*n_loc + 1):(r*n_loc), ]
}

get_xs_all_years <- function(lme_df, s){
  # s in 1:n_loc
  xs_list <- lapply(1:n_rep, function(r){
    blk <- get_block(lme_df, r)
    as.numeric(blk[s, -(1:2)])  # 年份欄
  })
  xs_kelvin <- unlist(xs_list)
  xs_kelvin - 273.15
}

build_Fx_inv_local <- function(x_s){
  x_s <- as.numeric(x_s)
  x_s <- x_s[is.finite(x_s)]

  bw <- npudistbw(dat = x_s)

  Fx_hat <- function(q){
    fit <- npudist(bws = bw, edat = data.frame(x = q))
    as.numeric(fitted(fit))
  }

  Fx_inv <- function(u){
    u <- pmin(pmax(u, 1e-8), 1 - 1e-8)
    xmin <- min(x_s) - 5*sd(x_s)
    xmax <- max(x_s) + 5*sd(x_s)
    sapply(u, function(uu){
      uniroot(function(q) Fx_hat(q) - uu,
              interval = c(xmin, xmax), tol = 1e-6)$root
    })
  }

  Fx_inv
}

# ---- F_Y,s(0) ----
FY_hat_1loc <- function(y, yhat, nu){
  yhat <- as.numeric(yhat)
  nu <- pmax(as.numeric(nu), 1e-8)
  stopifnot(length(yhat) == length(nu))
  mean(pnorm((y - yhat) / nu))
}

lme_base <- get_block(lme_df, 1) %>%
  transmute(lat = lati, lon = long)

reaches_base <- tempe_all1 %>%
  transmute(lon = long, lat = lat)

coord_key_lme <- lme_base %>% mutate(key = paste0(round(lon,3), "_", round(lat,3)))
coord_key_rea <- reaches_base %>% mutate(key = paste0(round(lon,3), "_", round(lat,3)))

idx_map <- match(coord_key_lme$key, coord_key_rea$key)

y_cols <- 3:ncol(tempe_all1)
nu_cols <- 3:ncol(nu_reaches1)
T_common <- min(length(y_cols), length(nu_cols))
y_cols <- y_cols[1:T_common]
nu_cols <- nu_cols[1:T_common]

g0_121 <- sapply(1:n_loc, function(s){
  # local LME sample
  x_s <- get_xs_all_years(lme_df, s)
  Fx_inv_s <- build_Fx_inv_local(x_s)

  # reaches row index for this location
  i_rea <- idx_map[s]

  yhat_s <- as.numeric(tempe_all1[i_rea, y_cols])
  nu_s   <- as.numeric(nu_reaches[i_rea, nu_cols])

  u0 <- FY_hat_1loc(y = 0, yhat = yhat_s, nu = nu_s)
  Fx_inv_s(u0)
})

plot_df <- lme_base %>%
  mutate(temp0_c = g0_121)


#jpeg("./Figure7e.png",width=5,height=4 , res = 300, units = "in")
    print(
      ggplot(plot_df,aes(lon,lat)) +
      geom_point(aes(colour=temp0_c),cex=8.99, shape=15) +
      coord_map(xlim=c(98,124.5),ylim=c(18,42.5)) + 
      scale_color_gradientn(colors = c("blue", "cyan", "green", "yellow", "red"), limits=c(-10,25),na.value="transparent", guide="colourbar") +
      #scale_colour_gradientn(colours=rev(brewer.pal(n=9,name='RdBu')),
      #                       limits=c(-6,25),na.value="transparent",
      #                       guide="colourbar") +
      borders(database="world",xlim=c(76,126),ylim=c(18,45),fill=NA,colour="grey50") +
      theme(text=element_text(size=15),legend.title=element_blank(),
        legend.position="right")
    )
#    dev.off()






  
