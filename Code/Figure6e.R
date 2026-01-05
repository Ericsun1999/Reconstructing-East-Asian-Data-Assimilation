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

