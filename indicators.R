# TITLE OF THE PAPER GOES HERE

# Description:
# This script contains the code used to fit the zero-dose indicators and
# construct the high-resolution prediction and uncertainty surfaces.


# ---
# PACKAGES

library(dplyr)
library(INLA)
library(sf)
library(raster)
library(gtools)


# ---
# DATA

dat <- read.csv()

# The dataset must contain columns titled LAT and LON corresponding to the 
# survey cluster coordinates. The dataset must contain columns titled count 
# and total corresponding to the counts of events and total number surveyed
# within the survey clusters.

finla1 <- formula(y ~ int + .)

# Define a formula() object that contains the covariates to be fitted in 
# the model later within INLA. The formula starts with int which denotes 
# the intercept. The preceding . denotes the covariates to be included.



# ---
# INLA

# Fitting a spatial model in INLA requires constructing a mesh.

mesh <- inla.mesh.2d(
  loc = cbind(dat$LON, dat$LAT),
  max.edge = c(0.3,5),
  cutoff = 0.3)

# Fitting a spatial model in INLA with the SPDE approach requires defining 
# the SPDE object. Supply r0 as an argument within the SPDE object. r0 is 
# 5% of the extent of the survey cluster locations.

r0 <- 0.05*pmax(
  abs(max(dat$LON) - min(dat$LON)),
  abs(max(dat$LAT) - min(dat$LAT)))

spde <- inla.spde2.pcmatern(
  alpha = 3/2, mesh = mesh, 
  prior.range = c(r0, 0.01),
  prior.sigma = c(1, 0.05))

s_index <- inla.spde.make.index("s", spde$n.spde)
r_index <- 1:nrow(dat)

# Fitting a spatial model in INLA requires creating an A object.

A_est <- inla.spde.make.A(
  loc = as.matrix(cbind(dat$LON, dat$LAT)),
  mesh = mesh)

# The INLA stack object stacks the pre-requisites and data together in a 
# format that the INLA function can use to fit the model.

stack_est <- inla.stack(
  tag = "est", 
  A = list(A_est, 1, 1),
  data = list(y = dat$count, n = dat$total),
  effects = list(s = s_index, r = r_index, data.frame(intercept = 1, dat)))

stack_full <- inla.stack(stack_est)

# The previously defined formula in the finla1 object is appended to include 
# the SPDE component and an iid component. Additionally, priors are specified
# for the model components.

hp <- list(theta = list(prior = "loggamma", param = c(2,1)))
finla2 <- update(finla1, . ~ . -1 + f(s, model = spde) + f(r, model = "iid", hyper = hp))

# Fit the model with the INLA function

imod <- inla(
  formula = finla,
  family = "binomial",
  data = inla.stack.data(stack_full),
  Ntrials = stack_full$data$data$n,
  control.family = list(link = "logit"),
  control.predictor = list(compute = T, link = 1, A = inla.stack.A(stack_full)),
  control.compute = list(dic = T, waic = T, config = T),
  verbose = F)


# ---
# PREDICTION

# Draw posterior samples of the estimated parameters from the 
# INLA  model object with the inla.posterior.sample() function.

nsamp <- 1000
contents <- imod$misc$configs$contents
postsamp <- inla.posterior.sample(nsamp, imod)

id_s <- which(contents$tag == "s")
id_r <- which(contents$tag == "r")
id_x <- which(contents$tag == "int")

index_S <- contents$start[id_s]-1 + (1:contents$length[id_s])
index_R <- contents$start[id_r]-1 + (1:contents$length[id_r])
index_X <- contents$start[id_x]-1 + (1:ncol(imod$model.matrix))

# Construct a data frame with columns from values of the covariate 
# rasters for prediction. Also include the prediction coordinates.

p_dat <- 1
for(ii in 2:length(imod$names.fixed)){
  
  r_vals <- scale(raster::values(raster(list.files(
    path = "rasters/",
    pattern = imod$names.fixed[ii],
    full.names = T))))
  
  p_dat <- cbind(p_dat, r_vals)
  
}
colnames(p_dat) <- imod$names.fixed
p_coord <- coordinates(raster(list.files(
  path = "rasters/", 
  pattern = imod$names.fixed[2], 
  full.names = T)))
p_dat <- cbind(p_coord, p_dat)

# To reduce computation burden, only predict for no NA 
# values from the constructed prediction data frame.

id <- apply(p_dat, 1, function(x) any(is.na(x)))
id_missing <- which(id == T)
id_nonmiss <- which(id == F)

p_dat_nonmiss <- p_dat[id_nonmiss,]
p_coord_nonmiss <- p_dat_nonmiss[,1:2]

# Create an INLA A object for prediction.

A_pred = inla.spde.make.A(
  loc = as.matrix(p_coord_nonmiss),
  mesh = mesh)

p_mat <- as.matrix(p_dat_nonmiss[,3:ncol(p_dat_nonmiss)])

xLat = matrix(0, nrow = length(postsamp[[1]]$latent), ncol = nsamp)
xHyp = matrix(0, nrow = length(postsamp[[1]]$hyperpar), ncol = nsamp)
sIID = matrix(0, nrow = nrow(p_mat), ncol = nsamp)

for(j in 1:nsamp){
  xLat[,j] = postsamp[[j]]$latent
  xHyp[,j] = postsamp[[j]]$hyperpar
} 

for(j in 1:nsamp){
  post_prec = xHyp[3,j]
  post_sigma = post_prec^(-0.5)
  sIID[,j] = rnorm(nrow(p_mat), sd = post_sigma)
} 

xS = xLat[index_S,]
xX = xLat[index_X,]

# Prediction calculations

inv_lpred <- gtools::inv.logit(as.matrix(A_pred %*% xS + p_mat %*% xX + sIID))

as.data.frame(t(apply(inv_lpred, 1, function(x) 
  c(mean(x), sd(x))))) %>%
  setNames(c("mean", "sd")) -> p_inla

# Write the prediction (mean) and uncertainty (sd) outputs to a raster.

tmp <- raster(list.files("rasters/", pattern =  imod$names.fixed[2], full.names = T)) 
  
values(tmp)[id_nonmiss] <- p_inla$mean
values(tmp)[id_missing] <- NA
writeRaster(tmp, filename = "mean.tif")
  
tmp <- raster(list.files("rasters/", pattern =  imod$names.fixed[2], full.names = T)) 

values(tmp)[id_nonmiss] <- p_inla$sd
values(tmp)[id_missing] <- NA
writeRaster(tmp, filename = "sd.tif")
