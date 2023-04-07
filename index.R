# A zero-dose vulnerability index for equity assessment and 
# spatial prioritization in low anad middle-income countries


# Description:
# This script contains the code used to calculate the vulnerability index
# from this research.


# ----
# PACKAGES

library(stringr)
library(dplyr)
library(purrr)
library(plyr)
library(maptools)
library(raster)
library(rgdal)
library(sf)
library(psych)


# ---
# FUNCTIONS

# fresc() is a function used to rescale the values of a raster to fall 
# between tmin and tmax which is set to as 0 and 100 respectively by 
# default

fresc <- function(xx, tmin = 0, tmax = 100){
  
  # tmin = new min value
  # tmax = new max value
  
  summ <- summary(as.vector(xx))
  a <- c(summ[1])
  b <- c(summ[6])
  resc <- ((xx-a)/(b-a))*(tmax-tmin) + tmin
  return(resc)
  
}


# ---
# DATA

rl <- list()

# Create a list of raster objects. The rasters objects are the high
# resolution prediction surfaces of the zero dose indicators.

rl_s <- lapply(rl, function(xx) fresc(xx))

# Rescale the rasters objects within the list with the fresc() 
# function defined and described above.  

inds <- rl_s


# ---
# VULNERABILITY INDEX 

# Method 1: Direct average 
# This methods involves taking a simple average across the zero-dose 
# indicators within the list object inds.

VI1 <- Reduce("+", inds)*(1/length(inds))
VI1_s <- fresc(VI1)


# Method 2: Group based average
# This method involves first grouping the zero-dose indicators into 
# categories then taking averages within those categories before 
# averaging across the categories.

VH <- inds$TT + inds$anc + inds$sba + inds$pnc + inds$contraception + 
  inds$malaria + inds$traveltime + inds$malnutrition + inds$insurance

VSC <- inds$meduc + inds$decision + inds$hwealth

VGD <- inds$mreligion + inds$ethnicity + inds$urban_rural + 
  inds$conflict + inds$mage + inds$birthq + inds$hsize

# The zero-dose indicators are grouped into the categories: health,
# socio-economic and demographic.

VH <- (1/9)*VH
VSC <- (1/3)*VSC
VGD <- (1/7)*VGD

VI2 <- (1/3)*(VH + VSC + VGD)
VI2_s <- fresc(VI2)


# Method 3: Regression-based weighting
# This method involves weighting the zero-dose indicators based on a 
# binomial regression. This method is split into two: Regression-based
# weighting 1 is the weighting based from the binomial regression and
# Regression-based weighting 2 further divides the weightings into
# quartiles.

cl <- list()
cl <- cl %>% purrr::reduce(inner_join, by = "DHSCLUST")

# Create a list object that contains the cluster level datasets of the 
# zero-dose indicators. Then merge them together as a data fram by the 
# cluster ID.

cl %>% 
  dplyr::select(-contains("DTP3")) %>%
  dplyr::select(-contains("MCV1")) %>%
  dplyr::select(-c("DHSCLUST", "DTP1_prop")) %>% 
  dplyr::rename(TotVax = DTP1_count, TotChild = DTP1_total) %>% 
  dplyr::relocate(TotChild, TotVax) -> dat_DTP1

cl %>% 
  dplyr::select(-contains("DTP1")) %>%
  dplyr::select(-contains("MCV1")) %>% 
  dplyr::select(-c("DHSCLUST", "DTP3_prop")) %>% 
  dplyr::rename(TotVax = DTP3_count, TotChild = DTP3_total) %>% 
  dplyr::relocate(TotChild, TotVax) -> dat_DTP3

cl %>% 
  dplyr::select(-contains("DTP")) %>%
  dplyr::select(-c("DHSCLUST", "MCV1_prop")) %>%
  dplyr::rename(TotVax = MCV1_count, TotChild = MCV1_total) %>% 
  dplyr::relocate(TotChild, TotVax) -> dat_MCV1

# Create three data frame objects where the they contain either DTP1,
# DTP3 or MCV1 and with the other other zero-dose indicators.

source("covrank.R")

# Apply the covrank() function defined and described within the 
# covrank_function.R script inside the covrank folder.

singles_dtp1 <- covrank(dat_DTP1)
singles_dtp3 <- covrank(dat_DTP3)
singles_mcv1 <- covrank(dat_MCV1)

# Regression-based weighting 1

singles_dtp1$rank1 <- nrow(singles_dtp1):1
singles_dtp3$rank1 <- nrow(singles_dtp3):1
singles_mcv1$rank1 <- nrow(singles_mcv1):1

# Regression-based weighting 2

singles_dtp1 <- within(singles_dtp1, rank2 = as.integer(cut(
  singles_dtp1$pr2, quantile(singles_dtp1$pr2, probs = 0:4/4), 
  include.lowest = T))) 

singles_dtp3 <- within(singles_dtp3, rank2 = as.integer(cut(
  singles_dtp3$pr2, quantile(singles_dtp3$pr2, probs = 0:4/4), 
  include.lowest = T))) 

singles_mcv1 <- within(singles_mcv1, rank2 = as.integer(cut(
  singles_mcv1$pr2, quantile(singles_mcv1$pr2, probs = 0:4/4), 
  include.lowest = T))) 

#' Compile a data frame called m to see the weightings.

rank_dtp1 <- singles_dtp1 %>% dplyr::select(cov = model, rank1, rank2)
rank_dtp3 <- singles_dtp3 %>% dplyr::select(cov = model, rank1, rank2)
rank_mcv1 <- singles_mcv1 %>% dplyr::select(cov = model, rank1, rank2)

m <- list(rank_dtp1, rank_dtp3, rank_mcv1) %>% purrr::reduce(left_join, by = "cov")

m$mean_rank1 <- round((m$rank1.x + m$rank1.y + m$rank1)*(1/3),2)
m$mean_rank2 <- round((m$rank2.x + m$rank2.y + m$rank2)*(1/3),2)

m <- m %>% dplyr::select(cov, mean_rank1, mean_rank2)
m$cov <- stringr::str_remove(m$cov, "_prop")

# Create vulnerability index with regression-based weighting 1. 
# The numbers are used as example in this demonstration code.

VI3a <- (
  inds$anc*18.67 +
    inds$sba*17.33 +
    inds$TT*15.67 +
    inds$contraception*17.67 +
    inds$meduc*15.00 +
    inds$hwealth*14.67 +
    inds$pnc*12.67 +
    inds$hsize*10.33 +
    inds$traveltime*10.00 +
    inds$mreligion*11.67 +
    inds$ethnicity*9.33 +
    inds$insurance*7.33 +
    inds$ malnutrition*5.33 +
    inds$decision*8.00 +
    inds$malaria*6.00 +
    inds$urban_rural*3.33 +
    inds$birthq*2.33 +
    inds$mage*3.67 +
    inds$conflict*1.00)/sum(m$mean_rank1)

# Create vulnerability index with regression-based weighting 2. 
# The numbers are used as example in this demonstration code.

VI3b <- (
  inds$anc*4 +
    inds$sba*4 +
    inds$TT*3.67 +
    inds$contraception*4 +
    inds$meduc*4 +
    inds$hwealth*3.33 +
    inds$pnc*3 +
    inds$hsize*2.33 +
    inds$traveltime*2.67 +
    inds$mreligion*2.67 +
    inds$ethnicity*2.00 +
    inds$insurance*2 +
    inds$ malnutrition*1.33 +
    inds$decision*2.33 +
    inds$malaria*1.67 +
    inds$urban_rural*1 +
    inds$birthq*1 +
    inds$mage*1 +
    inds$conflict*1.00)/sum(m$mean_rank2)

VI3a_s <- fresc(VI3a)
VI3b_s <- fresc(VI3b)


# Method 4. Factor analysis
# This method involves weighting the zero-dose indicators based on the 
# results of a factor analysis.

cl %>% 
  dplyr::select(-contains("DTP")) %>% 
  dplyr::select(-contains("MCV")) %>%
  dplyr::select(-DHSCLUST) -> covars

# Calculate the principal components.

covars_fa <- apply(covars, 2, function(x) fresc(x))
pr <- prcomp(covars_fa, scale = T)

pr_var <- pr$sdev^2
pve <- pr_var/sum(pr_var)

# Create a scree plot.

qplot(c(1:19), pve) + 
  geom_line() + 
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  scale_x_continuous(breaks = 1:19) + 
  ylim(0,0.5) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Perform factor analysis based on the results of the scree plot.

cfa <- psych::fa(covars_fa, 3, fm = "pa", covar = F, rotate = "varimax")
ll <- apply(as.matrix(cfa$loadings), 2, function(x) x^2)
ll <- apply(ll, 2, function(x) x/sum(x))
colSums(ll)

# Keep factors whose squared and normalized loadings > 0.1. If factors 
# appear on both groups retain the larger one.

as.data.frame(round(ll, 4)) %>% 
  filter(PA1 > 0.10 | PA2 > 0.10 | PA3 > 0.1) %>% 
  round(.,2)

# Assign the weights to the zero-dose indicators based on the results
# of the factor analysis.

G1 <- (inds$hwealth*0.17 + inds$meduc*0.12 + inds$sba*0.16)/(0.17 + 0.12 + 0.16)
G2 <- (inds$ethnicity*0.19 + inds$mreligion*0.34 + inds$contraception*0.17)/(0.19 + 0.34 + 0.17)
G3 <- (inds$anc*0.25 + inds$pnc*0.23 + inds$TT*0.13)/(0.25 + 0.23 + 0.13)

VIfa <- (G1*0.55 + G2*0.23 + G3*0.22)/(0.55 + 0.23 + 0.22)
VIfa_s <- fresc(VIfa)
