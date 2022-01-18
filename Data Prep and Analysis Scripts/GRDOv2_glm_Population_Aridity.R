#############################################
## Blaszczak et al. - Global extent, patterns, and drivers of hypoxia in rivers
## Analysis - watershed population density and aridity index glms
## Code author: J.R. Blaszczak
#############################################
## Load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate","boot",
         "tidyverse", "readxl","ggExtra", "data.table"), require, character.only=T) #any FALSE need to be installed
theme_set(theme_bw())

##########################
## Import Hypoxia Data
##########################
## Set wd again
#setwd("../Distribution, frequency, and global extent of hypoxia in rivers/Data")
## Import sites
dat <- fread("GRDO_GEE_HA_NHD.csv")

## Add in hypoxia (<2 mg/L) as a binomial
dat$Hyp_01 <- 0
dat[which(dat$Hyp_pr_sub2 > 0),]$Hyp_01 <- 1

#############################################
## Population density
############################################
df_pop <- na.omit(dat[,c("Hyp_01","NHD_PopDen2010Ws")])

test <- glm(Hyp_01 ~ NHD_PopDen2010Ws, df_pop, family = "binomial")
summary(test)

hist(df_pop$NHD_PopDen2010Ws)
pop_dat <- seq(from=quantile(df_pop$NHD_PopDen2010Ws, na.rm = T, probs = 0.025),
               to=quantile(df_pop$NHD_PopDen2010Ws, na.rm = T, probs = 0.975), by=10)
min(pop_dat); max(pop_dat) #0.23 - 1060

# predictions
pred <- rep(NA, times=length(pop_dat))
for (j in 1:length(pop_dat)){
  pred[j] <- pop_dat[j]*test$coefficients[2] + test$coefficients[1]
}

#inv.logit
z <- exp(pred)/(1+exp(pred))

plot(pop_dat, z)
max(z); max(pop_dat)
min(z); min(pop_dat)
max(z) - min(z)#6.4%

####################
## Aridity
#################
df_arid <- na.omit(dat[,c("Hyp_01","ari_ix_uav")])

arid.test <- glm(Hyp_01 ~ ari_ix_uav, df_arid, family = "binomial")
summary(arid.test)

hist(df_arid$ari_ix_uav)
arid_dat <- seq(from=quantile(df_arid$ari_ix_uav, na.rm = T, probs = 0.025),
                to=quantile(df_arid$ari_ix_uav, na.rm = T, probs = 0.975), by=1)
min(arid_dat);max(arid_dat) #34-208

# predictions
pred <- rep(NA, times=length(arid_dat))
for (j in 1:length(arid_dat)){
  pred[j] <- arid_dat[j]*arid.test$coefficients[2] + arid.test$coefficients[1]
}

#inv.logit
a <- exp(pred)/(1+exp(pred))

plot(arid_dat, a)
max(a); max(arid_dat)
min(a); min(arid_dat)
max(a) - min(a) #3.5%



