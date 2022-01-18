#####################################################################################
## Blaszczak et al. - Global extent, patterns, and drivers of hypoxia in rivers
## Figure SI - Date Range of data figure
## Code author: J.R. Blaszczak
######################################################################################

## Load packages
lapply(c("plyr","dplyr","ggplot2","ggspatial","cowplot",
         "tidyverse","sf","sp","rnaturalearth","rnaturalearthdata","raster",
         "ggmap","data.table","RColorBrewer","wesanderson", "ggExtra"), 
       require, character.only=T) #any FALSE need to be installed
theme_set(theme_bw())

##########################
## Import Hypoxia Data
##########################
## Setwd
setwd("../Distribution, frequency, and global extent of hypoxia in rivers/Data")
## Import sites
dat <- fread("GRDO_GEE_HA_NHD.csv")
head(dat)

##########################
## Subset and summarize
#########################

sub <- dat[,c("DO_mgL_mean","Start_time","End_time","Country")]
sub$Start_Yr <- year(as.POSIXct(as.character(sub$Start_time), format="%m/%d/%Y %H:%M"))
sub$End_Yr <- year(as.POSIXct(as.character(sub$End_time), format="%m/%d/%Y %H:%M"))

country_se <- sub %>%
  group_by(Country)%>%
  summarise_at(.vars = c("Start_Yr","End_Yr"), .funs = c(min, max))
colnames(country_se) <- c("Country","Start_Yr_min","End_Yr_min","Start_Yr_max","End_Yr_max")

c <- country_se[,c("Country","Start_Yr_min","End_Yr_max")]

## reorder based on start date
library(forcats)
c %>%
  mutate(name = fct_reorder(Country, Start_Yr_min)) %>%
  ggplot()+
  geom_segment( aes(x=name, xend=name, y=Start_Yr_min, yend=End_Yr_max), color="black") +
  geom_point( aes(x=name, y=Start_Yr_min), color="blue", size=3 ) +
  geom_point( aes(x=name, y=End_Yr_max), color="purple", size=3 ) +
  coord_flip()+
  theme(axis.text.y = element_text(size=8))+
  labs(y="Year",x="Country")+
  scale_y_continuous(breaks=c(1900,1925,1950,1975,2000,2020))











