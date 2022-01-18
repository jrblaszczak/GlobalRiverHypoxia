#####################################################################################
## Blaszczak et al. - Global extent, patterns, and drivers of hypoxia in rivers
## Figure 1 - Global map of hypoxia in river basins & data distribution
## Code author: J.R. Blaszczak
######################################################################################

## Load packages
lapply(c("plyr","dplyr","ggplot2","ggspatial","cowplot",
         "tidyverse","sf","sp","rnaturalearth","rnaturalearthdata","raster",
         "ggmap","data.table","RColorBrewer","wesanderson", "ggExtra"), 
       require, character.only=T) #any FALSE need to be installed
theme_set(theme_bw())

##################
## Import Basins
##################
## Download BasinATLAS_v10_shp for level 4 here: https://figshare.com/articles/dataset/HydroATLAS_version_1_0/9890531?file=20087237
## Reset wd
setwd("../Distribution, frequency, and global extent of hypoxia in rivers/Map Data/BasinATLAS_v10_lev04")
## Import
B <- st_read("BasinATLAS_v10_lev04.shp")
head(B)
st_crs(B)

##########################
## Import Hypoxia Data
##########################
## Reset wd again
setwd("../../Data")
## Import sites
dat <- fread("GRDO_GEE_HA_NHD.csv")

t <- fread("Temporal_Data_YearlySum.csv")
t <- t[which(t$Data_Subset == "Global_sum"),]

##########################
## Identify River Basin
##########################
#https://mgimond.github.io/Spatial/point-pattern-analysis-in-r.html

# Convert pointsDF to a SpatialPoints object 
pointsDF <- dat[,c("Lon_WGS84","Lat_WGS84")]
pointsSP <- SpatialPoints(pointsDF, proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs"))

# Convert B to SpatialPolygon
B.sp <- as(B, "Spatial")
B_sp <- geometry(B.sp)
projection(B_sp); projection(pointsSP)

# Use 'over' to get _indices_ of the Polygons object containing each point 
RiverBasin_indices <- over(pointsSP, B_sp)

# Return the RiverBasin of the Polygons object containing each point
RiverBasinNames <- sapply(B_sp@polygons, function(x) x@ID)
dat$HYBAS_ID <- B.sp$HYBAS_ID[RiverBasin_indices]

## Aggregate the number of points per Basin
## First summarise the data
dat_byBasin <- dat %>%
  group_by(HYBAS_ID) %>%
  tally()
colnames(dat_byBasin) <- c("HYBAS_ID", "Num_sites")

###############
## Map Basins
###############
## Join B and dat
Bdat <- left_join(B, dat_byBasin, by="HYBAS_ID") 

## Hypoxia prevalence - First summarize the data
dat_byBasin_sub2 <- dat[which(dat$Hyp_pr_sub2 > 0),] %>%
  group_by(HYBAS_ID) %>%
  tally()
colnames(dat_byBasin_sub2) <- c("HYBAS_ID", "Num_sub2_sites")

## Then join the two
Bdat <- left_join(Bdat, dat_byBasin_sub2, by="HYBAS_ID")
Bdat$percent_sub2 <- (Bdat$Num_sub2_sites/Bdat$Num_sites)*100
#Make sure that sites without any measurements are NA
nrow(Bdat[is.na(Bdat$Num_sub2_sites) == FALSE,]) #228
nrow(Bdat[is.na(Bdat$Num_sites) == FALSE,]) #405
Bdat$Num_sub2_sites <- ifelse(Bdat$Num_sites>=0 & is.na(Bdat$Num_sub2_sites)==TRUE,
                              yes=0, no= Bdat$Num_sub2_sites)
nrow(Bdat[is.na(Bdat$Num_sub2_sites) == FALSE,]) #Check - 405, problem solved

## Visualize distribution of sites
ggplot(Bdat, aes(percent_sub2))+geom_histogram(alpha=0.5, bins = 30)
ggplot(Bdat, aes(Num_sub2_sites))+geom_histogram(alpha=0.5, bins = 30)

## World map of hypoxic sites -- FIGURE USED IN MANUSCRIPT
# set color palette
pal <- wes_palette("Zissou1", type = "continuous")

## World map of hypoxic sites
world <- ggplot(data = Bdat) +
  geom_sf(aes(fill=Num_sub2_sites+0.1))+
  coord_sf(xlim = c(-180,180), ylim = c(-60,85), expand = FALSE)+
  xlab("Longitude") + ylab("Latitude")+
  scale_fill_gradientn("Number of Locations within\n Basin with Detected Hypoxia\n (<2 mg/L DO)",
                       colours=pal, na.value = "gray30",
                       breaks=c(0+0.1, 10+0.1, 100+0.1, 10000+0.1),
                       labels=c("0", "10", "100","10,000"),
                       trans="log", limits = c(0+0.1,10000+0.1))+
  theme(panel.background = element_rect(fill = "white", color="black"),
        panel.grid = element_line(color = "gray85", linetype = "dashed", size = 0.2),
        legend.position = "top",
        axis.text = element_text(size=15), axis.title = element_text(size=17))
world+theme(legend.position = "none")
plot_grid(get_legend(world))

## How many river basins had no measurements
nrow(Bdat[is.na(Bdat$Num_sites) == FALSE,]) #405
nrow(Bdat[is.na(Bdat$Num_sites) == TRUE,]) #937
nrow(Bdat[is.na(Bdat$Num_sites) == TRUE,])/(nrow(Bdat[is.na(Bdat$Num_sites) == TRUE,])+nrow(Bdat[is.na(Bdat$Num_sites) == FALSE,]))

## SI figure - world map of proportion of hypoxia relative to total measurements
world_SI <- ggplot(data = Bdat) +
  geom_sf(aes(fill=percent_sub2+0.01)) +
  coord_sf(xlim = c(-180,180), ylim = c(-60,85), expand = FALSE)+
  xlab("Longitude") + ylab("Latitude")+
  scale_fill_gradientn("% of Locations within\n Basin with Detected Hypoxia\n (<2 mg/L DO)",
                       colours=pal, na.value = "gray30",
                       breaks=c(0+0.01, 1+0.01, 10+0.01, 100+0.01),
                       labels=c("0%", "1%", "10%","100%"),
                       trans="log", limits = c(0+0.01,100+0.01))+
  theme(panel.background = element_rect(fill = "white", color="black"),
        panel.grid.major = element_line(color = "gray85", linetype = "dashed", size = 0.5),
        legend.position = "top",
        axis.text = element_text(size=15), axis.title = element_text(size=17))
world_SI+theme(legend.position = "none")
plot_grid(get_legend(world_SI))




#####################################################################
## Bottom Panels - Distribution of hypoxia across mean/sd
#####################################################################
## Order sites
dat <- dat[order(dat$DO_mgL_mean),]
dat$Order <- 1:nrow(dat)

## separating SD, mean, and proportion of hypoxic - USED IN MANUSCRIPT
nrow(dat[which(dat$n_time >1),]); nrow(dat[which(dat$n_time >1),])/nrow(dat) #94,645; 0.75
nrow(dat[which(dat$n_time >=3),]); nrow(dat[which(dat$n_time >=3),])/nrow(dat) #84,704; 0.67

## Define loess span for the sd of DO and the proportion of hypoxic measurements 
loessMod10_sd1 <- loess((15-DO_mgL_sd)/15 ~ DO_mgL_mean, data=dat[which(dat$n_time >=3),], span=0.01)
smoothed10_sd1 <- predict(loessMod10_sd1) 
loessMod10_pr1 <- loess(Hyp_pr_sub2 ~ DO_mgL_mean, data=dat[which(dat$n_time >=3),], span=0.01)
smoothed10_pr1 <- predict(loessMod10_pr1) 
## color palette
pal <- wes_palette("Zissou1", type = "continuous")

## loess with predefined span
mean_vs_sd_plot <- ggplot(dat[which(dat$n_time >=3),], aes(DO_mgL_mean,(15-DO_mgL_sd)/15))+
  geom_point(color="gray48",alpha=0.5, size=0.7)+
  geom_line(aes(y=smoothed10_sd1), color="gray25",size=1.1)+
  scale_y_continuous(limits=c(0,1), sec.axis = sec_axis(~rev(.*15)))+ #, name=expression("Within-site SD of [DO] (mg/L)")))+
  geom_point(aes(DO_mgL_mean, Hyp_pr_sub2,color=Hyp_pr_sub2+0.01),size=0.7)+
  scale_x_continuous(limits=c(0,15))+
  labs(x="Within-site Mean [DO] (mg/L)",y="Proportion of [DO] < 2 (mg/L)")+
  theme(axis.title = element_text(size=16),
        axis.text.x = element_text(size=14, hjust=1, angle=45),
        axis.text.y = element_text(size=14),
        axis.title.y.left = element_text(color="black"),
        panel.background = element_rect(fill="white", colour = "black"),
        #legend.position = "none",
        panel.grid.major = element_line(color="gray86"))+
  guides(trans="log", color=guide_legend(title="Proportion of\n[DO] < 2 (mg/L)"))+
  scale_color_gradientn(colours=pal, breaks=c(0.01+0.01, 0.1+0.01, 0.5+0.01, 1+0.01),
                        labels=c(0.01, 0.1, 0.5, 1),
                        trans="log", limits = c(0+0.01,1+0.01))+
  geom_line(aes(y=smoothed10_pr1), color="black",size=1.1)

## Add marginal distribution plots
ggExtra::ggMarginal(mean_vs_sd_plot+theme(legend.position = "none"), type="density", size=5,
                    xparams = list(fill="black", color="black"),
                    yparams = list(fill="gray75", color="gray25"))

plot_grid(get_legend(mean_vs_sd_plot+theme(legend.position = "left")))

## calcs
max(dat[which(dat$n_time >=3),]$DO_mgL_sd, na.rm=T) #14.5
nrow(dat[which(dat$n_time >=3 & dat$DO_mgL_mean >=15),]); 62/nrow(dat[which(dat$n_time >=3),]) #62; <0.1%


#####################################
## Temporal versus spatial variation
#####################################
df2 <-dat[which(dat$n_time >=3),]
df2 <- na.omit(df2[,c("DO_mgL_mean","DO_mgL_sd")])
df2$DO_mgL_CV <- df2$DO_mgL_sd/df2$DO_mgL_mean

## spatial variation
sd(df2$DO_mgL_mean)
mean(df2$DO_mgL_mean)
sd(df2$DO_mgL_mean)/mean(df2$DO_mgL_mean) # 0.27

## temporal variation
hist(df2$DO_mgL_CV)
sd(df2$DO_mgL_CV, na.rm = TRUE)
mean(df2$DO_mgL_CV, na.rm= TRUE) # 0.25

