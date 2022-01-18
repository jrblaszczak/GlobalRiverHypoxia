##################################################################################
## Blaszczak et al. - Global extent, patterns, and drivers of hypoxia in rivers
## Figure S - Sensor data trends through time hypoxia
## Code author: J.R. Blaszczak
#################################################################################
## Load packages
lapply(c("plyr","dplyr","ggplot2","cowplot",
         "lubridate","tidyverse","data.table",
         "stringr","zoo","wesanderson","trend",
         "rgdal","sp","sf","maptools","ggmap","maps","mapdata",
         "openintro"), require, character.only=T) ## if false package needs to be installed
theme_set(theme_bw())

####################################################
## Compile temporal data for sensors
####################################################
## Setwd to original data files
setwd("../../Organized GRDO Database for Publication/GRDO 2. Template ts with sumstat scripts/Year Count Files 2021")

YC <- ldply(list.files(pattern = "csv"), function(filename) {
  d <- read.csv(filename)
  d$file <- filename
  return(d)
})

## Fix column names
colnames(YC)[which(colnames(YC) == ".id")] <- "Year"

## Split file name
YC$file <- str_remove(YC$file, pattern = "YearCount_")
YC$file <- str_remove(YC$file, pattern = "_TS.csv")
YC$file <- str_remove(YC$file, pattern = "_formatted")
YC$string <- substr(YC$file,1,4)

## Specify separate dataframes
PC_sensor <- YC[which(YC$string %in% c("nwis")),]

## Calculate the percentage of meas. per year < 2 mg/L
PC_sensor$pct_sub2 <- PC_sensor$Count_sub2/PC_sensor$Count_Tot

################################################
## Split and calculate Sen slope per location
################################################

## Calculate sens.slope per sensor location
l <- split(PC_sensor, PC_sensor$file)

## only calculate sens.slope if three years
list.lengths <- ldply(lapply(l, function(x) nrow(x)), data.frame)
colnames(list.lengths) <- c("SiteID","years")
long_sites <- list.lengths[which(list.lengths$years >= 3),]
long_l <- l[names(l) %in% levels(as.factor(long_sites$SiteID))]

## calc sens.slope
senslopes_PC <- lapply(long_l, function(x) sens.slope(x$pct_sub2, conf.level = 0.95))
senslopes_PC_p <- ldply(lapply(senslopes_PC, function(x) x$p.value), data.frame); colnames(senslopes_PC_p) <- c("SiteID","p.value")
senslopes_PC_slope <- ldply(lapply(senslopes_PC, function(x) x$estimates), data.frame); colnames(senslopes_PC_slope) <- c("SiteID","Sens.slope")

ss <- left_join(senslopes_PC_p, senslopes_PC_slope, by="SiteID")

## summary stats
nrow(ss[is.na(ss$p.value),])/nrow(ss) ## 56% no trend because p-value = NA
nrow(ss[which(ss$Sens.slope == 0),])/nrow(ss) ## 77% no trend because slope = 0
nrow(ss[which(ss$p.value <= 0.05),]) ## only 10 sites with a p-value < 0.05

## visulize distribution of slopes
ss_trend <- ss[which(ss$p.value <= 0.1),]
ggplot(ss_trend, aes(Sens.slope))+
  geom_histogram()

############################
## Import metadata and map
############################
## reset wd
setwd("../../../Extra scripts for revisions/Distribution, frequency, and global extent of hypoxia in rivers/Data")

## import site data
site_info <- fread("GRDO_GEE_HA_NHD.csv")
head(site_info)

## merge mean duration data with site_info
trend_map <- left_join(ss, site_info, by="SiteID")

## Categorize Sens.slopes
trend_map$ss_cat <- "No trend" #if a p-value > 0, set the category to 0
trend_map[which(trend_map$Sens.slope == 0),]$ss_cat <- "No trend"
trend_map[which(trend_map$p.value <= 0.1 & trend_map$Sens.slope > 0),]$ss_cat <- "Positive"
trend_map[which(trend_map$p.value <= 0.1 & trend_map$Sens.slope < 0),]$ss_cat <- "Negative"
nrow(trend_map[which(trend_map$ss_cat == "No trend"),])/nrow(trend_map) ## 96% of sites zero trend
View(trend_map[which(trend_map$ss_cat %in% c("Positive","Negative")),])

nrow(trend_map[which(trend_map$ss_cat == "Positive"),])
nrow(trend_map[which(trend_map$ss_cat == "Negative"),])

##map
ss_map <- ggmap(get_stamenmap(bbox=c(-125, 25, -66, 50), zoom = 5, 
                    maptype='toner'))+
  geom_point(data = trend_map, aes(x = Longitude, y = Latitude, 
                                 fill=ss_cat, size=ss_cat), shape=21, alpha=0.9)+
  scale_fill_manual("Sen's Slope", values = c("Positive" = "red",
                                               "No trend" = "grey",
                                               "Negative" = "blue"))+
  scale_size_manual("Sen's Slope", values = c("Positive" = 5,
                                              "No trend" = 2,
                                              "Negative" = 5))+
  theme(legend.position = "right")+
  labs(x="Longitude", y="Latitude")
ss_map

## trends for sites with a positive or negative slope
t <- trend_map[which(trend_map$ss_cat %in% c("Positive","Negative")),]
t_l <- ldply(l[names(l) %in% levels(as.factor(t$SiteID))], data.frame)
colnames(t_l)[1] <- "SiteID"
t_l <- merge(t_l, trend_map, by="SiteID")

## plot trends only for those with p < 0.1 and non-zero slope
ss_lines <- ggplot(t_l, aes(Year, pct_sub2, group=SiteID, color=ss_cat))+
  geom_line()+ geom_point()+
  theme(legend.position = "right")+
  scale_color_manual("Sen's Slope", values = c("Positive" = "red",
                                              "Negative" = "blue"))+
  labs(x="Year",y="% Hypoxic Observations Per Year")+
  scale_y_continuous(labels = function(x) paste0(x * 100, '%'))
ss_lines

## plot together
plot_grid(
  ss_map, ss_lines, ncol=1
)


ggplot(t_l[which(t_l$pct_sub2 > 0.20),], aes(Year, pct_sub2, color=SiteID))+
  geom_line()+ geom_point()+
  theme(legend.position = "right")

## #1 Top line: Deep Creek
## #2 Top line: Hillsborough River




