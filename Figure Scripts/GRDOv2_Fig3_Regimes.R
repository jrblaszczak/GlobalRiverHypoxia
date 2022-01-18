#####################################################################################
## Blaszczak et al. - Global extent, patterns, and drivers of hypoxia in rivers
## Figure 3 - Compile hypoxic events and figure with maps and distributions
## Code author: J.R. Blaszczak
######################################################################################

## Load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "tidyverse","readxl","rMR","data.table", "filesstrings"), require, character.only=T)

## set wd
setwd("../../Organized GRDO Database for Publication/GRDO 2. Template ts with sumstat scripts/Hypoxic Events PC")

## Move empty files with no events
info <- file.info(list.files(pattern = ".csv"))
empty <- rownames(info[info$size == 65, ])
file.move(empty, destinations = "./No hypoxic event", overwrite = TRUE)

## Import files
PC_event <- ldply(list.files(pattern = "csv"), function(filename) {
  d <- read.csv(filename, header=T)
  d$file <- filename
  return(d)
})

head(PC_event)

## remove negative durations caused by data gaps
## remove individual point measurements
PC <- PC_event[which(PC_event$event_dur_min > 0),]

## threshold as factor
PC$DO_thresh <- as.factor(PC$DO_thresh)
PC$DO_Thresh <- revalue(PC$DO_thresh, 
                        c("1" = "[DO] < 1 (mg/L)",
                          "2" = "[DO] < 2 (mg/L)",
                          "3" = "[DO] < 3 (mg/L)",
                          "4" = "[DO] < 4 (mg/L)",
                          "5" = "[DO] < 5 (mg/L)"))

## minutes to hours
PC$event_dur_hour <- PC$event_dur_min/60

## Quick stats for DO < 2mg/L
PC2 <- PC[which(PC$DO_thresh == 2),]
mean(PC2$event_dur_hour) ## 16 hours
median(PC2$event_dur_hour) ## 3 hours
max(PC2$event_dur_hour)/(24*30) ## 2432 hours, 101 days, 3.4 months
PC2[which(PC2$event_dur_hour == max(PC2$event_dur_hour)),] ## nwis_02266495, 

nrow(PC2[which((PC2$event_dur_hour/(24)) >= 7),])/nrow(PC2) #2%

####################
## visualize
####################
theme_set(theme_bw())

## Color scheme
colfunc <- colorRampPalette(c("red", "gray45"))
redpal <- colfunc(5)

ggplot(PC, aes(event_dur_min, fill=DO_Thresh))+
  geom_histogram(position = "identity", alpha=0.3, color="black", binwidth = 0.1)+
  scale_x_continuous(trans="log10", breaks = c(5, 15, 60, 60*6, 60*24,
                                               60*24*7, 60*24*30, 60*24*30*6),
                     labels = c("5 min","15 min", "1 hour", "6 hours",
                                "1 day","1 week", "1 month","6 months"))+
  scale_fill_manual(name="Hypoxia Threshold", values = c("[DO] < 1 (mg/L)"=redpal[1],
                                                         "[DO] < 2 (mg/L)"=redpal[2],
                                                         "[DO] < 3 (mg/L)"=redpal[3],
                                                         "[DO] < 4 (mg/L)"=redpal[4],
                                                         "[DO] < 5 (mg/L)"=redpal[5]))+
  theme(panel.grid.major.y = element_line(color="gray85"),
        axis.title = element_text(size=18),
        title = element_text(size=18),
        axis.text.x = element_text(size=14, angle=45, hjust = 1),
        legend.position = "none")+
  labs(x="Event duration", y="Number of events")


## Groups
# <=1 hour, <=1 day, <=1 week, <= 1 month
# changing anything less than one hour to a group

quantiles<-c(60, 60*6, 60*12, 60*24,60*24*7, 60*24*30, 60*24*30*6)
PC$quant <- factor(findInterval(PC$event_dur_min,quantiles)+1)
PC$quant_val <- revalue(PC$quant, c("1" = "\u2264 1 hr",
                                    "2" = "1 hr to 6 hrs",
                                    "3" = "6 hrs to 12 hrs",
                                    "4" = "12 hrs to 1 day",
                                    "5" = "1 day to 1 week",
                                    "6" = "1 week to 1 month",
                                    "7" = "1 month to 6 months"))
## Plot

PC$DO_Thresh <- factor(PC$DO_Thresh, levels = c("5" = "[DO] < 5 (mg/L)",
                                                "4" = "[DO] < 4 (mg/L)",
                                                "3" = "[DO] < 3 (mg/L)",
                                                "2" = "[DO] < 2 (mg/L)",
                                                "1" = "[DO] < 1 (mg/L)"))


## Event duration
event_dur_hist <- ggplot(PC, aes(quant_val, fill=DO_Thresh))+
  geom_bar(alpha=0.4, color="black", position="identity")+
  scale_fill_manual(name="Hypoxia Threshold", values = c("[DO] < 1 (mg/L)"=redpal[1],
                                                         "[DO] < 2 (mg/L)"=redpal[2],
                                                         "[DO] < 3 (mg/L)"=redpal[3],
                                                         "[DO] < 4 (mg/L)"=redpal[4],
                                                         "[DO] < 5 (mg/L)"=redpal[5]))+
  theme(panel.grid.major.y = element_line(color="gray85"),
        axis.title = element_text(size=12),
        axis.text.x = element_text(size=12, angle=35, hjust = 1),
        axis.text.y = element_text(size=12),
        legend.position = "none")+
  labs(x="Event duration", y="Number of events")
event_dur_hist

##############################################
## What time of day is the onset or end?
############################################
head(PC)
PC$start_hour <- hour(as.POSIXct(PC$start_date, format="%Y-%m-%d %H:%M:%S"))
PC$end_hour <- hour(as.POSIXct(PC$end_date, format="%Y-%m-%d %H:%M:%S"))

PC_sub2 <- PC[which(PC$DO_thresh == 2),]
quant_val_sub_list <- c("1 hr to 6 hrs", "6 hrs to 12 hrs", "12 hrs to 1 day")

time_of_day <- ggplot(PC_sub2[which(PC_sub2$quant_val %in% quant_val_sub_list),], aes(start_hour))+
  geom_bar(fill="#010D26", alpha=0.4, color="black")+
  geom_bar(aes(end_hour), fill="#4CBFBB", alpha=0.5, color="black")+
  labs(x="Hour of Day", y="Number of Events")+
  facet_wrap(~quant_val, nrow=1)+
  scale_x_continuous(breaks=c(0,4,8,12,16,20,24))+
  theme(panel.grid.major.y = element_line(color="gray85"),
        axis.title = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        legend.position = "none",
        strip.background = element_rect(fill="white"))
time_of_day  

## All for SI
ggplot(PC_sub2, aes(start_hour))+
  geom_bar(fill="#010D26", alpha=0.4, color="black")+
  geom_bar(aes(end_hour), fill="#4CBFBB", alpha=0.5, color="black")+
  labs(x="Hour of Day", y="Number of Events")+
  facet_wrap(~quant_val, nrow=4)+
  scale_x_continuous(breaks=c(0,4,8,12,16,20,24))+
  theme(panel.grid.major.y = element_line(color="gray85"),
        axis.title = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        legend.position = "none",
        strip.background = element_rect(fill="white"))


#################################################
## Mean duration of hypoxia map
##################################################
##  reload and load more packages
lapply(c("plyr","dplyr","ggplot2","cowplot",
         "lubridate","tidyverse", "reshape2",
         "plotrix", "data.table","ggmap","maps","mapdata",
         "ggsn","wesanderson"), require, character.only=T)

## reset wd
setwd("../../../Extra scripts for revisions/Distribution, frequency, and global extent of hypoxia in rivers/Data")

## import site data
site_info <- fread("GRDO_GEE_HA_NHD.csv")
head(site_info)

## calculate mean duration of hypoxic events per site
mean_dur <- PC %>%
  group_by(SiteID, DO_Thresh) %>%
  summarize_at(.vars = "event_dur_min", .funs = mean)
mean_dur_sub2 <- mean_dur[which(mean_dur$DO_Thresh == "[DO] < 2 (mg/L)"),]

## merge mean duration data with site_info
dur_map <- left_join(mean_dur_sub2, site_info, by="SiteID")
max(dur_map$event_dur_min)/(60*24) # max mean duration is 14 days

## Which sites have the longest mean duration
ordered_dm <- dur_map[order(dur_map$event_dur_min),]
max_mean_duration <- tail(ordered_dm, n=10)
max_mean_duration[,c("SiteID","event_dur_min","US_state")] ## 6/10 in Florida

min_mean_duration <- head(ordered_dm, n=10)
min_mean_duration[,c("SiteID","event_dur_min","US_state")] ## spread around


# map
regime_map <- ggmap(get_stamenmap(bbox=c(-125, 25, -66, 50), zoom = 5, 
                    maptype='toner'))+
  geom_point(data = dur_map, aes(x = Longitude, y = Latitude, 
                                 fill=event_dur_min, size=event_dur_min),
             shape=21)+
  theme(legend.position = "right")+
  labs(x="Longitude", y="Latitude")+
  scale_fill_gradient("Mean Event Duration",
                      low = "yellow", high = "red",
                       breaks=c(60, 60*24, 60*24*14),
                       labels=c("1 hr", "1 day", "2 weeks"),
                       limits = c(30,60*24*14), trans="log")+
  scale_size_continuous("Mean Event Duration",
                        breaks = c(60, 60*24, 60*24*14),
                        labels=c("1 hr", "1 day", "2 weeks"),
                        limits = c(30,60*24*14))
regime_map
regime_map + theme(legend.position = "none")
get_legend(regime_map)

top_row <- plot_grid(event_dur_hist+theme(axis.title = element_blank()),
                     regime_map,#+theme(legend.position = "none"),
                      nrow=1)
top_row
plot_grid(top_row, time_of_day, nrow=2)


###########################################################################
## Relationships between mean hypoxic event duration and site attributes
###########################################################################
head(dur_map)


ggplot(dur_map, aes(WaterTempC_mean, event_dur_min))+
  geom_point()+
  scale_y_continuous(trans="log")
cor.test(log(dur_map$event_dur_min), dur_map$WaterTempC_mean) ## no correlation with mean water temperature
cor.test(log(dur_map$event_dur_min), dur_map$NHD_SLOPE)














