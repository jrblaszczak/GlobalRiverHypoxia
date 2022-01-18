##################################################################################
## Blaszczak et al. - Global extent, patterns, and drivers of hypoxia in rivers
## Figure 2 - Global and US EPA region trends
## Code author: J.R. Blaszczak
#################################################################################
## Load packages
lapply(c("plyr","dplyr","ggplot2","cowplot",
         "lubridate","tidyverse","data.table",
         "stringr","zoo","wesanderson","trend",
         "rgdal","sp","sf","maptools",
         "openintro"), require, character.only=T) ## if false package needs to be installed
theme_set(theme_bw())

########################
## Import Data
######################
## set wd to GRDO 6. Temporal trends
setwd("../Distribution, frequency, and global extent of hypoxia in rivers/Data")
## Import
temporal <- fread("Temporal_Data_YearlySum.csv")
## Split
global_df <- temporal[which(temporal$Data_Subset == "Global_sum"),]
global_nosensor_df <- temporal[which(temporal$Data_Subset == "Global_nosensor_sum"),]
wqp_df <- temporal[which(temporal$Data_Subset == "WQP_sum"),]

## Sensor diff
cbind(global_df$Year, global_df$Count_Tot - global_nosensor_df$Count_Tot) ## sensor data begins in 2005

###################
## Prep data
###################
## Proportion through time
global_df$Pct_sub1 <- global_df$Count_sub1/global_df$Count_Tot
global_df$Pct_sub2 <- global_df$Count_sub2/global_df$Count_Tot
global_df$Pct_sub3 <- global_df$Count_sub3/global_df$Count_Tot
global_df$Pct_sub4 <- global_df$Count_sub4/global_df$Count_Tot
global_df$Pct_sub5 <- global_df$Count_sub5/global_df$Count_Tot
global_nosensor_df$Pct_sub2 <- global_nosensor_df$Count_sub2/global_nosensor_df$Count_Tot
global_nosensor_df$Pct_sub5 <- global_nosensor_df$Count_sub5/global_nosensor_df$Count_Tot

global_plot <- gather(global_df[,c("Year","Pct_sub1", "Pct_sub2",
                                   "Pct_sub3", "Pct_sub4",
                                   "Pct_sub5")],
                      Meas_Thresh, Pct,
                      Pct_sub1:Pct_sub5, factor_key = TRUE)

global_plot$Meas_Thresh <- revalue(global_plot$Meas_Thresh, 
                                   c("Pct_sub1" = "% [DO] < 1 (mg/L)",
                                     "Pct_sub2" = "% [DO] < 2 (mg/L)",
                                     "Pct_sub3" = "% [DO] < 3 (mg/L)",
                                     "Pct_sub4" = "% [DO] < 4 (mg/L)",
                                     "Pct_sub5" = "% [DO] < 5 (mg/L)"))
global_plot$Pct_Thresh <- global_plot$Meas_Thresh


###################
## Global Figures
###################
## Color scheme
colfunc <- colorRampPalette(c("red", "black"))
redpal <- colfunc(7)

## Count of all measurements
Global_count_plot <- ggplot(global_df, aes(Year, log(Count_Tot)))+
  geom_rect(xmin=2005, xmax=2020, ymin=0, ymax=Inf, fill="grey90", alpha=0.1)+
  geom_vline(xintercept = 2005, color="black")+
  geom_line(size=1.5, alpha=0.8)+
  scale_x_continuous(limits=c(1950,2018), breaks = c(seq(from = 1950, to = 2020, by = 10)),
                     labels=c(seq(from = 1950, to = 2020, by = 10)),
                     expand = c(0,2))+
  scale_y_continuous(breaks = c(log(100),log(1e4),
                                log(1e6),log(1e8)),
                     labels=c(100,1e4,1e6,1e8),
                     limits=c(log(100),log(1e8)))+
  labs(x="Year",y="Obs. Per Year")+
  theme(panel.grid = element_line(color="gray85"),
        axis.text = element_text(size=12, angle=45, hjust = 1))
Global_count_plot ## removed years before 1950

#### Plot
Global_pct_plot <- ggplot(global_plot, aes(Year, Pct, color=Meas_Thresh))+
  geom_rect(xmin=2005, xmax=2020, ymin=0, ymax=Inf, fill="grey90", alpha=0.1)+
  #geom_area(data=global_nosensor_df, aes(Year, Pct_sub5), color=NA,
  #          alpha=0.4, fill="black")+
  geom_area(data=global_nosensor_df, aes(Year, Pct_sub2), color=NA,
               alpha=0.5, fill=redpal[3])+
  geom_vline(xintercept = 2005, color="black")+
  geom_line(size=1.25)+
  scale_x_continuous(limits=c(1950,2018), breaks = c(seq(from = 1950, to = 2020, by = 10)),
                     labels=c(seq(from = 1950, to = 2020, by = 10)),
                     expand = c(0,2))+
  scale_y_continuous(labels = function(x) paste0(x * 100, '%'))+
  scale_color_manual(name="% Threshold", values = c("% [DO] < 1 (mg/L)"=redpal[1],
                                                    "% [DO] < 2 (mg/L)"=redpal[3],
                                                    "% [DO] < 3 (mg/L)"=redpal[4],
                                                    "% [DO] < 4 (mg/L)"=redpal[5],
                                                    "% [DO] < 5 (mg/L)"="black"))+
  theme(panel.background = element_rect(fill="white"),
      panel.grid.major.y = element_line(color="gray85"),
        axis.text = element_text(size=12, angle=45, hjust = 1))+
  labs(x="Year",y="Percent Hypoxic Obs.")+
  coord_cartesian(ylim = c(0,0.18))
Global_pct_plot

combined <- plot_grid(
  Global_count_plot+theme(legend.position = "none"),
  Global_pct_plot+theme(legend.position = "none"),
  ncol=1, align="hv",rel_heights = c(0.5,0.8))
combined


#plot_grid(get_legend(Global_pct_plot))

## what is the max percentage of observations that are hypoxic without sensor
max(global_nosensor_df[which(global_nosensor_df$Year > 1950),]$Pct_sub2*100) #5.6
tail(global_nosensor_df,n=15)

####################
## By WQP state
####################
wqp_df$file <- str_remove(wqp_df$file,"USGSNWIS_ ")
colnames(wqp_df)[which(colnames(wqp_df) == "file")] <- "State"

## Calculate proportion of sub2
wqp_df$P_sub2 <- wqp_df$Count_sub2/wqp_df$Count_Tot
## Subset to post 1950
wqp_post1950 <- wqp_df[which(wqp_df$Year >= 1950 & wqp_df$Year <2019),]

## Add upper and low CI
wqp_post1950$P_sub2_lower <- wqp_post1950$P_sub2 - 1.96*sqrt((wqp_post1950$P_sub2*(1-wqp_post1950$P_sub2))/wqp_post1950$Count_Tot)
wqp_post1950[which(wqp_post1950$P_sub2_lower < 0),]$P_sub2_lower <- 0
wqp_post1950$P_sub2_upper <- wqp_post1950$P_sub2 + 1.96*sqrt((wqp_post1950$P_sub2*(1-wqp_post1950$P_sub2))/wqp_post1950$Count_Tot)

## Create an average across all states
wqp_avg <- wqp_post1950[,c("Year","State","Count_sub2","Count_Tot")] %>% 
  group_by(Year) %>%
  summarise_at(.vars=c("Count_sub2","Count_Tot"), .funs= sum)
## Calculate proportion of sub2
wqp_avg$P_sub2 <- wqp_avg$Count_sub2/wqp_avg$Count_Tot

## Select states to highlight
#View(wqp_post1950[order(-wqp_post1950$P_sub2),])
mean_psub2 <- wqp_post1950 %>%
  group_by(State) %>%
  summarise_at(.vars=c("P_sub2"), .funs= mean)
mean_psub2 <- mean_psub2[order(-mean_psub2$P_sub2),] ## FL = 10.3%; MS 7.0%; LA 5.3%
## Create color palette for top 3
colfunc2 <- colorRampPalette(c("red", "mediumpurple4"))
#colfunc2(3); plot(rep(1,3),col=colfunc2(3),pch=19,cex=3)
redpal2 <- colfunc2(3)
## Alternate plot structure
wqp_post1950$plot_color <- ifelse(wqp_post1950$State %in% c("FL ","MS ","LA "),
                                  yes=wqp_post1950$State, no="Other")
wqp_post1950 <- within(wqp_post1950, plot_color <- fct_rev(factor(plot_color,
                                                                  levels=c("Other","LA ",
                                                                           "MS ","FL "))))
wqp_post1950 <- wqp_post1950[order(fct_rev(factor(wqp_post1950$plot_color))),]


## Figure with moving average
ggplot(data=wqp_post1950,aes(Year, P_sub2))+
  geom_line(data=wqp_post1950, aes(Year, P_sub2, group=State), alpha=0.7)+
  geom_ribbon(data=wqp_post1950, aes(ymin=P_sub2_lower,ymax=P_sub2_upper, group=State),
              alpha=0.4, show.legend = FALSE, fill="grey65", color=NA)+
  geom_line(data=subset(wqp_post1950, plot_color %in% c("FL ","MS ","LA ")),
            aes(Year, P_sub2, group=State, color=plot_color, size=plot_color))+
  geom_line(data=wqp_avg, aes(y=rollmean(P_sub2, 5, na.pad=TRUE)), size=1.5)+
  scale_y_continuous(labels = function(x) paste0(x * 100, '%'), limits = c(0,0.5))+
  scale_x_continuous(limits=c(1950,2020), breaks = c(seq(from = 1950, to = 2020, by = 10)),
                     labels=c(seq(from = 1950, to = 2020, by = 10)),
                     expand = c(0,2))+
  scale_color_manual(name="States",values=c("FL "=redpal2[1],
                                            "MS "=redpal2[2],
                                            "LA "=redpal2[3],
                                            "Other"="gray80"),
                     labels=c("Florida","Mississippi","Louisiana","Other"))+
  scale_fill_manual(name="States",values=c("FL "=redpal2[1],
                                           "MS "=redpal2[2],
                                           "LA "=redpal2[3],
                                           "Other"="gray85"),
                    labels=c("Florida","Mississippi","Louisiana","Other"))+
  scale_size_manual(name="States",values= c("FL "=1,"MS "=1,"LA "=1,"Other"= 1),
                    labels=c("Florida","Mississippi","Louisiana","Other"))+
  labs(x="Year",y="Percent DO Measurements < 2 mg/L")+
  theme(panel.grid.major.y = element_line(color="gray80"), axis.title = element_text(size=16),
        axis.text = element_text(size=14))+
  ggtitle("Water Quality Data Portal Trends \n by US State, D.C., & Puerto Rico (1950-2018)")

#########################
## States by EPA Region
########################
R1_NE <- c("CT","ME", "MA", "NH", "RI", "VT") ## New England
R2_NJNY_PR <- c("NJ","NY","PR") ## New Jersey-New York-PR
R3_MA <- c("DE", "DC", "MD", "PA", "VA", "WV") ## Mid-Atlantic
R4_SE <- c("AL", "FL", "GA", "KY", "MS", "NC", "SC", "TN") ## Southeast
R5_UMW <- c("IL", "IN", "MI", "MN", "OH", "WI") ## Upper Midwest
R6_LMW <- c("AR", "LA", "NM", "OK", "TX") ## Lower Midwest
R7_MW <- c("IA", "KS", "MO", "NE") ## Midwest
R8_RM <- c("CO", "MT", "ND", "SD", "UT", "WY") ## Rocky Mountains
R9_WC <- c("AZ", "CA", "HI", "NV") ## West
R10_PNW <- c("AK", "ID", "OR", "WA") ## Pacific NW

EPA_Regions <- list(R1_NE, R2_NJNY_PR, R3_MA,R4_SE,R5_UMW,R6_LMW,R7_MW,R8_RM,R9_WC,R10_PNW)
names(EPA_Regions) <- c("Region 1", "Region 2", "Region 3",
                        "Region 4","Region 5","Region 6",
                        "Region 7","Region 8","Region 9","Region 10")
EPA_Reg_df <- ldply(EPA_Regions, data.frame)
colnames(EPA_Reg_df) <- c("EPA_Region","State")

## Merge
wqp_post1950$State <- trimws(wqp_post1950$State, which = c("right"))
wqp_post1950 <- merge(wqp_post1950, EPA_Reg_df, by="State", all.x=TRUE)

## Summarize data by EPA region
EPA_psub2 <- wqp_post1950 %>%
  group_by(EPA_Region, Year) %>%
  summarise_at(.vars=c("Count_sub2","Count_Tot"), .funs= sum)
EPA_psub2 <- EPA_psub2[which(EPA_psub2$Year >= 1970),]
EPA_psub2$P_sub2 <- EPA_psub2$Count_sub2/EPA_psub2$Count_Tot

## Rank the regions ahead of plotting
EPA_psub2 <- within(EPA_psub2, EPA_Region <- factor(EPA_Region, levels=c("Region 1", "Region 2", "Region 3",
                                                                         "Region 4","Region 5","Region 6",
                                                                         "Region 7","Region 8","Region 9","Region 10")))

ggplot(EPA_psub2, aes(Year, P_sub2+0.01, group=EPA_Region, color=EPA_Region))+
  geom_point()+
  facet_grid(.~EPA_Region)+
  geom_smooth(method="loess")+
  scale_y_continuous(trans= "log", labels = function(x) paste0(x * 100, '%'), expand = c(0,0))+ #, limits=c(0, NA))+
  scale_x_continuous(limits=c(1970,2020), breaks = c(seq(from = 1970, to = 2020, by = 15)),
                     labels=c(seq(from = 1970, to = 2020, by = 15)),
                     expand = c(0,0))+
  theme(axis.text.x = element_text(size=14, hjust=1, angle=45),
        axis.title = element_text(size=14),
        legend.position = "none")+
  labs(x="Year",y="Percent DO Measurements < 2 mg/L")


## Examine difference trend
EPA_psub2 <- EPA_psub2 %>% mutate(diff_psub2 = P_sub2 - lag(P_sub2)) # dplyr has its own version of lag
ggplot(EPA_psub2, aes(Year, diff_psub2))+geom_point()+
  scale_y_continuous(limits=c(-0.15,0.15))+
  facet_grid(.~EPA_Region)
EPA_split <- split(EPA_psub2, EPA_psub2$EPA_Region)
lapply(EPA_split, function(x) summary(lm(x$diff_psub2 ~ x$Year))) ## no trend

## examine sen slope
senslopes_EPA <- lapply(EPA_split, function(x) sens.slope(x$P_sub2, conf.level = 0.95))
senslopes_EPA_p <- ldply(lapply(senslopes_EPA, function(x) x$p.value), data.frame); colnames(senslopes_EPA_p) <- c("Region","p.value")
senslopes_EPA_slope <- ldply(lapply(senslopes_EPA, function(x) x$estimates), data.frame); colnames(senslopes_EPA_slope) <- c("Region","Sens.slope")

sensslopes_output <- left_join(senslopes_EPA_p, senslopes_EPA_slope, by="Region")
write.csv(sensslopes_output, "Sensslope_output.csv")

########################
## EPA region map
########################
## Download necessary mapping data at these two links:
## https://rpubs.com/technocrat/thematic-alaska-hawaii
## https://rpubs.com/technocrat/thematic-walkthrough


## Reset the wd
setwd("./Map Data")

# Remove some but not all territories
remove.territories = function(.df) {
  subset(.df, 
         .df$id != "AS" &
           .df$id != "MP" &
           .df$id != "GU" & 
           #.df$id != "PR" &
           .df$id != "VI" 
  )
}

save("remove.territories", file = "helpers/remove.territories")

# Set theme
plain_theme = theme(axis.text=element_blank()) + 
  theme(panel.background = element_blank(), 
        panel.grid = element_blank(), 
        axis.ticks = element_blank())

no_ylab = ylab("") 
no_xlab = xlab("")

poly = coord_map("polyconic")

# Load the basemap
layer <- "cb_2014_us_state_5m"
cb5 = readOGR(dsn="data", layer)
save(cb5, file ="data/cb5")
us = cb5
summary(us)

# Transform geographical coordinates to Lambert Azimuth Equal Area projection
us_aea = spTransform(us, CRS("+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs"))
us_aea@data$id = rownames(us_aea@data)

# Move Alaska (scaled down), Hawaii, and Puerto Rico
alaska = us_aea[us_aea$STATEFP=="02",]
alaska = elide(alaska, rotate=-50)
alaska = elide(alaska, scale=max(apply(bbox(alaska), 1, diff)) / 2.3)
alaska = elide(alaska, shift=c(-2100000, -2500000))
proj4string(alaska) = proj4string(us_aea)

hawaii = us_aea[us_aea$STATEFP=="15",]
hawaii = elide(hawaii, rotate=-35)
hawaii = elide(hawaii, shift=c(5400000, -1400000))
proj4string(hawaii) = proj4string(us_aea)

puertorico = us_aea[us_aea$STATEFP=="72",]
puertorico = elide(puertorico, shift=c(-1200000, -100000))
proj4string(puertorico) = proj4string(us_aea)

# Remove AK and HA from base map and substitute transforms
us_aea = us_aea[!us_aea$STATEFP %in% c("02", "15","72"),]
us_aea = rbind(us_aea, alaska, hawaii, puertorico)

# initial plot
us51 <- fortify(us_aea, region="STUSPS")
us51 = remove.territories(us51)
#check
p = ggplot(data=us51) + 
  geom_map(map=us51, aes(x=long, y=lat, map_id=id, group=group),
           fill="white", color="dark grey", size=0.15) + 
  plain_theme
p


## Combine with Data
EPA_Reg_df$State_LongName <- abbr2state(EPA_Reg_df$State)
EPA_Reg_df[which(EPA_Reg_df$State == "PR"),]$State_LongName <- "Puerto Rico"
EPA_Reg_df$id <- EPA_Reg_df$State

map_us51 <- merge(us51, EPA_Reg_df,by="id")

# order groups by region
map_us51 <- within(map_us51, EPA_Region <- factor(EPA_Region, levels=c("Region 1", "Region 2", "Region 3",
                                                                         "Region 4","Region 5","Region 6",
                                                                         "Region 7","Region 8","Region 9","Region 10")))
#plot
pal <- wes_palette("Royal1", 10, type = "continuous")
pal
ggplot(data = map_us51, aes(x=long, y=lat, group = group))+
  geom_polygon(aes(x=long, y=lat, group=group, fill=EPA_Region), color = "black", size = 0.3)+
  scale_fill_manual("", values = pal)+
  plain_theme+no_xlab+no_ylab+
  theme(legend.position = "top", legend.text = element_text(size=14))


# plot individual region trends
EPA_psub2_list <- split(EPA_psub2, EPA_psub2$EPA_Region)

EPA_plot <- function(x){
  ggplot(x, aes(Year, P_sub2, group=EPA_Region))+
    geom_point(aes(fill=EPA_Region), color="black", pch=21, size=2)+
    facet_grid(.~EPA_Region)+
    geom_smooth(method="loess", aes(color=EPA_Region))+
    scale_fill_manual("EPA Regions", values = c("Region 1" = pal[1], "Region 2" = pal[2],
                                                "Region 3" = pal[3], "Region 4" = pal[4],
                                                "Region 5" = pal[5], "Region 6" = pal[6],
                                                "Region 7" = pal[7], "Region 8" = pal[8],
                                                "Region 9" = pal[9],"Region 10" = pal[10]))+
    scale_color_manual("EPA Regions", values = c("Region 1" = pal[1], "Region 2" = pal[2],
                                                 "Region 3" = pal[3], "Region 4" = pal[4],
                                                 "Region 5" = pal[5], "Region 6" = pal[6],
                                                 "Region 7" = pal[7], "Region 8" = pal[8],
                                                 "Region 9" = pal[9],"Region 10" = pal[10]))+
    scale_y_continuous(trans="log10", labels = function(x) paste0(x * 100, '%'), expand = c(0.01,0), limits=c(0.0001, 0.17))+
    scale_x_continuous(limits=c(1970,2018), breaks = c(seq(from = 1970, to = 2018, by = 15)),
                     labels=c(seq(from = 1970, to = 2018, by = 15)),
                     expand = c(0,0))+
    theme(axis.text.x = element_text(size=14, hjust=1, angle=45),
          axis.text.y = element_text(size=14),
          strip.text.x = element_text(size=14),
          axis.title = element_text(size=14), axis.title.x = element_blank(),
        legend.position = "none")+
    labs(y="Measurements < 2 mg/L")
}

# test
lapply(EPA_psub2_list[1], function(x) EPA_plot(x))

## Next save individual pdfs
plots <- lapply(EPA_psub2_list, function(x) EPA_plot(x))
## Set wd to EPA region pdfs
lapply(names(plots), function(x) ggsave(filename = paste(x,".pdf",sep=""), plot=plots[[x]],
                                        width = 3, height = 3, units = "in", device = "pdf"))






