#####################################################################################
## Blaszczak et al. - Global extent, patterns, and drivers of hypoxia in rivers
## Hypoxic event duration calculation
## Code author: J.R. Blaszczak
######################################################################################

## Load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "tidyverse","readxl","rMR","data.table"), require, character.only=T)

######################################################
## Calculate hypoxia event duration - only do once
######################################################

## (1) Set working directory to formatted files from Appling et al. 2018
# See ScienceBase data release for data formatting script
# See https://www.sciencebase.gov/catalog/item/59bff507e4b091459a5e0982 for data

Hypoxia_Duration <- function(workdir, num){
  
  setwd(workdir)
  formatted_file <- list.files(pattern=".csv")[num]

  dat <- fread(formatted_file, header=T)
  ## Change to numeric if needed
  dat$DO_mgL <- as.numeric(dat$DO_mgL)
  dat$DO_psat <- as.numeric(dat$DO_psat)
  dat$WaterTempC <- as.numeric(dat$WaterTempC)
  
  ###################################
  ## Format & QA/QC Required Data ###
  ###################################
  
  ## First adjust time and time zone
  format_time <- function(xdat){
    xdat <- as.data.frame(xdat)
    
    xdat <- xdat %>% separate(DateTime, c("Date", "Time"),sep = " ")
    xdat$Time[is.na(xdat$Time) == TRUE] <- paste("12:00:00")
    xdat$DateTime <- lubridate::ymd_hms(paste(xdat$Date, xdat$Time))
    
    ## Remove any dates before 1900 or after 2019
    xdat <- subset(xdat, DateTime > "1900-01-01 00:00:00" & DateTime < "2020-01-01 00:00:00")
    
    return(xdat)
    
  }
  
  
  dat <- format_time(dat)
  
  
  ## Next convert any DO measurements in %sat without a corresponding mg/L value to units mg/L
  converted_DOmgL <- function(xdat){
    xdat <- as.data.frame(xdat)
    xdat$DO_mgL <- ifelse(is.na(xdat$DO_mgL) == TRUE,
                          no = xdat$DO_mgL,
                          yes = DO.unit.convert(xdat$DO_psat, DO.units.in = "pct",
                                                DO.units.out = "mg/L", bar.units.in = "mmHg",
                                                bar.press = 760,temp.C = xdat$WaterTempC,
                                                salinity = 0))
    
    return(xdat)
  }
  
  dat <- converted_DOmgL(dat)
  
  
  ## Next QA/QC DO and WaterTemp Data
  qaqc <- function(xdat){ 
    xdat <- as.data.frame(xdat)
    ## QA/QC Dissolved Oxygen Data
    xdat$DO_mgL[xdat$DO_mgL < -1] <- NA ## This falls outside of what might be considered calibration error & is most likely sensor failure
    xdat$DO_mgL[xdat$DO_mgL < 0] <- 0 ## Since we've probably removed sensor failure, values between -1 and 0 are likely calibration error
    xdat$DO_mgL[xdat$DO_mgL > 35] <- NA ## This is an arbitrarily high value to remove complete outliers (can be changed)
    
    ## QA/QC Dissolved Oxygen percent saturation data
    xdat$DO_psat[xdat$DO_psat < 0] <- NA
    xdat$DO_psat[xdat$DO_psat > 250] <- NA
    
    ## QA/QC Temperature Data
    xdat$WaterTempC[xdat$WaterTempC < 0 ] <- 0 ## Less than 0 C, freezing and most likely calibration error
    xdat$WaterTempC[xdat$WaterTempC > 50 ] <- NA ## Arbitrarily high value (can be changed)
    
    return(xdat)
  }
  
  dat <- qaqc(dat)
  
  
  ## Extract events and duration of events
  # First split data to remove long time gaps
  d <- dat
  d$diff_time <- NA
  d$diff_time[1] <- 0
  
  for(i in 2:nrow(d)){
    d$diff_time[i] = difftime(time1 = d$DateTime[i], time2 = d$DateTime[(i-1)], units="min")
  }
  
  d$diff_time <- as.character(as.numeric(d$diff_time))
  d$seq <- NA
  d$seq[1] <- 1

  for(i in 2:nrow(d)){
    if(d$diff_time[i] %in% c("5","10","15","30","45","60")){
      d$seq[i] = d$seq[(i-1)]
    } else{
      d$seq[i] = d$seq[(i-1)]+1
    }
  }
  
  l <- split(d, as.factor(d$seq))
  
  # Second calculate event duration
  events_calc <- function(z, t) {
    
    zz <- z %>%
      # add id for different periods/events
      mutate(below_DO = DO_mgL < t, id = rleid(below_DO)) %>%
      # keep only periods with hypoxia
      filter(below_DO) %>%
      # calculate duration
      group_by(id) %>%
      summarize(event_dur_min = difftime(time1 = last(DateTime), first(DateTime), units = "min"),
                start_date = first(DateTime),
                end_date = last(DateTime))
    
    zz[nrow(zz)+1,] <- NA
    
    return(zz)
  }
  
  event_sub1 <- ldply(lapply(l, function(x) events_calc(x, 1)), data.frame);event_sub1$DO_thresh <- 1
  event_sub2 <- ldply(lapply(l, function(x) events_calc(x, 2)), data.frame);event_sub2$DO_thresh <- 2
  event_sub3 <- ldply(lapply(l, function(x) events_calc(x, 3)), data.frame);event_sub3$DO_thresh <- 3
  event_sub4 <- ldply(lapply(l, function(x) events_calc(x, 4)), data.frame);event_sub4$DO_thresh <- 4
  event_sub5 <- ldply(lapply(l, function(x) events_calc(x, 5)), data.frame);event_sub5$DO_thresh <- 5
  
  events <- rbind(event_sub1, event_sub2, event_sub3, event_sub4, event_sub5)
  
  ## subset
  events_df <- events[,c("event_dur_min","start_date","end_date","DO_thresh")]
  events_df$SiteID <- dat$SiteID[1]
  events_df <- na.omit(events_df)
 
  
  ## Export
  setwd("../../Hypoxic Events PC")
  
  write.csv(events_df, paste("hypoxic_events_",formatted_file, sep = ""))
  return(events_df)
  
}


WD <- "C:/Users/jblaszczak/Dropbox/FLBS/Research/Hypoxia/Organized GRDO Database for Publication/GRDO 2. Template ts with sumstat scripts/Compiled formatted TS/PC formatted TS Raw"
#WD <- "C:/Users/Joanna/Dropbox (Duke Bio_Ea)/FLBS/Research/Hypoxia/Organized GRDO Database for Publication/GRDO 2. Template ts with sumstat scripts/Compiled formatted TS/PC formatted TS Raw"


Hypoxia_Duration(WD, 2)


for(i in (1:length(list.files(pattern = ".csv")))){
  Hypoxia_Duration(WD,i)
}
  
  
