##############################################
## Code for Processing GRDO Timeseries Data ##
##############################################
## Data should be in the format of the provided template
## Column names include:
## Required - SiteID, DateTime, UTC_TimeZone, DO_mgL, DO_psat, WaterTempC
## Optional - Discharge, Discharge Units, Does stage = channel depth?,
##            Depth, Depth Units, Stage above arbitrary datum, Stage Units, 
##            Baro Pressure, Baro Units, Salinity PSU

## List of parameters this code extracts from Required Data:
## (1) Time: Start Date, End Date, Timestep, Number of data points
## (2) For DO_mgL,DO_psat & WaterTempC: Mean, SD, Quantiles (0 (Min), 0.05, 0.1, 0.25, 0.5 (Median),
##                 0.75, 0.9, 0.95, 1 (Max)), 
## (3) Probability of Hypoxia of threshold (mg/L) 2, 3, 4, & 5 (add psat?)

## List of parameters this code extracts from Optional Data:
## Discharge: Summary stats, mean Q when hypoxic, Q antecedent to hypoxia(?)
## Depth: Summary stats, mean depth when hypoxic, depth antecendent to hypoxia
##        (Stage used when Stage = Depth)
## Baro used to convert DO_mgL to DO_psat when available
## Salinity: Summary Stats


## Load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","streamMetabolizer",
         "lubridate","tidyverse", "readxl"), require, character.only=T)

###########################
## Import Multiple Files ##
###########################

## (1) Set working directory to formatted files
getwd() # use bottom right panel to find folder & set directory

## Import and compile excel files
dat <- ldply(list.files(pattern = "xlsx"), function(filename) {
  d <- read_excel(filename, sheet=2)
  d$file <- filename
  return(d)
})
## Check
head(dat)

## Turn into list by splitting by filename
dat.list <- split(dat, dat$file)

###################################
## Format & QA/QC Required Data ###
###################################

## First adjust time and time zone
format_time <- function(xdat, time_format){
  xdat <- as.data.frame(xdat)
  ## Turn DateTime or Date into POSIXct
  xdat$DateTime <- as.POSIXct(as.character(xdat$DateTime), format=time_format)
  
  ## Specify the time zone
  xdat$GMT_tz <- paste("Etc/GMT",as.numeric(as.character(xdat$UTC_TimeZone)), sep="")
  xdat$DateTime <- force_tz(xdat$DateTime, tzone = xdat$GMT_tz[1]) ## need to work on a loop approach
  return(xdat)
  
}

dat.list <- lapply(dat.list, function(x) format_time(x, "%Y-%m-%d %H:%M:%S"))

## Next QA/QC DO and WaterTemp Data
qaqc <- function(xdat){ 
  xdat <- as.data.frame(xdat)
  ## QA/QC Dissolved Oxygen Data
  xdat$DO_mgL[xdat$DO_mgL < -1] <- NA ## This falls outside of what might be considered calibration error & is most likely sensor failure
  xdat$DO_mgL[xdat$DO_mgL < 0] <- 0 ## Since we've probably removed sensor failure, values between -1 and 0 are likely calibration error
  xdat$DO_mgL[xdat$DO_mgL > 35] <- NA ## This is an arbitrarily high value to remove complete outliers (can be changed)
  
  ## QA/QC Temperature Data
  xdat$WaterTempC[xdat$WaterTempC < 0 ] <- 0 ## Less than 0 C, freezing and most likely calibration error
  xdat$WaterTempC[xdat$WaterTempC > 50 ] <- 0 ## Arbitrarily high value (can be changed)

  return(xdat)
}

dat.list <- lapply(dat.list, function(x) qaqc(x))

## Calculate DO_psat using streamMetabolizer if not available
## Need barometric pressure in millibars, water temp in degrees C,
## salinity in PSU, & model selection = "garcia-benson"

## Use 1.33322 for mmHg to mbar

DOpsat_calc <- function(xdat, baro_unit_conversion){
  xdat <- as.data.frame(xdat)
  ## 
  xdat$Baro_mbar <- baro_unit_conversion*xdat$`Baro Pressure`
  xdat$DO_sat <- calc_DO_sat(xdat$WaterTempC, xdat$Baro_mbar, xdat$Salinity_PSU, model = "garcia-benson")
  xdat$DO_psat = 100 * (xdat$DO_mgL/ xdat$DO_sat)
  
  return(xdat)
}

dat.list <- lapply(dat.list, function(x) DOpsat_calc(x, 1.33322))
head(dat.list[2])

#####################
## Summary Metrics ##
#####################
## Recombine dat.list
df <- ldply(dat.list, data.frame)

## TIME
## Want: Start Date, End Date, Timestep, Number of data points
time_df <- df %>%
  group_by(.id) %>%
  summarise(Start_time = head(DateTime, n=1),
            End_time = tail(DateTime, n=1),
            tdiff_sec = DateTime[2] - DateTime[1],
            n_time = length(DateTime))

## DO_mgL, DO_psat, WaterTempC Summary Stats
## Want: Mean, SD, Quantiles (0 (Min), 0.05, 0.1, 0.25, 0.5 (Median),0.75, 0.9, 0.95, 1 (Max))

ss_df <- df %>%
  group_by(.id) %>%
  summarise_at(.vars= c("DO_mgL","DO_psat","WaterTempC"), .funs = c(mean, sd, min, max))

quant_DOmgL_df <- df %>%
  group_by(.id) %>% 
  do(data.frame(t(quantile(.$DO_mgL, probs = c(0.05, 0.1, 0.25, 0.50, 0.75, 0.9, 0.95)))))

quant_DOpsat_df <- df %>%
  group_by(.id) %>% 
  do(data.frame(t(quantile(.$DO_psat, probs = c(0.05, 0.1, 0.25, 0.50, 0.75, 0.9, 0.95)))))
                                                                  
quant_Temp_df <- df %>%
  group_by(.id) %>% 
  do(data.frame(t(quantile(.$WaterTempC, probs = c(0.05, 0.1, 0.25, 0.50, 0.75, 0.9, 0.95)))))


## Probability that DO < 2 mg/L
Hypox_Pr <- function(x){
  pr_sub2 <- nrow(x[which(x$DO_mgL < 2),])/nrow(x)
  pr_sub5 <- nrow(x[which(x$DO_mgL < 5),])/nrow(x)
  return(pr_sub2, pr_sub5)
}





