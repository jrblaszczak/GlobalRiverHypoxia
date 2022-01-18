#############################################
## Blaszczak et al. - Global extent, patterns, and drivers of hypoxia in rivers
## Figure 3 - Logistic regression of controls on river hypoxia
## Code author: J.R. Blaszczak
#############################################
## Load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","forcats","tidyr","tidyverse",
         "boot","ggridges","data.table"), require, character.only=T) #any FALSE need to be installed
theme_set(theme_bw())

##########################
## Import Hypoxia Data
##########################
## Set wd again
setwd("../Distribution, frequency, and global extent of hypoxia in rivers/Data")
## Import sites
dat <- fread("GRDO_GEE_HA_NHD.csv")
## Subset
df <- dat[,c("DB_ID","DB_Source","SiteID","DO_mgL_mean","DO_mgL_sd",
             "WaterTempC_mean","DOmgL_X10.","DOmgL_X50.",
             "DOmgL_X90.", "Hyp_pr_sub2", "Hyp_pr_sub3",
             "Hyp_pr_sub4","Hyp_pr_sub5","Lat_WGS84",
             "Lon_WGS84","Country","US_state","LU_category",
             "glc_cl_cmj","ORD_STRA","NHD_SLOPE")]
## Hypoxia as a binomial
df$Hyp_sub2_01 <- 0
df[which(df$Hyp_pr_sub2 > 0),]$Hyp_sub2_01 <- 1

#########################
## Data Preparation
#########################
# GEE
df %>%
  group_by(LU_category)%>%
  count()
## Remove NA (159), barren (206), and water (2347) before analysis
df_GEE <- df[-which(df$LU_category %in% c("barren","water")),]
df_GEE <- df_GEE[!is.na(df_GEE$LU_category),]
df_GEE$LU_category <- factor(df_GEE$LU_category)
levels(as.factor(df_GEE$LU_category))
nrow(df_GEE[is.na(df_GEE$LU_category),])

## Subset
sub_GEE <- df_GEE[,c("WaterTempC_mean","NHD_SLOPE","ORD_STRA", "Hyp_sub2_01","LU_category")]
sapply(sub_GEE, class)
sub_GEE <- na.omit(sub_GEE)

## normalize slope data
sub_GEE$NHD_SLOPE_norm <- log10(sub_GEE$NHD_SLOPE)
hist(sub_GEE$NHD_SLOPE_norm)

## standardize data
sub_GEE$WaterTempC_mean_std <- (sub_GEE$WaterTempC_mean - mean(sub_GEE$WaterTempC_mean))/sd(sub_GEE$WaterTempC_mean)
sub_GEE$NHD_SLOPE_std <- (sub_GEE$NHD_SLOPE_norm - mean(sub_GEE$NHD_SLOPE_norm))/sd(sub_GEE$NHD_SLOPE_norm)
sub_GEE$ORD_STRA_std <- (sub_GEE$ORD_STRA - mean(sub_GEE$ORD_STRA))/sd(sub_GEE$ORD_STRA)
par(mfrow=c(1,3)); hist(sub_GEE$WaterTempC_mean_std); hist(sub_GEE$NHD_SLOPE_std);hist(sub_GEE$ORD_STRA_std)
par(mfrow=c(1,1))

# pool stream orders into "Low", "Mid" and "High"
sub_GEE$ORD_STRA_pool <- fct_collapse(as.factor(sub_GEE$ORD_STRA), LOW = c("1","2","3"), MID = c("4","5","6"), HIGH = c("7","8","9"))

mean(sub_GEE$NHD_SLOPE)

######################################
## Bootstrapping two separate models
######################################

## bootstrap covariates
n<-5000
parmat<-matrix(NA,n,5)
for (i in 1:n){
  
  sampdata<- sub_GEE[sample(1:nrow(sub_GEE), nrow(sub_GEE), replace=TRUE),]
  
  riv_test <- glm(Hyp_sub2_01 ~ WaterTempC_mean_std + NHD_SLOPE_std + ORD_STRA_pool + 0, sampdata, family="binomial")
  parmat[i,]<-riv_test$coefficients
  
}

## bootstrap land use
n<-5000
parmat_GEE<-matrix(NA,n,5)
for (i in 1:n){
  
  sampdata_GEE<- df_GEE[sample(1:nrow(df_GEE), nrow(df_GEE), replace=TRUE),]
  
  LU_test_GEE <- glm(Hyp_sub2_01 ~ LU_category+0, sampdata_GEE, family="binomial")
  parmat_GEE[i,]<-LU_test_GEE$coefficients
  
}

## Save output
colnames(parmat) <- names(riv_test$coefficients)
colnames(parmat_GEE) <- names(LU_test_GEE$coefficients)
output_list <- list(parmat, parmat_GEE)
saveRDS(output_list, "../../Revised Figure Scripts/glm_Bootstraps_2021_10_17.rds")
## OR read in saved output
output_list <- readRDS("../../Revised Figure Scripts/glm_Bootstraps_2021_10_17.rds")

## split
parmat <- output_list[[1]]
parmat_GEE <- output_list[[2]]

########################################
## Visualize coefficient distributions
#######################################

## Convert to probability for histograms
logitp_GEE <- inv.logit(parmat_GEE)
pm_GEE <- as.data.frame(logitp_GEE)
pml_GEE_plot <- gather(pm_GEE)
levels(as.factor(pml_GEE_plot$key))

## Revalue
sub_GEE %>%
  group_by(LU_category)%>%
  count()
pml_GEE_plot <- gather(pm_GEE)
pml_GEE_plot$key <- revalue(pml_GEE_plot$key, c("LU_categoryforested"="Forested \n(11,025)", 
                                                "LU_categorygrassland"="Grassland \n(33,822)", 
                                                "LU_categoryurban"="Urban \n(9,322)",
                                                "LU_categorywetland"="Wetland \n(957)",
                                                "LU_categoryagricultural" = "Agriculture \n(11,344)"))

## Order by mean
pml_GEE_plot <- within(pml_GEE_plot, key <- fct_rev(factor(key, levels=c("Wetland \n(957)",
                                                                         "Urban \n(9,322)",
                                                                         "Grassland \n(33,822)",
                                                                         "Agriculture \n(11,344)",
                                                                         "Forested \n(11,025)"))))

GEE_dist_plot <- ggplot(pml_GEE_plot, aes(x = key, y = value, color= key)) +
  theme_bw()+
  #geom_jitter(size=0.5)+
  geom_boxplot()+#color="black")+
  scale_y_continuous(expand = c(0, 0), limits=c(0,0.28), breaks=c(0.05, 0.10, 0.15,0.20,0.25))+
  scale_color_manual(values=c("Forested \n(11,025)"="#285939", 
                              "Grassland \n(33,822)"="#869C80", 
                              "Urban \n(9,322)"="#666E73",
                              "Wetland \n(957)"="#93C2CC",
                              "Agriculture \n(11,344)"="#D9B814"))+
  theme(text = element_text(size=16), legend.position = "none",
        axis.title = element_blank(),
        #axis.text.x = element_text(angle=45, hjust=1),
        axis.text.y = element_text(size=14))


GEE_dist_plot

#####################
## River covariates
#####################
## Median of the bootstrapped estimates
parmatrix <- as.data.frame(parmat)

## Set range of data for data sources - standardized
temp_range_std <- (seq(from=0, to=30, by=0.5) - mean(sub_GEE$WaterTempC_mean))/sd(sub_GEE$WaterTempC_mean)
slope_std_flat <- mean(sub_GEE$NHD_SLOPE_std) - sd(sub_GEE$NHD_SLOPE_std)
slope_std_mid <- mean(sub_GEE$NHD_SLOPE_std)
slope_std_steep <- mean(sub_GEE$NHD_SLOPE_std) + sd(sub_GEE$NHD_SLOPE_std)

## Function to predict across all rows grouped by order pool
order_pool_cat <- c("ORD_STRA_poolLOW","ORD_STRA_poolMID","ORD_STRA_poolHIGH")

## Prediction function
pred_fxn <- function(x, order_pool, temp_range_std, slope_pool_std, slope_name){

  ## subset x
  sub_df <- x[,c("WaterTempC_mean_std","NHD_SLOPE_std",order_pool)]
  
  ## create matrix to receive predictions
  simmat<-matrix(NA, nrow = length(temp_range_std), ncol = length(sub_df$WaterTempC_mean_std))
  # Simulate
  for (i in 1:ncol(simmat)){
    
    #coefficients
    temp <- sub_df[i,]$WaterTempC_mean_std
    order <- sub_df[i,c(order_pool)]
    slope <- sub_df[i,]$NHD_SLOPE_std
    
    # predictions
    pred <- rep(NA, times=length(temp_range_std))
    for (j in 1:length(temp_range_std)){
      pred[j] <- order + temp*temp_range_std[j] + slope*slope_pool_std
    }
    
    #inv.logit
    z <- exp(pred)/(1+exp(pred))
    
    simmat[,i] <- z
  }
  
  # For every value extract median and CI
  median_simmat <- apply(simmat, 1, function(x) median(x))
  lower_simmat <- apply(simmat, 1, function(x) quantile(x, probs = 0.025))
  upper_simmat <- apply(simmat, 1, function(x) quantile(x, probs = 0.975))
  
  df_sim <- as.data.frame(cbind(temp_range_std, median_simmat,lower_simmat, upper_simmat))
  colnames(df_sim) <- c("temp_range_std","Median","Lower_CI","Upper_CI")
  
  # Unstandardized temp
  df_sim$temp_range <- df_sim$temp_range_std*sd(sub_GEE$WaterTempC_mean) + mean(sub_GEE$WaterTempC_mean)
  
  # Identify input combination
  df_sim$Order_cat <- order_pool
  df_sim$Slope_cat <- slope_name
  
  return(df_sim)
  
}

low_flat <- pred_fxn(parmatrix, "ORD_STRA_poolLOW", temp_range_std, slope_std_flat,"Flat")
low_mid <- pred_fxn(parmatrix, "ORD_STRA_poolLOW", temp_range_std, slope_std_mid,"Medium")
low_steep <- pred_fxn(parmatrix, "ORD_STRA_poolLOW", temp_range_std, slope_std_steep,"Steep")

mid_flat <- pred_fxn(parmatrix, "ORD_STRA_poolMID", temp_range_std, slope_std_flat,"Flat")
mid_mid <- pred_fxn(parmatrix, "ORD_STRA_poolMID", temp_range_std, slope_std_mid,"Medium")
mid_steep <- pred_fxn(parmatrix, "ORD_STRA_poolMID", temp_range_std, slope_std_steep,"Steep")

high_flat <- pred_fxn(parmatrix, "ORD_STRA_poolHIGH", temp_range_std, slope_std_flat,"Flat")
high_mid <- pred_fxn(parmatrix, "ORD_STRA_poolHIGH", temp_range_std, slope_std_mid,"Medium")
high_steep <- pred_fxn(parmatrix, "ORD_STRA_poolHIGH", temp_range_std, slope_std_steep,"Steep")

## rbind
preds <- rbind(low_flat, low_mid, low_steep,
               mid_flat, mid_mid, mid_steep,
               high_flat, high_mid, high_steep)

## Increase with water temperature
preds$groups <- paste(preds$Order_cat,preds$Slope_cat, sep="_")
preds_list <- split(preds, preds$groups)
diff_preds <- lapply(preds_list, function(x) cbind(x$temp_range_std[2:61], diff(x$Median)[1:60]))

## effect of strahler order
low_to_high_medslope <- (low_mid$Median - high_mid$Median)/low_mid$Median
mean(low_to_high_medslope)

## max difference
max(low_flat$Median) - min(low_flat$Median)
(max(low_flat$Upper_CI) - min(low_flat$Lower_CI))-(max(low_flat$Median) - min(low_flat$Median))

#####################
## Visualize!
####################

## rename the order categories
preds$RiverSize <- revalue(preds$Order_cat, c("ORD_STRA_poolLOW"="Small",
                                              "ORD_STRA_poolMID"="Medium",
                                              "ORD_STRA_poolHIGH"="Large"))
preds$RiverSize <- factor(preds$RiverSize, levels = c("Small","Medium","Large"))

# plot
Prox_plot <- ggplot(preds, aes(temp_range, Median, color=RiverSize))+
  geom_ribbon(aes(ymin = Lower_CI, ymax=Upper_CI, group=groups, fill=RiverSize), alpha=0.2, color=NA)+
  geom_line(aes(linetype=Slope_cat), size=1.75)+
  labs(x="Water Temperature (°C)", y="Probability of Riverine Hypoxia")+
  scale_color_manual("River Size", values=c("Small"="#FC4E07","Medium"="#E7B800","Large"="#2E9FDF"))+
  scale_fill_manual("River Size",guide=FALSE, values=c("Small"="#FC4E07","Medium"="#E7B800","Large"="#2E9FDF"))+
  scale_linetype_manual("Slope", values=c("Flat"="solid","Medium"="dashed","Steep"="dotted"))+
  theme_bw()+theme(axis.title = element_text(size=16),
                   axis.text.y = element_text(size=14),
                   axis.text.x = element_text(size=14))+
  scale_y_continuous(limits=c(0,0.28), expand=c(0,0), breaks=c(0.05, 0.10, 0.15,0.20,0.25))+
  scale_x_continuous(expand=c(0,0.1))
  

Prox_plot+theme(legend.position = "none")


##################################################
## COMBINE PLOTS
##################################################
## combined
plot_grid(Prox_plot+theme(legend.position = "none"),
          GEE_dist_plot+theme(axis.text.y=element_text(size=16)),
          ncol = 2, align="h", rel_widths = c(1,1))

## legend
plot_grid(get_legend(Prox_plot))



##################################################
## HydroAtlas plot - no longer included in the manuscript
##################################################

## Subset again to get the same data
sub_HA <- df_GEE[,c("WaterTempC_mean","NHD_SLOPE","ORD_STRA", "Hyp_sub2_01","glc_cl_cmj")]
sapply(sub_HA, class)
sub_HA <- na.omit(sub_HA)

## Consolidate glc HydroAtlas data
## This category corresponds to dominant land cover class
## that represents the spatial majority of the local
## catchment, HydroATLAS (unitless) (Bartholome & Belward 2005)
levels(as.factor(sub_HA$glc_cl_cmj))

## Category: Tree Cover, not regularly flooded (1,2,3,4,5,6)
## Category: Tree Cover, regularly flooded (7,8)
## Category: Shrubland, Burnt Tree Cover (9,10,11,12,18)
## Category: Herbaceous (13,14)
## Category: Wetlands (15)
## Category: Croplands (16,17)
## Category: Bare Areas (19)
## Category: Water Bodies (20)
## Category: Snow and Ice (21)
## Category: Urban (22)

## Assign grouped categories
sub_HA$GLC_Cat <- NA
sub_HA[which(sub_HA$glc_cl_cmj %in% c(1,2,3,4,5,6)),]$GLC_Cat <- "Tree Cover, Not Regularly Flooded"
sub_HA[which(sub_HA$glc_cl_cmj %in% c(9,10,11,12,18)),]$GLC_Cat <- "Shrubland, Burnt Tree Cover"
sub_HA[which(sub_HA$glc_cl_cmj %in% c(13,14)),]$GLC_Cat <- "Herbaceous"
sub_HA[which(sub_HA$glc_cl_cmj %in% c(7,8,15)),]$GLC_Cat <- "Wetlands and Regularly Flooded Tree Cover"
sub_HA[which(sub_HA$glc_cl_cmj %in% c(16,17)),]$GLC_Cat <- "Croplands"
sub_HA[which(sub_HA$glc_cl_cmj %in% c(19)),]$GLC_Cat <- "Bare Areas"
sub_HA[which(sub_HA$glc_cl_cmj %in% c(20)),]$GLC_Cat <- "Water Bodies"
sub_HA[which(sub_HA$glc_cl_cmj %in% c(21)),]$GLC_Cat <- "Snow and Ice"
sub_HA[which(sub_HA$glc_cl_cmj %in% c(22)),]$GLC_Cat <- "Urban"
sub_HA[which(sub_HA$glc_cl_cmj %in% c(23)),]$GLC_Cat <- NA
nrow(sub_HA[is.na(sub_HA$GLC_Cat),]) ## 46,011

## Summarize how many sites per grouped category
# GLC2000
sub_HA %>%
  group_by(GLC_Cat)%>%
  count()
## Remove snow and ice (22) and water bodies (298) before analysis
df_GLC <- sub_HA[-which(sub_HA$GLC_Cat %in% c("Snow and Ice","Water Bodies")),]
df_GLC <- df_GLC[!is.na(df_GLC$GLC_Cat),]
df_GLC$GLC_Cat <- factor(df_GLC$GLC_Cat)
levels(as.factor(df_GLC$GLC_Cat))
nrow(df_GLC[is.na(df_GLC$GLC_Cat),])


## bootstrap HydroAtlas
n<-5000
parmat_GLC<-matrix(NA,n,6)
for (i in 1:n){
  
  sampdata_GLC<- df_GLC[sample(1:nrow(df_GLC), nrow(df_GLC), replace=TRUE),]
  
  LU_test_GLC <- glm(Hyp_sub2_01 ~ GLC_Cat+0, sampdata_GLC, family="binomial")
  parmat_GLC[i,]<-LU_test_GLC$coefficients
  
}

#Save output
colnames(parmat_GLC) <- names(LU_test_GLC$coefficients)
saveRDS(parmat_GLC, "glm_GLC_Bootstrap_2021_01_08.rds")


## Convert to probability
logitp_GLC <- inv.logit(parmat_GLC)
pm_GLC <- as.data.frame(logitp_GLC)

## Visualize distributions & add number of sites to each category
df_GLC %>%
  group_by(GLC_Cat)%>%
  count()

pml_GLC_plot <- gather(pm_GLC)
pml_GLC_plot$key <- revalue(pml_GLC_plot$key, c("GLC_CatCroplands"="Croplands (22,022)", 
                                                "GLC_CatHerbaceous"="Herbaceous (5,441)", 
                                                "GLC_CatShrubland, Burnt Tree Cover"="Shrubland (3,844)",
                                                "GLC_CatTree Cover, Not Regularly Flooded"="Tree Cover (30,473)",
                                                "GLC_CatUrban"="Urban (4,712)",
                                                "GLC_CatWetlands and Regularly Flooded Tree Cover" = "Wetlands (193)"))
levels(as.factor(pml_GLC_plot$key))

## Order by mean
mean_GLC <- apply(logitp_GLC,2,mean)
mean_GLC
pml_GLC_plot <- within(pml_GLC_plot, key <- fct_rev(factor(key, levels=c("Wetlands (193)",
                                                                         "Urban (4,712)",
                                                                         "Herbaceous (5,441)",
                                                                         "Croplands (22,022)",
                                                                         "Tree Cover (30,473)",
                                                                         "Shrubland (3,844)"))))

GLC_dist_plot <- ggplot(pml_GLC_plot, aes(x = value, y = key)) +
  geom_density_ridges(scale = 1, rel_min_height = 0.01, 
                      quantile_lines = TRUE, fill = "#9986A5", alpha=0.5) +
  theme_bw()+
  scale_x_continuous(expand = c(0, 0), limits=c(0,0.45)) +
  theme(text = element_text(size=16),axis.title.y = element_blank(), legend.position = "none",
        axis.text.x = element_text(angle=45, hjust=1),axis.title.x = element_text(size=16),
        title = element_text(size=12))+
  labs(x = "Probability of Riverine Hypoxia")+
  ggtitle("HydroAtlas Dominant Local Catchment Land Cover (n = 66,685)")+
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))
GLC_dist_plot 





