#####################################################################################
## Blaszczak et al. - Global extent, patterns, and drivers of hypoxia in rivers
## Figure 4A - Boosted regression tree figure following tuning of tree
## Code author: J.R. Blaszczak
######################################################################################

## Boosted regression tree
## Tutorial: https://rspatial.org/raster/sdm/9_sdm_brt.html
## More information:
## https://towardsdatascience.com/understanding-gradient-boosting-machines-9be756fe76ab
## https://bradleyboehmke.github.io/HOML/gbm.html

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "tidyverse","dismo","gbm","rsample","stringr",
         "caret","data.table"), require, character.only=T)

# Functions:
# gbm.step - fits a gbm model, however need to determine optimal lr and tc
# gbm.simply - backwards elimination of variables to drop those without predictive performance
# gmb.plot  - plots partial dependence
# gmb.interactions - test whether interactions have been detected

################################
## Import & adjust hypoxia data
################################
df <- fread("../Distribution, frequency, and global extent of hypoxia in rivers/Data/GRDO_GEE_HA_NHD.csv", header=T)
names(df)

## Hypoxia as a binomial when <2
df$Hyp_sub2_01 <- 0
df[which(df$Hyp_pr_sub2 > 0),]$Hyp_sub2_01 <- 1

## Remove rows with maximum water temperatures of greater than 35 C
nrow(df[which(df$WaterTempC_mean >= 35),])
df <- df[which(df$WaterTempC_max <= 35),]

## for non-categorical variables, identify highly correlated variables
## Use HydroAtlas for global data predictions
vars <- df %>% dplyr::select(c(slope_calc, #streambed slope calculated
                               UPLAND_SKM, #total upstream area HydroAtlas
                               ORD_STRA, #Strahler stream order HydroAtlas
                               dis_m3_pyr, #Annual average discharge HydroAtlas
                               ele_mt_cav, #Average elevation for upstream watershed HydroAtlas
                               slp_dg_uav, #Average terrain slope for upstream watershed HydroAtlas
                               tmp_dc_uyr, #Annual average air temperature for upstream watershed HydroAtlas
                               pre_mm_uyr, #Annual average precip for upstream watershed HydroAtlas
                               ari_ix_uav, #Average global aridity index for upstream watershed HydroAtlas,
                               wet_pc_ug2, #Wetland pct cover for upstream watershed HydroAtlas
                               for_pc_use, #Forest pct cover for upstream watershed HydroAtlas
                               crp_pc_use, #Crop pct cover for upstream watershed HydroAtlas
                               urb_pc_use, #Urban pct cover for upstream watershed HydroAtlas
                               ppd_pk_uav, #Population density for upstream watershed HydroAtlas
                               rdd_mk_uav, #Road density for upstream watershed HydroAtlas
                               hdi_ix_cav, #Human development index for upstream watershed HydroAtlas
                               WaterTempC_mean,
                               WaterTempC_max)) #Maximum water temperature

check_cor <- 
  cor(vars, use = "pairwise.complete.obs", method = "spearman") #spearman is non-parametric

cor_high <-
  check_cor %>% 
  reshape2::melt() %>% 
  subset(Var1 != Var2 & is.finite(value)) %>% 
  subset(abs(value) >= 0.75) %>% # > 0.75
  dplyr::distinct() %>% 
  dplyr::arrange(Var1, Var2)

## Make choices between highly correlated variables
## Retain UPLAND_SKM over dis_m3_pyr (0.90), ORD_STRA (0.95)
## Retain urb_pc_use over ppd_pk_uav (0.85), rdd_mk_uav (0.82)

dat <- dplyr::select(df, names(vars))
dat <- dat %>% dplyr::select(-c("dis_m3_pyr","ORD_STRA",
                                "ppd_pk_uav","rdd_mk_uav"))
dat <- cbind(df[,c("DB_ID","DB_Source","Hyp_sub2_01")], dat)


## Prep data for model tuning:
#first, randomize the order of the data set
set.seed(42)
rows <- sample(nrow(dat))
dat <- dat[rows, ]
## Split data into training (now) and testing (later) sets
data_split <- initial_split(dat, prop = .75)
data_train <- training(data_split)
data_test <- testing(data_split)
head(data_train); head(data_test)

##############################################
## Figure: plot relative variable importance
###############################################
## Import data from tuned boosted regression tree model
hyp.mod <- readRDS("hyp.lr1.br7_maxwaterT_2021_10_14.rds")

var.info <- read.csv("BRT_pred_var_info.csv", header=T)
rel.imp <- summary(hyp.mod)
var.rel.imp <- merge(rel.imp, var.info, by="var")
names(var.rel.imp)

ggplot(data = var.rel.imp, aes(fct_reorder(short.name, rel.inf), rel.inf))+
  geom_col(aes(fill=var.category), color="black")+
  coord_flip()+theme_bw()+
  scale_fill_manual("Attributes", values = c("Topography"="#444F52",
                                                       "Land Cover"="#748C70",
                                                       "Climate"="#F2CE1B",
                                                       "Stream Attr."="#32738C"))+
  scale_y_continuous(limits=c(0,45),labels = function(x) paste0(x*1, "%"))+
  theme(axis.title.y = element_blank(), axis.title.x = element_text(size=14),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.position = "right",legend.text = element_text(size=12),
        title = element_text(size=14))+
  labs(y="Variable Relative Influence")
  #guides(fill = guide_legend(nrow = 2))




######################################################
## Use final model to predict results for test set
#####################################################
#gbm.plot(hyp.mod, n.plots = 13, write.title = FALSE)
#gbm.plot.fits(hyp.mod, v = 4)
print(hyp.mod)

find.int <- gbm.interactions(hyp.mod)
find.int$interactions
find.int$rank.list

par(mfrow=c(1,1))
gbm.perspec(hyp.mod, 7, 13, x.label = "Aridity Index",
            y.label = "Max Water Temp", z.label = "Fitted Values")

## prediction
preds <- predict.gbm(hyp.mod, data_test,
                     n.trees = hyp.mod$gbm.call$best.trees,
                     type="response")

preds.binary <- as.factor(ifelse(preds>0.7, 1, 0))
data_test$Hyp_sub2_01 <- as.factor(data_test$Hyp_sub2_01)
confusionMatrix(preds.binary, data_test$Hyp_sub2_01)




























