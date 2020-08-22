###############################################################################
####################### VALE INSTITUTE OF TECHNOLOGY ##########################
###############################################################################

#####################  LANDSCAPE GENOMICS TUTORIAL   ##########################
################   MLPE - STEP 03: RUNNING MLPE FOR IBR   #####################

### Script prepared by Carolina S. Carvalho, Jeronymo Dalapicolla, Luciana C. Resende-Moreira, Jamille C. Veiga, and Rodolfo Jaffé ###

#------------------------------------------------------------------------------
#                               PRE-ANALYSIS 
#------------------------------------------------------------------------------
##AIM = SELECTING BEST MODELS AND SHOW IMPORTANCE OF VARIABLES

##Load packages
library(corMLPE)
library(nlme)
library(MuMIn)
library(fmsb)
library(GeNetIt)
library(tidyverse)
library(ggradar)

#load this function:
max.r <- function(x){
  if(class(x)[length(class(x))] == "lm"){
    corm <- summary(x, correlation=TRUE)$correlation}
  else if(class(x) =="lmerMod"){
    corm <- cov2cor(vcov(x))}
  else if(class(x) =="lmerModLmerTest"){
    corm <- cov2cor(vcov(x))}
  else if(class(x) =="glmerMod"){
    corm <- cov2cor(vcov(x))}
  else if(class(x)=="gls"){
    corm <- summary(x)$corBeta} 
  else if(class(x)=="lme"){
    corm <- summary(x)$corFixed}
  else { print("Error: Invalid model class")}
  corm <- as.matrix(corm)
  if (length(corm)==1){
    corm <- 0
    max(abs(corm))
  } else if (length(corm)==4){
    cormf <- corm[2:nrow(corm),2:ncol(corm)]
    cormf <- 0
    max(abs(cormf))
  } else {
    cormf <- corm[2:nrow(corm),2:ncol(corm)]
    diag(cormf) <- 0
    max(abs(cormf))
  }
}


#------------------------------------------------------------------------------
#                            1. Loading Files
#------------------------------------------------------------------------------
#Load MLPE table:
mlpe_simonsi = read.csv("simonsi/MLPE_table_simonsi.csv", row.names = 1)
head(mlpe_simonsi)

#maximum of variable in MLPE. 15 observations to each predictor:
max_var = length(mlpe_steerei[,1])/15
max_var #11


#creating formula for MLPE:
#all variables:
vars_mlpe_st = names(mlpe_steerei)[5:length(mlpe_steerei)]
vars_mlpe_st


#Using PCA distance:
vars_mlpe_st_pca = vars_mlpe_st[1:7] #Only Wetlands_L. See Univariate results in step #2:
mlpe_formula_pca_st = as.formula(paste("PCAgen ~ ",paste(vars_mlpe_st_pca, collapse = " + "),sep = ""))
mlpe_formula_pca_st


#Using REL distance:
vars_mlpe_st_rel = vars_mlpe_st[c(1:6,10)] #Only Wetlands_VH. See Univariate results in step #2:
mlpe_formula_rel_st = as.formula(paste("RELgen ~ ",paste(vars_mlpe_st_rel, collapse = " + "),sep = ""))
mlpe_formula_rel_st



#------------------------------------------------------------------------------
#                            2. Running Full Model
#------------------------------------------------------------------------------
##Full model for REL:
Full_model_st_rel = gls(mlpe_formula_rel_st,
                 correlation = corMLPE(form = ~ from_ID + to_ID),
                 data = mlpe_steerei,
                 method = "ML")

##Run dredge function:
#specifying the number of predictor variables (max_var) and including the max.r function
Allmodels = dredge(Full_model_st_rel, rank = "AIC", m.lim=c(0, 11), extra= c(max.r))
length(Allmodels[,1]) # 128models

##Retrieve non-collinear models (max.r <=0.7)
NCM = get.models(Allmodels, subset = max.r<=0.7)
##Final model selection table
BM = model.sel(NCM)
length(BM[,1]) #63 models
write.csv(BM, row.names=TRUE, file="steerei/MLPE_REL/Allmodels_woCor_REL_sorted_steerei.csv")

## Selecting Best models - AIC <=2.
length(BM[BM$delta<=2, ][,1]) #9 models
write.csv(BM[BM$delta<=2, ], "steerei/MLPE_REL/BestModels_REL_steerei.csv")

## Calculate variable importance in all uncorrelated and best models:
impor_uncor = importance(BM)
write.csv(impor_uncor, "steerei/MLPE_REL/ImportanceVariables_Uncor_REL_steerei.csv")

impor_best = importance(BM[BM$delta<=2, ])
write.csv(impor_best, "steerei/MLPE_REL/ImportanceVariables_Best_REL_steerei.csv")

## Model-averaged coeficients
MA = model.avg(BM, subset = delta<=2)
confint(MA) ## Confidence Intervals for model-averaged coeficients
as.data.frame(MA$coefficients) ## Mean 

write.csv(as.data.frame(confint(MA)), row.names=TRUE, file="steerei/MLPE_REL/BestModels_CI_coefficients_REL_steerei.csv")
write.csv(as.data.frame(MA$coefficients), row.names=TRUE, file="steerei/MLPE_REL/BestModels_Mean_coefficients_REL_steerei.csv")


#------------------------------------------------------------------------------
#                            3. Running Best Model
#------------------------------------------------------------------------------
##Best model for REL:
BM[BM$delta<=2, ][1,]
names(mlpe_steerei)

##Run Best model 
Adjust_bestmodel = gls(RELgen ~ eucl_dist +PET +temperature,
                       correlation = corMLPE(form = ~ from_ID + to_ID),
                       data = mlpe_steerei,
                       method = "REML") #for better estimation of coefficients


#Information for Best Model:
model_st = as.data.frame(capture.output(Adjust_bestmodel))
summ_st = as.data.frame(capture.output(summary(Adjust_bestmodel)))
write.csv(model_st, file="steerei/MLPE_REL/AjustedBestModels1_steerei.csv", row.names = TRUE)
write.csv(summ_st, file="steerei/MLPE_REL/AjustedBestModels2_steerei.csv", row.names = TRUE)

## Model-averaged coeficients using REML - 9 models:
best_models = BM[BM$delta<=2, ]
names(mlpe_steerei)

best_models[2]
bestmodel_2nd = gls(RELgen ~ wetlands_VH +PET +temperature,
                    correlation = corMLPE(form = ~ from_ID + to_ID), data = mlpe_steerei, method = "REML")

best_models[3]
bestmodel_3rd = gls(RELgen ~ wetlands_VH +PET +temperature +precipitation +riverdistance,
                    correlation = corMLPE(form = ~ from_ID + to_ID), data = mlpe_steerei, method = "REML")

best_models[4]
bestmodel_4th = gls(RELgen ~ wetlands_VH +PET +temperature +riverdistance,
                    correlation = corMLPE(form = ~ from_ID + to_ID), data = mlpe_steerei, method = "REML")

best_models[5]
bestmodel_5th = gls(RELgen ~ PET +temperature +riverdistance +eucl_dist,
                    correlation = corMLPE(form = ~ from_ID + to_ID), data = mlpe_steerei, method = "REML")

best_models[6]
bestmodel_6th = gls(RELgen ~ PET +temperature +topography,
                    correlation = corMLPE(form = ~ from_ID + to_ID), data = mlpe_steerei, method = "REML")

best_models[7]
bestmodel_7th = gls(RELgen ~ PET +temperature +riverdistance +topography,
                    correlation = corMLPE(form = ~ from_ID + to_ID), data = mlpe_steerei, method = "REML")

best_models[8]
bestmodel_8th = gls(RELgen ~ PET +temperature +precipitation +wetlands_VH,
                    correlation = corMLPE(form = ~ from_ID + to_ID), data = mlpe_steerei, method = "REML")

best_models[9]
bestmodel_9th = gls(RELgen ~ PET +temperature +precipitation +eucl_dist,
                    correlation = corMLPE(form = ~ from_ID + to_ID), data = mlpe_steerei, method = "REML")

#Save Results:
MA_REML = model.avg(Adjust_bestmodel, bestmodel_2nd, bestmodel_3rd, bestmodel_4th, bestmodel_5th, bestmodel_6th, bestmodel_7th, bestmodel_8th, bestmodel_9th)
confint(MA_REML) ## Confidence Intervals for model-averaged coeficients
as.data.frame(MA_REML$coefficients) ## Mean 

write.csv(as.data.frame(confint(MA_REML)), row.names=TRUE, file="steerei/MLPE_REL/BestModels_CI_REML_coefficients_REL_steerei.csv")
write.csv(as.data.frame(MA_REML$coefficients), row.names=TRUE, file="steerei/MLPE_REL/BestModels_Mean_REML_coefficients_REL_steerei.csv")


#------------------------------------------------------------------------------
#              4. Spatial Autocorrelation in Residuals from Best Model
#------------------------------------------------------------------------------
#From Jaffé et al. 2019

## Load distances dataframe
df = mlpe_steerei

#rename cols:
colnames(df)[1:2] = c("ID1", "ID2")
colnames(df)[4] =c("Rij")
colnames(df)[5] =c("RD.Geo.Ind")
colnames(df)

## Sort data
df = df[order(df$ID1, df$ID2),]

## Identify which samples are from same location based on pairwise geographic distance
ulab = unique(c(as.character(df$ID1), as.character(df$ID2)))
dis = matrix(0, length(ulab), length(ulab))
rownames(dis) = colnames(dis) = ulab
for(i in 1:nrow(df)){
  dis[df$ID1[i],df$ID2[i]] = dis[df$ID2[i],df$ID1[i]] = df$RD.Geo.Ind[i]
  location = cutree(hclust(as.dist(dis)),h=0) #assign individuals to unique locations
}

location

### Model IBD using MLPE and NMLPE
m1 <- gls(Rij ~ RD.Geo.Ind, correlation = corMLPE(form = ~ID1+ID2), data = df) ## MLPE model
m2 <- gls(Rij ~ RD.Geo.Ind, correlation = corNMLPE2(form = ~ID1+ID2, clusters = location), data = df) ##NMLPE model

acf(resid(m1,type='normalized')) ## Spatial dependence pattern
acf(resid(m2,type='normalized')) ## The strong pattern disappears

summary(m1)
summary(m2)

#Testing significance. Models are the same.
anova(m1,m2)

#Save graphs
pdf("ACF_regular_steerei.pdf", onefile = T)
acf(resid(m1,type='normalized'))
dev.off()

pdf("ACF_correted_steerei.pdf", onefile = T)
acf(resid(m2,type='normalized'))
dev.off()


#------------------------------------------------------------------------------
#                                 5. Radar Graphs
#------------------------------------------------------------------------------
#df for ploting results for importance variables:
df = rbind(impor_uncor, impor_best)
df = as.data.frame(df)
group = c("Uncorrelated Models", "Best Models")
df = cbind(group,df)


plot_st = 
ggradar(df,
        base.size = 2,
        values.radar = c("0%", "50%", "100%"),
        axis.labels = c( "PET", "Temperature",  "River", "Wetlands", "Precipitation", "Euclidian", "Topography"),
        grid.min = 0,
        grid.mid = 0.5,
        grid.max = 1,
        grid.label.size = 4,
        axis.label.size = 5,
        axis.line.colour = "black",
        group.line.width = 1.5,
        group.point.size = 5,
        group.colours = c("red", "blue"),
        background.circle.colour = "grey",
        background.circle.transparency = 0.1,
        gridline.min.linetype = "dashed",
        gridline.mid.linetype = "dashed",
        gridline.max.linetype = "dashed",
        gridline.min.colour = "black",
        gridline.mid.colour = "black",
        gridline.max.colour = "black",
        legend.title = "",
        legend.text.size = 12,
        legend.position = "bottom"
        )

pdf("steerei/Radar_Variables_steerei.pdf", onefile = T)
plot_st
dev.off()




##----------------------------------------------------------
##SIMONSI
mlpe_simonsi = read.csv("simonsi/MLPE_table_simonsi.csv", row.names = 1)
head(mlpe_simonsi)

#maximum of variable in MLPE. 15 observations to each predictor:
max_var = length(mlpe_simonsi[,1])/15
max_var #9


#creating formula for MLPE:
#all variables:
vars_mlpe_si = names(mlpe_simonsi)[5:length(mlpe_simonsi)]
vars_mlpe_si


#Using PCA distance:
vars_mlpe_si_pca = vars_mlpe_si[1:7] #Only Wetlands_L. See Univariate results in step #2:
mlpe_formula_pca_si = as.formula(paste("PCAgen ~ ",paste(vars_mlpe_si_pca, collapse = " + "),sep = ""))
mlpe_formula_pca_si


#Using REL distance:
vars_mlpe_si_rel = vars_mlpe_si[c(1:6,10)] #Only Wetlands_VH. See Univariate results in step #2:
mlpe_formula_rel_si = as.formula(paste("RELgen ~ ",paste(vars_mlpe_si_rel, collapse = " + "),sep = ""))
mlpe_formula_rel_si




#------------------------------------------------------------------------------
#                            2. Running Full Model
#------------------------------------------------------------------------------
##Full model for REL:
Full_model_si_rel = gls(mlpe_formula_rel_si,
                        correlation = corMLPE(form = ~ from_ID + to_ID),
                        data = mlpe_simonsi,
                        method = "ML")

##Run dredge function:
#specifying the number of predictor variables (max_var) and including the max.r function
Allmodels = dredge(Full_model_si_rel, rank = "AIC", m.lim=c(0, 9), extra= c(max.r))
length(Allmodels[,1]) # 128models

##Retrieve non-collinear models (max.r <=0.7)
NCM = get.models(Allmodels, subset = max.r<=0.7)
##Final model selection table
BM = model.sel(NCM)
length(BM[,1]) #64 models
write.csv(BM, row.names=TRUE, file="simonsi/MLPE_REL/Allmodels_woCor_REL_sorted_simonsi.csv")

## Selecting Best models - AIC <=2.
length(BM[BM$delta<=2, ][,1]) #9 models
write.csv(BM[BM$delta<=2, ], "simonsi/MLPE_REL/BestModels_REL_simonsi.csv")

## Calculate variable importance in all uncorrelated and best models:
impor_uncor = importance(BM)
write.csv(impor_uncor, "simonsi/MLPE_REL/ImportanceVariables_Uncor_REL_simonsi.csv")

impor_best = importance(BM[BM$delta<=2, ])
write.csv(impor_best, "simonsi/MLPE_REL/ImportanceVariables_Best_REL_simonsi.csv")

## Model-averaged coeficients
MA = model.avg(BM, subset = delta<=2)
confint(MA) ## Confidence Intervals for model-averaged coeficients
as.data.frame(MA$coefficients) ## Mean 

write.csv(as.data.frame(confint(MA)), row.names=TRUE, file="simonsi/MLPE_REL/BestModels_CI_coefficients_REL_simonsi.csv")
write.csv(as.data.frame(MA$coefficients), row.names=TRUE, file="simonsi/MLPE_REL/BestModels_Mean_coefficients_REL_simonsi.csv")


#------------------------------------------------------------------------------
#                            3. Running Best Model
#------------------------------------------------------------------------------
##Best model for REL:
BM[BM$delta<=2, ][1,]
names(mlpe_simonsi)

##Run Best model 
Adjust_bestmodel = gls(RELgen ~ PET +wetlands_VH,
                       correlation = corMLPE(form = ~ from_ID + to_ID),
                       data = mlpe_simonsi,
                       method = "REML") #for better estimation of coefficients


#Information for Best Model:
model_st = as.data.frame(capture.output(Adjust_bestmodel))
summ_st = as.data.frame(capture.output(summary(Adjust_bestmodel)))
write.csv(model_st, file="simonsi/MLPE_REL/AjustedBestModels1_simonsi.csv", row.names = TRUE)
write.csv(summ_st, file="simonsi/MLPE_REL/AjustedBestModels2_simonsi.csv", row.names = TRUE)

## Model-averaged coeficients using REML - 2 models:
best_models = BM[BM$delta<=2, ]
names(mlpe_simonsi)

best_models[2]
bestmodel_2nd = gls(RELgen ~ wetlands_VH +PET +riverdistance,
                    correlation = corMLPE(form = ~ from_ID + to_ID), data = mlpe_simonsi, method = "REML")

#Save Results:
MA_REML = model.avg(Adjust_bestmodel, bestmodel_2nd)
confint(MA_REML) ## Confidence Intervals for model-averaged coeficients

as.data.frame(MA_REML$coefficients) ## Mean 

write.csv(as.data.frame(confint(MA_REML)), row.names=TRUE, file="simonsi/MLPE_REL/BestModels_CI_REML_coefficients_REL_simonsi.csv")
write.csv(as.data.frame(MA_REML$coefficients), row.names=TRUE, file="simonsi/MLPE_REL/BestModels_Mean_REML_coefficients_REL_simonsi.csv")


#------------------------------------------------------------------------------
#              4. Spatial Autocorrelation in Residuals from Best Model
#------------------------------------------------------------------------------
#From Jaffé et al. 2019

## Load distances dataframe
df = mlpe_simonsi

#rename cols:
colnames(df)[1:2] = c("ID1", "ID2")
colnames(df)[4] =c("Rij")
colnames(df)[5] =c("RD.Geo.Ind")
colnames(df)

## Sort data
df = df[order(df$ID1, df$ID2),]

## Identify which samples are from same location based on pairwise geographic distance
ulab = unique(c(as.character(df$ID1), as.character(df$ID2)))
dis = matrix(0, length(ulab), length(ulab))
rownames(dis) = colnames(dis) = ulab
for(i in 1:nrow(df)){
  dis[df$ID1[i],df$ID2[i]] = dis[df$ID2[i],df$ID1[i]] = df$RD.Geo.Ind[i]
  location = cutree(hclust(as.dist(dis)),h=0) #assign individuals to unique locations
}

location

### Model IBD using MLPE and NMLPE
m1 <- gls(Rij ~ RD.Geo.Ind, correlation = corMLPE(form = ~ID1+ID2), data = df) ## MLPE model
m2 <- gls(Rij ~ RD.Geo.Ind, correlation = corNMLPE2(form = ~ID1+ID2, clusters = location), data = df) ##NMLPE model

acf(resid(m1,type='normalized')) ## Spatial dependence pattern
acf(resid(m2,type='normalized')) ## The strong pattern disappears

summary(m1)
summary(m2)

#Testing significance. p>0.5 Models are the same.
anova(m1,m2)

#Save graphs
pdf("ACF_regular_simonsi.pdf", onefile = T)
acf(resid(m1,type='normalized'))
dev.off()

pdf("ACF_correted_simonsi.pdf", onefile = T)
acf(resid(m2,type='normalized'))
dev.off()


#------------------------------------------------------------------------------
#                             5. Radar Graphs
#------------------------------------------------------------------------------
#df for ploting results for importance variables:
impor_best2 = impor_best

x = as.data.frame(impor_best2)

x[nrow(x) + 1,] = 0
x[nrow(x) + 1,] = 0
x[nrow(x) + 1,] = 0
x[nrow(x) + 1,] = 0
row.names(x)[4:7] = c("precipitation", "temperature", "topography", "eucl_dist")


df = rbind(impor_uncor, t(x))
df = as.data.frame(df)
group = c("Uncorrelated Models", "Best Models")
df = cbind(group,df)
df

#order as another species:
names(df)
df = df[c(1,2,6,4,3,5,8,7)]

plot_si = 
  ggradar(df,
          base.size = 2,
          values.radar = c("0%", "50%", "100%"),
          axis.labels = c( "PET", "Temperature",  "River", "Wetlands", "Precipitation", "Euclidian", "Topography"),
          grid.min = 0,
          grid.mid = 0.5,
          grid.max = 1,
          grid.label.size = 4,
          axis.label.size = 5,
          axis.line.colour = "black",
          group.line.width = 1.5,
          group.point.size = 5,
          group.colours = c("red", "blue"),
          background.circle.colour = "grey",
          background.circle.transparency = 0.1,
          gridline.min.linetype = "dashed",
          gridline.mid.linetype = "dashed",
          gridline.max.linetype = "dashed",
          gridline.min.colour = "black",
          gridline.mid.colour = "black",
          gridline.max.colour = "black",
          legend.title = "",
          legend.text.size = 12,
          legend.position = "bottom"
  )

pdf("simonsi/Radar_Variables_simonsi.pdf", onefile = T)
plot_si
dev.off()


#END  