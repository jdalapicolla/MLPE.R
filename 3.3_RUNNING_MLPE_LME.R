#####################  LANDSCAPE GENOMICS TUTORIAL   ##################
########################   STEP 03: MLPE MODELS  ######################

### Script prepared by Jeronymo Dalapicolla, Jamille C. Veiga, Carolina S. Carvalho, Luciana C. Resende-Moreira, and Rodolfo Jaffé ###


##### PRE-ANALYSIS -------

#AIM = RUNNING MIXED MODELS TO IBR AND CALCULATE IMPORTANCE OF VARIABLES

### Load packages
library(corMLPE)
library(nlme)
library(MuMIn)
library(fmsb)
library(GeNetIt)
library(tidyverse)
library(r2glmm)
library(snow)
library(parallel)
library(GGally)
library(reshape2)
library(gridExtra)
library(RColorBrewer)
library(ggradar) #devtools::install_github("ricardo-bion/ggradar", dependencies = TRUE)
library(scales)

### Load auxiliary functions:
## Check max correlation
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

# Calculate decorrelated model residuals for GLS objects:
decorltd_res_gls = function(object){
  val <- object[["residuals"]]
  fit <- fitted(object) 
  rawfit <- fit
  raw <- resid(object) + fit
  fit <- fit/attr(val, "std")
  val <- val/attr(val, "std")
  cSt <- object$modelStruct$corStruct
  val <- corMLPE:::recalc.corMLPE(cSt, list(Xy = as.matrix(val)))$Xy
  return(data.frame(pr=fit+val,fit=fit,raw=raw,rawfit=rawfit))
}

# Calculate decorrelated model residuals for LME objects
decorltd_res_lme = function(object){
  object$fitted
  object$residuals
  val <- object[["residuals"]]
  fit <- as.data.frame(object$fitted)$fixed
  rawfit <- fit
  raw <- resid(object) + fit
  fit <- fit/attr(val, "std")
  val <- val/attr(val, "std")
  cSt <- object$modelStruct$corStruct
  val <- corMLPE:::recalc.corMLPE(cSt, list(Xy = as.matrix(val)))$Xy
  return(data.frame(pr=fit+val,fit=fit,raw=raw,rawfit=rawfit))
}


##### ANALYSIS


##### 1. LOAD FILES --------------
#A. Load MLPE table:
mlpe_steerei = read.csv("Metafiles/MLPE_table_steerei.csv", row.names = 1)
head(mlpe_steerei)

#B. Scale predictors/variables for use in models
mlpe_steerei[, 5:(length(mlpe_steerei)-1)] = as.data.frame(scale(mlpe_steerei[, 5:(length(mlpe_steerei)-1)]))
head(mlpe_steerei)

#C. Creating formula for MLPE:
#select variables for models:
vars_mlpe_st = names(mlpe_steerei)[5:(length(mlpe_steerei)-1)]
vars_mlpe_st

#use REL distance plus varaibles of your choice keeping one for geographic distance (Null Model):
vars_mlpe_st_rel = vars_mlpe_st[c(1:5)]
mlpe_formula_rel_st = as.formula(paste("RELgen ~ ",paste(vars_mlpe_st_rel, collapse = " + "),sep = ""))
mlpe_formula_rel_st #5 variables




##### 2. CREATE A CLUSTER TO RUN MLPE IN PARALLEL ------
#A. Make cluster
cluster = makeCluster(4, type = "SOCK")  ## also need snow installed

#B. Use clusterExport to send data to global environment (aka ‘workspace’) of each node. 
clusterExport(cluster,"mlpe_steerei")

#C. Load packages in the cluster
clusterEvalQ(cluster,
             c(library(nlme), library(MuMIn), library(corMLPE)))



##### 3. RUNNING FULL MODEL -----------
#A. Full model:
Fullmodel = nlme::lme(mlpe_formula_rel_st,
                        random = ~1|POP,
                        correlation = corMLPE(form = ~ from_ID + to_ID),
                        data = mlpe_steerei,
                        method = "ML") #For model selection is ML

#B. Check model:
summary(Fullmodel)
RES <- residuals(Fullmodel, type="normalized")
FIT <- fitted(Fullmodel)
plot(FIT, RES) ; abline(0,0, col="red") #OK
acf(RES) ##there small or none spatial autocorrelation in residuals



##### 4. RUN MODEL SELECTION WITH PARALLELIZED DREDGE ---- 
#A. Check maximum correlation
max.r(Fullmodel) ##  0.9529138

#B. Specify the number of predictor variables and including the max.r function. 15 observations to each predictor:
max_var = length(mlpe_steerei[,1])/15
max_var #11
nrow(mlpe_steerei) ## 171
options(na.action = na.fail)

#C. Run pdredge:
start_time = Sys.time()
Allmodels = MuMIn::pdredge(Fullmodel, rank = "AIC", 
                            m.lim=c(0, 11), extra = c(max.r), cluster)
end_time = Sys.time()

#D. Compute run time
start_time1 = start_time
end_time1 = end_time
diff1 = end_time - start_time
diff1 ##2.72156 secs

#E. Number of models
nrow(Allmodels) ## 32

#F. Save and load All Models:
save(Allmodels, file="./Results/steerei/FullModel/Allmodels_steerei.RData")
load(file = "./Results/steerei/FullModel/Allmodels_steerei.RData")

#G. Retrieve Not Colinearity Models, with max.r <=0.6. Get Not colinearity Models using cluster too!
NCM = get.models(Allmodels, subset = max.r <= 0.6, cluster)

#H. Number of Not colinear models
length(NCM) #7

#I. Save and load Not Colinear models
save(NCM, file="./Results/steerei/FullModel/NCM_steerei.RData")
load(file="./Results/steerei/FullModel/NCM_steerei.RData")

#J. Select Best Models using AIC
BM = model.sel(NCM, rank=AIC)
nrow(BM) #7
#Check best models
best_models = BM[BM$delta <= 2, ]
best_models

#K. Save and load Model selection
save(BM, file="./Results/steerei/FullModel/BM_steerei.RData")
load(file="./Results/steerei/FullModel/BM_steerei.RData")


#L. Save model selection as table/dataframe
df_model_sel = as.data.frame(BM)
head(df_model_sel)
write.csv(df_model_sel, "./Results/steerei/FullModel/ModelsSelection_AIC_steerei.csv", row.names = F)


#M. Save best models tables
as.data.frame(best_models)
write.csv(as.data.frame(best_models), "./Results/steerei/FullModel/BestModels_AIC_steerei.csv", row.names = F)


#####5. REFIT BEST MODELS USING REML AND MODEL VALIDATION ----
#A. TopM1 ----
TopM1 = get.models(BM, subset = 1)[[1]]
summary(TopM1)
intervals(TopM1)
TopM1

#B. Refit TopM1
Refit_TopM1 <- nlme::lme(RELgen ~ riverdistance + wetlands_L,
                         random = ~1|POP,
                         correlation = corMLPE(form = ~ from_ID + to_ID),
                         data = mlpe_steerei,
                         method = "REML") #get accurated estimates
summary(Refit_TopM1)

#C. Check autocorrelation of residuals
RES <- residuals(Refit_TopM1, type="normalized")
FIT <- fitted(Refit_TopM1)
plot(FIT, RES) ; abline(0,0, col="red") #ok
acf(RES) #ok

#D. Plot predictors x residuals
plot(mlpe_steerei$riverdistance, RES) ; abline(0,0, col="red") #ok
plot(mlpe_steerei$wetlands_L, RES) ; abline(0,0, col="red") #ok

#E. Save ACF graph:
acf(resid(Refit_TopM1,type='normalized')) ## Spatial dependence pattern
#Save graphs
pdf("Results/steerei/Figures/refit_ACF_LME_steerei_TopM1.pdf", onefile = T)
acf(resid(Refit_TopM1,type='normalized'))
dev.off()

#F. Summary for Best Models:
summ_st = as.data.frame(capture.output(summary(Refit_TopM1)))
write.csv(summ_st, file="Results/steerei/FullModel/refit_Summary_AjustedBestModels_steerei.csv", row.names = TRUE)

#G. Save as Rdata:
save(Refit_TopM1, file="./Results/steerei/FullModel/refit_REML_TopM1_steerei.RData")



##### 8. COEFFICIENT OF DETERMINATION R² -----
#Not applicable in GLS models, only LME objects
#ML models and REML models show differences < 0.01 or <1%. I chose the REML for comaparison with Estimates.

#A. RΣ2, the proportion of generalized variance explained by the fixed predictors - Jaeger et al. (2017)
r2_steerei = r2beta(Refit_TopM1, method = 'sgv', partial = 'TRUE')
r2_steerei
#               Effect   Rsq upper.CL lower.CL
#1         Model 0.668    0.732    0.599
#2 riverdistance 0.565    0.645    0.479
#3    wetlands_L 0.099    0.194    0.031

#B. Nakagawa and Schielzeth (2013)
#marginal R² statistic = variance explained by fixed effects
#conditional R² statistic variance explained by the entire model (fixed & random effects)
r2beta(Refit_TopM1, method = 'nsj', partial = 'TRUE')

#C. Based on Nakagawa et al. (2017)
MuMIn::r.squaredGLMM(Refit_TopM1)
#      R2m      R2c
#  0.06557444 0.8566096

#D. Save the results
write.csv(as.data.frame(r2_steerei), row.names=TRUE, file="Results/steerei/FullModel/refit_BestModel_PartialR2_REML_steerei.csv")



#### 9. PLOT RADAR CHART FOR IMPORTANCE VARIABLES ----
#A. Load if it is necessary:
load(file="./Results/steerei/FullModel/BM_steerei.RData")

#A. Calculate variable importance using best models. Need to be more than 1 model. Not used in this case, because the results will be the same than using uncorrelated models.
impor_best = importance(BM[BM$delta <= 2, ])
impor_best
write.csv(impor_best, "./Results/steerei/FullModel/ImportanceVariables_Best_steerei.csv")

#B. Create a dataframe for ploting results:
impor_best

impor_best2 = impor_best
df = as.data.frame(impor_best2)
df[nrow(df) + 1,] = 0
df[nrow(df) + 1,] = 0
df[nrow(df) + 1,] = 0
row.names(df)[5:7] = c("wetlands_VH", "topography", "eucl_dist")
df = as.data.frame(t(df))
group = c("Best Models")
df = cbind(group,df)
df

#C. Plot radar:
plot_radar = 
  ggradar(df,
          base.size = 2,
          values.radar = c("0%", "50%", "100%"),
          axis.labels = c( "PET", "Precipitation", "River", "Temperature", "Habitat", "Topography", "Euclidean"),
          grid.min = 0,
          grid.mid = 0.5,
          grid.max = 1,
          grid.label.size = 4,
          axis.label.size = 5,
          axis.line.colour = "black",
          group.line.width = 1.5,
          group.point.size = 5,
          group.colours = c("blue"),
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

#D. Verify radar:
plot_radar

#E. Save
pdf("./Results/steerei/Figures/RelativePredictorWeight_MLPE_AllVars_steerei.pdf")
plot_radar
dev.off()


#F. Calculate variable importance using uncorrlebest models. Need to be more than 1 model:
impor_best = importance(BM)
impor_best
write.csv(impor_best, "./Results/steerei/FullModel/ImportanceVariables_Best_steerei.csv")

#G. Create a dataframe for ploting results:
impor_best

impor_best2 = impor_best
df = as.data.frame(impor_best2)
df = as.data.frame(t(df))
group = c("Uncorrelated Models")
df = cbind(group,df)
df

#C. Reorder like simonsi:
#steerei order = c( "Habitat", "Productivity", "River", "Euclidean", "Topography")
names(df)
df = df[,c(1,3,6,2,5,4)]

#H. Plot radar:
plot_radar = 
  ggradar(df,
          base.size = 2,
          values.radar = c("0%", "50%", "100%"),
          axis.labels = c( "Habitat", "Productivity", "River", "Euclidean", "Topography"),
          grid.min = 0,
          grid.mid = 0.5,
          grid.max = 1,
          grid.label.size = 4,
          axis.label.size = 5,
          axis.line.colour = "black",
          group.line.width = 1.5,
          group.point.size = 5,
          group.colours = c("blue"),
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

#I. Verify radar:
plot_radar

#J. Save
pdf("./Results/steerei/Figures/RelativePredictorWeight_MLPE_AllVars_steerei_uncorrelated.pdf")
plot_radar
dev.off()




#### 10. PLOT UNIVARIATES PREDICTORS USING BEST MODEL ----
#A. LOAD UNSCALED DATA
DT_unscaled =  read.csv("Metafiles/MLPE_table_steerei.csv", row.names = 1)
names(DT_unscaled)

#B. Verify variables in Best Model
TopM1
#riverdistance + wetlands_L

#C. Build model with River distance ----
M1 = nlme::lme(RELgen ~ riverdistance,
                random = ~ 1|POP,
                correlation = corMLPE(form = ~ from_ID + to_ID),
                data = DT_unscaled, method = "REML")

# Decorrelate residulas
dec_resids = decorltd_res_lme(M1)
#dec_resids$dist_interval = DT_unscaled$dist_interval
dec_resids$geoDist = DT_unscaled$eucl_dist

# Change covar in dataframe to plot
cov = DT_unscaled$riverdistance
df2plot = as_tibble(data.frame(dec_resids, covar = cov))
names(df2plot)


## Plot relationship
plotA = ggplot(df2plot, 
                aes(x = covar, y = pr.fixed)) + 
  geom_point(aes(y=pr.fixed), alpha=0.3, size = 4) +
  geom_line(aes(y=fit), size=1) +
  geom_rug(sides = "b", alpha = 0.02) +
  ylab("Relatedness (decorrelated)") + 
  annotate("text", x=c(5), y=c(15),label=c(""), size=7) +
  xlab("River Distance") + 
  theme_bw() + 
  theme(axis.title.y = element_text(size=15, color = "black", face = "bold"),
        axis.title.x = element_text(size=15, color = "black", face = "bold")) +
  theme(axis.text = element_text(face = "bold", color = "black", size = 15))

#Verify
plotA

pdf("./Results/steerei/Figures/REL_deco_RiverDistance_steerei.pdf")
plotA
dev.off()




#D. Build model with Habitat ----
M2 = nlme::lme(RELgen ~ wetlands_L,
                random = ~ 1|POP,
                correlation = corMLPE(form = ~ from_ID + to_ID),
                data = DT_unscaled, method = "REML")

# Decorrelate residulas
dec_resids = decorltd_res_lme(M2)
#dec_resids$dist_interval = DT_unscaled$dist_interval
dec_resids$geoDist = DT_unscaled$eucl_dist

# Change covar in dataframe to plot
cov = DT_unscaled$wetlands_L
df2plot = as_tibble(data.frame(dec_resids, covar = cov))
names(df2plot)


# Plot relationship
plotB = ggplot(df2plot, 
                aes(x = covar, y = pr.fixed)) + 
  geom_point(aes(y=pr.fixed), alpha=0.3, size = 4) +
  geom_line(aes(y=fit), size=1) +
  geom_rug(sides = "b", alpha = 0.02) +
  ylab("Relatedness (decorrelated)") + 
  annotate("text", x=c(5), y=c(15),label=c(""), size=7) +
  xlab("Habitat Resistance") + 
  theme_bw() + 
  theme(axis.title.y = element_text(size=15, color = "black", face = "bold"),
        axis.title.x = element_text(size=15, color = "black", face = "bold")) +
  theme(axis.text = element_text(face = "bold", color = "black", size = 15))

#Verify
plotB

pdf("./Results/steerei/Figures/REL_deco_Habitat_steerei.pdf")
plotB
dev.off()


#E. Build model with Eucliadean Distance ----
M4 = nlme::lme(RELgen ~ eucl_dist,
                random = ~ 1|POP,
                correlation = corMLPE(form = ~ from_ID + to_ID),
                data = DT_unscaled, method = "REML")

## Decorrelate residulas
dec_resids = decorltd_res_lme(M4)
#dec_resids$dist_interval = DT_unscaled$dist_interval
dec_resids$geoDist = DT_unscaled$eucl_dist

## Change covar in dataframe to plot
cov = DT_unscaled$eucl_dist
df2plot = as_tibble(data.frame(dec_resids, covar = cov))
names(df2plot)


## Plot relationship
plotD = ggplot(df2plot, 
                aes(x = covar, y = pr.fixed)) + 
  geom_point(aes(y=pr.fixed), alpha=0.3, size = 4) +
  geom_line(aes(y=fit), size=1) +
  geom_rug(sides = "b", alpha = 0.02) +
  ylab("Relatedness (decorrelated)") + 
  annotate("text", x=c(5), y=c(15),label=c(""), size=7) +
  xlab("Euclidean Distance") + 
  theme_bw() + 
  theme(axis.title.y = element_text(size=15, color = "black", face = "bold"),
        axis.title.x = element_text(size=15, color = "black", face = "bold")) +
  theme(axis.text = element_text(face = "bold", color = "black", size = 15))

plotD

pdf("./Results/steerei/Figures/REL_deco_EuclideanDis_steerei.pdf")
plotD
dev.off()

#END
  
