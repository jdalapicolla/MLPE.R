#####################  LANDSCAPE GENOMICS TUTORIAL   ##################
###################   STEP 03: MLPE MODELS WITH LME  ##################

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
library(ggradar)
library(r2glmm)
library(snow)
library(parallel)
library(GGally)
library(reshape2)
library(gridExtra)
library(RColorBrewer)
library(ggradar)
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
mlpe_steerei[, 5:15] = as.data.frame(scale(mlpe_steerei[, 5:15]))
head(mlpe_steerei)

#C. Creating formula for MLPE:
#select variables for models:
vars_mlpe_st = names(mlpe_steerei)[5:15]
vars_mlpe_st
#use REL distance plus varaibles of your choice keeping one for geographic distance (Null Model):
vars_mlpe_st_rel = vars_mlpe_st[c(1:6,10)]
mlpe_formula_rel_st = as.formula(paste("RELgen ~ ",paste(vars_mlpe_st_rel, collapse = " + "),sep = ""))
mlpe_formula_rel_st #7 variables




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
acf(RES) #OK

#Don't need any other action!


##### 4. RUN MODEL SELECTION WITH PARALLELIZED DREDGE ---- 
#A. Check maximum correlation
max.r(Fullmodel) ## 0.7920

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
diff1 ##12.09021 secs

#E. Number of models
nrow(Allmodels) ## 128

#F. Save and load All Models:
save(Allmodels, file="./Results/steerei/FullModel/Allmodels_steerei.RData")
load(file = "./Results/steerei/FullModel/Allmodels_steerei.RData")

#G. Retrieve Not Collinear Models, with max.r <=0.6. Get Not Collinear Models using cluster too!
NCM = get.models(Allmodels, subset = max.r <= 0.6, cluster)

#H. Number of Not Collinear models
length(NCM) #25

#I. Save and load Not Collinear models
save(NCM, file="./Results/steerei/FullModel/NCM_steerei.RData")
load(file="./Results/steerei/FullModel/NCM_steerei.RData")

#J. Select Best Models using AIC
BM = model.sel(NCM, rank=AIC)
nrow(BM) #25
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
Refit_TopM1 <- nlme::lme(RELgen ~ PET + precipitation + riverdistance,
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
plot(mlpe_steerei$PET, RES) ; abline(0,0, col="red") #ok
plot(mlpe_steerei$precipitation, RES) ; abline(0,0, col="red") #ok
plot(mlpe_steerei$riverdistance, RES) ; abline(0,0, col="red") #ok

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



#H. TopM2 ----
TopM2 <- get.models(BM, subset = 2)[[1]]
summary(TopM2)
intervals(TopM2)
TopM2

#I. Refit TopM2
Refit_TopM2 <- nlme::lme(RELgen ~ PET + precipitation + riverdistance + temperature,
                         random = ~1|POP,
                         correlation = corMLPE(form = ~ from_ID + to_ID),
                         data = mlpe_steerei,
                         method = "REML")
summary(Refit_TopM2)

#J. Check autocorrelation of residuals
RES <- residuals(Refit_TopM2, type="normalized")
FIT <- fitted(Refit_TopM2)
plot(FIT, RES) ; abline(0,0, col="red") #ok
acf(RES) #ok

#K. Plot predictors x residuals
plot(mlpe_steerei$PET, RES) ; abline(0,0, col="red") #ok
plot(mlpe_steerei$precipitation, RES) ; abline(0,0, col="red") #ok
plot(mlpe_steerei$riverdistance, RES) ; abline(0,0, col="red") #ok
plot(mlpe_steerei$temperature, RES) ; abline(0,0, col="red") #ok

#L. Save ACF graph:
acf(resid(Refit_TopM2,type='normalized')) ## Spatial dependence pattern
pdf("Results/steerei/Figures/refit_ACF_LME_steerei_TopM2.pdf", onefile = T)
acf(resid(Refit_TopM2,type='normalized'))
dev.off()

#M. Summary for Best Models:
summ_st2 = as.data.frame(capture.output(summary(Refit_TopM2)))
write.csv(summ_st2, file="Results/steerei/FullModel/refit_Summary_AjustedBestModels_2nd_steerei.csv", row.names = TRUE)

#N. Save as Rdata:
save(Refit_TopM2, file="./Results/steerei/FullModel/refit_REML_TopM2_steerei.RData")

load(file="./Results/steerei/FullModel/refit_REML_TopM2_steerei.RData")


##### 6. CALCULATE AND PLOT THE ESTIMATES -----
#A. Using best models with REML
MA_REML = model.avg(Refit_TopM1, Refit_TopM2)
confint(MA_REML) ## Confidence Intervals for model-averaged coeficients
as.data.frame(MA_REML$coefficients) ## Mean

#B. Save as dataframe to plot estimates
df1 = as.data.frame(MA_REML$coefficients[1, ])
rownames(df1)
rownames(df1) <- c("Intercept",
                   "PET",
                   "Precipitation",
                   "River Distance",
                   "Temperature")

df1$term <- rownames(df1)
colnames(df1) <- c("coef", "term")
df2 <- bind_cols(df1$term, df1$coef, as.data.frame(confint(MA_REML)[,1]), as.data.frame(confint(MA_REML)[,2]))
colnames(df2) <- c("term", "estimate", "conf.low", "conf.high" )
rownames(df2) <- NULL
head(df2)

#C. Save and load if necessary:
write.csv(df2, "./Results/steerei/FullModel/refit_model_avg_df_steerei.csv", row.names = F)

#D. Load model average dataframe
model_avg_estimates = read.csv("./Results/steerei/FullModel/refit_model_avg_df_steerei.csv", h=T)
head(model_avg_estimates)

#E. Plot
p1 = ggcoef(model_avg_estimates,
            conf.int = TRUE,
            conf.level = 0.95,
            size = 4.5,
            vline = TRUE,
            vline_intercept = "auto",
            vline_color = "gray50",
            vline_linetype = "dashed",
            vline_size = 1,
            errorbar_color = "black",
            errorbar_height = 0,
            errorbar_linetype = "solid",
            errorbar_size = 0.6,
            sort = "ascending",
            mapping = aes_string(y = "term", x = "estimate")) +
  labs(x = "Estimates", y = "Terms") +
  theme_bw() +
  theme(axis.text=element_text(size=14),
      axis.title=element_text(size=14,face="bold"))

#F. Verify Plot 
p1


#G. Plot confidence intervals of estimates
pdf("./Results/steerei/Figures/Estimates_avg_df_steerei.pdf", onefile = T)
p1
dev.off()


##### 7. TEST SIGNIFICANCE OF PREDICTORS - LRT ----
#A. Models with ML because REML do note provide AIC values
LRT = drop1(TopM1, test="Chisq")
LRT

#B. Save Results
save(LRT, file = "./Results/steerei/FullModel/LRT_TopM1_steerei.RData")
load(file = "./Results/steerei/FullModel/LRT_TopM1_steerei.RData")
write.csv(LRT, file="./Results/steerei/FullModel/LRT_Model1_steerei.csv", row.names = TRUE)

#C. TopM2
LRT = drop1(TopM2, test="Chisq")
LRT

#D. Save Results
save(LRT, file = "./Results/steerei/FullModel/LRT_TopM2_steerei.RData")
load(file = "./Results/steerei/FullModel/LRT_TopM2_steerei.RData")
write.csv(LRT, file="./Results/steerei/FullModel/LRT_Model2_steerei.csv", row.names = TRUE)

#E. LTR for Nested Models:
compare_bestmodels = anova.lme(TopM1, TopM2)
compare_bestmodels
#F. Save LTR results:
write.csv(as.data.frame(compare_bestmodels), file="Results/steerei/FullModel/LTR_BestModels_steerei.csv", row.names = TRUE)


##### 8. COEFFICIENT OF DETERMINATION R² -----
#Not applicable in GLS models, only LME objects
#ML models and REML models show differences < 0.01 or <1%. I chose the REML for comaparison with Estimates.

#A. RΣ2, the proportion of generalized variance explained by the fixed predictors - Jaeger et al. (2017)
r2_steerei = r2beta(Refit_TopM1, method = 'sgv', partial = 'TRUE')
r2_steerei

#B. Nakagawa and Schielzeth (2013)
#marginal R² statistic = variance explained by fixed effects
#conditional R² statistic variance explained by the entire model (fixed & random effects)
r2beta(Refit_TopM1, method = 'nsj', partial = 'TRUE')

#C. Based on Nakagawa et al. (2017)
MuMIn::r.squaredGLMM(Refit_TopM1)
#      R2m      R2c
# 0.023964 0.818  671

#D. Save the results
write.csv(as.data.frame(r2_steerei), row.names=TRUE, file="Results/steerei/FullModel/refit_BestModel_PartialR2_REML_steerei.csv")



#### 9. PLOT RADAR CHART FOR IMPORTANCE VARIABLES ----
#A. Load if it is necessary:
load(file="./Results/steerei/FullModel/BM_steerei.RData")

#A. Calculate variable importance using best models. Need to be more than 1 model:
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



#### 10. PLOT UNIVARIATES PREDICTORS USING BEST MODEL ----
#A. LOAD UNSCALED DATA
DT_unscaled =  read.csv("Metafiles/MLPE_table_steerei.csv", row.names = 1)
names(DT_unscaled)

#B. Verify variables in Best Model
TopM1
#PET + precipitation + riverdistance

#C. Build model with PET ----
M1 = nlme::lme(RELgen ~ PET,
                random = ~ 1|POP,
                correlation = corMLPE(form = ~ from_ID + to_ID),
                data = DT_unscaled, method = "REML")

# Decorrelate residulas
dec_resids = decorltd_res_lme(M1)
dec_resids$dist_interval = DT_unscaled$dist_interval
dec_resids$geoDist = DT_unscaled$geoDist

# Change covar in dataframe to plot
cov = DT_unscaled$PET
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
  xlab("PET Dissimilarity") + 
  theme_bw() + 
  theme(axis.title.y = element_text(size=15, color = "black", face = "bold"),
        axis.title.x = element_text(size=15, color = "black", face = "bold")) +
  theme(axis.text = element_text(face = "bold", color = "black", size = 15))

#Verify
plotA

pdf("./Results/steerei/Figures/REL_deco_PET_steerei.pdf")
plotA
dev.off()




#D. Build model with Precipitation ----
M2 = nlme::lme(RELgen ~ precipitation,
                random = ~ 1|POP,
                correlation = corMLPE(form = ~ from_ID + to_ID),
                data = DT_unscaled, method = "REML")

# Decorrelate residulas
dec_resids = decorltd_res_lme(M2)
dec_resids$dist_interval = DT_unscaled$dist_interval
dec_resids$geoDist = DT_unscaled$geoDist

# Change covar in dataframe to plot
cov = DT_unscaled$precipitation
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
  xlab("Precipitation Dissimilarity") + 
  theme_bw() + 
  theme(axis.title.y = element_text(size=15, color = "black", face = "bold"),
        axis.title.x = element_text(size=15, color = "black", face = "bold")) +
  theme(axis.text = element_text(face = "bold", color = "black", size = 15))

#Verify
plotB

pdf("./Results/steerei/Figures/REL_deco_Precipi_steerei.pdf")
plotB
dev.off()



#E. Build model with River Distance ----
M3 = nlme::lme(RELgen ~ riverdistance,
                random = ~ 1|POP,
                correlation = corMLPE(form = ~ from_ID + to_ID),
                data = DT_unscaled, method = "REML")

# Decorrelate residulas
dec_resids = decorltd_res_lme(M3)
dec_resids$dist_interval = DT_unscaled$dist_interval
dec_resids$geoDist = DT_unscaled$geoDist

# Change covar in dataframe to plot
cov = DT_unscaled$riverdistance
df2plot = as_tibble(data.frame(dec_resids, covar = cov))
names(df2plot)


## Plot relationship
plotC = ggplot(df2plot, 
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
plotC

pdf("./Results/steerei/Figures/REL_deco_RiverDis_steerei.pdf")
plotC
dev.off()


#F. Build model with Eucliadean Distance ----
M4 = nlme::lme(RELgen ~ eucl_dist,
                random = ~ 1|POP,
                correlation = corMLPE(form = ~ from_ID + to_ID),
                data = DT_unscaled, method = "REML")

## Decorrelate residulas
dec_resids = decorltd_res_lme(M4)
dec_resids$dist_interval = DT_unscaled$dist_interval
dec_resids$geoDist = DT_unscaled$geoDist

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


#G. Build model with Temperature
M5 = nlme::lme(RELgen ~ temperature,
                random = ~ 1|POP,
                correlation = corMLPE(form = ~ from_ID + to_ID),
                data = DT_unscaled, method = "REML")

## Decorrelate residulas
dec_resids = decorltd_res_lme(M5)
dec_resids$dist_interval = DT_unscaled$dist_interval
dec_resids$geoDist = DT_unscaled$geoDist

## Change covar in dataframe to plot
cov = DT_unscaled$temperature
df2plot = as_tibble(data.frame(dec_resids, covar = cov))
names(df2plot)


## Plot relationship
plotE = ggplot(df2plot, 
                aes(x = covar, y = pr.fixed)) + 
  geom_point(aes(y=pr.fixed), alpha=0.3, size = 4) +
  geom_line(aes(y=fit), size=1) +
  geom_rug(sides = "b", alpha = 0.02) +
  ylab("Relatedness (decorrelated)") + 
  annotate("text", x=c(5), y=c(15),label=c(""), size=7) +
  xlab("Temperature Dissimilarity") + 
  theme_bw() + 
  theme(axis.title.y = element_text(size=15, color = "black", face = "bold"),
        axis.title.x = element_text(size=15, color = "black", face = "bold")) +
  theme(axis.text = element_text(face = "bold", color = "black", size = 15))

plotE

pdf("./Results/steerei/Figures/REL_deco_Temperature_steerei.pdf")
plotE
dev.off()

#END
  
