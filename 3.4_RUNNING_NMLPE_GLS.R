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
library(ggradar)
library(r2glmm)
library(snow)
library(parallel)
library(GGally)
library(reshape2)
library(gridExtra)
library(RColorBrewer)
library(scales)

### Load auxiliary functions:
## Check max correlation
max.r = function(x){
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

# For LME objects
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
mlpe_simonsi = read.csv("Metafiles/MLPE_table_simonsi.csv", row.names = 1)
head(mlpe_simonsi)

#B. Scale predictors/variables for use in models
mlpe_simonsi[, 5:(length(mlpe_simonsi)-1)] = as.data.frame(scale(mlpe_simonsi[, 5:(length(mlpe_simonsi)-1)]))
head(mlpe_simonsi)

#C. Creating formula for MLPE:
#select variables for models:
vars_mlpe_st = names(mlpe_simonsi)[5:(length(mlpe_simonsi)-1)]
vars_mlpe_st

#use REL distance plus variables of your choice keeping the geographic distance (Null Model):
vars_mlpe_st_rel = vars_mlpe_st[c(1:5)]
mlpe_formula_rel_st = as.formula(paste("RELgen ~ ",paste(vars_mlpe_st_rel, collapse = " + "),sep = ""))
mlpe_formula_rel_st #5 variables



##### 2. RUNNING FULL MODEL -----------
#A. Full model:
Fullmodel = nlme::lme(mlpe_formula_rel_st,
                      random = ~1|POP,
                      correlation = corMLPE(form = ~ from_ID + to_ID),
                      data = mlpe_simonsi,
                      method = "ML") #For model selection is ML

#B. Check model:
summary(Fullmodel)
RES <- residuals(Fullmodel, type="normalized")
FIT <- fitted(Fullmodel)
plot(FIT, RES) ; abline(0,0, col="red") #OK
acf(RES) #there is no spatial autocorrelation in residuals


##3. TESTING THE GLS MODEL BECAUSE WE HAVE ONLY ONE POPULATION ----
#A. Compare LME with 1 population to GLS:
#LME
Fullmodel

#GLS
Fullmodel_2 = nlme::gls(mlpe_formula_rel_st,
                        correlation = corMLPE(form = ~ from_ID + to_ID),
                        data = mlpe_simonsi,
                        method = "ML")

## Check model
summary(Fullmodel_2)
RES <- residuals(Fullmodel_2, type="normalized")
FIT <- fitted(Fullmodel_2)
plot(FIT, RES) ; abline(0,0, col="red")
acf(RES) ##there is no spatial autocorrelation in residuals

#B. Comparing LME and GLS results by LTR:
compare_bestmodels = anova.lme(Fullmodel, Fullmodel_2)
compare_bestmodels #models are not different (p = 0.999)

#C. Compare Importance of Variables. Maximum of 9 variables (15 observation by predictor)
#LME:
Allmodels = dredge(Fullmodel, rank = "AIC", m.lim=c(0, 9), extra= c(max.r))
NCM = get.models(Allmodels, subset = max.r<=0.6)
BM = model.sel(NCM)
impor_uncor_lme = importance(BM)

#GLS
Allmodels = dredge(Fullmodel_2, rank = "AIC", m.lim=c(0, 9), extra= c(max.r))
NCM = get.models(Allmodels, subset = max.r<=0.6)
BM = model.sel(NCM)
impor_uncor_gls = importance(BM)

#Verify
impor_uncor_lme
impor_uncor_gls
# variation less than 1%. So we keep the lme to compare to P. steerei that needs lme! 


#4. REFIT TO NESTED MLPE IF YOU HAVE SPATIAL STRUCTURE IN RESIDUALS
#A. Sort data so we can see acf better
#mlpe_simonsi2 = mlpe_simonsi[order(mlpe_simonsi$from_ID, mlpe_simonsi$to_ID),]

#B. Identify which individuals are from same location based on pairwise geographic distance
#ulab = unique(c(as.character(mlpe_simonsi2$from_ID), as.character(mlpe_simonsi2$to_ID)))
#dis = matrix(0, length(ulab), length(ulab))
#colnames(dis) = ulab
#rownames(dis) = ulab
#dis[1:10, 1:10]

#for(i in 1:nrow(mlpe_simonsi2)){
 # dis[mlpe_simonsi2$from_ID[i], mlpe_simonsi2$to_ID[i]] = 
  #  dis[mlpe_simonsi2$to_ID[i],mlpe_simonsi2$from_ID[i]] = 
  #  mlpe_simonsi2$eucl_dist[i]
  #location = cutree(hclust(as.dist(dis)), h=0) #assign individuals to unique locations
#}

#C. Plot dendogram
#location
#plot(hclust(as.dist(dis)), h=0)


#D. Run nested MLPE --- 
#Model = nlme::gls(mlpe_formula_rel_st,
 #                  correlation = corNMLPE2(form = ~ from_ID + to_ID, 
  #                                         clusters = location),
   #               data = mlpe_simonsi2, method = "ML")

#E. Check residuals  
#RES <- residuals(Model, type="normalized")
#FIT <- fitted(Model)
#plot(FIT, RES) ; abline(0,0, col="red") #OK
#acf(RES) #small spatial correlation




##### 4. CREATE A CLUSTER TO RUN MLPE IN PARALLEL ------
#A. Choose the full model to continue the analysis: LME (Fullmodel), GLS (Fullmodel_2), or Nested MLPE (Model). In this case, LME (if you choose gls and other some fuctions habe to be adapted such as the r² estimates)
Fullmodel = Fullmodel

#A. Make cluster
cluster = makeCluster(4, type = "SOCK")  ## also need snow installed

#B. Use clusterExport to send data to global environment (aka ‘workspace’) of each node. 
clusterExport(cluster, c("mlpe_simonsi")) #add other objects if you need

#C. Load packages in the cluster
clusterEvalQ(cluster,
             c(library(nlme), library(MuMIn), library(corMLPE)))




##### 5. RUN MODEL SELECTION WITH PARALLELIZED DREDGE ---- 
#A. Check maximum correlation
max.r(Fullmodel) ## 0.9966734

#B. Specify the number of predictor variables and including the max.r function. 15 observations to each predictor:
max_var = length(mlpe_simonsi[,1])/15
max_var #9
nrow(mlpe_simonsi) ## 136
options(na.action = na.fail)

#C. Run pdredge:
start_time = Sys.time()
Allmodels = MuMIn::pdredge(Fullmodel, rank = "AIC", 
                           m.lim=c(0, 9), extra = c(max.r), cluster)
end_time = Sys.time()

#D. Compute run time
start_time1 = start_time
end_time1 = end_time
diff1 = end_time - start_time
diff1 ##7.350864 secs

#E. Number of models
nrow(Allmodels) ## 32

#F. Save and load All Models:
save(Allmodels, file="./Results/simonsi/FullModel/Allmodels_simonsi.RData")
load(file = "./Results/simonsi/FullModel/Allmodels_simonsi.RData")

#G. Retrieve Not Collinear Models, with max.r <=0.6. Get Not Collinear Models using cluster too!
NCM = get.models(Allmodels, subset = max.r <= 0.6, cluster)

#H. Number of Not Collinear models
length(NCM) #16

#I. Save and load Not Collinear models
save(NCM, file="./Results/simonsi/FullModel/NCM_simonsi.RData")
load(file="./Results/simonsi/FullModel/NCM_simonsi.RData")

#J. Select Best Models using AIC
BM = model.sel(NCM, rank=AIC)
nrow(BM) #16
#Check best models
best_models = BM[BM$delta <= 2, ]
best_models

#K. Save and load Model selection
save(BM, file="./Results/simonsi/FullModel/BM_simonsi.RData")
load(file="./Results/simonsi/FullModel/BM_simonsi.RData")


#L. Save model selection as table/dataframe
df_model_sel = as.data.frame(BM)
head(df_model_sel)
write.csv(df_model_sel, "./Results/simonsi/FullModel/ModelsSelection_AIC_simonsi.csv", row.names = T)


#M. Save best models tables
as.data.frame(best_models)
write.csv(as.data.frame(best_models), "./Results/simonsi/FullModel/BestModels_AIC_simonsi.csv", row.names = F)





##### 6. REFIT BEST MODELS USING REML AND MODEL VALIDATION ----
#A. TopM1 ----
TopM1 = get.models(BM, subset = 1)[[1]]
summary(TopM1)
intervals(TopM1, which = "fixed")
TopM1

#B. Refit TopM1
Refit_TopM1 = nlme::lme(RELgen ~ productivity_dist + wetlands_L,
                        random = ~1|POP,
                        correlation = corMLPE(form = ~ from_ID + to_ID),
                        data = mlpe_simonsi,
                        method = "REML") #best estimates
summary(Refit_TopM1)

#C. Check autocorrelation of residuals
RES <- residuals(Refit_TopM1, type="normalized")
FIT <- fitted(Refit_TopM1)
plot(FIT, RES) ; abline(0,0, col="red") #ok
acf(RES) #ok

#D. Plot predictors x residuals
plot(mlpe_simonsi$productivity_dist, RES) ; abline(0,0, col="red") #ok
plot(mlpe_simonsi$wetlands_L, RES) ; abline(0,0, col="red") #ok

#E. Save ACF graph:
acf(resid(Refit_TopM1,type='normalized')) ## Spatial dependence pattern
#Save graphs
pdf("Results/simonsi/Figures/refit_ACF_LME_simonsi_TopM1.pdf", onefile = T)
acf(resid(Refit_TopM1,type='normalized'))
dev.off()

#F. Summary for Best Models:
summ_st = as.data.frame(capture.output(summary(Refit_TopM1)))
write.csv(summ_st, file="Results/simonsi/FullModel/refit_Summary_AjustedBestModels_simonsi.csv", row.names = TRUE)

#G. Save as Rdata:
save(Refit_TopM1, file="./Results/simonsi/FullModel/refit_REML_TopM1_simonsi.RData")



#H. TopM2 ----
TopM2 <- get.models(BM, subset = 2)[[1]]
summary(TopM2)
intervals(TopM2, which = "fixed")
TopM2

#I. Refit TopM2
Refit_TopM2 = nlme::lme(RELgen ~ wetlands_L,
                        random = ~1|POP,
                        correlation = corMLPE(form = ~ from_ID + to_ID),
                        data = mlpe_simonsi,
                        method = "REML") #best estimates
summary(Refit_TopM2)

#J. Check autocorrelation of residuals
RES <- residuals(Refit_TopM2, type="normalized")
FIT <- fitted(Refit_TopM2)
plot(FIT, RES) ; abline(0,0, col="red") #ok
acf(RES) #ok

#K. Plot predictors x residuals
plot(mlpe_simonsi$wetlands_L, RES) ; abline(0,0, col="red") #ok

#L. Save ACF graph:
acf(resid(Refit_TopM2,type='normalized')) ## Spatial dependence pattern
pdf("Results/simonsi/Figures/refit_ACF_LME_simonsi_TopM2.pdf", onefile = T)
acf(resid(Refit_TopM2,type='normalized'))
dev.off()

#M. Summary for Best Models:
summ_st2 = as.data.frame(capture.output(summary(Refit_TopM2)))
write.csv(summ_st2, file="Results/simonsi/FullModel/refit_Summary_AjustedBestModels_2nd_simonsi.csv", row.names = TRUE)

#N. Save as Rdata:
save(Refit_TopM2, file="./Results/simonsi/FullModel/refit_REML_TopM2_simonsi.RData")


#O. TopM3 ----
TopM3 = get.models(BM, subset = 3)[[1]]
summary(TopM3)
intervals(TopM3, which = "fixed")
TopM3

#P. Refit TopM3
Refit_TopM3 = nlme::lme(RELgen ~ wetlands_L + riverdistance,
                                     random = ~1|POP,
                                     correlation = corMLPE(form = ~ from_ID + to_ID),
                                     data = mlpe_simonsi,
                                     method = "REML") #best estimates
summary(Refit_TopM3)

#Q. Check autocorrelation of residuals
RES <- residuals(Refit_TopM3, type="normalized")
FIT <- fitted(Refit_TopM3)
plot(FIT, RES) ; abline(0,0, col="red") #ok
acf(RES) #ok

#R. Plot predictors x residuals
plot(mlpe_simonsi$wetlands_L, RES) ; abline(0,0, col="red") #ok
plot(mlpe_simonsi$riverdistance, RES) ; abline(0,0, col="red") #ok


#S. Save ACF graph:
acf(resid(Refit_TopM3,type='normalized')) ## Spatial dependence pattern
pdf("Results/simonsi/Figures/refit_ACF_LME_simonsi_TopM3.pdf", onefile = T)
acf(resid(Refit_TopM3,type='normalized'))
dev.off()

#T. Summary for Best Models:
summ_st3 = as.data.frame(capture.output(summary(Refit_TopM3)))
write.csv(summ_st3, file="Results/simonsi/FullModel/refit_Summary_AjustedBestModels_3rd_simonsi.csv", row.names = TRUE)

#U. Save as Rdata:
save(Refit_TopM3, file="./Results/simonsi/FullModel/refit_REML_TopM3_simonsi.RData")



#V. TopM4 ----
TopM4 = get.models(BM, subset = 4)[[1]]
summary(TopM4)
intervals(TopM4, which = "fixed")
TopM4


#V. Refit TopM4
Refit_TopM4 = nlme::lme(RELgen ~ wetlands_L + riverdistance + productivity_dist,
                        random = ~1|POP,
                        correlation = corMLPE(form = ~ from_ID + to_ID),
                        data = mlpe_simonsi,
                        method = "REML") #best estimates
summary(Refit_TopM4)

#J. Check autocorrelation of residuals
RES <- residuals(Refit_TopM4, type="normalized")
FIT <- fitted(Refit_TopM4)
plot(FIT, RES) ; abline(0,0, col="red") #ok
acf(RES) #ok

#K. Plot predictors x residuals
plot(mlpe_simonsi$wetlands_L, RES) ; abline(0,0, col="red") #ok
plot(mlpe_simonsi$riverdistance, RES) ; abline(0,0, col="red") #ok
plot(mlpe_simonsi$productivity_dist, RES) ; abline(0,0, col="red") #ok


#L. Save ACF graph:
acf(resid(Refit_TopM4,type='normalized')) ## Spatial dependence pattern
pdf("Results/simonsi/Figures/refit_ACF_LME_simonsi_TopM4.pdf", onefile = T)
acf(resid(Refit_TopM4,type='normalized'))
dev.off()

#M. Summary for Best Models:
summ_st4 = as.data.frame(capture.output(summary(Refit_TopM4)))
write.csv(summ_st4, file="Results/simonsi/FullModel/refit_Summary_AjustedBestModels_4TH_simonsi.csv", row.names = TRUE)

#N. Save as Rdata:
save(Refit_TopM4, file="./Results/simonsi/FullModel/refit_REML_TopM4_simonsi.RData")



##### 7. CALCULATE AND PLOT THE ESTIMATES -----
#A. Using best models with REML
MA_REML = model.avg(Refit_TopM1, Refit_TopM2, Refit_TopM3, Refit_TopM4)
confint(MA_REML) ## Confidence Intervals for model-averaged coeficients
as.data.frame(MA_REML$coefficients) ## Mean

#B. Save as dataframe to plot estimates
df1 = as.data.frame(MA_REML$coefficients[1, ])
rownames(df1)
rownames(df1) <- c("Intercept",
                   "Habitat",
                   "Productivity",
                   "River Distance")

df1$term <- rownames(df1)
colnames(df1) <- c("coef", "term")
df2 <- bind_cols(df1$term, df1$coef, as.data.frame(confint(MA_REML)[,1]), as.data.frame(confint(MA_REML)[,2]))
colnames(df2) <- c("term", "estimate", "conf.low", "conf.high" )
rownames(df2) <- NULL
head(df2)

#C. Save and load if necessary:
write.csv(df2, "./Results/simonsi/FullModel/refit_model_avg_df_simonsi.csv", row.names = F)

#D. Load model average dataframe
model_avg_estimates = read.csv("./Results/simonsi/FullModel/refit_model_avg_df_simonsi.csv", h=T)
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
pdf("./Results/simonsi/Figures/Estimates_avg_df_simonsi.pdf", onefile = T)
p1
dev.off()



##### 8. TEST SIGNIFICANCE OF PREDICTORS - LRT ----
#A. Models with ML because REML do note provide AIC values
LRT = drop1(TopM1, test="Chisq")
LRT

#B. Save Results
save(LRT, file = "./Results/simonsi/FullModel/LRT_TopM1_simonsi.RData")
load(file = "./Results/simonsi/FullModel/LRT_TopM1_simonsi.RData")
write.csv(LRT, file="./Results/simonsi/FullModel/LRT_Model1_simonsi.csv", row.names = TRUE)

#C. TopM2
LRT = drop1(TopM2, test="Chisq")
LRT

#D. Save Results
save(LRT, file = "./Results/simonsi/FullModel/LRT_TopM2_simonsi.RData")
load(file = "./Results/simonsi/FullModel/LRT_TopM2_simonsi.RData")
write.csv(LRT, file="./Results/simonsi/FullModel/LRT_Model2_simonsi.csv", row.names = TRUE)

#E. TopM3
LRT = drop1(TopM3, test="Chisq")
LRT

#F. Save Results
save(LRT, file = "./Results/simonsi/FullModel/LRT_TopM3_simonsi.RData")
load(file = "./Results/simonsi/FullModel/LRT_TopM3_simonsi.RData")
write.csv(LRT, file="./Results/simonsi/FullModel/LRT_Model3_simonsi.csv", row.names = TRUE)


#G. TopM4
LRT = drop1(TopM4, test="Chisq")
LRT

#F. Save Results
save(LRT, file = "./Results/simonsi/FullModel/LRT_TopM4_simonsi.RData")
load(file = "./Results/simonsi/FullModel/LRT_TopM4_simonsi.RData")
write.csv(LRT, file="./Results/simonsi/FullModel/LRT_Model4_simonsi.csv", row.names = TRUE)



#E. LTR for Nested Models:
compare_bestmodels = anova.lme(TopM1, TopM2)
compare_bestmodels # L.Ratio = 2.08432 p = 0.1488
compare_bestmodels2 = anova.lme(TopM3, TopM2)
compare_bestmodels2 #L.Ratio = 0.2794834 p = 0.597
compare_bestmodels2 = anova.lme(TopM4, TopM1)
compare_bestmodels2 #L.Ratio = 0.04543168 p = 0.8312
compare_bestmodels2 = anova.lme(TopM4, TopM2)
compare_bestmodels2 #L.Ratio = 2.129751 p = 0.3448
compare_bestmodels2 = anova.lme(TopM4, TopM3)
compare_bestmodels2 #L.Ratio = 1.850268 p = 0.1738


#F. Save LTR results:
#write.csv(as.data.frame(compare_bestmodels), file="Results/simonsi/FullModel/LTR_BestModels_simonsi.csv", row.names = TRUE)
#write.csv(as.data.frame(compare_bestmodels2), file="Results/simonsi/FullModel/LTR_BestModels_simonsi.csv", row.names = TRUE)




##### 9. COEFFICIENT OF DETERMINATION R² -----
#Not applicable in GLS models, only LME objects
#ML models and REML models show differences < 0.01 or <1%. I chose the REML for comparison with Estimates.
Refit_TopM1_lme = Refit_TopM2

#A. RΣ2, the proportion of generalized variance explained by the fixed predictors - Jaeger et al. (2017)
r2_simonsi = r2beta(Refit_TopM1_lme, method = 'sgv', partial = 'TRUE')
r2_simonsi

#              Effect   Rsq upper.CL lower.CL
#      Model    0.349    0.467    0.234
#  wetlands_L   0.349    0.467    0.234


#B. Nakagawa and Schielzeth (2013)
#marginal R² statistic = variance explained by fixed effects
#conditional R² statistic variance explained by the entire model (fixed & random effects)
r2beta(Refit_TopM1_lme, method = 'nsj', partial = 'TRUE')
#              Effect   Rsq upper.CL lower.CL
#      Model    0.338    0.457    0.223
#  wetlands_L   0.338    0.457    0.223

#C. Based on Nakagawa et al. (2017)
MuMIn::r.squaredGLMM(Refit_TopM1_lme)
#      R2m      R2c
# 0.3380082 0.3709024

#D. Save the results
write.csv(as.data.frame(r2_simonsi), row.names=TRUE, file="Results/simonsi/FullModel/refit_BestModel_PartialR2_REML_simonsi.csv")



#### 10. PLOT RADAR CHART FOR IMPORTANCE VARIABLES ----
#A. Load if it is necessary. Need to be ML models, for AIC calculation:
load(file="./Results/simonsi/FullModel/BM_simonsi.RData")

#A. Calculate variable importance using best models. Need to be more than 1 model:
impor_best = importance(BM[BM$delta <= 2, ])
impor_best
write.csv(impor_best, "./Results/simonsi/FullModel/ImportanceVariables_Best_simonsi.csv")

#B. Create a dataframe for ploting results:
impor_best
impor_best2 = impor_best
df = as.data.frame(impor_best2)
df[nrow(df) + 1,] = 0
df[nrow(df) + 1,] = 0
row.names(df)[4:5] = c("topography", "eucl_dist")
df = as.data.frame(t(df))
group = c("Best Models")
df = cbind(group,df)
df
#C. Reorder like steerei:
# order = c( "Habitat", "Productivity", "River", "Euclidean", "Topography")
#df = df[,c(1,2,5,6,3,4,7,8)]

#C. Plot radar:
plot_radar = 
  ggradar(df,
          base.size = 2,
          values.radar = c("0%", "50%", "100%"),
          axis.labels = c( "Habitat", "Productivity", "River",  "Euclidean", "Topography"),
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
pdf("./Results/simonsi/Figures/RelativePredictorWeight_MLPE_AllVars_simonsi.pdf")
plot_radar
dev.off()




#A. Calculate variable importance using all uncorrelated models. Need to be more than 1 model:
impor_best = importance(BM)
impor_best
write.csv(impor_best, "./Results/simonsi/FullModel/ImportanceVariables_Uncorrelated_simonsi.csv")

#B. Create a dataframe for ploting results:
impor_best
impor_best2 = impor_best
df = as.data.frame(impor_best2)
#df[nrow(df) + 1,] = 0
#df[nrow(df) + 1,] = 0
#row.names(df)[4:5] = c("topography", "eucl_dist")
df = as.data.frame(t(df))
group = c("Uncorrelated Models")
df = cbind(group,df)
df
#C. Reorder like steerei:
#steerei order = c( "Habitat", "Productivity", "River",  "Euclidean", "Topography")
#df = df[,c(1,2,5,6,3,4,7,8)]

#C. Plot radar:
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

#D. Verify radar:
plot_radar

#E. Save
pdf("./Results/simonsi/Figures/RelativePredictorWeight_MLPE_AllVars_simonsi_uncorrelated.pdf")
plot_radar
dev.off()






#### 11. PLOT UNIVARIATES PREDICTORS USING BEST MODEL ----
#A. LOAD UNSCALED DATA
DT_unscaled =  read.csv("Metafiles/MLPE_table_simonsi.csv", row.names = 1)
names(DT_unscaled)

#B. Verify variables in Best Model
TopM2
# wetlands_L

#C. Build model with Habitat ----
M1 = nlme::lme(RELgen ~ wetlands_L,
               random = ~1|POP,
               correlation = corMLPE(form = ~ from_ID + to_ID),
               data = DT_unscaled, method = "REML")

# Decorrelated residuals
dec_resids = decorltd_res_lme(M1)
#dec_resids$dist_interval = DT_unscaled$dist_interval
dec_resids$geoDist = DT_unscaled$eucl_dist

# Change covar in dataframe to plot
cov = DT_unscaled$wetlands_L
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
  xlab("Habitat Resistance") + 
  theme_bw() + 
  theme(axis.title.y = element_text(size=15, color = "black", face = "bold"),
        axis.title.x = element_text(size=15, color = "black", face = "bold")) +
  theme(axis.text = element_text(face = "bold", color = "black", size = 15))

#Verify
plotA

pdf("./Results/simonsi/Figures/REL_deco_Habitat_simonsi.pdf")
plotA
dev.off()




#D. Build model with Productivity ----
M2 = nlme::lme(RELgen ~ productivity_dist,
               random = ~ 1|POP,
               correlation = corMLPE(form = ~ from_ID + to_ID),
               data = DT_unscaled, method = "REML")

# Decorrelate residulas
dec_resids = decorltd_res_lme(M2)
#dec_resids$dist_interval = DT_unscaled$dist_interval
dec_resids$geoDist = DT_unscaled$eucl_dist

# Change covar in dataframe to plot
cov = DT_unscaled$productivity_dist
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
  xlab("Productivity Resistance") + 
  theme_bw() + 
  theme(axis.title.y = element_text(size=15, color = "black", face = "bold"),
        axis.title.x = element_text(size=15, color = "black", face = "bold")) +
  theme(axis.text = element_text(face = "bold", color = "black", size = 15))

#Verify
plotB

pdf("./Results/simonsi/Figures/REL_deco_Productivity_simonsi.pdf")
plotB
dev.off()




#E. Build model with River Distance ----
M2 = nlme::lme(RELgen ~ riverdistance,
               random = ~ 1|POP,
               correlation = corMLPE(form = ~ from_ID + to_ID),
               data = DT_unscaled, method = "REML")

# Decorrelate residulas
dec_resids = decorltd_res_lme(M2)
#dec_resids$dist_interval = DT_unscaled$dist_interval
dec_resids$geoDist = DT_unscaled$eucl_dist

# Change covar in dataframe to plot
cov = DT_unscaled$riverdistance
df2plot = as_tibble(data.frame(dec_resids, covar = cov))
names(df2plot)


# Plot relationship
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

pdf("./Results/simonsi/Figures/REL_deco_River_simonsi.pdf")
plotC
dev.off()


#F. Build model with Eucliadean Distance ----
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

pdf("./Results/simonsi/Figures/REL_deco_EuclideanDis_simonsi.pdf")
plotD
dev.off()

#END
