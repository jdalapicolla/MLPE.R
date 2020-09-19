#####################  LANDSCAPE GENOMICS TUTORIAL   ##########################
##################   MLPE - STEP 01: DISTANCE MATRICES   ######################

### Script prepared by Jeronymo Dalapicolla, Jamille C. Veiga, Carolina S. Carvalho, Luciana C. Resende-Moreira, and Rodolfo Jaffé ###


##### PRE-ANALYSIS -------------------------------------------------------------

#AIM = SELECTING BEST MODELS AND SHOW IMPORTANCE OF VARIABLES

##Load packages
library(corMLPE)
library(nlme)
library(MuMIn)
library(fmsb)
library(GeNetIt)
library(tidyverse)
library(ggradar)
library(r2glmm)

#load functions:
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

pr.norm.singlecov = function(object){
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




##################################################### STEEREI - SEASONAL FLOODPLAIN FOREST

##### 1. LOAD FILES ---------------------------------------------
#A. Load MLPE table:
mlpe_steerei = read.csv("Metafiles/MLPE_table_steerei.csv", row.names = 1)
head(mlpe_steerei)

#B. Maximum of variable in MLPE. 15 observations to each predictor:
max_var = length(mlpe_steerei[,1])/15
max_var #11

#C. Creating formula for MLPE:
#all variables:
vars_mlpe_st = names(mlpe_steerei)[5:(length(mlpe_steerei)-1)]
vars_mlpe_st

#B Use REL distance because in PCA river distance was not linear:
vars_mlpe_st_rel = vars_mlpe_st[c(1:6,10)] #Only Wetlands_VH. See Univariate results in step #2. Keep Null or Euclidean distance for Null hyphotesis.
mlpe_formula_rel_st = as.formula(paste("RELgen ~ ",paste(vars_mlpe_st_rel, collapse = " + "),sep = ""))
mlpe_formula_rel_st







##### 2. RUNNING FULL MODEL ---------------------------------------------------
#A. Full model for REL:
Full_model_st_rel = lme(mlpe_formula_rel_st,
                        random = ~1|POP,
                        correlation = corMLPE(form = ~ from_ID + to_ID|POP),
                        data = mlpe_steerei,
                        method = "ML")

#B. Run dredge function:
#specifying the number of predictor variables (max_var) and including the max.r function
Allmodels = dredge(Full_model_st_rel, rank = "AIC", m.lim=c(0, 11), extra= c(max.r))
length(Allmodels[,1]) # 128models

#C. Retrieve non-collinear models (max.r <=0.7) - Uncorrelated models
NCM = get.models(Allmodels, subset = max.r<=0.7)
##Final model selection table
BM = model.sel(NCM)
length(BM[,1]) #51 models
write.csv(BM, row.names=TRUE, file="steerei/MLPE_REL/Allmodels_woCor_REL_sorted_steerei.csv")

#D. Selecting Best models - AIC <=2. - Best Models
best_models = BM[BM$delta<=2, ]
length(best_models[,1]) #2 models
write.csv(best_models, "steerei/MLPE_REL/BestModels_REL_steerei.csv")

#E. Calculate variable importance in all uncorrelated and best models. Need to be more than 1 model:
impor_uncor = importance(BM)
impor_uncor
write.csv(impor_uncor, "steerei/MLPE_REL/ImportanceVariables_Uncor_REL_steerei.csv")

impor_best = importance(best_models)
impor_best
write.csv(impor_best, "steerei/MLPE_REL/ImportanceVariables_Best_REL_steerei.csv")


#D. Model-averaged coeficients
MA = model.avg(BM, subset = delta<=2)
confint(MA) ## Confidence Intervals for model-averaged coeficients
as.data.frame(MA$coefficients) ## Mean 

write.csv(as.data.frame(confint(MA)), row.names=TRUE, file="steerei/MLPE_REL/BestModels_CI_coefficients_REL_steerei.csv")
write.csv(as.data.frame(MA$coefficients), row.names=TRUE, file="steerei/MLPE_REL/BestModels_Mean_coefficients_REL_steerei.csv")






##### 3. SIGNIFICANCE OF BEST MODELS ---------------------------------------------------
#A. Best model for REL:

best_models[1,]
names(mlpe_steerei)

#B. Run Best model 
bestmodel = lme(RELgen ~ PET +precipitation +riverdistance,
                       random = ~1|POP,
                       correlation = corMLPE(form = ~ from_ID + to_ID|POP),
                       data = mlpe_steerei,
                       method = "ML")


#C. 2nd Best Model - 2 models:
best_models[2,]
names(mlpe_steerei)

bestmodel_2nd =lme(RELgen ~ PET +precipitation +riverdistance +wetlands_VH,
                   random = ~1|POP,
                   correlation = corMLPE(form = ~ from_ID + to_ID|POP),
                   data = mlpe_steerei,
                   method = "ML")

#D. LTR for Nested Models:
compare_bestmodels = anova.lme(bestmodel, bestmodel_2nd)
compare_bestmodels
#Save LTR results:
write.csv(as.data.frame(compare_bestmodels), file="LTR_BestModels_steerei.csv", row.names = TRUE)








##### 4. COEFFICIENTS OF BEST MODELS ---------------------------------------------------
#A. Best model with REML:
best_models[1,]
names(mlpe_steerei)

#B. Run Best model 
Ajusted_bestmodel = lme(RELgen ~ PET +precipitation +riverdistance,
                random = ~1|POP,
                correlation = corMLPE(form = ~ from_ID + to_ID|POP),
                data = mlpe_steerei,
                method = "REML")


#C. 2nd Best Model:
best_models[2,]
names(mlpe_steerei)

Ajusted_bestmodel_2nd = lme(RELgen ~ PET +precipitation +riverdistance +wetlands_VH,
                   random = ~1|POP,
                   correlation = corMLPE(form = ~ from_ID + to_ID|POP),
                   data = mlpe_steerei,
                   method = "REML")

#D. Testing Spatial Dependence pattern
acf(resid(Ajusted_bestmodel,type='normalized')) ## Spatial dependence pattern
#Save graphs
pdf("ACF_regular_steerei.pdf", onefile = T)
acf(resid(Ajusted_bestmodel,type='normalized'))
dev.off()


#Summary for Best Models:
summ_st = as.data.frame(capture.output(summary(Ajusted_bestmodel)))
write.csv(summ_st, file="Summary_AjustedBestModels_steerei.csv", row.names = TRUE)

summ_st2 = as.data.frame(capture.output(summary(Ajusted_bestmodel_2nd)))
write.csv(summ_st, file="Summary_AjustedBestModels_2nd_steerei.csv", row.names = TRUE)

#Estimates:
MA_REML = model.avg(Ajusted_bestmodel, Ajusted_bestmodel_2nd)
confint(MA_REML) ## Confidence Intervals for model-averaged coeficients
as.data.frame(MA_REML$coefficients) ## Mean 

write.csv(as.data.frame(confint(MA_REML)), row.names=TRUE, file="BestModels_CI_REML_coefficients_REL_steerei.csv")
write.csv(as.data.frame(MA_REML$coefficients), row.names=TRUE, file="BestModels_Mean_REML_coefficients_REL_steerei.csv")


#Partial R2:
r2_steerei = r2beta(Ajusted_bestmodel, partial = 'TRUE')
r2_steerei

write.csv(as.data.frame(r2_steerei), row.names=TRUE, file="BestModel_PartialR2_REML_steerei.csv")



######### 5. ESTIMATES GRAPHS -----------------------------------------------------
#A, Calculated decorrelated Relatedness:
Ajusted_bestmodel
goo1 = pr.norm.singlecov(Ajusted_bestmodel)

#A. PET:
df = data.frame(goo1, covar=mlpe_steerei$PET)

PET_graph =
ggplot(df, aes(x = covar)) +
  geom_point(aes(y=pr.fixed), alpha=0.50, color= 'black') +
  geom_point(aes(y=pr.POP), alpha=0.50, color='red') +
  geom_smooth(aes(y=fit), method='lm', color="black") +
  ylab("Relatedness (Decorrelated)") +
  xlab("PET Dissimilarity") +
  theme_bw() +
  theme(axis.title.y = element_text(size=15, color = "black", face = "bold"),
        axis.title.x = element_text(size=15, color = "black", face = "bold"))

## Save
pdf("./PET_Estimates_steerei.pdf")
PET_graph
dev.off()


#A. PET:
df = data.frame(goo1, covar=mlpe_steerei$PET)

PET_graph =
  ggplot(df, aes(x = covar)) +
  geom_point(aes(y=pr.fixed), alpha=0.50, color= 'black') +
 # geom_point(aes(y=pr.POP), alpha=0.50, color='red') +
  geom_smooth(aes(y=fit), method='lm', color="black") +
  ylab("Relatedness (Decorrelated)") +
  xlab("PET Dissimilarity") +
  theme_bw() +
  theme(axis.title.y = element_text(size=15, color = "black", face = "bold"),
        axis.title.x = element_text(size=15, color = "black", face = "bold"))

PET_graph

## Save
pdf("./PET_Estimates_steerei.pdf")
PET_graph
dev.off()



#B. Precipitation:
df = data.frame(goo1, covar=mlpe_steerei$precipitation)

Preci_graph =
  ggplot(df, aes(x = covar)) +
  geom_point(aes(y=pr.fixed), alpha=0.50, color= 'black') +
  # geom_point(aes(y=pr.POP), alpha=0.50, color='red') +
  geom_smooth(aes(y=fit), method='lm', color="black") +
  ylab("Relatedness (Decorrelated)") +
  xlab("Precipitation Dissimilarity") +
  theme_bw() +
  theme(axis.title.y = element_text(size=15, color = "black", face = "bold"),
        axis.title.x = element_text(size=15, color = "black", face = "bold"))

Preci_graph

## Save
pdf("./Preci_Estimates_steerei.pdf")
Preci_graph
dev.off()


#C. River Distance:
df = data.frame(goo1, covar=mlpe_steerei$riverdistance)

River_graph =
  ggplot(df, aes(x = covar)) +
  geom_point(aes(y=pr.fixed), alpha=0.50, color= 'black') +
  # geom_point(aes(y=pr.POP), alpha=0.50, color='red') +
  geom_smooth(aes(y=fit), method='lm', color="black") +
  ylab("Relatedness (Decorrelated)") +
  xlab("River Distance") +
  theme_bw() +
  theme(axis.title.y = element_text(size=15, color = "black", face = "bold"),
        axis.title.x = element_text(size=15, color = "black", face = "bold"))

River_graph

## Save
pdf("./River_Estimates_steerei.pdf")
River_graph
dev.off()



#D. Euclidean Distance:
df = data.frame(goo1, covar=mlpe_steerei$eucl_dist)

Eucli_graph =
  ggplot(df, aes(x = covar)) +
  geom_point(aes(y=pr.fixed), alpha=0.50, color= 'black') +
  # geom_point(aes(y=pr.POP), alpha=0.50, color='red') +
  geom_smooth(aes(y=fit), method='lm', color="black") +
  ylab("Relatedness (Decorrelated)") +
  xlab("Euclidean (Geographic) Distance") +
  theme_bw() +
  theme(axis.title.y = element_text(size=15, color = "black", face = "bold"),
        axis.title.x = element_text(size=15, color = "black", face = "bold"))

Eucli_graph

## Save
pdf("./Eucli_Estimates_steerei.pdf")
Eucli_graph
dev.off()




######### 6. RADAR GRAPHS -----------------------------------------------------
#df for ploting results for importance variables:
impor_uncor
impor_best

impor_best2 = impor_best
x = as.data.frame(impor_best2)
x[nrow(x) + 1,] = 0
x[nrow(x) + 1,] = 0
x[nrow(x) + 1,] = 0
row.names(x)[5:7] = c("temperature", "topography", "eucl_dist")
x = x[c(3,1,2,4,6,7,5),]

df = rbind(impor_uncor, t(x))
df = as.data.frame(df)
group = c("Uncorrelated Models", "Best Models")
df = cbind(group,df)
df

plot_st = 
ggradar(df,
        base.size = 2,
        values.radar = c("0%", "50%", "100%"),
        axis.labels = c( "River", "PET", "Precipitation", "Wetlands", "Topography", "Euclidean", "Temperature"),
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





########################################################### SIMONSI - NON-FLOOED-FOREST

##### 1b. LOAD FILES ---------------------------------------------
#A. Load MLPE table:
mlpe_simonsi = read.csv("Metafiles/MLPE_table_simonsi.csv", row.names = 1)
head(mlpe_simonsi)

#B. Maximum of variable in MLPE. 15   observations to each predictor:
max_var = length(mlpe_simonsi[,1])/15
max_var #9


#C. Creating formula for MLPE:
#all variables:
vars_mlpe_si = names(mlpe_simonsi)[5:(length(mlpe_simonsi)-1)]
vars_mlpe_si


#D. Using REL distance:
vars_mlpe_si_rel = vars_mlpe_si[c(1:6,10)] #Only Wetlands_VH. See Univariate results in step #2:
mlpe_formula_rel_si = as.formula(paste("RELgen ~ ",paste(vars_mlpe_si_rel, collapse = " + "),sep = ""))
mlpe_formula_rel_si





##### 2b. RUNNING FULL MODEL ---------------------------------------------------
#A. Full model for REL:
Full_model_si_rel = lme(mlpe_formula_rel_si,
                        random = ~1|POP,
                        correlation = corMLPE(form = ~ from_ID + to_ID|POP),
                        data = mlpe_simonsi,
                        method = "ML")

#B. Run dredge function:
#specifying the number of predictor variables (max_var) and including the max.r function
Allmodels = dredge(Full_model_si_rel, rank = "AICc", m.lim=c(0, 9), extra= c(max.r))
length(Allmodels[,1]) # 128models

#C. Retrieve non-collinear models (max.r <=0.7)
NCM = get.models(Allmodels, subset = max.r<=0.7)
##Final model selection table
BM = model.sel(NCM)
length(BM[,1]) #64 models
write.csv(BM, row.names=TRUE, file="simonsi/MLPE_REL/Allmodels_woCor_REL_sorted_simonsi.csv")

#D. Selecting Best models - AIC <=2.
best_models = BM[BM$delta<=2, ]
length(best_models[,1]) #2 models
write.csv(best_models, "simonsi/MLPE_REL/BestModels_REL_simonsi.csv")

#E. Calculate variable importance in all uncorrelated and best models:
impor_uncor = importance(BM)
write.csv(impor_uncor, "simonsi/MLPE_REL/ImportanceVariables_Uncor_REL_simonsi.csv")

impor_best = importance(BM[BM$delta<=2, ])
write.csv(impor_best, "simonsi/MLPE_REL/ImportanceVariables_Best_REL_simonsi.csv")

#F. Model-averaged coeficients
MA = model.avg(BM, subset = delta<=2)
confint(MA) ## Confidence Intervals for model-averaged coeficients
as.data.frame(MA$coefficients) ## Mean 

write.csv(as.data.frame(confint(MA)), row.names=TRUE, file="simonsi/MLPE_REL/BestModels_CI_coefficients_REL_simonsi.csv")
write.csv(as.data.frame(MA$coefficients), row.names=TRUE, file="simonsi/MLPE_REL/BestModels_Mean_coefficients_REL_simonsi.csv")






##### 3b. SIGNIFICANCE OF BEST MODELS ---------------------------------------------------
#A. Best model for REL:
best_models[1,]
names(mlpe_simonsi)

#B. Run Best model 
bestmodel = lme(RELgen ~ PET +wetlands_VH,
                random = ~1|POP,
                correlation = corMLPE(form = ~ from_ID + to_ID|POP),
                data = mlpe_simonsi,
                method = "ML")


#C. 2nd Best Model - 2 models:
best_models[2,]
names(mlpe_simonsi)

bestmodel_2nd =lme(RELgen ~ PET +wetlands_VH +riverdistance ,
                   random = ~1|POP,
                   correlation = corMLPE(form = ~ from_ID + to_ID|POP),
                   data = mlpe_simonsi,
                   method = "ML")

#D. LTR for Nested Models:
compare_bestmodels = anova.lme(bestmodel, bestmodel_2nd)
compare_bestmodels
#Save LTR results:
write.csv(as.data.frame(compare_bestmodels), file="LTR_BestModels_simonsi.csv", row.names = TRUE)








##### 4b. COEFFICIENTS OF BEST MODELS ---------------------------------------------------
#A. Best model with REML:
best_models[1,]
names(mlpe_simonsi)

#B. Run Best model 
Ajusted_bestmodel =lme(RELgen ~ PET +wetlands_VH,
                      random = ~1|POP,
                      correlation = corMLPE(form = ~ from_ID + to_ID|POP),
                      data = mlpe_simonsi,
                      method = "REML") #for better estimation of coefficients

#B. Run 2nd Best model 
best_models[2]
Ajusted_bestmodel_2nd = lme(RELgen ~ wetlands_VH +PET +riverdistance,
                    random = ~1|POP,
                    correlation = corMLPE(form = ~ from_ID + to_ID|POP),
                    data = mlpe_simonsi,
                    method = "REML")


#D. Testing Spatial Dependence pattern
acf(resid(Ajusted_bestmodel,type='normalized')) ## Spatial dependence pattern
#Save graphs
pdf("ACF_regular_simonsi.pdf", onefile = T)
acf(resid(Ajusted_bestmodel,type='normalized'))
dev.off()


#Summary for Best Models:
summ_st = as.data.frame(capture.output(summary(Ajusted_bestmodel)))
write.csv(summ_st, file="Summary_AjustedBestModels_simonsi.csv", row.names = TRUE)

summ_st2 = as.data.frame(capture.output(summary(Ajusted_bestmodel_2nd)))
write.csv(summ_st, file="Summary_AjustedBestModels_2nd_simonsi.csv", row.names = TRUE)

#Save Results:
MA_REML = model.avg(Ajusted_bestmodel, Ajusted_bestmodel_2nd)
confint(MA_REML) ## Confidence Intervals for model-averaged coeficients
as.data.frame(MA_REML$coefficients) ## Mean 

write.csv(as.data.frame(confint(MA_REML)), row.names=TRUE, file="BestModels_CI_REML_coefficients_REL_simonsi.csv")
write.csv(as.data.frame(MA_REML$coefficients), row.names=TRUE, file="BestModels_Mean_REML_coefficients_REL_simonsi.csv")



#Partial R2:
r2_simonsi = r2beta(Ajusted_bestmodel, partial = 'TRUE')
r2_simonsi

write.csv(as.data.frame(r2_simonsi), row.names=TRUE, file="BestModel_PartialR2_REML_simons.csv")









######### 5b. ESTIMATES GRAPHS -----------------------------------------------------
#A, Calculated decorrelated Relatedness:
Ajusted_bestmodel
goo1 = pr.norm.singlecov(Ajusted_bestmodel)

#A. PET:
df = data.frame(goo1, covar=mlpe_simonsi$PET)

PET_graph =
  ggplot(df, aes(x = covar)) +
  geom_point(aes(y=pr.fixed), alpha=0.50, color= 'black') +
  #geom_point(aes(y=pr.POP), alpha=0.50, color='red') +
  geom_smooth(aes(y=fit), method='lm', color="black") +
  ylab("Relatedness (Decorrelated)") +
  xlab("PET Dissimilarity") +
  theme_bw() +
  theme(axis.title.y = element_text(size=15, color = "black", face = "bold"),
        axis.title.x = element_text(size=15, color = "black", face = "bold"))

PET_graph

## Save
pdf("./PET_Estimates_simonsi.pdf")
PET_graph
dev.off()


#B. Habitat:
df = data.frame(goo1, covar=mlpe_simonsi$wetlands_VH)

HAB_graph =
  ggplot(df, aes(x = covar)) +
  geom_point(aes(y=pr.fixed), alpha=0.50, color= 'black') +
  # geom_point(aes(y=pr.POP), alpha=0.50, color='red') +
  geom_smooth(aes(y=fit), method='lm', color="black") +
  ylab("Relatedness (Decorrelated)") +
  xlab("Habitat Resistance") +
  theme_bw() +
  theme(axis.title.y = element_text(size=15, color = "black", face = "bold"),
        axis.title.x = element_text(size=15, color = "black", face = "bold"))

HAB_graph

## Save
pdf("./HAB_Estimates_simonsi.pdf")
HAB_graph
dev.off()

#C. Euclidean Distance:
df = data.frame(goo1, covar=mlpe_simonsi$eucl_dist)

Eucli_graph =
  ggplot(df, aes(x = covar)) +
  geom_point(aes(y=pr.fixed), alpha=0.50, color= 'black') +
  # geom_point(aes(y=pr.POP), alpha=0.50, color='red') +
  geom_smooth(aes(y=fit), method='lm', color="black") +
  ylab("Relatedness (Decorrelated)") +
  xlab("Euclidean (Geographic) Distance") +
  theme_bw() +
  theme(axis.title.y = element_text(size=15, color = "black", face = "bold"),
        axis.title.x = element_text(size=15, color = "black", face = "bold"))

Eucli_graph

## Save
pdf("./Eucli_Estimates_simonsi.pdf")
Eucli_graph
dev.off()




######### 6b. RADAR GRAPHS -----------------------------------------------------
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
df = df[c(1,4,2,5,3,7,8,6)]

plot_si = 
  ggradar(df,
          base.size = 2,
          values.radar = c("0%", "50%", "100%"),
          axis.labels = c( "River", "PET", "Precipitation", "Wetlands", "Topography", "Euclidean", "Temperature"),
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
