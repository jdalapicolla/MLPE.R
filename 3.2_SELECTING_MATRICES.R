#####################  LANDSCAPE GENOMICS TUTORIAL   ##########################
##################   MLPE - STEP 01: DISTANCE MATRICES   ######################

### Script prepared by Jeronymo Dalapicolla, Jamille C. Veiga, Carolina S. Carvalho, Luciana C. Resende-Moreira, and Rodolfo Jaffé ###


##### PRE-ANALYSIS -------------------------------------------------------------

#AIM = ORGANAZING AND SELECTING VARIABLES FOR MODELS:

#Load packages:
library(tidyverse)
library(ggplot2)
library(plyr)
library(GeNetIt)
library(reshape2)
library(corrplot)
library(r2vcftools)
library(MuMIn)
library(corMLPE) #devtools::install_github("nspope/corMLPE")
library(nlme)
library(snow)
library(parallel)


#Turn off scientific notation
options(scipen=999)

#Load auxiliary functions:
VCFsummary <- function(snps){
  Nid <- nrow(snps@meta)
  Nsnps <- length(snps@site_id)
  cat(paste(Nid, "individuals and", Nsnps, "SNPs."), sep="\n")
}


##### 1. REARRANGING DISTANCE MATRICES FILES ---------------------------------------------
#Calculated in Pipeline for Genetic Structure #1.4 https://github.com/jdalapicolla/LanGen_pipeline_version2:


########################################################### STEEREI - NON-FLOOED-FOREST:

#A. Load landscape distance matrices and reorder:
euclidean = read.csv("Distances/eucl_dist_steerei.csv", row.names = 1)
head(euclidean)
euclidean = euclidean[order(euclidean$X1, euclidean$X2), ]
head(euclidean)

productivity = read.csv("Distances/productivity_dist_steerei.csv", row.names = 1)
head(productivity)
productivity = productivity[order(productivity$X1, productivity$X2), ]
head(productivity)

river = read.csv("Distances/river_dist_steerei.csv", row.names = 1)
head(river)
river = river[order(river$X1, river$X2), ]
head(river)

topo = read.csv("Distances/topo_dist_steerei.csv", row.names = 1)
head(topo)
topo = topo[order(topo$X1, topo$X2), ]
head(topo)

wet_L = read.csv("Distances/wetlands_dist_L_steerei.csv", row.names = 1)
head(wet_L)
wet_L = wet_L[order(wet_L$X1, wet_L$X2), ]
head(wet_L)

wet_M = read.csv("Distances/wetlands_dist_M_steerei.csv", row.names = 1)
head(wet_M)
wet_M = wet_M[order(wet_M$X1, wet_M$X2), ]
head(wet_M)

wet_H = read.csv("Distances/wetlands_dist_H_steerei.csv", row.names = 1)
head(wet_H)
wet_H = wet_H[order(wet_H$X1, wet_H$X2), ]
head(wet_H)

wet_VH = read.csv("Distances/wetlands_dist_VH_steerei.csv", row.names = 1)
head(wet_VH)
wet_VH = wet_VH[order(wet_VH$X1, wet_VH$X2), ]
head(wet_VH)


#B. Load genetic distance matrices - Relatedness.
REL = read.csv("Distances/genetic_dist_Yang_Rel_steerei.csv", row.names = 1)
head(REL)
length(REL[,3]) #190

#C. Remove comparisons between same individuals, change the colnames, and invert columns of individuals because Relatedness and PCA function use start their estimation by cols and not by rows
REL = REL[REL$INDV1 != REL$INDV2,]
length(REL[,3]) #171
ind1 = REL$INDV2
ind2 = REL$INDV1
colnames(REL) = c("X1", "X2","RELgen")
REL$X1 = ind1
REL$X2 = ind2
head(REL)

#D.Reordering according to the first two cols from other distance matrices
REL = REL[order(REL$X1, REL$X2),]
head(REL)

#E. Load genetic distance matrices - PCA-distance.
pca = as.matrix(read.csv("Distances/genetic_dist_PCA_BSR_steerei.csv", row.names = 1))
pca[upper.tri(pca, diag = T)] = NA
head(pca)

#F. Convert in data frame and reorder columns
pca = pca %>%
  melt %>%
  na.omit %>%
  arrange(., Var1) %>%
  setNames(c("X1", "X2","PCAgen"))
head(pca)
length(pca[,3])

pca = pca[order(as.character(pca$X1), as.character(pca$X2)), ]
head(pca)

#G. Verify individuals order
#of colunm #1: 2 by turn
identical(as.character(REL$X1), as.character(pca$X1))
identical(as.character(euclidean$X1), as.character(river$X1))
identical(as.character(river$X1), as.character(productivity$X1))
identical(as.character(productivity$X1), as.character(topo$X1))
identical(as.character(wet_L$X1), as.character(topo$X1))
identical(as.character(wet_L$X1), as.character(wet_M$X1))
identical(as.character(wet_H$X1), as.character(wet_M$X1))
identical(as.character(wet_H$X1), as.character(wet_VH$X1))

#colunm #2:
identical(as.character(REL$X2), as.character(pca$X2))
identical(as.character(euclidean$X2), as.character(river$X2))
identical(as.character(river$X2), as.character(productivity$X2))
identical(as.character(productivity$X2), as.character(topo$X2))
identical(as.character(wet_L$X2), as.character(topo$X2))
identical(as.character(wet_L$X2), as.character(wet_M$X2))
identical(as.character(wet_H$X2), as.character(wet_M$X2))
identical(as.character(wet_H$X2), as.character(wet_VH$X2))



##### 2. MERGING ALL DATA IN A SINGLE FILE ---------------------------------------------
#A. United in a single data frame
mlpe_table = cbind(REL, pca[3], euclidean[3], productivity[3], river[3], topo[3], wet_L[3], wet_M[3], wet_H[3], wet_VH[3])
#verify
head(mlpe_table)
colnames(mlpe_table)[1:2] = c("from_ID", "to_ID")
head(mlpe_table)


#B. Creating a structure population column in the table for classifying the distances in between individuals from same or different populations:
#Load neutral .vcf file with population geographical information and genetic clusters ID:
snps_neutral = vcfLink("vcf/steerei_filtered_neutral_LEA_DAPC_TESS.vcf", overwriteID=T)
VCFsummary(snps_neutral) #19 individuals and 13971 SNPs.
names(snps_neutral@meta) #verify col names in metafile

#extract popualtion and sample ID info
pop1 = snps_neutral@meta$sample_name[snps_neutral@meta$PopID_snmf==1]
pop2 = snps_neutral@meta$sample_name[snps_neutral@meta$PopID_snmf==2]
pop3 = snps_neutral@meta$sample_name[snps_neutral@meta$PopID_snmf==3]

#create a df to represent distances
teste = as.data.frame(mlpe_table[,1:2])
head(teste)
teste[] = lapply(teste, as.character)
head(teste)

#replace individual names by numbers
teste$from_ID[teste$from_ID %in% pop1] = 1
teste$to_ID[teste$to_ID %in% pop1] = 1

teste$from_ID[teste$from_ID %in% pop2] = 2
teste$to_ID[teste$to_ID %in% pop2] = 2

teste$from_ID[teste$from_ID %in% pop3] = 3
teste$to_ID[teste$to_ID %in% pop3] = 3
#verify
head(teste)

#classify distances between same or different populations
#same genetic cluster
teste$POP_dist[teste$from_ID == teste$to_ID] = 1
#genetic cluster 1 and 2
teste$POP_dist[teste$from_ID == 1 & teste$to_ID == 2 | teste$from_ID == 2 & teste$to_ID == 1] = 2
#genetic cluster 1 and 3
teste$POP_dist[teste$from_ID == 1 & teste$to_ID == 3 | teste$from_ID == 3 & teste$to_ID == 1] = 3
#genetic cluster 2 and 3
teste$POP_dist[teste$from_ID == 2 & teste$to_ID == 3 | teste$from_ID == 3 & teste$to_ID == 2] = 4
head(teste)
teste

#add to mlpe table:
mlpe_table$POP = teste$POP_dist
head(mlpe_table)

#save df 
write.csv(mlpe_table, "Metafiles/MLPE_table_steerei.csv")




##### 3. COVERAGE DEPTH BY SAMPLE ---------------------------------------------
#A. Estimating coverage depth by sample: 
coverage_ind = c()
for (p in 1:length(snps_neutral@sample_id)){
  beta = Query(Subset(snps_neutral, samples = p), type="site-mean-depth")
  coverage_ind[p] = mean(beta$MEAN_DEPTH, na.rm = T)}
#verify
coverage_ind
#save as metafile
snps@meta$coverage = coverage_ind

#B. save coverage depth by sample:
write.csv(as.data.frame(cbind(snps_neutral@meta$ind_ID,coverage_ind)), "./Metafiles/Depth_coverage_bysamples_steerei.csv")





##### 4. VARIABLE RANGES -------------------------------------------------------------
#Variable ranges for all pairwise resistance distances used to run MLPE regression models.

#Load MLPE table if necessary:
mlpe_table = read.csv("Metafiles/MLPE_table_steerei.csv", row.names = 1)
head(mlpe_table)

#Preparing input for ggplot:
df = mlpe_table[5:(length(mlpe_table)-1)] #explanatory variables start in col #5
df = melt(df)
head(df)

#D. Plot the graph
pdf("boxplot_variables_mlpe_steerei.pdf", onefile = F)
ggplot(data=df, mapping = aes(x= value, fill = variable)) +
  geom_boxplot() +
  theme_bw()
dev.off()







#####5. VARIABLE LINEARITY ---------------------------------------------------------
## Assess linearity between response and predictor variables. The graphs show if the slope (trend line) is continually changing; if so, it is not a constant! These variables do not show a linear relationship. With a linear relationship, the slope never changes.

#A. Load MLPE table if necessary:
mlpe_table = read.csv("Metafiles/MLPE_table_steerei.csv", row.names = 1)
head(mlpe_table)
mlpe_table[, 5:(length(mlpe_table)-1)] = as.data.frame(scale(mlpe_table[, 5:(length(mlpe_table)-1)]))
head(mlpe_table)

names = names(mlpe_table[, 5:(length(mlpe_table)-1)]) #explanatory variables start in col #5
names


#Using Relatedness:
for(i in 1:length(names)){
  pdf(paste0("Results/steerei/Linearity/Relatedness_ScatterSmooth_steerei_",names[i],".pdf"), onefile = T)
  scatter.smooth(x = mlpe_table[, names[i]], y=mlpe_table[, "RELgen"], ylab = "Relatedness", xlab=colnames(mlpe_table[names[i]]), lpars = list(col = "red", lwd = 3))
  dev.off()
}

#Using PCA distance:
for(i in 1:length(names)){
  pdf(paste0("Results/steerei/Linearity/PCA_ScatterSmooth_steerei_",names[i],".pdf"), onefile = T)
  scatter.smooth(x = mlpe_table[, names[i]], y=mlpe_table[, "PCAgen"], ylab = "PCA-Distance", xlab=colnames(mlpe_table[names[i]]), lpars = list(col = "red", lwd = 3))
  dev.off()
}

##NOT LINEAR:
#River Distance using PCA distance - Use Relatedness




#####6. UNIVARIATE MODELS ---------------------------------------------------------
#Univariate Models to choose one Habitat Distance matrix
#A. Load MLPE table if necessary:
mlpe_table = read.csv("Metafiles/MLPE_table_steerei.csv", row.names = 1)
head(mlpe_table)
mlpe_table[, 5:(length(mlpe_table)-1)] = as.data.frame(scale(mlpe_table[, 5:(length(mlpe_table)-1)]))
head(mlpe_table)
 

#A. Using Relatedness and controlling by population structure
## BUILD FULL MODEL
form = as.formula(paste("RELgen ~ ", 
                         paste(names(mlpe_table)[9:12], collapse = " + "), sep = ""))

Fullmodel <- lme(form,
                 random = ~ 1|POP,
                 correlation = corMLPE(form = ~ from_ID + to_ID),
                 data = mlpe_table, method = "ML")


## CHECK RESIDUALS
RES <- residuals(Fullmodel, type="normalized")
FIT <- fitted(Fullmodel)
plot(FIT, RES) ; abline(0,0, col="red")
acf(RES)

## ALL VARIABLE PLOTS
names(mlpe_table)[9:12]
plot(mlpe_table$wetlands_L, RES) ; abline(0,0, col="red")
plot(mlpe_table$wetlands_M, RES) ; abline(0,0, col="red")
plot(mlpe_table$wetlands_H, RES) ; abline(0,0, col="red")
plot(mlpe_table$wetlands_VH, RES) ; abline(0,0, col="red") 


## RUN PDREDGE 

## Set cluster
## Run paralleled dredge
## by Ben Bolker in: https://stackoverflow.com/questions/30261473/mumin-pdredge-error-with-glmer

## Make cluster
cluster = makeCluster(4, type = "SOCK")  ## also need snow installed

## clusterExport assigns the values on the master R process of the variables named in varlist to variables of the same names in the global environment (aka ‘workspace’) of each node. 
## The environment on the master from which variables are exported defaults to the global environment.
clusterExport(cluster,"mlpe_table")

## clusterEvalQ evaluates a literal expression on each cluster node. 
## It is a parallel version of evalq, and is a convenience function invoking clusterCall
clusterEvalQ(cluster,
             c(library(nlme), library(MuMIn), library(corMLPE)))

## Specify the number of predictor variables and including the max.r function
nrow(mlpe_table) ## 171
options(na.action = na.fail)

## Run pdredge
Allmodels <- MuMIn::pdredge(Fullmodel, rank = "AIC", 
                            m.lim=c(0, 1), cluster)

## SAVE ALL MODELS
save(Allmodels, file="./Results/steerei/Univariates/Allmodels_univariates_steerei.RData")

### LOAD ALL MODELS
load(file="./Results/steerei/Univariates/Allmodels_univariates_steerei.RData")

## BEST MODELS
nrow(Allmodels) ## 5
BM = model.sel(Allmodels, rank=AIC)
df = as.data.frame(BM[BM$delta <=2, ])
df 
write.csv(as.data.frame(BM), "Results/steerei/Univariates/Unimodels_habitats_steerei.csv")

#REL-distance: Low Resistance





##### 7. CORRELATIONS -------------------------------------------------------------
#perform Correlogram showing the correlation between all resistance distances included as predictors in MLPE regression models. Pearson’s correlation coefficients (r). 

#A. Load MLPE table if necessary:
mlpe_table = read.csv("Metafiles/MLPE_table_steerei.csv", row.names = 1)
head(mlpe_table)

## Check the correlation among environmental variables
ALLvars = as.data.frame(mlpe_table[5:(length(mlpe_table)-1)])
head(ALLvars)
Allvars_cor = cor(ALLvars)
head(Allvars_cor)

write.csv(Allvars_cor, "Results/steerei/Correlation/Correlation_steerei.csv")

## Corplot
pdf("Results/steerei/Correlation/Correlogram_steerei.pdf", onefile = T)
corrplot(Allvars_cor, 
         method="color", 
         type = "upper", 
         outline = FALSE, 
         diag = F,
         tl.col = "black")
dev.off()






##### 1b. REARRANGING DISTANCE MATRICES FILES ---------------------------------------------
#Calculated in Pipeline for Genetic Structure #1.4 https://github.com/jdalapicolla/LanGen_pipeline_version2:

## Clean Global Environment 
rm(list=ls())

#Turn off scientific notation
options(scipen=999)

#Load auxiliary functions:
VCFsummary <- function(snps){
  Nid <- nrow(snps@meta)
  Nsnps <- length(snps@site_id)
  cat(paste(Nid, "individuals and", Nsnps, "SNPs."), sep="\n")
}


################################################## SIMONSI - SEASONAL FLOODPLAIN FORESTS:

#A. Load landscape distance matrices and reorder:
euclidean = read.csv("Distances/eucl_dist_simonsi.csv", row.names = 1)
head(euclidean)
euclidean = euclidean[order(euclidean$X1, euclidean$X2), ]
head(euclidean)

productivity = read.csv("Distances/productivity_dist_simonsi.csv", row.names = 1)
head(productivity)
productivity = productivity[order(productivity$X1, productivity$X2), ]
head(productivity)

river = read.csv("Distances/river_dist_simonsi.csv", row.names = 1)
head(river)
river = river[order(river$X1, river$X2), ]
head(river)

topo = read.csv("Distances/topo_dist_simonsi.csv", row.names = 1)
head(topo)
topo = topo[order(topo$X1, topo$X2), ]
head(topo)

wet_L = read.csv("Distances/wetlands_dist_L_simonsi.csv", row.names = 1)
head(wet_L)
wet_L = wet_L[order(wet_L$X1, wet_L$X2), ]
head(wet_L)

wet_M = read.csv("Distances/wetlands_dist_M_simonsi.csv", row.names = 1)
head(wet_M)
wet_M = wet_M[order(wet_M$X1, wet_M$X2), ]
head(wet_M)

wet_H = read.csv("Distances/wetlands_dist_H_simonsi.csv", row.names = 1)
head(wet_H)
wet_H = wet_H[order(wet_H$X1, wet_H$X2), ]
head(wet_H)

wet_VH = read.csv("Distances/wetlands_dist_VH_simonsi.csv", row.names = 1)
head(wet_VH)
wet_VH = wet_VH[order(wet_VH$X1, wet_VH$X2), ]
head(wet_VH)


#B. Load genetic distance matrices - Relatedness.
REL = read.csv("Distances/genetic_dist_Yang_Rel_simonsi.csv", row.names = 1)
head(REL)
length(REL[,3]) #190

#C. Remove comparisons between same individuals, change the colnames, and invert columns of individuals because Relatedness and PCA function use start their estimation by cols and not by rows
REL = REL[REL$INDV1 != REL$INDV2,]
length(REL[,3]) #171
ind1 = REL$INDV2
ind2 = REL$INDV1
colnames(REL) = c("X1", "X2","RELgen")
REL$X1 = ind1
REL$X2 = ind2
head(REL)

#D.Reordering according to the first two cols from other distance matrices
REL = REL[order(REL$X1, REL$X2),]
head(REL)

#E. Load genetic distance matrices - PCA-distance.
pca = as.matrix(read.csv("Distances/genetic_dist_PCA_BSR_simonsi.csv", row.names = 1))
pca[upper.tri(pca, diag = T)] = NA
head(pca)

#F. Convert in data frame and reorder columns
pca = pca %>%
  melt %>%
  na.omit %>%
  arrange(., Var1) %>%
  setNames(c("X1", "X2","PCAgen"))
head(pca)
length(pca[,3])
pca = pca[order(as.character(pca$X1), as.character(pca$X2)), ]
head(pca)

#G. Verify individuals order
#of colunm #1: 2 by turn
identical(as.character(REL$X1), as.character(pca$X1))
identical(as.character(euclidean$X1), as.character(river$X1))
identical(as.character(river$X1), as.character(productivity$X1))
identical(as.character(productivity$X1), as.character(topo$X1))
identical(as.character(wet_L$X1), as.character(topo$X1))
identical(as.character(wet_L$X1), as.character(wet_M$X1))
identical(as.character(wet_H$X1), as.character(wet_M$X1))
identical(as.character(wet_H$X1), as.character(wet_VH$X1))

#colunm #2:
identical(as.character(REL$X2), as.character(pca$X2))
identical(as.character(euclidean$X2), as.character(river$X2))
identical(as.character(river$X2), as.character(productivity$X2))
identical(as.character(productivity$X2), as.character(topo$X2))
identical(as.character(wet_L$X2), as.character(topo$X2))
identical(as.character(wet_L$X2), as.character(wet_M$X2))
identical(as.character(wet_H$X2), as.character(wet_M$X2))
identical(as.character(wet_H$X2), as.character(wet_VH$X2))



##### 2b. MERGING ALL DATA IN A SINGLE FILE ---------------------------------------------
#A. United in a single data frame
mlpe_table = cbind(REL, pca[3], euclidean[3], productivity[3], river[3], topo[3], wet_L[3], wet_M[3], wet_H[3], wet_VH[3])
#verify
head(mlpe_table)
colnames(mlpe_table)[1:2] = c("from_ID", "to_ID")
head(mlpe_table)


#B. Creating a structure population column in the table for classifying the distances
#Just a single population
#add to mlpe table:
mlpe_table$POP = 1
head(mlpe_table)

#save df 
write.csv(mlpe_table, "Metafiles/MLPE_table_simonsi.csv")




##### 3b. COVERAGE DEPTH BY SAMPLE ---------------------------------------------
#A. Creating a structure population column in the table for classifying the distances in between individuals from same or different populations:
#Load neutral .vcf file with population geographical information and genetic clusters ID:
snps_neutral = vcfLink("vcf/simonsi_filtered_neutral_LEA_DAPC_TESS.vcf", overwriteID=T)
VCFsummary(snps_neutral) #17 individuals and 12784 SNPs.
names(snps_neutral@meta) #verify col names in metafiles

#B. Estimating coverage depth by sample: 
coverage_ind = c()
for (p in 1:length(snps_neutral@sample_id)){
  beta = Query(Subset(snps_neutral, samples = p), type="site-mean-depth")
  coverage_ind[p] = mean(beta$MEAN_DEPTH, na.rm = T)}
#verify
coverage_ind
#save as metafile
snps@meta$coverage = coverage_ind

#C. save coverage depth by sample:
write.csv(as.data.frame(cbind(snps_neutral@meta$ind_ID,coverage_ind)), "./Metafiles/Depth_coverage_bysamples_simonsi.csv")





##### 4b. VARIABLE RANGES -------------------------------------------------------------
#Variable ranges for all pairwise resistance distances used to run MLPE regression models.

#Load MLPE table if necessary:
mlpe_table = read.csv("Metafiles/MLPE_table_simonsi.csv", row.names = 1)
head(mlpe_table)

#Preparing input for ggplot:
df = mlpe_table[5:(length(mlpe_table)-1)] #explanatory variables start in col #5
df = melt(df)
head(df)

#D. Plot the graph
pdf("boxplot_variables_mlpe_simonsi.pdf", onefile = F)
ggplot(data=df, mapping = aes(x= value, fill = variable)) +
  geom_boxplot() +
  theme_bw()
dev.off()







#####5b. VARIABLE LINEARITY ---------------------------------------------------------
## Assess linearity between response and predictor variables. The graphs show if the slope (trend line) is continually changing; if so, it is not a constant! These variables do not show a linear relationship. With a linear relationship, the slope never changes.

#A. Load MLPE table if necessary:
mlpe_table = read.csv("Metafiles/MLPE_table_simonsi.csv", row.names = 1)
head(mlpe_table)
mlpe_table[, 5:(length(mlpe_table)-1)] = as.data.frame(scale(mlpe_table[, 5:(length(mlpe_table)-1)]))
head(mlpe_table)

names = names(mlpe_table[, 5:(length(mlpe_table)-1)]) #explanatory variables start in col #5
names


#Using Relatedness:
for(i in 1:length(names)){
  pdf(paste0("Results/simonsi/Linearity/Relatedness_ScatterSmooth_simonsi_",names[i],".pdf"), onefile = T)
  scatter.smooth(x = mlpe_table[, names[i]], y=mlpe_table[, "RELgen"], ylab = "Relatedness", xlab=colnames(mlpe_table[names[i]]), lpars = list(col = "red", lwd = 3))
  dev.off()
}

#Using PCA distance:
for(i in 1:length(names)){
  pdf(paste0("Results/simonsi/Linearity/PCA_ScatterSmooth_simonsi_",names[i],".pdf"), onefile = T)
  scatter.smooth(x = mlpe_table[, names[i]], y=mlpe_table[, "PCAgen"], ylab = "PCA-Distance", xlab=colnames(mlpe_table[names[i]]), lpars = list(col = "red", lwd = 3))
  dev.off()
}

##NOT LINEAR:
#River Distance using PCA distance - Use Relatedness




#####6b. UNIVARIATE MODELS ---------------------------------------------------------
#Univariate Models to choose one Habitat Distance matrix
#A. Load MLPE table if necessary:
mlpe_table = read.csv("Metafiles/MLPE_table_simonsi.csv", row.names = 1)
head(mlpe_table)
mlpe_table[, 5:(length(mlpe_table)-1)] = as.data.frame(scale(mlpe_table[, 5:(length(mlpe_table)-1)]))
head(mlpe_table)


#A. Using Relatedness and controlling by population structure
## BUILD FULL MODEL
form = as.formula(paste("RELgen ~ ", 
                        paste(names(mlpe_table)[9:12], collapse = " + "), sep = ""))

Fullmodel <- lme(form,
                 random = ~ 1|POP,
                 correlation = corMLPE(form = ~ from_ID + to_ID),
                 data = mlpe_table, method = "ML")


## CHECK RESIDUALS
RES <- residuals(Fullmodel, type="normalized")
FIT <- fitted(Fullmodel)
plot(FIT, RES) ; abline(0,0, col="red")
acf(RES)

## ALL VARIABLE PLOTS
names(mlpe_table)[9:12]
plot(mlpe_table$wetlands_L, RES) ; abline(0,0, col="red")
plot(mlpe_table$wetlands_M, RES) ; abline(0,0, col="red")
plot(mlpe_table$wetlands_H, RES) ; abline(0,0, col="red")
plot(mlpe_table$wetlands_VH, RES) ; abline(0,0, col="red") 


## RUN PDREDGE 

## Set cluster
## Run paralleled dredge
## by Ben Bolker in: https://stackoverflow.com/questions/30261473/mumin-pdredge-error-with-glmer

## Make cluster
cluster = makeCluster(4, type = "SOCK")  ## also need snow installed

## clusterExport assigns the values on the master R process of the variables named in varlist to variables of the same names in the global environment (aka ‘workspace’) of each node. 
## The environment on the master from which variables are exported defaults to the global environment.
clusterExport(cluster,"mlpe_table")

## clusterEvalQ evaluates a literal expression on each cluster node. 
## It is a parallel version of evalq, and is a convenience function invoking clusterCall
clusterEvalQ(cluster,
             c(library(nlme), library(MuMIn), library(corMLPE)))

## Specify the number of predictor variables and including the max.r function
nrow(mlpe_table) ## 171
options(na.action = na.fail)

## Run pdredge
Allmodels <- MuMIn::pdredge(Fullmodel, rank = "AIC", 
                            m.lim=c(0, 1), cluster)

## SAVE ALL MODELS
save(Allmodels, file="./Results/simonsi/Univariates/Allmodels_univariates_simonsi.RData")

### LOAD ALL MODELS
load(file="./Results/simonsi/Univariates/Allmodels_univariates_simonsi.RData")

## BEST MODELS
nrow(Allmodels) ## 5
BM = model.sel(Allmodels, rank=AIC)
df = as.data.frame(BM[BM$delta <=2, ])
df 
write.csv(as.data.frame(BM), "Results/simonsi/Univariates/Unimodels_habitats_simonsi.csv")

#REL-distance: Low Resistance





##### 8. CORRELATIONS -------------------------------------------------------------
#perform Correlogram showing the correlation between all resistance distances included as predictors in MLPE regression models. Pearson’s correlation coefficients (r). 

#A. Load MLPE table if necessary:
mlpe_table = read.csv("Metafiles/MLPE_table_simonsi.csv", row.names = 1)
head(mlpe_table)

## Check the correlation among environmental variables
ALLvars = as.data.frame(mlpe_table[5:(length(mlpe_table)-1)])
head(ALLvars)
Allvars_cor = cor(ALLvars)
head(Allvars_cor)

write.csv(Allvars_cor, "Results/simonsi/Correlation/Correlation_simonsi.csv")

## Corplot
pdf("Results/simonsi/Correlation/Correlogram_simonsi.pdf", onefile = T)
corrplot(Allvars_cor, 
         method="color", 
         type = "upper", 
         outline = FALSE, 
         diag = F,
         tl.col = "black")
dev.off()


#END;
