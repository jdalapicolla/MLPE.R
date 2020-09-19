#####################  LANDSCAPE GENOMICS TUTORIAL   ##########################
##################   MLPE - STEP 01: DISTANCE MATRICES   ######################

### Script prepared by Jeronymo Dalapicolla, Jamille C. Veiga, Carolina S. Carvalho, Luciana C. Resende-Moreira, and Rodolfo Jaffé ###


##### PRE-ANALYSIS -------------------------------------------------------------

#AIM = ORGANAZING AND SELECTING VARIABLES FOR MODELS:

#Load packages:
library(tidyverse)
library(reshape2)
library(corrplot)
library(r2vcftools)
library(corMLPE)
library(nlme)


#Turn off scientific notation
options(scipen=999)

#Load auxiliary functions:
VCFsummary <- function(snps){
  Nid <- nrow(snps@meta)
  Nsnps <- length(snps@site_id)
  cat(paste(Nid, "individuals and", Nsnps, "SNPs."), sep="\n")
}


##### 1. REARRANGING GENETIC MATRICES FILES ---------------------------------------------
#Calculated in Pipeline for Genetic Structure #Step 3:


########################################################### STEEREI - NON-FLOOED-FOREST:
#A. Load genetic distance matrices.
REL = read.csv("Distances/genetic_dist_Yang_Rel_steerei.csv", row.names = 1)
head(REL)
length(REL[,3]) #190

#B. Format Relatedness data frame in the same order than others dataframe:
#1st) Matrix
n_ind = nrow(unique(REL[1])) #number of individuals
n_ind
mt_all = matrix(0, n_ind, n_ind)
mt_all[lower.tri(mt_all, diag=T)] = as.vector(REL$RELATEDNESS_AJK_Yang)
mt_all = mt_all + t(mt_all) 
diag(mt_all) = diag(mt_all)/2 
mt_all
#verify the matrix
class(mt_all)
dim(mt_all)
nrow(mt_all)
ncol(mt_all)
mt_all[1:10, 1:10]
#Names in rows and cols
rownames(mt_all) = REL$INDV2[1:19]
colnames(mt_all) = REL$INDV2[1:19]
mt_all[1:10, 1:10]

#removing upper diagonal:
mt_all[upper.tri(mt_all, diag = T)] = NA
mt_all

#D. Converting distance matrix into data frame for analyses:
dmat_st = mt_all %>%
  melt %>%
  na.omit %>%
  arrange(., Var1) %>%
  setNames(c("X1", "X2","RELgen"))

#E. Verify for NA or 0 values
head(dmat_st)
summary(dmat_st)
length(dmat_st[,3]) #171


#F. Save rearranged table:
REL = dmat_st
write.csv(dmat_st, "Distances/RELgen_organized_steerei.csv")


#G. PCA distance:
pca = as.matrix(read.csv("Distances/genetic_dist_PCA_BSR_steerei.csv", row.names = 1))
head(pca)
pca[upper.tri(pca, diag = T)] = NA

#H. Converting distance matrix into data frame for analyses:
pca = pca %>%
  melt %>%
  na.omit %>%
  arrange(., Var1) %>%
  setNames(c("X1", "X2","PCAgen"))

#I. Verify for NA or 0 values
head(pca)
summary(pca)
length(pca[,3]) #171

#J. Save rearranged table:
write.csv(pca, "Distances/PCAgen_organized_steerei.csv")








##### 2. LOAD RESISTANCE MATRICES FILES ---------------------------------------------
#A. Load landscape distance matrices:
euclidian = read.csv("Distances/eucl_dist_steerei.csv", row.names = 1)
head(euclidian)

pet = read.csv("Distances/PET_dist_steerei.csv", row.names = 1)
head(pet)

preci = read.csv("Distances/preci_dist_steerei.csv", row.names = 1)
head(preci)

river = read.csv("Distances/river_dist_steerei.csv", row.names = 1)
head(river)

temp = read.csv("Distances/temp_dist_steerei.csv", row.names = 1)
head(temp)

topo = read.csv("Distances/topo_dist_steerei.csv", row.names = 1)
head(topo)

wet_NULL = read.csv("Distances/wetlands_dist_NULL_steerei.csv", row.names = 1)
head(wet_NULL)

wet_L = read.csv("Distances/wetlands_dist_L_steerei.csv", row.names = 1)
head(wet_L)

wet_M = read.csv("Distances/wetlands_dist_M_steerei.csv", row.names = 1)
head(wet_M)

wet_H = read.csv("Distances/wetlands_dist_H_steerei.csv", row.names = 1)
head(wet_H)

wet_VH = read.csv("Distances/wetlands_dist_VH_steerei.csv", row.names = 1)
head(wet_VH)


#B. Verify individuals order
#of colunm #1: 2 by turn
identical(as.character(REL$X1), as.character(pca$X1))
identical(as.character(euclidian$X1), as.character(pca$X1))
identical(as.character(euclidian$X1), as.character(pet$X1))
identical(as.character(preci$X1), as.character(pet$X1))
identical(as.character(preci$X1), as.character(river$X1))
identical(as.character(temp$X1), as.character(river$X1))
identical(as.character(temp$X1), as.character(topo$X1))
identical(as.character(wet_L$X1), as.character(topo$X1))
identical(as.character(wet_L$X1), as.character(wet_M$X1))
identical(as.character(wet_H$X1), as.character(wet_M$X1))
identical(as.character(wet_H$X1), as.character(wet_VH$X1))
identical(as.character(wet_NULL$X1), as.character(wet_VH$X1))

#colunm #2:
identical(as.character(REL$X2), as.character(pca$X2))
identical(as.character(euclidian$X2), as.character(pca$X2))
identical(as.character(euclidian$X2), as.character(pet$X2))
identical(as.character(preci$X2), as.character(pet$X2))
identical(as.character(preci$X2), as.character(river$X2))
identical(as.character(temp$X2), as.character(river$X2))
identical(as.character(temp$X2), as.character(topo$X2))
identical(as.character(wet_L$X2), as.character(topo$X2))
identical(as.character(wet_L$X2), as.character(wet_M$X2))
identical(as.character(wet_H$X2), as.character(wet_M$X2))
identical(as.character(wet_H$X2), as.character(wet_VH$X2))
identical(as.character(wet_NULL$X2), as.character(wet_VH$X2))







##### 3. MERGING ALL DATA IN A SINGLE FILE ---------------------------------------------
#A. United in a single data frame
mlpe_table = cbind(REL, pca[3], euclidian[3], pet[3], preci[3], river[3], temp[3], topo[3], wet_L[3], wet_M[3], wet_H[3], wet_VH[3], wet_NULL[3])
#verify
head(mlpe_table)
colnames(mlpe_table)[1:2] = c("from_ID", "to_ID")
head(mlpe_table)


#B. Creating a struture population column in the table for classifing the distances in between individuals from same or different populations:
#Load neutral .vcf file with populationgeographical information and genetic clusters ID:
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
teste$POP_dist[teste$from_ID == teste$to_ID] =1
teste$POP_dist[teste$from_ID != teste$to_ID] =2
head(teste)

#add to mlpe table:
mlpe_table$POP = teste$POP_dist
head(mlpe_table)

#save df 
write.csv(mlpe_table, "Metafiles/MLPE_table_steerei.csv")








##### 4. COVERAGE DEPTH BY SAMPLE ---------------------------------------------
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





##### 5. VARIABLE RANGES -------------------------------------------------------------
#Variable ranges for all pairwise resistance distances used to run MLPE regression models.
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







#####6. VARIABLE LINEARITY ---------------------------------------------------------
## Assess linearity between response and predictor variables. The graphs show if the slope (trend line) is continually changing; if so, it is not a constant! These variables do not show a linear relationship. With a linear relationship, the slope never changes.

names = names(mlpe_table[, 5:ncol(mlpe_table)]) #explanatory variables start in col #5

#Using Relatedness:
for(i in 1:length(names)){
  pdf(paste0("Linearity/Relatedness_ScatterSmooth_steerei_",names[i],".pdf"), onefile = T)
  scatter.smooth(x = mlpe_table[, names[i]], y=mlpe_table[, "RELgen"], ylab = "Relatedness", xlab=colnames(mlpe_table[names[i]]), lpars = list(col = "red", lwd = 3))
  dev.off()
}

#Using PCA distance:
for(i in 1:length(names)){
  pdf(paste0("Linearity/PCA_ScatterSmooth_steerei_",names[i],".pdf"), onefile = T)
  scatter.smooth(x = mlpe_table[, names[i]], y=mlpe_table[, "PCAgen"], ylab = "PCA-Distance", xlab=colnames(mlpe_table[names[i]]), lpars = list(col = "red", lwd = 3))
  dev.off()
}

##NOT LINEAR:
#River Distance using PCA distance - Use Relatedness








#####7. UNIVARIATE MODELS ---------------------------------------------------------
#Univariates Models to choose one Habitat Distance matrix

#A. Using Relatedness and controlling by population structure
mod_H_r = lme(RELgen ~ wetlands_H, random = ~1|POP, correlation = corMLPE(form = ~ from_ID + to_ID|POP), data = mlpe_table, method = "ML")

mod_L_r = lme(RELgen ~ wetlands_L, random = ~1|POP, correlation = corMLPE(form = ~ from_ID + to_ID|POP), data = mlpe_table, method = "ML")

mod_M_r = lme(RELgen ~ wetlands_M, random = ~1|POP, correlation = corMLPE(form = ~ from_ID + to_ID|POP), data = mlpe_table, method = "ML")

mod_VH_r = lme(RELgen ~ wetlands_VH, random = ~1|POP, correlation = corMLPE(form = ~ from_ID + to_ID|POP), data = mlpe_table, method = "ML")


#B. Comparing models:
uni_models_wet_RELdist_st = MuMIn::model.sel(mod_L_r, mod_M_r, mod_H_r, mod_VH_r, rank= "AICc")

#C, Save Selected models:
uni_models_wet_RELdist_st
write.csv(uni_models_wet_RELdist_st, "Distances/uni_models_wet_RELdist_steerei.csv")

#REL-distance: Very-High Resistance








##### 8. CORRELATIONS -------------------------------------------------------------
#perform Correlogram showing the correlation between all resistance distances included as predictors in MLPE regression models. Pearson’s correlation coefficients (r). 

## Check the correlation among enviromental variables
ALLvars = as.data.frame(mlpe_table[5:(length(mlpe_table)-1)])
head(ALLvars)
Allvars_cor = cor(ALLvars)
head(Allvars_cor)

write.csv(Allvars_cor, "Metafiles/Correlation_steerei.csv")

pdf("Correlogram_steerei.pdf", onefile = T)
corrplot(Allvars_cor, method="pie")
dev.off()






##### 1b. REARRANGING GENETIC MATRICES FILES ---------------------------------------------
#Calculated in Pipeline for Genetic Structure #Step 3:


################################################## SIMONSI - SEASONAL FLOODPLAIN FORESTS:
#Load genetic distance matrices:
REL = read.csv("Distances/genetic_dist_Yang_Rel_simonsi.csv", row.names = 1)
head(REL)
length(REL[,3]) #153

# putting Relatedness data frame in the same order than other dataframe: 1st) Matrix
n_ind = nrow(unique(REL[1])) #number of individuals
n_ind
mt_all = matrix(0, n_ind, n_ind)
mt_all[lower.tri(mt_all, diag=T)] = as.vector(REL$RELATEDNESS_AJK_Yang)
mt_all = mt_all + t(mt_all) 
diag(mt_all) = diag(mt_all)/2 
mt_all
#verify the matrix
class(mt_all)
dim(mt_all)
nrow(mt_all)
ncol(mt_all)
mt_all[1:10, 1:10]
#Names in rows and cols
rownames(mt_all) = REL$INDV2[1:17]
colnames(mt_all) = REL$INDV2[1:17]
mt_all[1:10, 1:10]

#removing upper diagonal:
mt_all[upper.tri(mt_all, diag = T)] = NA
mt_all

#D. Converting distance matrix into data frame for analyses:
dmat_st = mt_all %>%
  melt %>%
  na.omit %>%
  arrange(., Var1) %>%
  setNames(c("X1", "X2","RELgen"))

#E. Verify for NA or 0 values
head(dmat_st)
summary(dmat_st)
length(dmat_st[,3]) #136

#F. Save rearranged table:
REL = dmat_st
write.csv(dmat_st, "Distances/RELgen_organized_steerei.csv")

#G. PCA distance
pca = as.matrix(read.csv("Distances/genetic_dist_PCA_BSR_simonsi.csv", row.names = 1))
head(pca)
pca[upper.tri(pca, diag = T)] = NA

#H. Converting distance matrix into data frame for analyses:
pca = pca %>%
  melt %>%
  na.omit %>%
  arrange(., Var1) %>%
  setNames(c("X1", "X2","PCAgen"))

#I. Verify for NA or 0 values
head(pca)
summary(pca)
length(pca[,3]) #136







##### 2b. LOAD RESISTANCE MATRICES FILES ---------------------------------------------
#A. Load landscape distance matrices:
euclidian = read.csv("Distances/eucl_dist_simonsi.csv", row.names = 1)
head(euclidian)

pet = read.csv("Distances/PET_dist_simonsi.csv", row.names = 1)
head(pet)

preci = read.csv("Distances/preci_dist_simonsi.csv", row.names = 1)
head(preci)

river = read.csv("Distances/river_dist_simonsi.csv", row.names = 1)
head(river)

temp = read.csv("Distances/temp_dist_simonsi.csv", row.names = 1)
head(temp)

topo = read.csv("Distances/topo_dist_simonsi.csv", row.names = 1)
head(topo)

wet_NULL = read.csv("Distances/wetlands_dist_NULL_simonsi.csv", row.names = 1)
head(wet_NULL)

wet_L = read.csv("Distances/wetlands_dist_L_simonsi.csv", row.names = 1)
head(wet_L)

wet_M = read.csv("Distances/wetlands_dist_M_simonsi.csv", row.names = 1)
head(wet_M)

wet_H = read.csv("Distances/wetlands_dist_H_simonsi.csv", row.names = 1)
head(wet_H)

wet_VH = read.csv("Distances/wetlands_dist_VH_simonsi.csv", row.names = 1)
head(wet_VH)


#B. Verify order
#of colunm #1: 2 by turn
identical(as.character(REL$X1), as.character(pca$X1))
identical(as.character(euclidian$X1), as.character(pca$X1))
identical(as.character(euclidian$X1), as.character(pet$X1))
identical(as.character(preci$X1), as.character(pet$X1))
identical(as.character(preci$X1), as.character(river$X1))
identical(as.character(temp$X1), as.character(river$X1))
identical(as.character(temp$X1), as.character(topo$X1))
identical(as.character(wet_L$X1), as.character(topo$X1))
identical(as.character(wet_L$X1), as.character(wet_M$X1))
identical(as.character(wet_H$X1), as.character(wet_M$X1))
identical(as.character(wet_H$X1), as.character(wet_VH$X1))
identical(as.character(wet_NULL$X1), as.character(wet_VH$X1))

#colunm #2:
identical(as.character(REL$X2), as.character(pca$X2))
identical(as.character(euclidian$X2), as.character(pca$X2))
identical(as.character(euclidian$X2), as.character(pet$X2))
identical(as.character(preci$X2), as.character(pet$X2))
identical(as.character(preci$X2), as.character(river$X2))
identical(as.character(temp$X2), as.character(river$X2))
identical(as.character(temp$X2), as.character(topo$X2))
identical(as.character(wet_L$X2), as.character(topo$X2))
identical(as.character(wet_L$X2), as.character(wet_M$X2))
identical(as.character(wet_H$X2), as.character(wet_M$X2))
identical(as.character(wet_H$X2), as.character(wet_VH$X2))
identical(as.character(wet_NULL$X2), as.character(wet_VH$X2))








##### 3b. MERGING ALL DATA IN A SINGLE FILE ---------------------------------------------
#A. United in a single data frame
mlpe_table = cbind(REL, pca[3], euclidian[3], pet[3], preci[3], river[3], temp[3], topo[3], wet_L[3], wet_M[3], wet_H[3], wet_VH[3], wet_NULL[3])
head(mlpe_table)
colnames(mlpe_table)[1:2] = c("from_ID", "to_ID")
head(mlpe_table)

#B. Creating a column with population structure. P. simonsi is panmitic
mlpe_table$POP = 1
head(mlpe_table)

write.csv(mlpe_table, "Metafiles/MLPE_table_simonsi.csv")








##### 4b. COVERAGE DEPTH BY SAMPLE ---------------------------------------------
#A. Load neutral .vcf file:
snps_neutral = vcfLink("vcf/simonsi_filtered_neutral_LEA_DAPC_TESS.vcf", overwriteID=T)
VCFsummary(snps_neutral) #17 individuals and 12784 SNPs.
names(snps_neutral@meta) #verify col names in metafile


#B. Estimating coverage depth by sample:
coverage_ind = c()
for (p in 1:length(snps_neutral@sample_id)){
  beta = Query(Subset(snps_neutral, samples = p), type="site-mean-depth")
  coverage_ind[p] = mean(beta$MEAN_DEPTH, na.rm = T)}
#verify
coverage_ind
#save as metafile
snps_neutral@meta$coverage = coverage_ind


#C. save coverage depth by sample:
write.csv(as.data.frame(cbind(snps_neutral@meta$ind_ID,coverage_ind)), "./Metafiles/Depth_coverage_bysamples_simonsi.csv")








##### 5b. VARIABLE RANGES -------------------------------------------------------------
#Variable ranges for all pairwise resistance distances used to run MLPE regression models.
#Preparing input for ggplot:
df = mlpe_table[5:(length(mlpe_table)-1)]
df = melt(df)
head(df)

#D. Plot the graph
pdf("boxplot_variables_mlpe_simonsi.pdf", onefile = F)
ggplot(data=df, mapping = aes(x= value, fill = variable)) +
  geom_boxplot() +
  theme_bw()
dev.off()












##### 6b. VARIABLE LINEARITY -------------------------------------------------------------
## Assess linearity between response and predictor variables. The graphs show if the slope (trend line) is continually changing; if so, it is not a constant! These variables do not show a linear relationship. With a linear relationship, the slope never changes.

names = names(mlpe_table[, 5:(ncol(mlpe_table)-1)]) #explanatory variables start in col #5

#Using Relatedness:
for(i in 1:length(names)){
  pdf(paste0("Linearity/Relatedness_ScatterSmooth_simonsi_",names[i],".pdf"), onefile = T)
  scatter.smooth(x = mlpe_table[, names[i]], y=mlpe_table[, "RELgen"], ylab = "Relatedness", xlab=colnames(mlpe_table[names[i]]), lpars = list(col = "red", lwd = 3))
  dev.off()
}

#Using PCA distance:
for(i in 1:length(names)){
  pdf(paste0("Linearity/PCA_ScatterSmooth_simonsi_",names[i],".pdf"), onefile = T)
  scatter.smooth(x = mlpe_table[, names[i]], y=mlpe_table[, "PCAgen"], ylab = "PCA-Distance", xlab=colnames(mlpe_table[names[i]]), lpars = list(col = "red", lwd = 3))
  dev.off()
}

##NOT LINEAR:
#River Distance using PCA distance - Use Relatedness








###### 7b. UNIVARIATE MODELS --------------------------------------------------------------
#Univariates Models to choose one Habitat Distance

#A. Performing models:
mod_H_r = lme(RELgen ~ wetlands_H, random = ~1|POP, correlation = corMLPE(form = ~ from_ID + to_ID|POP), data = mlpe_table, method = "ML")
mod_L_r = lme(RELgen ~ wetlands_L, random = ~1|POP, correlation = corMLPE(form = ~ from_ID + to_ID|POP), data = mlpe_table, method = "ML")
mod_M_r = lme(RELgen ~ wetlands_M, random = ~1|POP, correlation = corMLPE(form = ~ from_ID + to_ID|POP), data = mlpe_table, method = "ML")
mod_VH_r = lme(RELgen ~ wetlands_VH, random = ~1|POP, correlation = corMLPE(form = ~ from_ID + to_ID|POP), data = mlpe_table, method = "ML")

#B. Comparing Models:
uni_models_wet_RELdist_si = MuMIn::model.sel(mod_L_r, mod_M_r, mod_H_r, mod_VH_r, rank= "AICc")

#C. Save Results:
uni_models_wet_RELdist_si
write.csv(uni_models_wet_RELdist_si, "Distances/uni_models_wet_RELdist_simonsi.csv")

#REL-distance: Very-High Resistance




###### 8b. Correlations ----------------------------------------------------
#perform Correlogram showing the correlation between all resistance distances included as predictors in MLPE regression models. Pearson’s correlation coefficients (r). 

#A. Check the correlation among enviromental variables
ALLvars = as.data.frame(mlpe_table[5:length(mlpe_table)])
head(ALLvars)
Allvars_cor = cor(ALLvars)
head(Allvars_cor)

#B. Save results
write.csv(Allvars_cor, "Correlation_simonsi.csv")

pdf("Correlogram_simonsi.pdf", onefile = T)
corrplot(Allvars_cor, method="pie")
dev.off()



##END
