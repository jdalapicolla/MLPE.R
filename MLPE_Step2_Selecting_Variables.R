###############################################################################
####################### VALE INSTITUTE OF TECHNOLOGY ##########################
###############################################################################

#####################  LANDSCAPE GENOMICS TUTORIAL   ##########################
#############   MLPE - STEP 02: SELECTING DISTANCE MATRICES   #################

### Script prepared by Carolina S. Carvalho, Jeronymo Dalapicolla, Luciana C. Resende-Moreira, Jamille C. Veiga, and Rodolfo Jaffé ###

#------------------------------------------------------------------------------
#                               PRE-ANALYSIS 
#------------------------------------------------------------------------------

##AIM = ORGANAZING AND SELECTING VARIABLES FOR MODELS:

#Load packages:
library(GeNetIt)
library(tidyverse)
library(reshape2)
library(corrplot)


#------------------------------------------------------------------------------
#                            1. Loading Files
#------------------------------------------------------------------------------

#STEEREI - FLOOED-FOREST:
#Load genetic distance matrices. Calculated in Pipeline for Genetic Structure #Step 3:
REL = read.csv("Distances/genetic_dist_Yang_Rel_steerei.csv", row.names = 1)
head(REL)
length(REL[,3]) #190

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
rownames(mt_all) = REL$INDV2[1:19]
colnames(mt_all) = REL$INDV2[1:19]
mt_all[1:10, 1:10]

#removing upper diagonal:
mt_all[upper.tri(mt_all, diag = T)] = NA
mt_all

#convert distance matrix into data frame for analyses
dmat_st = GeNetIt::dmatrix.df(mt_all)
head(dmat_st)
colnames(dmat_st) = c("X1", "X2","RELgen")

#verify for NA or 0 values
summary(dmat_st)
length(dmat_st[,3]) #171

REL = dmat_st

pca = as.matrix(read.csv("Distances/genetic_dist_PCA_BSR_steerei.csv", row.names = 1))
head(pca)
pca[upper.tri(pca, diag = T)] = NA
pca= GeNetIt::dmatrix.df(pca)
head(pca)
colnames(pca) = c("X1", "X2","PCAgen")
#verify for NA or 0 values
summary(pca)
length(pca[,3]) #171

#----------------------------------------------------
##Load landscape distance matrices:
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

wet_L = read.csv("Distances/wetlands_dist_L_steerei.csv", row.names = 1)
head(wet_L)

wet_M = read.csv("Distances/wetlands_dist_M_steerei.csv", row.names = 1)
head(wet_M)

wet_H = read.csv("Distances/wetlands_dist_H_steerei.csv", row.names = 1)
head(wet_H)

wet_VH = read.csv("Distances/wetlands_dist_VH_steerei.csv", row.names = 1)
head(wet_VH)


#----------------------------------------------- Verify individuals order
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



#-------------------------------------------------------------------#
#                       2. Merge Distances                          #
#-------------------------------------------------------------------#
#United in a single data frame

mlpe_table = cbind(REL, pca[3], euclidian[3], pet[3], preci[3], river[3], temp[3], topo[3], wet_L[3], wet_M[3], wet_H[3], wet_VH[3])
head(mlpe_table)
colnames(mlpe_table)[1:2] = c("from_ID", "to_ID")
head(mlpe_table)

write.csv(mlpe_table, "Metafiles/MLPE_table_steerei.csv")



#-------------------------------------------------------------------#
#                       3. Variable Ranges                          #
#-------------------------------------------------------------------#
#Variable ranges for all pairwise resistance distances used to run MLPE regression models.

#Preparing input for ggplot:  
df = mlpe_table[5:length(mlpe_table)] #explanatory variables start in col #5
df = melt(df)
head(df)

#D. Plot the graph
pdf("boxplot_variables_mlpe_steerei.pdf", onefile = F)
ggplot(data=df, mapping = aes(x= value, fill = variable)) +
  geom_boxplot() +
  theme_bw()
dev.off()



#-------------------------------------------------------------------#
#                         4. Linearity                              #
#-------------------------------------------------------------------#
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



#-------------------------------------------------------------------#
#                     5. Univariate Models                          #
#-------------------------------------------------------------------#
#Univariates Models to choose one Wetlands Distance

##Using PCA distance - Not necessary - you will use only REL:
mod_H = lm(PCAgen ~ wetlands_H, data = mlpe_table)
mod_L = lm(PCAgen ~ wetlands_L, data = mlpe_table)
mod_M = lm(PCAgen ~ wetlands_M, data = mlpe_table)
mod_VH = lm(PCAgen ~ wetlands_VH, data = mlpe_table)
uni_models_wet_PCAdist_st = MuMIn::model.sel(mod_L, mod_M, mod_H, mod_VH, rank= "AICc")
uni_models_wet_PCAdist_st

write.csv(uni_models_wet_PCAdist_st, "Distances/uni_models_wet_PCAdist_steerei.csv")

#Using Relatedness
mod_H_r = lm(RELgen ~ wetlands_H, data = mlpe_table)
mod_L_r = lm(RELgen ~ wetlands_L, data = mlpe_table)
mod_M_r = lm(RELgen ~ wetlands_M, data = mlpe_table)
mod_VH_r = lm(RELgen ~ wetlands_VH, data = mlpe_table)
uni_models_wet_RELdist_st = MuMIn::model.sel(mod_L_r, mod_M_r, mod_H_r, mod_VH_r, rank= "AICc")
uni_models_wet_RELdist_st

write.csv(uni_models_wet_RELdist_st, "Distances/uni_models_wet_RELdist_steerei.csv")


#PCA-distance: Low Resistance
#REL-distance: Very-High Resistance



#-------------------------------------------------------------------#
#                       6. Correlations                             #
#-------------------------------------------------------------------#
#perform Correlogram showing the correlation between all resistance distances included as predictors in MLPE regression models. Pearson’s correlation coefficients (r). 

## Check the correlation among enviromental variables
ALLvars = as.data.frame(mlpe_table[5:length(mlpe_table)])
head(ALLvars)
Allvars_cor = cor(ALLvars)
head(Allvars_cor)

write.csv(Allvars_cor, "Correlation_steerei.csv")

pdf("Correlogram_steerei.pdf", onefile = T)
corrplot(Allvars_cor, method="pie")
dev.off()


#--------------------------------------------------------------------


#-------------------------------------------------------------------#
#                       1. Loading Files                            #
#-------------------------------------------------------------------#

#SIMONSI - DRY-FOREST:
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

#convert distance matrix into data frame for analyses
dmat_st = GeNetIt::dmatrix.df(mt_all)
head(dmat_st)
colnames(dmat_st) = c("X1", "X2","RELgen")

#verify for NA or 0 values
summary(dmat_st)
length(dmat_st[,3]) #136

REL = dmat_st

pca = as.matrix(read.csv("Distances/genetic_dist_PCA_BSR_simonsi.csv", row.names = 1))
head(pca)
pca[upper.tri(pca, diag = T)] = NA
pca= GeNetIt::dmatrix.df(pca)
head(pca)
colnames(pca) = c("X1", "X2","PCAgen")
#verify for NA or 0 values
summary(pca)
length(pca[,3]) #136



#----------------------------------------------------
##Load landscape distance matrices:
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

wet_L = read.csv("Distances/wetlands_dist_L_simonsi.csv", row.names = 1)
head(wet_L)

wet_M = read.csv("Distances/wetlands_dist_M_simonsi.csv", row.names = 1)
head(wet_M)

wet_H = read.csv("Distances/wetlands_dist_H_simonsi.csv", row.names = 1)
head(wet_H)

wet_VH = read.csv("Distances/wetlands_dist_VH_simonsi.csv", row.names = 1)
head(wet_VH)


#-----------------------------------------------verify order
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


#-------------------------------------------------------------------#
#                       2. Merge Distances                          #
#-------------------------------------------------------------------#
#United in a single data frame

mlpe_table = cbind(REL, pca[3], euclidian[3], pet[3], preci[3], river[3], temp[3], topo[3], wet_L[3], wet_M[3], wet_H[3], wet_VH[3])
head(mlpe_table)
colnames(mlpe_table)[1:2] = c("from_ID", "to_ID")
head(mlpe_table)

write.csv(mlpe_table, "Metafiles/MLPE_table_simonsi.csv")



#-------------------------------------------------------------------#
#                       3. Varaibles Range                          #
#-------------------------------------------------------------------#
#Variable ranges for all pairwise resistance distances used to run MLPE regression models.

df = mlpe_table[5:length(mlpe_table)]
df = melt(df)
head(df)

#D. Plot the graph
pdf("boxplot_variables_mlpe_simonsi.pdf", onefile = F)
ggplot(data=df, mapping = aes(x= value, fill = variable)) +
  geom_boxplot() +
  theme_bw()
dev.off()



#-------------------------------------------------------------------#
#                         4. Linearity                              #
#-------------------------------------------------------------------#
## Assess linearity between response and predictor variables. The graphs show if the slope (trend line) is continually changing; if so, it is not a constant! These variables do not show a linear relationship. With a linear relationship, the slope never changes.

names = names(mlpe_table[, 5:ncol(mlpe_table)]) #explanatory variables start in col #5

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



#-------------------------------------------------------------------#
#                     5. Univariate Models                          #
#-------------------------------------------------------------------#
#Univariates Models to choose one Wetlands Distance

##Using PCA distance - Not necessary - you will use only REL:
mod_H = lm(PCAgen ~ wetlands_H, data = mlpe_table)
mod_L = lm(PCAgen ~ wetlands_L, data = mlpe_table)
mod_M = lm(PCAgen ~ wetlands_M, data = mlpe_table)
mod_VH = lm(PCAgen ~ wetlands_VH, data = mlpe_table)
uni_models_wet_PCAdist_si = MuMIn::model.sel(mod_L, mod_M, mod_H, mod_VH, rank= "AICc")
uni_models_wet_PCAdist_si

write.csv(uni_models_wet_PCAdist_si, "Distances/uni_models_wet_PCAdist_simonsi.csv")

#Using Relatedness
mod_H_r = lm(RELgen ~ wetlands_H, data = mlpe_table)
mod_L_r = lm(RELgen ~ wetlands_L, data = mlpe_table)
mod_M_r = lm(RELgen ~ wetlands_M, data = mlpe_table)
mod_VH_r = lm(RELgen ~ wetlands_VH, data = mlpe_table)
uni_models_wet_RELdist_si = MuMIn::model.sel(mod_L_r, mod_M_r, mod_H_r, mod_VH_r, rank= "AICc")
uni_models_wet_RELdist_si

write.csv(uni_models_wet_RELdist_si, "Distances/uni_models_wet_RELdist_simonsi.csv")


#PCA-distance: Low Resistance
#REL-distance: Very-High Resistance



#-------------------------------------------------------------------#
#                       6. Correlations                             #
#-------------------------------------------------------------------#
#perform Correlogram showing the correlation between all resistance distances included as predictors in MLPE regression models. Pearson’s correlation coefficients (r). 

## Check the correlation among enviromental variables
ALLvars = as.data.frame(mlpe_table[5:length(mlpe_table)])
head(ALLvars)
Allvars_cor = cor(ALLvars)
head(Allvars_cor)

write.csv(Allvars_cor, "Correlation_simonsi.csv")

pdf("Correlogram_simonsi.pdf", onefile = T)
corrplot(Allvars_cor, method="pie")
dev.off()



##END
