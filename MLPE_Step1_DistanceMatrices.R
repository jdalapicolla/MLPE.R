#####################  LANDSCAPE GENOMICS TUTORIAL   ##########################
##################   MLPE - STEP 01: DISTANCE MATRICES   ######################

### Script prepared by Jeronymo Dalapicolla, Jamille C. Veiga, Carolina S. Carvalho, Luciana C. Resende-Moreira, and Rodolfo Jaffé ###


##### PRE-ANALYSIS -------------------------------------------------------------

#AIM = CREATING DIFFERENT DISTANCE MATRICES TO USE AS VARIABLES FOR IBR MODELS

#Load packages
library(r2vcftools)
library(usedist)
library(raster)
library(rgdal)
library(geosphere)
library(ade4)
library(vegan)
library(reshape2)
library(topoDistance)
library(riverdist)
library(adegenet)
library(tidyverse)
library(sf)
library(gdistance)

#Turn off scientific notation
options(scipen=999)






##### 1. LOADING FILES   ------------------------------------------------------
#A. Load CSV files with geographical information on samples:
steerei = read.csv("Metafiles/metafiles_steerei.csv", row.names = 1)
head(steerei)

simonsi = read.csv("Metafiles/metafiles_simonsi.csv", row.names = 1)
head(simonsi)

#B. Create a data frame with the geographical coordenates:
steerei_coord = steerei[,c(7,4:5)]
head(steerei_coord)

simonsi_coord = simonsi[,c(8,5:6)]
head(simonsi_coord)

#C. Set long and lat colunms
coordinates(steerei_coord) = steerei_coord[,c(2,3)]
projection(steerei_coord) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84")

coordinates(simonsi_coord) = simonsi_coord[,c(2,3)]
projection(simonsi_coord) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84")






###### 2. EUCLIDEAN (GEOGRAPHIC) DISTANCE --------------------------------------------
#STEEREI:
#A. Creating euclidean distance:
eucl_dist_steerei = distm(steerei_coord, fun=distGeo)
rownames(eucl_dist_steerei) = steerei[,7]
colnames(eucl_dist_steerei) = steerei[,7]
eucl_dist_steerei

#B. Converting in KM
eucl_dist_steerei = eucl_dist_steerei/100

#C. Removing upper diagonal:
eucl_dist_steerei[upper.tri(eucl_dist_steerei, diag = T)] = NA
eucl_dist_steerei

#D. Convert distance matrix into data frame for analyses
eucl_dist_steerei = eucl_dist_steerei %>%
  melt %>%
  na.omit %>%
  arrange(., Var1) %>%
  setNames(c("X1", "X2","eucl_dist"))

#E. Verify for NA or 0 values
head(eucl_dist_steerei)
summary(eucl_dist_steerei)
length(eucl_dist_steerei[,3]) #171

#F. Save results
write.csv(eucl_dist_steerei, "Distances/eucl_dist_steerei.csv")



#SIMONSI:
#A. Creating euclidean distance:
eucl_dist_simonsi = distm(simonsi_coord, fun=distGeo)
rownames(eucl_dist_simonsi) = simonsi[,8]
colnames(eucl_dist_simonsi) = simonsi[,8]
eucl_dist_simonsi

#B. Converting in KM
eucl_dist_simonsi = eucl_dist_simonsi/100

#C. Removing upper diagonal:
eucl_dist_simonsi[upper.tri(eucl_dist_simonsi, diag = T)] = NA
eucl_dist_simonsi

#D. Converting distance matrix into data frame for analyses:
eucl_dist_simonsi = eucl_dist_simonsi %>%
  melt %>%
  na.omit %>%
  arrange(., Var1) %>%
  setNames(c("X1", "X2","eucl_dist"))

#E. Verify for NA or 0 values
head(eucl_dist_steerei)
summary(eucl_dist_simonsi)
length(eucl_dist_simonsi[,3]) #136

#F. Save results
write.csv(eucl_dist_simonsi, "Distances/eucl_dist_simonsi.csv")






##### 3. TEMPERATURE DISSIMILARITY ------------------------------------------------
### Based on WorldClim 2.1 - Monthly Temperature
#A. Loading rasters
current.list = list.files(path="Rasters/WC2_tavg", pattern =".tif", full.names=TRUE)
Temp_WC2 = raster::stack(current.list)
names(Temp_WC2)
projection(Temp_WC2) = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +units=m"

#STEEREI:
#A. Extracting climatic data
Temp_values_WC2_st = as.data.frame(raster::extract(Temp_WC2, steerei_coord[,c(2,3)]))
class(Temp_values_WC2_st)
head(Temp_values_WC2_st)
summary(Temp_values_WC2_st)

#B. Run PCA and extract the first components
pca_temp_WC2_st = dudi.pca(Temp_values_WC2_st, center = TRUE, scale = TRUE, scannf = FALSE, nf=length(Temp_values_WC2_st))

#C. Saving PCA results
#% of PC variation
perc_pca = as.data.frame(round(pca_temp_WC2_st$eig,3))
soma = sum (perc_pca)
perc_pca[,2] = round(((perc_pca/soma)*100),3)
colnames(perc_pca) = c("Eigenvalues", "Contribution (%)")
perc_pca
write.csv(perc_pca, file = "./PCA/PCA_temp/contribution_pc_eig_temp_WC2_steerei.csv")

#row coordinates i.e. the principal components
pca_temp_WC2_st$li
write.csv(pca_temp_WC2_st$li, file = "./PCA/PCA_temp/PCompPCA_coord_temp_WC2_steerei.csv")

#column normed scores i.e. the principal axes
pca_temp_WC2_st$c1
write.csv(as.data.frame(pca_temp_WC2_st$c1), file = "./PCA/PCA_temp/AXcompPCA_coord_temp_WC2_steerei.csv")

#D. Selecting PC by Broken Stick: 1PCs
pca_temp_BS_st = prcomp(Temp_values_WC2_st, center=TRUE, scale=TRUE)
summary(pca_temp_BS_st) #verify if same PC% was recovered by dudi.pca function
screeplot(pca_temp_BS_st, bstick=TRUE, type="lines")
screeplot(pca_temp_BS_st, bstick=TRUE, type="barplot")
#Broken Stick Rule: Principal components should be retained as long as observed eigenvalues are higher than corresponding random broken stick components. Jackson 1993 & Legendre & Legendre 2012

summary(pca_WC2_BS)
#1PCs are 86.13% of variance
n_pcs = 1

#E. Calculate PCA-based distance based on Broken Stick Rule
PC_distLEA_prcomp_st = ecodist::distance(pca_temp_BS_st$x[,1:n_pcs], method = "euclidean")
#verify
head(PC_distLEA_prcomp_st)
class(PC_distLEA_prcomp_st)
#convert to matrix
t_prcompLEA_st = as.matrix(PC_distLEA_prcomp_st)
t_prcompLEA_st[1:10,1:10]
#row and col names:
rownames(t_prcompLEA_st) = steerei[,7]
colnames(t_prcompLEA_st) = steerei[,7]
t_prcompLEA_st[1:10,1:10]

#F. Removing upper diagonal:
t_prcompLEA_st[upper.tri(t_prcompLEA_st, diag = T)] = NA
t_prcompLEA_st

#G. Converting distance matrix into data frame for analyses:
temp_dist_st = t_prcompLEA_st %>%
  melt %>%
  na.omit %>%
  arrange(., Var1) %>%
  setNames(c("X1", "X2","temperature"))

#H. Verify for NA or 0 values
head(temp_dist_st)
summary(temp_dist_st)
length(temp_dist_st[,3]) #171

#I. Replace 0 to 0.001
temp_dist_st[,3][temp_dist_st[,3] == 0] = 0.001

#J. Verify for NA or 0 values
summary(temp_dist_st)
length(temp_dist_st[,3]) #171

#K. Save results
write.csv(temp_dist_st, "Distances/temp_dist_steerei.csv")




#SIMONSI:
#A. Extracting climatic data
Temp_values_WC2_si = as.data.frame(raster::extract(Temp_WC2, simonsi_coord[,c(2,3)]))
class(Temp_values_WC2_si)
head(Temp_values_WC2_si)
summary(Temp_values_WC2_si)

#B. Run PCA and extract the first components
pca_temp_WC2_si = dudi.pca(Temp_values_WC2_si, center = TRUE, scale = TRUE, scannf = FALSE, nf=length(Temp_values_WC2_si))

#C. Saving PCA results
#% of PC variation
perc_pca = as.data.frame(round(pca_temp_WC2_si$eig,3))
soma = sum (perc_pca)
perc_pca[,2] = round(((perc_pca/soma)*100),3)
colnames(perc_pca) = c("Eigenvalues", "Contribution (%)")
perc_pca
write.csv(perc_pca, file = "./PCA/PCA_temp/contribution_pc_eig_temp_WC2_simonsi.csv")
#row coordinates i.e. the principal components
pca_temp_WC2_si$li
write.csv(pca_temp_WC2_si$li, file = "./PCA/PCA_temp/PCompPCA_coord_temp_WC2_simonsi.csv")
#column normed scores i.e. the principal axes
pca_temp_WC2_si$c1
write.csv(as.data.frame(pca_temp_WC2_si$c1), file = "./PCA/PCA_temp/AXcompPCA_coord_temp_WC2_simonsi.csv")

#D. Selecting PC by Broken Stick: 1PCs
pca_temp_BS_si = prcomp(Temp_values_WC2_si, center=TRUE, scale=TRUE)
summary(pca_temp_BS_si) #verify if same PC% was recovered by dudi.pca function
screeplot(pca_temp_BS_si, bstick=TRUE, type="lines")
screeplot(pca_temp_BS_si, bstick=TRUE, type="barplot")
#Broken Stick Rule: Principal components should be retained as long as observed eigenvalues are higher than corresponding random broken stick components. Jackson 1993 & Legendre & Legendre 2012

summary(pca_temp_BS_si)
#1PCs are 94.97% of variance
n_pcs = 1

#E. Calculate PCA-based distance based on Broken Stick Rule
PC_distLEA_prcomp_si = ecodist::distance(pca_temp_BS_si$x[,1:n_pcs], method = "euclidean")
#verify
head(PC_distLEA_prcomp_si)
class(PC_distLEA_prcomp_si)
#convert to matrix
t_prcompLEA_si = as.matrix(PC_distLEA_prcomp_si)
t_prcompLEA_si[1:10,1:10]
#row and col names:
rownames(t_prcompLEA_si) = simonsi[,8]
colnames(t_prcompLEA_si) = simonsi[,8]
t_prcompLEA_si[1:10,1:10]

#F. Removing upper diagonal:
t_prcompLEA_si[upper.tri(t_prcompLEA_si, diag = T)] = NA
t_prcompLEA_si

#G. Converting distance matrix into data frame for analyses:
temp_dist_si = t_prcompLEA_si %>%
  melt %>%
  na.omit %>%
  arrange(., Var1) %>%
  setNames(c("X1", "X2","temperature"))

#H. Verify for NA or 0 values
head(temp_dist_si)
summary(temp_dist_si)
length(temp_dist_si[,3]) #136

#I. Replace 0 to 0.001
temp_dist_si[,3][temp_dist_si[,3] == 0] = 0.001

#J. Verify for NA or 0 values
summary(temp_dist_si)
length(temp_dist_si[,3]) #171

#K. Save results
write.csv(temp_dist_si, "Distances/temp_dist_simonsi.csv")









##### 4. PRECIPITATION DISSIMILARITY ------------------------------------------
### Based on WorldClim 2.1 - Monthly Precipitation
#A. Loading rasters
current.list = list.files(path="Rasters/precipitation", pattern =".tif", full.names=TRUE)
preci_WC2 = raster::stack(current.list)
names(preci_WC2)
projection(preci_WC2) = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +units=m"


#STEEREI
#A. Extracting climatic data
preci_values_WC2_st = as.data.frame(raster::extract(preci_WC2, steerei_coord[,c(2,3)]))
class(preci_values_WC2_st)
head(preci_values_WC2_st)
summary(preci_values_WC2_st)

#B. Run PCA and extract the first components
pca_preci_WC2_st = dudi.pca(preci_values_WC2_st, center = TRUE, scale = TRUE, scannf = FALSE, nf=length(preci_values_WC2_st))

#C. Saving PCA results
#% of PC variation
perc_pca = as.data.frame(round(pca_preci_WC2_st$eig,3))
soma = sum (perc_pca)
perc_pca[,2] = round(((perc_pca/soma)*100),3)
colnames(perc_pca) = c("Eigenvalues", "Contribution (%)")
perc_pca
write.csv(perc_pca, file = "./PCA/PCA_preci/contribution_pc_eig_preci_WC2_steerei.csv")

#row coordinates i.e. the principal components
pca_preci_WC2_st$li
write.csv(pca_preci_WC2_st$li, file = "./PCA/PCA_preci/PCompPCA_coord_preci_WC2_steerei.csv")

#column normed scores i.e. the principal axes
pca_preci_WC2_st$c1
write.csv(as.data.frame(pca_preci_WC2_st$c1), file = "./PCA/PCA_temp/AXcompPCA_coord_preci_WC2_steerei.csv")

#D. Selecting PC by Broken Stick: 1PCs
pca_preci_BS_st = prcomp(preci_values_WC2_st, center=TRUE, scale=TRUE)
summary(pca_preci_BS_st) #verify if same PC% was recovered by dudi.pca function
screeplot(pca_preci_BS_st, bstick=TRUE, type="lines")
screeplot(pca_preci_BS_st, bstick=TRUE, type="barplot")
#Broken Stick Rule: Principal components should be retained as long as observed eigenvalues are higher than corresponding random broken stick components. Jackson 1993 & Legendre & Legendre 2012

summary(pca_preci_BS_st)
#1PCs are 72.74% of variance
n_pcs = 1

#E. Calculate PCA-based distance based on Broken Stick Rule
PC_distLEA_prcomp_st = ecodist::distance(pca_preci_BS_st$x[,1:n_pcs], method = "euclidean")
#verify
head(PC_distLEA_prcomp_st)
class(PC_distLEA_prcomp_st)
#convert to matrix
t_prcompLEA_st = as.matrix(PC_distLEA_prcomp_st)
t_prcompLEA_st[1:10,1:10]
#row and col names:
rownames(t_prcompLEA_st) = steerei[,7]
colnames(t_prcompLEA_st) = steerei[,7]
t_prcompLEA_st[1:10,1:10]

#F. Removing upper diagonal:
t_prcompLEA_st[upper.tri(t_prcompLEA_st, diag = T)] = NA
t_prcompLEA_st

#G. Converting distance matrix into data frame for analyses:
preci_dist_st = t_prcompLEA_st %>%
  melt %>%
  na.omit %>%
  arrange(., Var1) %>%
  setNames(c("X1", "X2","precipitation"))

#H. Verify for NA or 0 values
head(preci_dist_st)
summary(preci_dist_st)
length(preci_dist_st[,3]) #171

#I. replace 0 to 0.001
preci_dist_st[,3][preci_dist_st[,3] == 0] = 0.00001

#J. verify for NA or 0 values
summary(preci_dist_st)
length(preci_dist_st[,3]) #171

#K. Save results
write.csv(preci_dist_st, "Distances/preci_dist_steerei.csv")



#SIMONSI
#A. Extracting climatic data
preci_values_WC2_si = as.data.frame(raster::extract(preci_WC2, simonsi_coord[,c(2,3)]))
class(preci_values_WC2_si)
head(preci_values_WC2_si)
summary(preci_values_WC2_si)

#B. Run PCA and extract the first components
pca_preci_WC2_si = dudi.pca(preci_values_WC2_si, center = TRUE, scale = TRUE, scannf = FALSE, nf=length(preci_values_WC2_si))

#C. Saving PCA results
#% of PC variation
perc_pca = as.data.frame(round(pca_preci_WC2_si$eig,3))
soma = sum (perc_pca)
perc_pca[,2] = round(((perc_pca/soma)*100),3)
colnames(perc_pca) = c("Eigenvalues", "Contribution (%)")
perc_pca
write.csv(perc_pca, file = "./PCA/PCA_preci/contribution_pc_eig_preci_WC2_simonsi.csv")

#row coordinates i.e. the principal components
pca_preci_WC2_si$li
write.csv(pca_preci_WC2_si$li, file = "./PCA/PCA_preci/PCompPCA_coord_preci_WC2_simonsi.csv")

#column normed scores i.e. the principal axes
pca_preci_WC2_si$c1
write.csv(as.data.frame(pca_preci_WC2_si$c1), file = "./PCA/PCA_temp/AXcompPCA_coord_preci_WC2_simonsi.csv")

#D. Selecting PC by Broken Stick: 2PCs
pca_preci_BS_si = prcomp(preci_values_WC2_si, center=TRUE, scale=TRUE)
summary(pca_preci_BS_si) #verify if same PC% was recovered by dudi.pca function
screeplot(pca_preci_BS_si, bstick=TRUE, type="lines")
screeplot(pca_preci_BS_si, bstick=TRUE, type="barplot")
#Broken Stick Rule: Principal components should be retained as long as observed eigenvalues are higher than corresponding random broken stick components. Jackson 1993 & Legendre & Legendre 2012

summary(pca_preci_BS_si)
#2PCs are 86.44% of variance
n_pcs = 2

#E. Calculate PCA-based distance based on Broken Stick Rule
PC_distLEA_prcomp_si = ecodist::distance(pca_preci_BS_si$x[,1:n_pcs], method = "euclidean")
#verify
head(PC_distLEA_prcomp_si)
class(PC_distLEA_prcomp_si)
#convert to matrix
t_prcompLEA_si = as.matrix(PC_distLEA_prcomp_si)
t_prcompLEA_si[1:10,1:10]
#row and col names:
rownames(t_prcompLEA_si) = simonsi[,8]
colnames(t_prcompLEA_si) = simonsi[,8]
t_prcompLEA_si[1:10,1:10]

#F. removing upper diagonal:
t_prcompLEA_si[upper.tri(t_prcompLEA_si, diag = T)] = NA
t_prcompLEA_si

#G. Converting distance matrix into data frame for analyses:
preci_dist_si = t_prcompLEA_si %>%
  melt %>%
  na.omit %>%
  arrange(., Var1) %>%
  setNames(c("X1", "X2","precipitation"))

#H. Verify for NA or 0 values
head(preci_dist_si)
summary(preci_dist_si)
length(preci_dist_si[,3]) #136

#I. Replace 0 to 0.001
preci_dist_si[,3][preci_dist_si[,3] == 0] = 0.00001

#J. Verify for NA or 0 values
summary(preci_dist_si)
length(preci_dist_si[,3]) #136

#K. Save results
write.csv(preci_dist_si, "Distances/preci_dist_simonsi.csv")









##### 5. PET DISSIMILARITY -----------------------------------------------------
### Based on 6 Environ variables - PET (Potential Evapotranspiration)
#A. Loading rasters
current.list = list.files(path="Rasters/PET", pattern =".tif", full.names=TRUE)
PET_WC2 = raster::stack(current.list)
names(PET_WC2)
projection(PET_WC2) = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +units=m"


##STEEREI:
#A. Extracting climatic data
PET_values_st = as.data.frame(raster::extract(PET_WC2, steerei_coord[,c(2,3)]))
class(PET_values_st)
head(PET_values_st)
summary(PET_values_st)

#B. Run PCA and extract the components
pca_PET_st = dudi.pca(PET_values_st, center = TRUE, scale = TRUE, scannf = FALSE, nf=length(PET_values_st))

#C. Saving PCA results
#% of PC variation
perc_pca = as.data.frame(round(pca_PET_st$eig,3))
soma = sum (perc_pca)
perc_pca[,2] = round(((perc_pca/soma)*100),3)
colnames(perc_pca) = c("Eigenvalues", "Contribution (%)")
perc_pca
write.csv(perc_pca, file = "./PCA/PCA_PET/contribution_pc_eig_PET_WC2_steerei.csv")

#row coordinates i.e. the principal components
pca_PET_st$li
write.csv(pca_PET_st$li, file = "./PCA/PCA_PET/PCompPCA_coord_PET_WC2_steerei.csv")

#column normed scores i.e. the principal axes
pca_PET_st$c1
write.csv(as.data.frame(pca_PET_st$c1), file = "./PCA/PCA_PET/AXcompPCA_coord_PET_WC2_steerei.csv")

#D. Selecting PC by Broken Stick: 3PCs
pca_PET_BS_st = prcomp(PET_values_st, center=TRUE, scale=TRUE)
summary(pca_PET_BS_st) #verify if same PC% was recovered by dudi.pca function
screeplot(pca_PET_BS_st, bstick=TRUE, type="lines")
screeplot(pca_PET_BS_st, bstick=TRUE, type="barplot")
#Broken Stick Rule: Principal components should be retained as long as observed eigenvalues are higher than corresponding random broken stick components. Jackson 1993 & Legendre & Legendre 2012

summary(pca_PET_BS_st)
#3PCs are 95.65% of variance
n_pcs = 3

#E. Calculate PCA-based distance based on Broken Stick Rule
PC_distLEA_prcomp_st = ecodist::distance(pca_PET_BS_st$x[,1:n_pcs], method = "euclidean")
#verify
head(PC_distLEA_prcomp_st)
class(PC_distLEA_prcomp_st)
#convert to matrix
t_prcompLEA_st = as.matrix(PC_distLEA_prcomp_st)
t_prcompLEA_st[1:10,1:10]
#row and col names:
rownames(t_prcompLEA_st) = steerei[,7]
colnames(t_prcompLEA_st) = steerei[,7]
t_prcompLEA_st[1:10,1:10]

#F. Removing upper diagonal:
t_prcompLEA_st[upper.tri(t_prcompLEA_st, diag = T)] = NA
t_prcompLEA_st

#G. Converting distance matrix into data frame for analyses:
PET_dist_st = t_prcompLEA_st %>%
  melt %>%
  na.omit %>%
  arrange(., Var1) %>%
  setNames(c("X1", "X2","PET"))

#H. Verify for NA or 0 values
head(PET_dist_st)
summary(PET_dist_st)
length(PET_dist_st[,3]) #171

#I. Replace 0 to 0.001
PET_dist_st[,3][PET_dist_st[,3] == 0] = 0.001

#J. Verify for NA or 0 values
summary(PET_dist_st)
length(PET_dist_st[,3]) #171

#K. Save results
write.csv(PET_dist_st, "Distances/PET_dist_steerei.csv")



##SIMONSI:
#A. Extracting climatic data
PET_values_si = as.data.frame(raster::extract(PET_WC2, simonsi_coord[,c(2,3)]))
class(PET_values_si)
head(PET_values_si)
summary(PET_values_si)

#B. Run PCA and extract the first components
pca_PET_si = dudi.pca(PET_values_si, center = TRUE, scale = TRUE, scannf = FALSE, nf=length(PET_values_si))

#C. Saving PCA results
#% of PC variation
perc_pca = as.data.frame(round(pca_PET_si$eig,3))
soma = sum (perc_pca)
perc_pca[,2] = round(((perc_pca/soma)*100),3)
colnames(perc_pca) = c("Eigenvalues", "Contribution (%)")
perc_pca
write.csv(perc_pca, file = "./PCA/PCA_PET/contribution_pc_eig_PET_WC2_simonsi.csv")

#row coordinates i.e. the principal components
pca_PET_si$li
write.csv(pca_PET_si$li, file = "./PCA/PCA_PET/PCompPCA_coord_PET_WC2_simonsi.csv")

#column normed scores i.e. the principal axes
pca_PET_si$c1
write.csv(as.data.frame(pca_PET_si$c1), file = "./PCA/PCA_PET/AXcompPCA_coord_PET_WC2_simonsi.csv")

#D. Selecting PC by Broken Stick: 2PCs
pca_PET_BS_si = prcomp(PET_values_si, center=TRUE, scale=TRUE)
summary(pca_PET_BS_si) #verify if same PC% was recovered by dudi.pca function
screeplot(pca_PET_BS_si, bstick=TRUE, type="lines")
screeplot(pca_PET_BS_si, bstick=TRUE, type="barplot")
#Broken Stick Rule: Principal components should be retained as long as observed eigenvalues are higher than corresponding random broken stick components. Jackson 1993 & Legendre & Legendre 2012

summary(pca_PET_BS_si)
#2PCs are 93.93% of variance
n_pcs = 2

#E. Calculate PCA-based distance based on Broken Stick Rule
PC_distLEA_prcomp_si = ecodist::distance(pca_PET_BS_si$x[,1:n_pcs], method = "euclidean")
#verify
head(PC_distLEA_prcomp_si)
class(PC_distLEA_prcomp_si)
#convert to matrix
t_prcompLEA_si = as.matrix(PC_distLEA_prcomp_si)
t_prcompLEA_si[1:10,1:10]
#row and col names:
rownames(t_prcompLEA_si) = simonsi[,8]
colnames(t_prcompLEA_si) = simonsi[,8]
t_prcompLEA_si[1:10,1:10]

#F. Removing upper diagonal:
t_prcompLEA_si[upper.tri(t_prcompLEA_si, diag = T)] = NA
t_prcompLEA_si

#G. Converting distance matrix into data frame for analyses:
PET_dist_si = t_prcompLEA_si %>%
  melt %>%
  na.omit %>%
  arrange(., Var1) %>%
  setNames(c("X1", "X2","PET"))

#H. Verify for NA or 0 values
head(PET_dist_si)
summary(PET_dist_si)
length(PET_dist_si[,3]) #136

#I. Replace 0 to 0.001
PET_dist_si[,3][PET_dist_si[,3] == 0] = 0.001

#J. Verify for NA or 0 values
summary(PET_dist_si)
length(PET_dist_si[,3]) #136

#K. Save results
write.csv(PET_dist_si, "Distances/PET_dist_simonsi.csv")









##### 6. TOPOGRAPHIC DISTANCE ----------------------------------------------------
### Based on SRTM elevation data from WorldClim 2.1
#A. Loading rasters
elevation = raster("Rasters/DEM/wc2.1_30s_elev.tif")

#B. Choose a extension for study area to reduce time for analyses:
ext = extent(-80, -60, -20, -0)
elevation = crop(elevation, ext)
projection(elevation) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84")

#C. Save the study area raster:
writeRaster(elevation, filename="Rasters/DEM/elevation_WA.tif", format="GTiff", overwrite=TRUE, options=c('TFW=YES'), bylayer=F)

#D. Load the reduced DEM raster
elevation = raster("Rasters/DEM/elevation_WA.tif")
projection(elevation) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84", doCheckCRSArgs = F)


#STEEREI:
#A. create a topographic distance object:
topo_dist_st = topoDist(elevation, steerei_coord[2:3], directions = 8, paths = FALSE, zweight = 1)
class(topo_dist_st)
topo_dist_st[1:10,1:10]

#B. Convert in KM
topo_dist_st = topo_dist_st/100
topo_dist_st[1:10,1:10]

#C. Rremoving upper diagonal:
topo_dist_st[upper.tri(topo_dist_st, diag = T)] = NA
topo_dist_st

#D. Converting distance matrix into data frame for analyses:
topo_dist_st = topo_dist_st %>%
  melt %>%
  na.omit %>%
  arrange(., Var1) %>%
  setNames(c("X1", "X2","topography"))

#E. Verify for NA or 0 values
head(topo_dist_st)
summary(topo_dist_st)
length(topo_dist_st[,3]) #171

#F. Replace 0 to 0.001
topo_dist_st[,3][topo_dist_st[,3] == 0] = 0.001

#G. Verify for NA or 0 values
summary(topo_dist_st)
length(topo_dist_st[,3]) #136

#H. Save results
write.csv(topo_dist_st, "Distances/topo_dist_steerei.csv")



#SIMONSI:
#A. create a topographic distance object:
topo_dist_si = topoDist(elevation, simonsi_coord[2:3], directions = 8, paths = FALSE, zweight = 1)
class(topo_dist_si)
topo_dist_si[1:10,1:10]

#B. Convert in KM
topo_dist_si = topo_dist_si/100
topo_dist_si[1:10,1:10]

#C. Removing upper diagonal:
topo_dist_si[upper.tri(topo_dist_si, diag = T)] = NA
topo_dist_si

#D. Converting distance matrix into data frame for analyses:
topo_dist_si = topo_dist_si %>%
  melt %>%
  na.omit %>%
  arrange(., Var1) %>%
  setNames(c("X1", "X2","topography"))

#E. Verify for NA or 0 values
head(topo_dist_si)
summary(topo_dist_si)
length(topo_dist_si[,3]) #136

#F. Replace 0 to 0.001
topo_dist_si[,3][topo_dist_si[,3] == 0] = 0.001

#G. Verify for NA or 0 values
summary(topo_dist_si)
length(topo_dist_si[,3]) #136

#H. Save results
write.csv(topo_dist_si, "Distances/topo_dist_simonsi.csv")








##### 7. RIVER DISTANCE --------------------------------------------------------------
### Based on https://www.naturalearthdata.com/downloads/ The The most detailed. Suitable for making zoomed-in maps of countries and regions. 1:10,000,000 / 1″ = 158 miles / 1 cm = 100 km 1:10m physical rivers: "ne_10m_rivers_lake_centerlines"

#A. Loading river shapefile:

#it is strongly recommended to simplify the river shapefile as much as possible before importing into R. In particular, a spatial dissolve will likely be very helpful, if GIS software is available. This will create a few long line segments instead of many, many short segments.
river = readOGR("Shapefiles/rivers_WA_d.shp") #shoud be dissolved and in projected in the case I used the South America Albers Equal Area.
plot(river)
class(river)
crs(river)

#B. South America Albers Equal Area:
albers = "+proj=aea +lat_1=-5 +lat_2=-42 +lat_0=-32 +lon_0=-60 +x_0=0 +y_0=0 +ellps=aust_SA +units=m +no_defs"

#C. Creating a Network based on river shapefile:
wa_steerei = line2network(river)

#D. Displaying
plot(wa_steerei)

#E. Cleaning network
#Recommend saving output to a .Rdata or .rda file.
#There a lot of questions to answering. cleanup () will identify all the segments and attempt to remove duplicate lines and dissolve segments that are split (uncessarily). When it asks whether to "Dissolve" or not, select yes.

#The "Insert vertices" question makes a more evenly distributed network so when we try snap points, there will be more vertices available in the network. The histogram provides a sense of the distribution of the segement sizes. I try to select something near the peak so most segments are a consistent size. I picked 500 meters for this example.

#Please identify segment number of river mouth and the vertix of the mouth.

#The cleanup function will then continue and ask if you want to remove additional segments, and check for braiding. Just walk through and determine what you want to keep or not. If there are braided segment, you’ll need to pick out the components you want to keep and which to delete. Here I remove as many of the braided segments as I can.

#You need to choose how the segment are connected, if they are connected by end of the segment or the begining (considering the mouth you have delimitated)

wa_steerei_fixed = cleanup(wa_steerei)

#F. Saving the cleaned topology a Rdata file:
save(wa_steerei_fixed, file = "Shapefiles/rivers_network_WA_fixed.Rdata")

#G. Topology
topologydots(rivers=wa_steerei_fixed)


#STEEREI
#A. Make coordinates data sf object (spatial): 
st_locs = st_as_sf(steerei, #data frame
                    coords = c("longitude", "latitude"), # for point data
                    remove = F, # don't remove these lat/lon cols from df
                    crs = 4326) # add projection (this is WGS84)

#B. Converting lat/lon coordinates to Albers (planar)
st_locs_converted = st_locs %>% 
  st_transform(crs = "+proj=aea +lat_1=-5 +lat_2=-42 +lat_0=-32 +lon_0=-60 +x_0=0 +y_0=0 +ellps=aust_SA +units=m +no_defs") # convert to Albers Equal Area

#C. Add COORDS in projected form
st_locs_converted$X = st_coordinates(st_locs_converted)[,1]
st_locs_converted$Y = st_coordinates(st_locs_converted)[,2]

#D. Snap points to the closest vertix/line
st_riv_points = xy2segvert(x=st_locs_converted$X, y=st_locs_converted$Y, rivers=wa_steerei_fixed)
head(st_riv_points)

#E. See the distribution of the distances:
hist(st_riv_points$snapdist, main="snapping distance (m)")

#F. Displaying point data in river location anda save it as pdf
pdf("./Network_WA_steerei.pdf", onefile = F)
zoomtoseg(seg=c(1:11), rivers=wa_steerei_fixed)
riverpoints(seg=st_riv_points$seg, vert=st_riv_points$vert, rivers=wa_steerei_fixed, pch=15, col="blue")
dev.off()

#G. Computing a a matrix between all observations
dmat_st = as.matrix(riverdistancemat(st_riv_points$seg, st_riv_points$vert, wa_steerei_fixed))
head(dmat_st)
dmat_st = dmat_st/1000 # convert to km
head(dmat_st)
#row and col names:
rownames(dmat_st) = steerei[,7]
colnames(dmat_st) = steerei[,7]
dmat_st[1:10,1:10]

#H. removing upper diagonal:
dmat_st[upper.tri(dmat_st, diag = T)] = NA
dmat_st

#I. Converting distance matrix into data frame for analyses:
dmat_st = dmat_st %>%
  melt %>%
  na.omit %>%
  arrange(., Var1) %>%
  setNames(c("X1", "X2","riverdistance"))

#J. Verify for NA or 0 values
head(dmat_st)
summary(dmat_st)
length(dmat_st[,3]) #171

#K. replace 0 to 0.001
dmat_st[,3][dmat_st[,3] == 0] = 0.001

#L. verify for NA or 0 values
summary(dmat_st)
length(dmat_st[,3]) #171

#M. Save results
write.csv(dmat_st, "Distances/river_dist_steerei.csv")



#SIMONSI
#A. Make coordinates data sf object (spatial): 
si_locs = st_as_sf(simonsi, #data frame
                   coords = c("longitude", "latitude"), # for point data
                   remove = F, # don't remove these lat/lon cols from df
                   crs = 4326) # add projection (this is WGS84)

#B. Converting lat/lon coordinates to Albers (planar)
si_locs_converted = si_locs %>% 
  st_transform(crs = "+proj=aea +lat_1=-5 +lat_2=-42 +lat_0=-32 +lon_0=-60 +x_0=0 +y_0=0 +ellps=aust_SA +units=m +no_defs") # convert to Albers Equal Area

#C. Add COORDS in projected form
si_locs_converted$X = st_coordinates(si_locs_converted)[,1]
si_locs_converted$Y = st_coordinates(si_locs_converted)[,2]

#D. Snap points to vertix/line
si_riv_points = xy2segvert(x=si_locs_converted$X, y=si_locs_converted$Y, rivers=wa_steerei_fixed)
head(si_riv_points)

#E. See the distribution of the distances:
hist(si_riv_points$snapdist, main="snapping distance (m)")

#F. Displaying point data in river location and save it as pdf
pdf("./Network_WA_simonsi.pdf", onefile = F)
zoomtoseg(seg=c(1:11), rivers=wa_steerei_fixed)
riverpoints(seg=si_riv_points$seg, vert=si_riv_points$vert, rivers=wa_steerei_fixed, pch=15, col="blue")
dev.off()

#G. Computing a a matrix between all observations
dmat_si = as.matrix(riverdistancemat(si_riv_points$seg, si_riv_points$vert, wa_steerei_fixed))
head(dmat_si)
dmat_si = dmat_si/1000 # convert to km
head(dmat_si)
#row and col names:
rownames(dmat_si) = simonsi[,8]
colnames(dmat_si) = simonsi[,8]
dmat_si[1:10,1:10]

#H. Removing upper diagonal:
dmat_si[upper.tri(dmat_si, diag = T)] = NA
dmat_si

#I. Converting distance matrix into data frame for analyses:
dmat_si = dmat_si %>%
  melt %>%
  na.omit %>%
  arrange(., Var1) %>%
  setNames(c("X1", "X2","riverdistance"))

#J. Verify for NA or 0 values
head(dmat_si)
summary(dmat_si)
length(dmat_si[,3]) #136

#K. Replace 0 to 0.001
dmat_si[,3][dmat_si[,3] == 0] = 0.001

#L. verify for NA or 0 values
summary(dmat_si)
length(dmat_si[,3]) #136

#M. Save results
write.csv(dmat_si, "Distances/river_dist_simonsi.csv")




##### 8. HABITAT RESISTANCE -------------------------------------------------
### Based on Tropical and Subtropical Wetlands Distribution version 2
#Gumbricht, T.; Román-Cuesta, R.M.; Verchot, L.V.; Herold, M.; Wittmann, F; Householder, E.; Herold, N.; Murdiyarso, D., 2017, "Tropical and Subtropical Wetlands Distribution version 2", https://doi.org/10.17528/CIFOR/DATA.00058, Center for International Forestry Research (CIFOR), V3, UNF:6:Bc9aFtBpam27aFOCMgW71Q== [fileUNF]

#A. Loading rasters
wetlands = raster("Rasters/Wetlands/TROP-SUBTROP_PeatV21_2016_CIFOR.tif")
plot(wetlands)

#B. Choose a extension for study area to reduce time for analyses:
ext = extent(-80, -60, -20, -0)
wetlands = crop(wetlands, ext)
projection(wetlands) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84")
plot(wetlands)

#D. clip the shoreline
southamerica = shapefile("Shapefiles/americadosul_d.shp") #shoud be dissolved and in projected in the case I used the South America Albers Equal Area.
plot(southamerica)
projection(southamerica) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84")
wetlands = crop(wetlands, southamerica)
plot(wetlands)

#E. Resampling wetlands layers. This step will take a long time.
#load the reduced DEM raster as example for resolution:
elevation = raster("Rasters/DEM/elevation_WA.tif")
projection(elevation) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84", doCheckCRSArgs = F)
#resample to same resolution than elevation data:
wetlands = resample(wetlands, elevation, method="bilinear", bylayer=F, progress='text', snap="out")
wetlands
plot(wetlands)

#E. Save the study area raster:
writeRaster(wetlands, filename="Rasters/Wetlands/Wetlands_WA.asc", format="ascii", overwrite=TRUE, options=c('TFW=YES'), bylayer=F)



###############################################STEEREI (SEASONAL FLOODPLAIN FORESTS):
#A. Load raster:
wetlands_st = raster("Rasters/Wetlands/Wetlands_WA.asc")
ext = extent(-74, -61, -14, -3) #reducing more
wetlands_st = crop(wetlands_st, ext)
wetlands_st
plot(wetlands_st)

#B. Correct original values 0 == dry forest / != 0 is wetlands. Observation: Resample function created different values for wetland, not only 1, we gonna corrected this
wetlands_st[wetlands_st != 0] = 1
wetlands_st
plot(wetlands_st)


#C. Create a data frame with the geographical coordenates. Need to remove the row and col names:
steerei_gdist = as.matrix(steerei[,c(4:5)])
rownames(steerei_gdist) = NULL
colnames(steerei_gdist) = NULL
head(steerei_gdist)


######################################## D. NULL HYPOTHESIS. ALL VALUES ARE 1.
NULL_wetlands_st = wetlands_st
#change the values
NULL_wetlands_st[NULL_wetlands_st == 0] = 1
NULL_wetlands_st[NULL_wetlands_st == 1] = 1
NULL_wetlands_st
plot(NULL_wetlands_st)

#Create a cost distance
NULL_wetlands_st_tr = transition(NULL_wetlands_st, transitionFunction=mean, 16)
NULL_wetlands_st_trC = geoCorrection(NULL_wetlands_st_tr, type="c", multpl=FALSE, scl=TRUE)
wetDIST_st_NULL = costDistance(NULL_wetlands_st_trC, steerei_gdist)

#Computing a matrix between all observations
dmat_st = as.matrix(wetDIST_st_NULL)
head(dmat_st)
#row and col names:
rownames(dmat_st) = steerei[,7]
colnames(dmat_st) = steerei[,7]
dmat_st[1:10,1:10]

#Removing upper diagonal:
dmat_st[upper.tri(dmat_st, diag = T)] = NA
dmat_st

#Converting distance matrix into data frame for analyses:
dmat_st = dmat_st %>%
  melt %>%
  na.omit %>%
  arrange(., Var1) %>%
  setNames(c("X1", "X2","wetlands_NULL"))

#Verify for NA or 0 values
head(dmat_st)
summary(dmat_st)
length(dmat_st[,3]) #171

#Replace 0 to 0.001
dmat_st[,3][dmat_st[,3] == 0] = 0.001

#Verify for NA or 0 values
summary(dmat_st)
length(dmat_st[,3]) #171

#Save results
write.csv(dmat_st, "Distances/wetlands_dist_NULL_steerei.csv")



###################### 0.1 - 0.9 - VERY HIGH RESISTENCE HYPOTHESIS | 1 = WETLANDS!
VH_wetlands_st = wetlands_st

#change values to represent wetlands promoting movements
VH_wetlands_st[VH_wetlands_st == 0] = 0.9
VH_wetlands_st[VH_wetlands_st == 1] = 0.1
VH_wetlands_st
plot(VH_wetlands_st)

#Create a cost-disntance
#transition matrix
VH_wetlands_st_tr = transition(VH_wetlands_st, transitionFunction=mean, 16)

#correct the transition matrix
#The argument scl is set to TRUE to scale the transition values to a reasonable range. If the transition values are too large, commute distance and randomized shortest path functions will not work well. No scaling should be done if the user wants to obtain absolute distance values as output.
VH_wetlands_st_trC = geoCorrection(VH_wetlands_st_tr, type="c", multpl=FALSE, scl=TRUE)

#Need to remove the row and col names
# Create a data frame with the geographical coordenates:
steerei_gdist = as.matrix(steerei[,c(4:5)])
rownames(steerei_gdist) = NULL
colnames(steerei_gdist) = NULL
head(steerei_gdist)

#calculate the cost distance:
wetDIST_st_VH = costDistance(VH_wetlands_st_trC, steerei_gdist)

#Computing a a matrix between all observations
dmat_st = as.matrix(wetDIST_st_VH)
head(dmat_st)
#row and col names:
rownames(dmat_st) = steerei[,7]
colnames(dmat_st) = steerei[,7]
dmat_st[1:10,1:10]

#removing upper diagonal:
dmat_st[upper.tri(dmat_st, diag = T)] = NA
dmat_st

#Converting distance matrix into data frame for analyses:
dmat_st = dmat_st %>%
  melt %>%
  na.omit %>%
  arrange(., Var1) %>%
  setNames(c("X1", "X2","wetlands_VH"))

#Verify for NA or 0 values
head(dmat_st)
summary(dmat_st)
length(dmat_st[,3]) #171

#replace 0 to 0.001
dmat_st[,3][dmat_st[,3] == 0] = 0.001

#verify for NA or 0 values
summary(dmat_st)
length(dmat_st[,3]) #171

#Save results
write.csv(dmat_st, "Distances/wetlands_dist_VH_steerei.csv")



############################################## 0.2 - 0.8 - HIGH RESISTENCE HYPOTHESIS
H_wetlands_st = wetlands_st
#change values to represent wetlands promoting movements
H_wetlands_st[H_wetlands_st == 0] = 0.8
H_wetlands_st[H_wetlands_st == 1] = 0.2
H_wetlands_st
plot(H_wetlands_st)

#Create a cost-distance
H_wetlands_st_tr = transition(H_wetlands_st, transitionFunction=mean, 16)
H_wetlands_st_trC = geoCorrection(H_wetlands_st_tr, type="c", multpl=FALSE, scl=TRUE)
wetDIST_st_H = costDistance(H_wetlands_st_trC, steerei_gdist)

#Computing a a matrix between all observations
dmat_st = as.matrix(wetDIST_st_H)
head(dmat_st)
#row and col names:
rownames(dmat_st) = steerei[,7]
colnames(dmat_st) = steerei[,7]
dmat_st[1:10,1:10]

#removing upper diagonal:
dmat_st[upper.tri(dmat_st, diag = T)] = NA
dmat_st

#Converting distance matrix into data frame for analyses:
dmat_st = dmat_st %>%
  melt %>%
  na.omit %>%
  arrange(., Var1) %>%
  setNames(c("X1", "X2","wetlands_H"))

#Verify for NA or 0 values
head(dmat_st)

#verify for NA or 0 values
summary(dmat_st)
length(dmat_st[,3]) #171

#replace 0 to 0.001
dmat_st[,3][dmat_st[,3] == 0] = 0.001

#verify for NA or 0 values
summary(dmat_st)
length(dmat_st[,3]) #171

#Save results
write.csv(dmat_st, "Distances/wetlands_dist_H_steerei.csv")



###################################### 0.3 - 0.7 - MODERATE RESISTENCE HYPOTHESIS
M_wetlands_st = wetlands_st
#change values to represent wetlands promoting movements
M_wetlands_st[M_wetlands_st == 0] = 0.7
M_wetlands_st[M_wetlands_st == 1] = 0.3
M_wetlands_st
plot(M_wetlands_st)

#Create a cost-distance
M_wetlands_st_tr = transition(M_wetlands_st, transitionFunction=mean, 16)
M_wetlands_st_trC = geoCorrection(M_wetlands_st_tr, type="c", multpl=FALSE, scl=TRUE)
wetDIST_st_M = costDistance(M_wetlands_st_trC, steerei_gdist)

#Computing a a matrix between all observations
dmat_st = as.matrix(wetDIST_st_M)
head(dmat_st)
#row and col names:
rownames(dmat_st) = steerei[,7]
colnames(dmat_st) = steerei[,7]
dmat_st[1:10,1:10]

#removing upper diagonal:
dmat_st[upper.tri(dmat_st, diag = T)] = NA
dmat_st

#Converting distance matrix into data frame for analyses:
dmat_st = dmat_st %>%
  melt %>%
  na.omit %>%
  arrange(., Var1) %>%
  setNames(c("X1", "X2","wetlands_M"))

#Verify for NA or 0 values
head(dmat_st)
summary(dmat_st)
length(dmat_st[,3]) #171

#replace 0 to 0.001
dmat_st[,3][dmat_st[,3] == 0] = 0.001

#verify for NA or 0 values
summary(dmat_st)
length(dmat_st[,3]) #171

#Save results
write.csv(dmat_st, "Distances/wetlands_dist_M_steerei.csv")



############################################### 0.4 - 0.6 - LOW RESISTENCE HYPOTHESIS
L_wetlands_st = wetlands_st
#change values to represent wetlands promoting movements
L_wetlands_st[L_wetlands_st == 0] = 0.6
L_wetlands_st[L_wetlands_st == 1] = 0.4
L_wetlands_st
plot(L_wetlands_st)

#Create a cost-distance
L_wetlands_st_tr = transition(L_wetlands_st, transitionFunction=mean, 16)
L_wetlands_st_trC = geoCorrection(L_wetlands_st_tr, type="c", multpl=FALSE, scl=TRUE)
wetDIST_st_L = costDistance(L_wetlands_st_trC, steerei_gdist)

#Computing a a matrix between all observations
dmat_st = as.matrix(wetDIST_st_L)
head(dmat_st)
#row and col names:
rownames(dmat_st) = steerei[,7]
colnames(dmat_st) = steerei[,7]
dmat_st[1:10,1:10]

#removing upper diagonal:
dmat_st[upper.tri(dmat_st, diag = T)] = NA
dmat_st

#Converting distance matrix into data frame for analyses:
dmat_st = dmat_st %>%
  melt %>%
  na.omit %>%
  arrange(., Var1) %>%
  setNames(c("X1", "X2","wetlands_L"))

#Verify for NA or 0 values
head(dmat_st)
summary(dmat_st)
length(dmat_st[,3]) #171

#replace 0 to 0.001
dmat_st[,3][dmat_st[,3] == 0] = 0.001

#verify for NA or 0 values
summary(dmat_st)
length(dmat_st[,3]) #171

#Save results
write.csv(dmat_st, "Distances/wetlands_dist_L_steerei.csv")





###################################################### SIMONSI (NON-FLOODED FOREST):
#A. Load raster:
wetlands_si = raster("Rasters/Wetlands/Wetlands_WA.asc")
ext = extent(-74, -61, -14, -3) #reducing more
wetlands_si = crop(wetlands_si, ext)
wetlands_si
plot(wetlands_si)

#B. Correcting values, original 0 == dry forest / != 0 is wetlands
wetlands_si[wetlands_si != 0] = 1
wetlands_si
plot(wetlands_si)


#C. Create a data frame with the geographical coordenates. Need to remove the row and col names:
simonsi_gdist = as.matrix(simonsi[,c(5:6)])
rownames(simonsi_gdist) = NULL
colnames(simonsi_gdist) = NULL
head(simonsi_gdist)



#################################################### D. NULL HYPOTHESIS. ALL VALUES ARE 1
NULL_wetlands_si = wetlands_si
#change the values
NULL_wetlands_si[NULL_wetlands_si == 0] = 1
NULL_wetlands_si[NULL_wetlands_si == 1] = 1
NULL_wetlands_si
plot(NULL_wetlands_si)

#Create a distance matrix:
NULL_wetlands_si_tr = transition(NULL_wetlands_si, transitionFunction=mean, 16)
NULL_wetlands_si_trC = geoCorrection(NULL_wetlands_si_tr, type="c", multpl=FALSE, scl=TRUE)
wetDIST_si_NULL = costDistance(NULL_wetlands_si_trC, simonsi_gdist)

#Computing a matrix between all observations
dmat_si = as.matrix(wetDIST_si_NULL)
head(dmat_si)
#row and col names:
rownames(dmat_si) = simonsi[,8]
colnames(dmat_si) = simonsi[,8]
dmat_si[1:10,1:10]

#Removing upper diagonal:
dmat_si[upper.tri(dmat_si, diag = T)] = NA
dmat_si

#Converting distance matrix into data frame for analyses:
dmat_si = dmat_si %>%
  melt %>%
  na.omit %>%
  arrange(., Var1) %>%
  setNames(c("X1", "X2","wetlands_NULL"))

#Verify for NA or 0 values
head(dmat_si)
summary(dmat_si)
length(dmat_si[,3]) #136

#replace 0 to 0.001
dmat_si[,3][dmat_si[,3] == 0] = 0.001

#verify for NA or 0 values
summary(dmat_si)
length(dmat_si[,3]) #136

#Save results
write.csv(dmat_si, "Distances/wetlands_dist_NULL_simonsi.csv")


########################################## 0.1 - 0.9 - VERY HIGH RESISTENCE | 1 = wetland
VH_wetlands_si = wetlands_si
#change values to represent wetlands preventing movements
VH_wetlands_si[VH_wetlands_si == 0] = 0.1
VH_wetlands_si[VH_wetlands_si == 1] = 0.9
VH_wetlands_si
plot(VH_wetlands_si)

#Create a cost-distance
VH_wetlands_si_tr = transition(VH_wetlands_si, transitionFunction=mean, 16)
VH_wetlands_si_trC = geoCorrection(VH_wetlands_si_tr, type="c", multpl=FALSE, scl=TRUE)
wetDIST_si_VH = costDistance(VH_wetlands_si_trC, simonsi_gdist)

#Computing a a matrix between all observations
dmat_si = as.matrix(wetDIST_si_VH)
head(dmat_si)
#row and col names:
rownames(dmat_si) = simonsi[,8]
colnames(dmat_si) = simonsi[,8]
dmat_si[1:10,1:10]

#removing upper diagonal:
dmat_si[upper.tri(dmat_si, diag = T)] = NA
dmat_si

#Converting distance matrix into data frame for analyses:
dmat_si = dmat_si %>%
  melt %>%
  na.omit %>%
  arrange(., Var1) %>%
  setNames(c("X1", "X2","wetlands_VH"))

#Verify for NA or 0 values
head(dmat_si)
summary(dmat_si)
length(dmat_si[,3]) #136

#replace 0 to 0.001
dmat_si[,3][dmat_si[,3] == 0] = 0.001

#verify for NA or 0 values
summary(dmat_si)
length(dmat_si[,3]) #136

#Save results
write.csv(dmat_si, "Distances/wetlands_dist_VH_simonsi.csv")



################################################# 0.2 - 0.8 - HIGH RESISTENCE | 1 = wetland
H_wetlands_si = wetlands_si
#change values to represent wetlands preventing movements
H_wetlands_si[H_wetlands_si == 0] = 0.2
H_wetlands_si[H_wetlands_si == 1] = 0.8
H_wetlands_si
plot(H_wetlands_si)

#Create a cost-distance
H_wetlands_si_tr = transition(H_wetlands_si, transitionFunction=mean, 16)
H_wetlands_si_trC = geoCorrection(H_wetlands_si_tr, type="c", multpl=FALSE, scl=TRUE)
wetDIST_si_H = costDistance(H_wetlands_si_trC, simonsi_gdist)

#Computing a a matrix between all observations
dmat_si = as.matrix(wetDIST_si_H)
head(dmat_si)
#row and col names:
rownames(dmat_si) = simonsi[,8]
colnames(dmat_si) = simonsi[,8]
dmat_si[1:10,1:10]

#removing upper diagonal:
dmat_si[upper.tri(dmat_si, diag = T)] = NA
dmat_si

#Converting distance matrix into data frame for analyses:
dmat_si = dmat_si %>%
  melt %>%
  na.omit %>%
  arrange(., Var1) %>%
  setNames(c("X1", "X2","wetlands_H"))

#Verify for NA or 0 values
head(dmat_si)
summary(dmat_si)
length(dmat_si[,3]) #136

#replace 0 to 0.001
dmat_si[,3][dmat_si[,3] == 0] = 0.001

#verify for NA or 0 values
summary(dmat_si)
length(dmat_si[,3]) #136

#Save results
write.csv(dmat_si, "Distances/wetlands_dist_H_simonsi.csv")




######################################### 0.3 - 0.7 - MODERATE RESISTENCE | 1 = wetland
M_wetlands_si = wetlands_si
#change values to represent wetlands preventing movements
M_wetlands_si[M_wetlands_si == 0] = 0.3
M_wetlands_si[M_wetlands_si == 1] = 0.7
M_wetlands_si
plot(M_wetlands_si)

#Create a cost-distance
M_wetlands_si_tr = transition(M_wetlands_si, transitionFunction=mean, 16)
M_wetlands_si_trC = geoCorrection(M_wetlands_si_tr, type="c", multpl=FALSE, scl=TRUE)
wetDIST_si_M = costDistance(M_wetlands_si_trC, simonsi_gdist)

#Computing a a matrix between all observations
dmat_si = as.matrix(wetDIST_si_M)
head(dmat_si)
#row and col names:
rownames(dmat_si) = simonsi[,8]
colnames(dmat_si) = simonsi[,8]
dmat_si[1:10,1:10]

#removing upper diagonal:
dmat_si[upper.tri(dmat_si, diag = T)] = NA
dmat_si

#Converting distance matrix into data frame for analyses:
dmat_si = dmat_si %>%
  melt %>%
  na.omit %>%
  arrange(., Var1) %>%
  setNames(c("X1", "X2","wetlands_M"))

#Verify for NA or 0 values
head(dmat_si)
summary(dmat_si)
length(dmat_si[,3]) #136

#replace 0 to 0.001
dmat_si[,3][dmat_si[,3] == 0] = 0.001

#verify for NA or 0 values
summary(dmat_si)
length(dmat_si[,3]) #171

#Save results
write.csv(dmat_si, "Distances/wetlands_dist_M_simonsi.csv")



################################################ 0.4 - 0.6 - LOW RESISTENCE | 1 = wetland
L_wetlands_si = wetlands_si
#change values to represent wetlands preventing movements
L_wetlands_si[L_wetlands_si == 0] = 0.4
L_wetlands_si[L_wetlands_si == 1] = 0.6
L_wetlands_si
plot(L_wetlands_si)

#Create a cost-distance
L_wetlands_si_tr = transition(L_wetlands_si, transitionFunction=mean, 16)
L_wetlands_si_trC = geoCorrection(L_wetlands_si_tr, type="c", multpl=FALSE, scl=TRUE)
wetDIST_si_L = costDistance(L_wetlands_si_trC, simonsi_gdist)

#Computing a a matrix between all observations
dmat_si = as.matrix(wetDIST_si_L)
head(dmat_si)
#row and col names:
rownames(dmat_si) = simonsi[,8]
colnames(dmat_si) = simonsi[,8]
dmat_si[1:10,1:10]

#removing upper diagonal:
dmat_si[upper.tri(dmat_si, diag = T)] = NA
dmat_si

#Converting distance matrix into data frame for analyses:
dmat_si = dmat_si %>%
  melt %>%
  na.omit %>%
  arrange(., Var1) %>%
  setNames(c("X1", "X2","wetlands_L"))

#Verify for NA or 0 values
head(dmat_si)
summary(dmat_si)
length(dmat_si[,3]) #136

#replace 0 to 0.001
dmat_si[,3][dmat_si[,3] == 0] = 0.001

#verify for NA or 0 values
summary(dmat_si)
length(dmat_si[,3]) #136

#Save results
write.csv(dmat_si, "Distances/wetlands_dist_L_simonsi.csv")


##END
