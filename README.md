# Isolation by Resistance using MLPE.R

Script in 3 steps to run Isolation-by-Resistance (IBR) using Maximum‐likelihood population‐effects (MLPE) mixed models:

  STEP 1: Creating Distances Matrices using:
      *Euclidian Distance (Geographic Distance)
      *Topographic Distance in TopoDistance R package
      *RiverNetwork Distance in riverdist R package
      *Environmental Distances based on PCA-distance of Continuous variables as WorldClim-Bioclimatic
      *Habitat Distance based on a Categorial raster LandUse/Biomes/etc.
 
 STEP 2: Organizing and choose Distances Matrices
 
 STEP 3: Running Mixed models to estimate Isolation-by-Resistance (IBR) between individuals
