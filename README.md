### LANDSCAPE GENOMICS PIPELINE IN R ###
### PART 3 - ISOLATION BY RESISTANCE USING MLPE MODELS ###

Scripts in 3 steps to run Isolation-by-Resistance (IBR) using Maximum‐likelihood population‐effects (MLPE) mixed models:

  STEP 1: Creating distances matrices to represent resistance values based on:
  1. Euclidian Distance (Geographic Distance)
  2. Topographic Distance in TopoDistance R package
  3. RiverNetwork Distance in riverdist R package
  4. Productivity Distance based on Species Distribution Models (SDM) based on PET, Temperature, and Precipitation. To see the scripts on SDM: https://github.com/jdalapicolla/SDM_biomod2
  5. Habitat Distance based on a categorial raster LandUse/Biomes/etc.
 
 STEP 2: Organizing and choose Distances Matrices
 
 STEP 3: Running Mixed models to estimate Isolation-by-Resistance (IBR) between individuals using lme()
 
 STEP 4: Running Mixed models to estimate Isolation-by-Resistance (IBR) between individuals using gls() and Nested MLPE
 
 
## PART 1 - GENETIC STRUCTURE AND GENETIC DIVERSITY https://github.com/jdalapicolla/LanGen_pipeline_version2
## PART 2 - ISOLATION BY DISTANCE AND FINE-SCALE SPATIAL GENETIC STRUCTURE https://github.com/jdalapicolla/IBD_models.R
## PART 3 - ISOLATION BY RESISTANCE USING MLPE MODELS - https://github.com/jdalapicolla/MLPE.R
## PART 4 - LOCAL ADAPTATION ANALYSES - IDENTIFICATION OF CANDIDATES LOCI UNDER SELECTION - https://github.com/jdalapicolla/LOCAL_ADAPTATION.R
