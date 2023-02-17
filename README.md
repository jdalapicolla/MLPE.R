# Part 3 - Isolation by Resistance using MLPE Models

This suit of scripts were developed based on [LanGen_pipeline](https://github.com/rojaff/LanGen_pipeline)

For a similar pipeline without using *VCFtools v.0.1.16*, check our [Populational Genomics Pipeline](https://github.com/jdalapicolla/PopGenPipe)

&nbsp;

## Scripts:
 __3.1_DISTANCE_MATRICES.R__: Creating distances matrices to represent resistance values based on:
  1. Euclidian Distance (Geographic Distance)
  2. Topographic Distance in TopoDistance R package
  3. RiverNetwork Distance in riverdist R package
  4. Productivity Distance based on Species Distribution Models (SDM). [To see the scripts on SDM](https://github.com/jdalapicolla/SDM_biomod2)
  5. Habitat Distance based on a categorial raster LandUse/Biomes/etc
  
 __3.2_SELECTING_MATRICES.R__: Organizing and choose Distances Matrices
 
 __3.3_RUNNING_MLPE_LME.R__: Running Mixed models to estimate Isolation-by-Resistance (IBR) between individuals using `lme()`
 
 __3.4_RUNNING_NMLPE_GLS.R__: Running Mixed models to estimate Isolation-by-Resistance (IBR) between individuals using `gls()` and Nested MLPE
 

&nbsp;

[INITIAL PAGE](https://github.com/jdalapicolla/LanGen_pipeline_version2)

[PART 1 - GENETIC STRUCTURE AND GENETIC DIVERSITY](https://github.com/jdalapicolla/LanGen_pipeline_version2/tree/master/PART1)

[PART 2 - ISOLATION BY DISTANCE AND FINE-SCALE SPATIAL GENETIC STRUCTURE](https://github.com/jdalapicolla/IBD_models.R)

[PART 4 - LOCAL ADAPTATION ANALYSES - IDENTIFICATION OF CANDIDATES LOCI UNDER SELECTION](https://github.com/jdalapicolla/LOCAL_ADAPTATION.R)
