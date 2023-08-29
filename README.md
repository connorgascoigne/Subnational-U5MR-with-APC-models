# Subnational Under-Five Mortality Rates with Age-Period-Cohort models

## Summary

This document provides the code to run the analysis for the manuscript Estimating subnational Under-Five Mortality Rates (U5MR) using a spatio-temporal Age-Period-Cohort model. Within the repository, there are seperate folders for the [Code](https://github.com/connorgascoigne/Subnational-U5MR-with-APC-models/tree/main/Code) and additional [Information on Kenya](https://github.com/connorgascoigne/Subnational-U5MR-with-APC-models/tree/main/Info) used within the code.

## Data 

- Kenyan U5MR data: we cannot share the 2014 Kenyan Demographic and Health Surveys (DHS) data. To access this, the user will need to sign up to the [DHS](https://dhsprogram.com/) with a suitable project that allows access to the Kenyan dataset specifically.
- Kenyan spatial polygons: these can be downloaded from the [Database of Global Administrative Areas (GADM)](https://gadm.org/)
- Population totals: the population totals for the proportional aggregation were calculated using the [World Population](https://www.worldpop.org/)

## [Code](https://github.com/connorgascoigne/Subnational-U5MR-with-APC-models/tree/main/Code)

0. Pre-analysis before running code for analysis 
  - Run the `createDirectoryStructure.R` folder to create all the folders where data will be stored
  - Download the shapefiles from GADM
  - Get access to the Kenyan DHS after registration with the DHS
1. Run `dataProcessing.R` to download the 2014 KDHS and organise it into the correct format
2. Defining the proportions for aggregation:
  - Run `adminWeights.R` to download the yearly-regional proportions from [worldpop](https://www.worldpop.org/)
  - Find the most recent census for Kenya and write a `.csv` file with the most recent urban and rural fractions for each region 
  - Run `urProportions.R` to assign the correct names to the urban/rural fractions in the `.csv` file
  - Run `urThreshold.R` to create the yearly-regional-urban/rural proportions using the yearly-regional proportions and the urban/rural fractions as described by [Wu and Wakefield](https://arxiv.org/abs/2209.10619) 
3. Run `dataExploration.R` to generate plots of the cluster locations and population proportions
4. Generate estimates:
  - Run `summerEstimates.R` to generate and save the direct and Fay-Heriot estimates for Kenyan U5MR
  - Run `apcEsimates.R` to generate and save the Age-Period, Age-Cohort and Age-Period-Cohort subnational model estimates
  - Run `crossValidation.R` to perform the cross-validation for the Age-Period, Age-Cohort, and Age-Period-Cohort models
5. Run `comparisonPlots.R` to make any plots and tables seen in the manuscript and supplementary material
