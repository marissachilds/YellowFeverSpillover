# YellowFeverSpillover
These are the scripts and data used for analysis in the manuscript "Mosquito and primate ecology predict human risk of yellow fever virus spillover in Brazil" by Marissa L. Childs*, Nicole Nova, Justine Colvin, Erin A. Mordecai. The manuscript in its pre-print form is available on bioRxiv (https://doi.org/10.1101/523704). The manuscript describes the mechanistic model of yellow fever spillover parameterized from published data and the findings from the analysis.

*marissac at stanford dot edu

## File organization 

The files are organized into code, output, and data directories. Most analyses rely on code or scripts in the "code" directory, and data either in the "data/raw" or data/cleaned" directory. The code/script will usually save its output either as an image, R data file, or csv in the "output" directory. Note that for each R script included, the working directory should be set to this top directory (i.e. to "YellowFeverSpillover"). 

In some cases, we have not included files due to data use restrictions or size constraints. These include the MODIS MOD44B data for all of Brazil from 2001 - 2016, the IUCN shapefiles used for primate ranges, and the GBIF species occurrences. 

Below are the doi's for the GBIF data and the locations to save these datasets. Once the datasets are downloaded and saved to the appropriate locations, run the script "code/SDM/clean_species_occurrences.R" to generate the cleaned dataset of species occurrences used in the SDM. This should be done before step 2 in the workflow below.  

https://doi.org/10.15468/dl.wvs9g2 saved as "data/raw/SDM/all_mosquitoes_SA.csv"

https://doi.org/10.15468/dl.ozsvnj saved as "data/raw/SDM/Hg.leucoceleanus.gbif.csv"

https://doi.org/10.15468/dl.gxbxtq saved as "data/raw/SDM/Hg.janthinomys.gbif.csv"

https://doi.org/10.15468/dl.1uo4ty saved as "data/raw/SDM/Sa.chloropterus.gbif.csv"

## Workflow

The following is the workflow for the analysis in this manuscript. Note that R refers to the R programming language and GEE refers to Google Earth Engine. 

1. Unzip the shapefiles in the data/raw directory. They are in the Primates directory and brazil_border_shapefiles directory. 
2. Run submodels of spillover components. These include the dispersal, EIP (infectiousness), species distribution model (SDM), seasonality, survival, and phenomenological primate dynamics models. Code for each of these submodels can be found in the "code" directory   
3. Input the parameters from the submodels into the GEE code (lines 89-105 of "code/GoogleEarthEngine/spillover_model") and upload the SDM map to a GEE asset. 
4. Run the GEE script to extract monthly average temperatures for each pixel ("code/GoogleEarthEngine/get_monthly_data") and save the output as a GEE asset. 
5. Run the spillover model for environmental, immunological, population-scaled and periodic risk by changing lines 203, 268, and 270 in the "code/GoogleEarthEngine/spillover_model" script. The risk estimates are saved as GEE image assets.
6. Extract the municipality maximum risk and mean risk for each risk metric in GEE using the script "code/GoogleEarthEngine/extract_mean_max_model_est". Extracting for all 4 risk metrics will require cycling through the different GEE image assets that are the risk estimates. They should be placed in the directory "output/municipality_model_estimates".
7. Also extract data for the BRT for each municipality-month using all of the scripts in "code/GoogleEarthEngine/BRT" directory. The extracted data should then be placed in "data/cleaned/BRT_data".
8. Now in R, combine the model estimates and case data using the script "code/other/combine_model_case.R". This will save an R data file "all_data.rds" in the directory "data/cleaned". It will also save a subset of the data that is from 2001-2016 to "data/cleaned/BRT_data/brt_data.rds".
9. The BRT data should then be split into test and training data ensuring spatial and temporal balance using the script "code/other/split_brt_data.R". The test and training data will be saved to the directory "data/cleaned/BRT_data".
10. Next run the boosted regression trees with the script "code/other/run_BRT.R" and fit to the training data. Currently this script runs the boosted regression tree for the optimal tree complexity (10) and learning rate (0.001) but changing the variables tc and lr in the script will adjust these parameters. The resulting BRT will be saved as an R data file in the directory "output/BRT_fits". Note that for size purposes, we have only included the best fitting BRT in this GitHub respository.
11. Now compare the BRT fits to find the one with the lowest residual deviance using "code/other/BRT_compare". 
12. Finally, run the script "code/other/model_eval" to calculate and produce results, tables, and figures 3-6 in the manuscript.  
