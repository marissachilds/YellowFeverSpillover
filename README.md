# YellowFeverSpillover
These are the scripts and data used for analysis in the manuscript "Mosquito and primate ecology predict human risk of yellow fever virus spillover in Brazil" by Marissa L. Childs*, Nicole Nova, Justine Colvin, Erin A. Mordecai. The manuscript in its pre-print form is available on bioRxiv (https://doi.org/10.1101/523704). The manuscript describes the mechanistic model of yellow fever spillover parameterized from published data and the findings from the analysis.

*marissac at stanford dot edu

## File organization 

only include the best BRT for space

For every R script, set the working directory to the top directory of this (i.e. to "YellowFeverSpillover")

## Workflow

The following is the workflow for the analysis in this manuscript. Note that R refers to the R programming language and GEE refers to Google Earth Engine. 

1. Unzip the shapefiles in the data/raw directory. They are in the Primates directory and brazil_border_shapefiles directory. 
2. Run submodels of spillover components. 

input the parameters and SDM map to google earth engine/google earth engine code
GEE: get monthly max temperatures
GEE: run spillover model  for env, imm env, pop, and phenom env risk by changing XYZ in the script
GEE: extract the muni maxes and means 
GEE: extract other BRT data 
R: combine model and case data 
R: split BRT data into test and training
R: run BRTs
R: compare BRTs and calculate AUC of the best with test data
for size, we've only included the best brt result, others available on request
R: run model evals, making figures 3-6

for all code, set wd to YellowFeverSpillover 

include link to GEE