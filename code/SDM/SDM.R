# Script fits a species distribution model to occurence data and environmental 
# covariates, saving the resulting SDM to a tif for use in the spillover model

# Setup ----
library(raster)
library(rJava)
library(maxnet)
library(dismo)
library(rgdal)
library(plyr)
library(dplyr)
library(magrittr)
library(pdp)

# Read in the occurrence data, and clean it if necessary ----
if(!file.exists("./data/cleaned/SDM/species_occurrences.csv")){
  print("missing cleaned occurrences data. cleaning occurences")
  source("./code/SDM/clean_species_occurrences.R")
}
occ.data = read.csv("./data/cleaned/SDM/species_occurrences.csv", 
                    header = T, stringsAsFactors = F)

# Identify the unique GPS locations with occurence data
occ.GPS <- dplyr::select(occ.data, c(decimallongitude, decimallatitude)) %>% 
  unique

# Read in environmental data ----
predictors <- alply(list.files("./data/cleaned/SDM/env_covariates", 
                               pattern = ".tif",
                               full.names = T), 1, function(file){
  rast <- raster(file)
  rast[] %<>% plyr::mapvalues(from = c(Inf, -Inf, NaN), to = c(NA, NA, NA))
  return(rast)
}) 
predictors %<>% stack()

# Calculate latitude and add it to the environmental predictors
lat_raster = raster::init(x = predictors@layers[[1]], fun = 'y')
lat_raster[] %<>% abs()
predictors = stack(predictors, lat_raster)
names(predictors) = c(gsub(".tif", "", 
                           list.files("./data/cleaned/SDM/env_covariates", 
                                      pattern = ".tif")),"AbsLatitude")
plot(predictors)

# Construct biased background points ----
# Read in all mosquito data from GBIF and exclude species of interest 
all_mosq <- read.csv(
  "./data/raw/SDM/species_occurrence_data/all_mosquitoes_SA.csv", 
  sep = "\t", header = T, stringsAsFactors = F)
all_mosq %<>% base::subset(!(species %in% c("Haemagogus janthinomys", 
                                        "Haemagogus leucocelaenus", 
                                        "Sabethes chloropterus")))

# construct a masked raster for the biases, and make each cell value the number of occurences in that cell
mosq_bias <- sum(predictors)
mosq_bias[!is.na(mosq_bias)] <- 0
pixels <- cellFromXY(mosq_bias, 
                     cbind(all_mosq$decimallongitude, 
                           all_mosq$decimallatitude)) %>% table
mosq_bias[as.numeric(names(pixels))] <- mosq_bias[as.numeric(names(pixels))] + 
  as.vector(pixels)

# Sample points from the biased background
a = 1000 # number of points to sample
bg = randomPoints(mosq_bias, n = a, p = occ.GPS, excludep = T, prob = T)

# Run the SDM ----
# Remove occurence points that have NAs for any environmental covariates
occ.GPS <- occ.GPS[-c(which(is.na(raster::extract(mosq_bias, occ.GPS)))),] 
# Identify the raster cells that the GPS points fall in and get the unique ones
occ.GPS.unique <- xyFromCell(predictors, 
                             cellFromXY(predictors, occ.GPS) %>% 
                               unique) %>% as.data.frame()

# extract covariates for background points and occurence points
predictors_df <- rbind(data.frame(raster::extract(predictors, occ.GPS.unique)), 
                        data.frame(raster::extract(predictors, bg))) 
predictors_df %<>% mutate(LandCover = as.factor(LandCover))

# fit a maxnet model, calibrating for regularization parameters and feature classes
source("code/SDM/SDM_functions.R")
maxnet_fit <- trainMaxNet(data=cbind(c(rep(1, nrow(occ.GPS.unique)), 
                                       rep(0, a)), predictors_df),
                          classes = "lqhp",
                          verbose=TRUE,
                          out = c("model", "tuning"))

# Plot responses ----
# Set up a dictionary to go from var.name to clean name for plotting
var.names.plotting = data.frame(
  var.names = c("AnnualPrecip", "DriestMonthPrecip", "Elev", "ForestCover", 
                "LandCover", "MaxAnnualLST", "MedAnnualLST", "MedianAnnualEVI",   
                "MinAnnualLST", "WettestMonthPrecip", "AbsLatitude"),
  plotting.names = c("Total Annual Precipitation", "Driest Month Precipitation", 
                     "Elevation", "Forest Cover (%)", "Land Cover", 
                     "Maximum Annual LST", "Median Annual LST", 
                     "Median Annual EVI", "Minimum Annual LST",
                     "Wettest Month Precipitation", "Absolute Latitude"),
  stringsAsFactors = F)

# j tracks the jth most important variable that we are looking at
pdf("./output/submodels/SDM_covariate_responses.pdf", width = 8, height = 10)
par(mfrow = c(4,3))
par(mar=c(5, 6, 1.5, 2.5) + 0.1)
par(oma = c(0, 0, 0, 4))
par(cex.axis = 1.2, cex.lab = 1.2)
for (j in 1:11){
  if( j %% par('mfrow')[2] == 1){left_ax = T
  } else{left_ax = F}
  if(j %% par('mfrow')[2] == 0){right_ax = T
  } else{right_ax = F}
  var.name = var.names.plotting[j,1]
  print(var.name)
  # Figure out which of the predictor names it was
  var.data = predictors_df[[var.name]]
  test_resp <- pdp::partial(maxnet_fit$model, pred.var = var.name,
                            quantiles = T, probs = c(1:99/100),
                            plot=F, 
                            train = cbind(c(rep(1, nrow(occ.GPS.unique)), 
                                                    rep(0, a)), predictors_df))
  if(is.factor(var.data)) {
    test_sums = table(var.data)
    test_hist <- graphics::barplot(test_sums, col = "grey", border = "white",
                                   las = 1, ylab = "",
                                   xlab = "")
    title(xlab = paste0(letters[j], ". ", 
                        mapvalues(var.name, 
                                  from = var.names.plotting$var.names, 
                                  to = var.names.plotting$plotting.names,
                                  warn_missing = F)),
          line = 3)
    title(ylab = ifelse(left_ax, "Density", ""),
          line = 4)
    par(new = TRUE)
    plot(test_hist, 
         test_resp$yhat,
         xlab = "", ylab = "", axes = FALSE, bty = "n",
         type = "l")
    axis(side=4, at = pretty(range(test_resp$yhat)), las = 1)
    mtext(ifelse(right_ax, "Marginal Effect", ""), side=4, line=4, cex = 0.75)
    box()
  } else {
    resp_x = test_resp[,1]
    test_hist <- graphics::hist(var.data,
                                plot = F)
    plot(test_hist, freq = F, col = "grey", border = "white", 
         las = 1, 
         ylab = "", xlab = "",
         main = NULL)
    title(ylab = ifelse(left_ax, "Density", ""),
          line = 4)
    title(xlab = paste0(letters[j], ". ", 
                        mapvalues(var.name, 
                                  from = var.names.plotting$var.names, 
                                  to = var.names.plotting$plotting.names,
                                  warn_missing = F)),
          line = 3)
    par(new = TRUE)
    plot(resp_x,
         test_resp$yhat,
         xlab = "", ylab = "", axes = FALSE, bty = "n",
         type = "l")
    axis(side=4, at = pretty(range(test_resp$yhat)), las = 1)
    mtext(ifelse(right_ax, "Marginal Effect", ""), side=4, line=4, cex = 0.75)
    box()
  }
}
dev.off()

# Predict over South America ----
maxnet_predict <- raster::predict(object = predictors, 
                                  model = maxnet_fit$model, type = "cloglog")
par(mfrow = c(1,1))
plot(maxnet_predict)
writeRaster(maxnet_predict, "./output/submodels/vector_distribution_SDM.tif", 
             format = "GTiff", overwrite = T)


# Calcuate AUC ----
library(pROC)
roc(c(rep(1, nrow(occ.GPS.unique)), rep(0, a)), as.vector(predict(maxnet_fit$model,predictors_df))) %>% auc

# Plot predictions ----
pdf("./output/submodels/SDM_data_comparison.pdf")
par(mfrow = c(1,1))
plot(maxnet_predict, axes = F, box = F, legend = T)
points(occ.GPS, pch = 16, cex = 0.75, col = "black")
dev.off()


