# This code splits the data for the boosted regression tree into test and training data
# and save the two datasets as r data files (.rds)

# Load the libraries and read in the BRT data
library(plyr)
library(dplyr)
library(raster)
library(magrittr)
library(pracma)
library(ggplot2)
library(stringr)
library(rgeos)
library(rgdal)
library(BalancedSampling)
brt_subset = readRDS("./data/cleaned/BRT_data/brt_data.rds")

# Split data into training and test data, so as to spatially and temporally -----
# balance the training and test data for spillover and nonspillover observations 

# Read in the shapefile of brazil municipalities to get the municipality centroids
shape = readOGR(dsn = 'data/raw/brazil_border_shapefiles/br_municipios')
cent = data.frame(shape@data$CD_GEOCMU, gCentroid(shape, byid = T)@coords)
colnames(cent) = c("code", "x", "y")
cent %<>% mutate(code = as.character(code))

# For spillover and non-spillover observations, join the data with the centroids
# calculated above, set the seed, and then make balanced subsets of 20% and 80% 
# for testing and training
training_test <- dlply(brt_subset, .(spillover), function(data_sub){
  all_sub <- left_join(data_sub, cent, by = c("CD_GEOCMU" = "code")) 
  test_size = ceiling(0.2*nrow(all_sub))
  test_p = rep(test_size/nrow(all_sub), nrow(all_sub))
  set.seed(23049172)
  test_inds = lcube(test_p, 
                    as.matrix(dplyr::select(all_sub, c(seq_ord, x, y))),
                    cbind(test_p))
  test_sub = all_sub[test_inds,]
  train_sub = all_sub[-test_inds,]
  return(list(test_sub, train_sub))
})

# Combine the training and testing data
test = rbind(training_test[[1]][[1]], training_test[[2]][[1]])
train = rbind(training_test[[1]][[2]], training_test[[2]][[2]])

# Visualize balancing
# plot(train$x, train$y)
# points(test$x, test$y, col = "red")
# plot(train$x, train$seq_ord)
# points(test$x, test$seq_ord, col = "red", pch = 19)
# plot(train$y, train$seq_ord)
# points(test$y, test$seq_ord, col = "red")

# Removed the centroids
test %<>% dplyr::select(-c(x, y))
train %<>% dplyr::select(-c(x, y))

# Save the test and training as .rds  files
saveRDS(test, "data/cleaned/BRT_data/brt_test.rds")
saveRDS(train, "data/cleaned/BRT_data/brt_train.rds")
