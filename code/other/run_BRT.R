# This code runs a boosted regression tree for the tree complexity (tc) 
# and learning rate (lr) specified. The boosted regression tree is fit by
# selecting optimal number of trees based on cross validation residual deviance
# with gbm.step, and the resulting gbm model is saved as an r data file (.rds)

# Load packages ----
library(plyr)
library(dplyr)
library(magrittr)
library(dismo)
library(gbm)

# Read in the training data
df_brt = readRDS("./data/BRT_data/brt_train.rds")

# Clean up the data to make sure its all in order
df_brt$region = as.factor(df_brt$region)

# Identify which columns of the data frame have covariates and response 
covars = which(
  colnames(df_brt) %in% c("max_EnvRisk", "month", "vax", "fireArea",
                          "fire_pct", "primates_mean", "primates_max", 
                          "Tair_mean", "precipitation", "pop_density","phenom",
			  "region", "lagged_max_EnvRisk", "lagged_fireArea")
)
spill = which(colnames(df_brt) == "spillover")

# Set the max number of trees, tree comlexity, and learning rate
ntrees_max = 5000
tc = 10
lr = 0.001

# Set the seed for repeatability
set.seed(50348561)

# run the brt, selecting optimal number of trees and saving
boost_tree <- gbm.step(data=df_brt, 
                         gbm.x = covars, gbm.y = spill, 
                         max.trees = ntrees_max, 
                         tolerance.method = "auto", n.trees = 500,
                         tolerance = 0.001, step.size = 500,
                         family = "bernoulli", tree.complexity = tc,
                         learning.rate = lr, bag.fraction = 0.5)
saveRDS(boost_tree, file = paste0("./output/BRT_fits/spillover_brt_tc", tc, "_lr", lr, ".rds"))
