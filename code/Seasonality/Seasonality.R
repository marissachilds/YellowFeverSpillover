# Here, we predict the relationship between rainfall and relative  mosquito 
# abundance to get a seasonal pattern of mosquitoes for the spillover model

# Load libraries
library(plyr)
library(magrittr)
library(ggplot2)

# Read in the seasonality data and get the needed columns
season_data_orig = read.csv(
  "./data/raw/Seasonality/seasonality_data_extraction.csv", 
  stringsAsFactors = F)
season_data = season_data_orig[c("NumberOfMosquitoes", "Rainfall", 
                                 "SequenceNumber", "Month", 
                                 "MosquitoSpecies", "Source", "Location")]

# Remove the data from sequences which have non adult mosquito data (1,2),
# relative humidity (4,5,6) or non-consecutive months (8, 9, 13, 14, 15)
season_data = season_data[!(season_data$SequenceNumber %in% 
                              c(1,2,4,5,6, 8, 9, 13, 14, 15)),]

# Operating one sequence at a time, convert mosquitoes and rainfall to 
# percentages for relative abundance and also add lagged rainfall
scaled_season_data <- alply(unique(season_data$SequenceNumber), 1, function(i){
  seq_data = season_data[season_data$SequenceNumber == i,]
  pct_mosquitoes = seq_data$NumberOfMosquitoes/max(seq_data$NumberOfMosquitoes)
  lagged_pct_mosquitoes = c(NA, pct_mosquitoes[-length(pct_mosquitoes)])
  pct_rainfall = seq_data$Rainfall/max(seq_data$Rainfall)
  rainfall = seq_data$Rainfall
  lagged_pct_rainfall = c(NA, pct_rainfall[-length(pct_rainfall)])
  lagged_rainfall = c(NA, rainfall[-length(rainfall)])
  new_data = data.frame(pct_mosquitoes, lagged_pct_mosquitoes, 
                        pct_rainfall, lagged_pct_rainfall, 
                        rainfall, lagged_rainfall,
                        month = seq_data$Month, 
                        species = seq_data$MosquitoSpecies,
                        i, seq_order = 1:nrow(seq_data), 
                        source = seq_data$Source,
                        location = seq_data$Location)
  return(new_data)
})

scaled_season_data = do.call(rbind, scaled_season_data)
scaled_season_data$i %<>% as.factor
scaled_season_data$month %<>% as.factor

# Look at the data just to see what we have
ggplot(scaled_season_data, mapping = aes(x = seq_order, y = pct_mosquitoes, 
                                         group = i, color = species)) + 
  geom_line(linetype = "solid") + 
  geom_line(aes(x = seq_order, y = pct_rainfall, group = i, color = species),
            linetype = "dashed") + 
  facet_grid(species~location , scales = "free_x") + 
  theme_bw()

# Model percent of maximum mosquitos from lagged and current percent rainfall
seasonal_mod <- glm(pct_mosquitoes ~ lagged_pct_rainfall + pct_rainfall, 
                    data = na.omit(scaled_season_data),
                    family = binomial)

# Save the coefficients to a csv for use in the spillover model
write.csv(data.frame(coefficients = c("season_int", 
                                      "season_lagged_rain", 
                                      "season_present_rain"),
                     values = seasonal_mod$coefficients),
          "./output/submodels/coefficient_estimates/seasonality_coefficients.csv")

# Save a summary of the model for use in the supplementary materials
write.csv(summary(seasonal_mod)$coefficients, 
          "./output/submodels/seasonality_glm.csv")

# Add a column to the data with the predictions
seasonal_predict <- na.omit(scaled_season_data)
seasonal_predict %<>% cbind(
  predict(seasonal_mod, 
          new_data = c(seasonal_predict$lagged_pct_rainfall, 
                       seasonal_predict$pct_rainfall),
          type = "response", se.fit = T))
seasonal_predict$errors = seasonal_predict$pct_mosquitoes - seasonal_predict$fit
library(stringr)
seasonal_predict$source_short = str_split(seasonal_predict$source %>% 
                                            as.character, 
                                          pattern = " ", n = 3) %>% 
  laply(function(list_strings){
  paste(list_strings[[1]], list_strings[[2]])
})

# Now plot the predicted and actual seasonal abundances and save it 
season_fig <- ggplot(seasonal_predict, 
                     mapping = aes(x = seq_order, 
                                   y = pct_mosquitoes, 
                                   group = i)) + 
  geom_line(aes(color = species), linetype = "solid") +
  geom_line(aes(x = seq_order, y = fit, group = i), color = "black",
            linetype = "dashed") +
  facet_wrap(species~source_short, scales = "free_x") + 
  ylab("Relative mosquito abundance") + xlab("Month") +
  theme_bw() + guides(color = F) 
pdf('output/submodels/seasonality_data_comparison.pdf', width = 7, height = 6)
print(season_fig)
dev.off()


