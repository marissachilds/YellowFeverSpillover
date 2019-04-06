# This code fits a sine function to the number of annual number of 
# municipality-months with spillover so approximate primate dynamics for use
# in the periodic disease risk metric and the boosted regression tree analysis.
# The coefficients for the sine function are saved to a csv

# Load libraries and source the clean_case_data script to get case data
library(plyr)
library(dplyr)
library(magrittr)
source("./code/other/clean_case_data.R")

# For each year, count the number of municipality-months with non-zero cases
year_spillovers = ddply(case.data, .(year), function(sub){
  sum(sub$cases > 0)
})
colnames(year_spillovers) = c("year", "spillovers")

# Clean up a bit and add that in 2011 and 2012 there were no cases
year_spillovers %<>% mutate(year = as.character(year) %>% as.numeric)
year_spillovers %<>% rbind(c(2011, 0)) %>% rbind(c(2012,0))
year_spillovers = year_spillovers[order(year_spillovers$year),]

# Run a linear model with sin and cosine of year as predictors
spill_period <- lm(spillovers ~ sin(2*pi*year/7)+cos(2*pi*year/7), 
                   data=year_spillovers)

# Lets just take a look at the results 
summary(spill_period)
par(mfrow = c(1,1))
plot(spillovers~year,data=year_spillovers)
lines(year_spillovers$year,spill_period$fitted,col=2)

# Extract the coefficients and plot the sine function compared to the data
sin_coef = spill_period$coefficients[2]
cos_coef = spill_period$coefficients[3]
years = seq(2001, 2018, by = 1/12)
plot(years, sin_coef*sin(2*pi*years/7) + cos_coef*cos(2*pi*years/7), type = "l")

# Using wolframalpha we calculate min and max to be sqrt(10919733629929/10)/100000 
# Lets rescale the coefs to make sure our periodic function goes from 0 to 1
sin_coef_rescaled = sin_coef/(2*9.781)
cos_coef_rescaled = cos_coef/(2*9.781)

# Save the scaled coefficients
write.csv(data.frame(
  coefficient = c("sin_coef", "cos_coef"),
  value = c(sin_coef_rescaled, cos_coef_rescaled)
), './output/submodels/coefficient_estimates/phenom_periodic_coefficients.csv')

# Can also save a pdf of the plot of comparing rescaled sine curve with case data
# pdf("./output/submodels/phenom_data_comparison.pdf", width = 9, height = 6)
# par(mfrow = c(1,1))
# par(oma = c(0,1,0,3))
# years = seq(2001, 2018, by = 1/12)
# plot(years, sin_coef_rescaled*sin(2*pi*years/7) + 
#        cos_coef_rescaled*cos(2*pi*years/7) + 1/2, type = "l",
#      ylab = "Periodic Reservoir Infection Prevalence", xlab = "Year")
# par(new = TRUE)
# plot(year_spillovers$spillovers ~ year_spillovers$year, col = "red", pch = 19,
#      xlab = "", ylab = "", xlim = c(2001, 2018), axes = FALSE, bty = "n")
# axis(side=4, at = pretty(range(year_spillovers$spillovers)), las = 1)
# mtext("Number of municipality-months\nwith spillover", side=4, line=3, cex = 1)
# dev.off()
