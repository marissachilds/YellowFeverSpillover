# This code models mosquito survival as a quadratic function of temperature. 
# It take in data for both field and lab observations, assumes each setting has
# a different max lifespan but the same critical thermal limits and calculates
# max lifespans in both and thermal limits. These coefficients are saved to a csv file

# Setup ----
library(plyr)
library(dplyr)
library(magrittr)
library(rstan)

# Read in data and format for model fitting ----
survival_data = read.csv("./data/cleaned/Survival/survival_data.csv", 
                         stringsAsFactors = F)

model_data = list(k = nrow(survival_data),
                  L = survival_data[,2], 
                  temp = survival_data[,1],
                  loc = (survival_data[,4] == "lab"))

# Define initializing function ----
init_fun <- function() { 
  list(max_L_lab = runif(1, 0, 200),
       max_L_field = runif(1, 0, 200),
       T0 = rnorm(1, 10, 2), # Relatively confined initial values for thermal limits so that it starts
       Tm = rnorm(1, 40, 2), # Relatively confined initial values for thermal limits so that it starts
       sigma = rexp(1, 0.1))
} 

# Fit model with stan ----
survival.stan = stan(file = "./code/Survival/survival_model.stan", 
                data = model_data, 
                iter = 6000, chains = 6,
                save_dso = F,
                init = init_fun,
                control = list(adapt_delta = 0.999))
print(survival.stan)
traceplot(survival.stan)

# Extract distributions from stanfits ----
survival.est <- rstan::extract(survival.stan) 
survival.est <- do.call("cbind", survival.est) %>% as.data.frame
param.dist <- aaply(as.matrix(survival.est), 2, function(param){
  return(c(mean = mean(param), sd = sd(param), 
           quantile(param, probs = c(0.025)),
           quantile(param, probs = c(0.5)),
           quantile(param, probs = c(0.975))
  ))
})

# Save median parameter estimate to csv ----
T0 = param.dist["T0", "50%"]
Tm = param.dist["Tm", "50%"]
c_lab = param.dist["c_lab", "50%"]
c_field = param.dist["c_field", "50%"]
max_L_lab = param.dist["max_L_lab", "50%"]
max_L_field = param.dist["max_L_field", "50%"]
sigma = param.dist["sigma", "50%"]
write.csv(data.frame(
  coefficients = c("s_T0", "s_Tm", "s_c"),
  values = c(T0, Tm, c_field)),
  "./output/submodels/coefficient_estimates/survival_coefficients.csv")

# Plot lifespan as a function of temperature
png("./output/submodels/survival_data_comparison.png")
temps = seq(0,50, by = 0.1)
plot(temps, ifelse(T0 < temps & temps < Tm, c_field*(temps-T0)*(temps-Tm), 0), type = "l",
     xlab = "temperatures (degrees C)",
     ylab = "lifespan (days)", lty = "solid",
     ylim = c(0,45))
points(temps, ifelse(T0 < temps & temps < Tm, c_lab*(temps-T0)*(temps-Tm), 0), type = "l",
    lty = "dashed")
points(survival_data$temperature_mean, survival_data$longevity_days,
       pch = ifelse(survival_data$location == "lab", 1, 16),
       col = "red")
dev.off()

# Compare found distributions with priors
# Distribution is in red, starting parameters in blue
par(mfrow = c(2,2))
hist(survival.est$max_L, freq = F, xlim = c(0,100))
points(seq(0,100, by = 0.01), dunif(seq(0,200, by = 0.01), 0, 100), type = "l", col = "red")
abline(v = laply(survival.stan@inits, function(list_elem){list_elem$max_L}), col = "blue")
# abline(v = param.dist["max_L", 3:5], col = "green")

hist(survival.est$T0, freq = F, xlim = c(-10, 20))
points(seq(-10,20, by = 0.01), dnorm(seq(-10,20, by = 0.01), 10, 5), type = "l", col = "red")
abline(v = laply(survival.stan@inits, function(list_elem){list_elem$T0}), col = "blue")
# abline(v = param.dist["T0", 3:5], col = "green")

hist(survival.est$Tm, freq = F, xlim = c(30, 60))
points(seq(30,60, by = 0.01), dnorm(seq(30,60, by = 0.01), 40, 5), type = "l", col = "red")
abline(v = laply(survival.stan@inits, function(list_elem){list_elem$Tm}), col = "blue")
# abline(v = param.dist["Tm", 3:5], col = "green")

hist(survival.est$sigma, freq = F, xlim = c(0, 20))
points(seq(0,20, by = 0.01), dexp(seq(0,20, by = 0.01), 0.1), type = "l", col = "red")
abline(v = laply(survival.stan@inits, function(list_elem){list_elem$sigma}), col = "blue")
# abline(v = param.dist["sigma", 3:5], col = "green")

# Calculate daily survival probability from expected lifespan
par(mfrow = c(1,1))
temps = seq(0,50, by = 0.1)
plot(temps, exp(-1/ifelse(T0 < temps & temps < Tm, c*(temps-T0)*(temps-Tm), 0)), type = "l",
     xlab = "temperatures (degrees C)",
     ylab = "probability of survival (per day)")

# Calculate probability of surviving 10 days
plot(temps, exp(-1/ifelse(T0 < temps & temps < Tm, c*(temps-T0)*(temps-Tm), 0))^10, type = "l",
     xlab = "temperatures (degrees C)",
     ylab = "probability",
     main = "10-day survival",
     ylim = c(0,1))

# Compare it to data from same Bates 1947 paper, but survival data from Table VII
pct_alive = c(53, 47, 43, 2)/100
alive_temps = c(25,  26.6,  30, 35)
points(alive_temps, pct_alive, col = "red", pch = 16)


