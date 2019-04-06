# This is code for estimating the portion of mosquitos that will be infection 
# a certain number of days after an infectious blood meal as a function of 
# temperature. Its a little slow to run the rstan model, maybe 5ish minutes per chain.

# Set up ----
library(plyr)
library(dplyr)
library(magrittr)
library(rstan)

# Read in data and set it up ----
EIP.orig = read.csv("./data/cleaned/EIP/EIP_cleaned_data.csv", header = T, stringsAsFactors = F)
# For transmission results repeated multiple times, repeat the results multiple times in the data frame
EIP = aaply(as.matrix(EIP.orig[,- which(colnames(EIP.orig) == "repetitions")]), 2,
             function(column){base::rep(column, times = as.vector(EIP.orig$repetitions))}) %>%
  t %>%
  as.data.frame(stringsAsFactors = F,
                colnames = setdiff(colnames(EIP.orig), "repetitions")) %>% 
  mutate(record = as.integer(record),
         refeeding1_dpi = as.numeric(refeeding1_dpi),
         refeeding1_mosquitoes = as.numeric(refeeding1_mosquitoes),
         refeeding2_dpi = as.numeric(refeeding2_dpi),
         refeeding2_mosquitoes = as.numeric(refeeding2_mosquitoes),
         minT = as.numeric(minT), 
         maxT = as.numeric(maxT),
         meanT = as.numeric(meanT),
         std_meanT = (meanT - mean(meanT))/sd(meanT))

model_data = list(k = nrow(EIP),
                  X1 = EIP$refeeding1_dpi,
                  N1 = EIP$refeeding1_mosquitoes, 
                  X2 = EIP$refeeding2_dpi,
                  N2 = EIP$refeeding2_mosquitoes,
                  temp = EIP$meanT)

# Define initializing function and run stanfit ----
init_fun <- function() { 
  list(M_max_inf = runif(1, 0.5, 1),
       M_T0 = rnorm(1, 10, 1),#rnorm(1, -3, 0.1), 
       M_Tm = rnorm(1, 50, 1), #rnorm(1, 3, 0.1),
       mu0 = rnorm(1, 1.5, 0.1), 
       muT = rnorm(1, 0, 0.1), 
       sigma0 = rnorm(1, -1, 0.1),
       sigmaT = rnorm(1, 0, 0.01))
  } 

EIP.stan = stan(file = "./code/EIP/EIP_model.stan", 
                data = model_data, 
                iter = 4000, chains = 4,
                save_dso = F,
                init = init_fun,
                control = list(adapt_delta = 0.999,
                               max_treedepth = 15))

# check the output to make sure things converged----
traceplot(EIP.stan, inc_warmup = T)
print(EIP.stan)

# Extract means from stanfits----
EIP.est <- rstan::extract(EIP.stan) 
EIP.est <- do.call("cbind", EIP.est) %>% as.data.frame
param.dist <- aaply(as.matrix(EIP.est), 2, function(param){
  return(c(mean = mean(param), sd = sd(param), 
           quantile(param, probs = c(0.025)),
           quantile(param, probs = c(0.5)),
           quantile(param, probs = c(0.975))
  ))
})

# Save median parameter estimates to csv----
M_T0 = param.dist["M_T0", "50%"]
M_Tm = param.dist["M_Tm", "50%"]
c = param.dist["c", "50%"]
muT = param.dist["muT", "50%"]
mu0 = param.dist["mu0", "50%"]
sigma0 = param.dist["sigma0", "50%"]
sigmaT = param.dist["sigmaT", "50%"]

write.csv(data.frame(
  coefficient = c("T0", "Tm", "c", "muT", "mu0", "sigma0", "sigmaT"),
  value = c(M_T0, M_Tm, c, muT, mu0, sigma0, sigmaT)
), './output/submodels/coefficient_estimates/EIP_coefficients.csv')

# Plot EIP for each temperature we have, and plot the EIP functions----
pdf('./output/submodels/EIP_data_comparison.pdf', height = 4.5, width = 9)
par(mfrow = c(2,5))
for(i in sort(unique(EIP$meanT))){
  print(exp(i*muT + mu0))
  dpis <- seq(0, 45, by = 0.1)
  plot(dpis, ifelse(M_T0 < i && i < M_Tm, c*(i-M_T0)*(i-M_Tm), 0)*plnorm(dpis, meanlog = i*muT + mu0, sdlog = exp(i*sigmaT + sigma0)), type = "l", 
       xlab = "days post infection",
       ylab = "Percent infected",
       # main = paste0("std_T = ", round(i,2)), 
       main = paste0("T = ", round(i,2)), 
       ylim = c(0,1),
       las = 1, cex = 1.25)
  sub.data = subset(EIP, meanT == i)
  daily.data = data.frame(rbind(cbind(sub.data$refeeding1_dpi, sub.data$refeeding1_mosquitoes, 0),
                                cbind(sub.data$refeeding2_dpi, sub.data$refeeding2_mosquitoes, 1))) 
  colnames(daily.data) = c("day", "no_mosquitoes", "transmitting")
  daily.data %<>% subset(no_mosquitoes > 0)
  points(daily.data$day, daily.data$transmitting, 
         pch = 16, col = "blue", cex = 2)
  # abline(v = exp(i*muT + mu0), col = "red")
}
dev.off()

# Plots of median, max, and variance as functions of temperature ----
brazil_temperature = read.csv("./data/cleaned/BRT_data/mean_temperature_200002_201612.csv",
                       stringsAsFactors = F) %>% 
  dplyr::select(-c("system.index", ".geo", "CD_GEOCMU", "NM_MUNICIP")) 
brazil_temperature = c(min(brazil_temperature, na.rm = T), max(brazil_temperature, na.rm = T)) 
brazil_temperature = brazil_temperature - 273
obs_range = c(min(EIP$meanT),max(EIP$meanT) ) 
pdf('./output/submodels/EIP_model.pdf')
par(mfrow = c(2,2))
temps = seq(0, 50, by = 0.1)
plot(temps, 
     ifelse(temps < M_T0 | temps > M_Tm, 0, c*(temps-M_T0)*(temps - M_Tm)), 
     type = "l", main = "Vector competence",
     ylab = "proportion",
     xlab = "Temperature",
     las = 1)
points(obs_range, c(0,0), col = "blue", type = "l", lwd = 3)
points(brazil_temperature, c(0,0)+0.02, col = "red", type = "l", lwd = 3)
plot(temps, 
     exp(temps*muT + mu0), 
     type = "l", main = "Time until 50% infectious",
     ylab = expression('EIP'[50]),
     xlab = "Temperature",
     las = 1)
points(obs_range, rep(min(exp(temps*muT + mu0)), 2), col = "blue", type = "l", lwd = 3)
points(brazil_temperature, rep(min(exp(temps*muT + mu0)), 2) + 0.4, col = "red", type = "l", lwd = 3)
plot(temps, 
     exp(sigma0 + sigmaT*temps), 
     type = "l", main = "Standard deviation",
     ylab = "Standard deviation of EIP",
     xlab = "Temperature",
     las = 1)
points(obs_range, rep(min(exp(temps*sigmaT + sigma0)), 2), col = "blue", type = "l", lwd = 3)
points(brazil_temperature, rep(min(exp(temps*sigmaT + sigma0)), 2) + 0.01, col = "red", type = "l", lwd = 3)
dev.off()

# Compare found distributions with priors
# Distribution is in red, starting parameters in blue
par(mfrow = c(3,3))
hist(EIP.est$M_max_inf, freq = F, xlim = c(0,1))
points(seq(0,1, by = 0.01), dunif(seq(0,1, by = 0.01), 0, 1), type = "l", col = "red")
abline(v = laply(EIP.stan@inits, function(list_elem){list_elem$M_max_inf}), col = "blue")
abline(v = param.dist["M_max_inf", 3:5], col = "green")

hist(EIP.est$M_T0, freq = F, xlim = c(0, 50)) 
points(seq(0,50, by = 0.01), dnorm(seq(0,50, by = 0.01), 10, 5), type = "l", col = "red")
abline(v = laply(EIP.stan@inits, function(list_elem){list_elem$M_T0}), col = "blue")
abline(v = param.dist["M_T0", 3:5], col = "green")

hist(EIP.est$M_Tm, freq = F, xlim = c(0, 50))
points(seq(0,50, by = 0.01), dnorm(seq(0,50, by = 0.01), 40, 5), type = "l", col = "red")
abline(v = laply(EIP.stan@inits, function(list_elem){list_elem$M_Tm}), col = "blue")
abline(v = param.dist["M_Tm", 3:5], col = "green")

hist(EIP.est$mu0, freq = F, xlim = c(0, 4))
points(seq(0,4, by = 0.01), dnorm(seq(0,4, by = 0.01), 2, 0.5), type = "l", col = "red")
abline(v = laply(EIP.stan@inits, function(list_elem){list_elem$mu0}), col = "blue")
abline(v = param.dist["mu0", 3:5], col = "green")

hist(EIP.est$muT, freq = F, xlim = c(-0.5, 0.5))
points(seq(-5,5, by = 0.01), dnorm(seq(-5,5, by = 0.01), 0, 0.05), type = "l", col = "red")
abline(v = laply(EIP.stan@inits, function(list_elem){list_elem$muT}), col = "blue")
abline(v = param.dist["muT", 3:5], col = "green")

hist(EIP.est$sigma0, freq = F, xlim = c(-4, 0))
points(seq(-4,0, by = 0.01), dnorm(seq(-4,0, by = 0.01), -1, 0.5), type = "l", col = "red")
abline(v = laply(EIP.stan@inits, function(list_elem){list_elem$sigma0}), col = "blue")
abline(v = param.dist["sigma0", 3:5], col = "green")

hist(EIP.est$sigmaT, freq = F, xlim = c(-0.5, 0.5))
points(seq(-5,5, by = 0.01), dnorm(seq(-5,5, by = 0.01), 0, 0.1), type = "l", col = "red")
abline(v = laply(EIP.stan@inits, function(list_elem){list_elem$sigmaT}), col = "blue")
abline(v = param.dist["sigmaT", 3:5], col = "green")

