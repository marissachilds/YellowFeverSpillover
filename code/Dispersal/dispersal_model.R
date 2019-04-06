# code for fitting a model to mosquito dispersal data. save csv with the 
# coefficient for use in the overall spillover model and save figures for
# use in the supplementary materials

# set up packages and set wd
library(plyr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(rstan)

# Read in data of capture locations
F1 <- read.csv("data/cleaned/dispersal/Preliminary_F1.csv", header = F, 
               col.names = c("R", "theta", "names"), 
               stringsAsFactors = F)
F1$names = trimws(F1$names)

# Define a function to convert from polar coordinates to x-y
polar_to_xy <- function(r, theta, degrees = T){
  if(degrees == T){
    theta = theta*pi/180
  }
  x <- r*cos(theta)
  y <- r*sin(theta)
  return(cbind(x,y))
}

# Convert capture locations to x-y cartesian
F1_xy <- polar_to_xy(F1$R, F1$theta) %>% as.data.frame 
E.x = F1_xy[F1_xy$names == "E", "x"]
E.y = F1_xy[F1_xy$names == "E", "y"]

# Calculate distance from release location
F1_xy %<>% cbind(names = F1$names) %>%
  mutate(names = as.character(names),
         distance_from_release = sqrt((x)^2+(y)^2))

# Plot the locations just to check
plot(F1_xy$x, F1_xy$y, pch = F1_xy$names)

# read in the data on number of mosquitoes caught at each location
F1mosquitoes_orig <- read.csv("./data/cleaned/Dispersal/dispersal_data_extraction.csv", 
                              stringsAsFactors = F)

# a little cleaning and renaming of the mosquito capture data, keeping only
# necessary variables, and then joining with the capture location data
F1mosquitoes <- transmute(F1mosquitoes_orig, 
                          station_name = do.call(rbind, 
                                                 strsplit(Capture.station, 
                                                          split = " "))[,2],
                          species = Species, 
                          number_caught = Number.caught,
                          number_per_time = Number.caught/
                            Time.spent.in.man.hours,
                          time_spent = Time.spent.in.man.hours,
                          standardized_catch = (3600/Time.spent.in.man.hours)*
                            number_caught) %>% 
  left_join(F1_xy, by = c("station_name" = "names")) 
F1mosquitoes <- F1mosquitoes[order(F1mosquitoes$distance_from_release),]  

# calculate total captures for each station
total_data = F1mosquitoes %>% select(c(station_name, number_caught, time_spent, 
                                       distance_from_release)) %>% 
  ddply(.(station_name, distance_from_release, time_spent), 
        summarise, number_caught = sum(number_caught))

# Set up data for use in the rstan model
model_data = list(n = nrow(total_data),
                  Y = total_data$number_caught,
                  r = total_data$distance_from_release,
                  t = total_data$time_spent)

# Define initial parameters
init_fun <- function() { 
  list(alpha = runif(1, -10, 10),
       beta = runif(1, 0, 1/2),
       k = runif(1, 0, 10))
} 

# Run the model with rstan, using the model defined elsewhere
dispersal.stan = stan(file = "./code/Dispersal/dispersal_model.stan", 
                data = model_data, 
                iter = 4000, chains = 4,
                save_dso = F,
                init = init_fun,
                control = list(adapt_delta = 0.99999,
                               max_treedepth = 20))

# Look at the draws to make sure they've converged
traceplot(dispersal.stan, inc_warmup = T, pars = c("alpha", "beta", "k"))

# Get the post-warmup draws from the stan object
dispersal.dist <- do.call("cbind", rstan::extract(dispersal.stan) ) %>% 
  as.data.frame

# For each parameter, calculate mean, sd, and quantiles
param.dist <- aaply(as.matrix(dispersal.dist), 2, function(param){
  return(c(mean = mean(param), sd = sd(param), 
           quantile(param, probs = c(0.025)),
           quantile(param, probs = c(0.5)),
           quantile(param, probs = c(0.975))
  ))
})

# Use the median as our parameters
alpha = param.dist["alpha", "50%"]
beta = param.dist["beta", "50%"]
k = param.dist["k", "50%"]

# Save beta, which is the only relevant parameter, the others are just scaling factors
write.csv(data.frame(coefficient = "beta", value = beta), 
          file = "./output/submodels/coefficient_estimates/dispersal_coefficients.csv")

# compare posterior with prior distributions to see how much they have changed
par(mfrow = c(2,2))
hist(dispersal.dist$alpha, freq = F, xlim = c(-10, 10))
points(seq(-50,50, by = 0.1), dunif(seq(-50,50, by = 0.1), -100, 100), 
       type = "l", col = "red")
abline(v = laply(dispersal.stan@inits, function(list_elem){list_elem$alpha}), 
       col = "blue")

hist(dispersal.dist$beta, freq = F, xlim = c(0,1/2))
points(seq(-50,50, by = 0.1), dunif(seq(-50,50, by = 0.1), 0, 100), 
       type = "l", col = "red")
abline(v = laply(dispersal.stan@inits, function(list_elem){list_elem$beta}), 
       col = "blue")

hist(dispersal.dist$k, freq = F, xlim = c(0,10))
points(seq(-50,50, by = 0.1), dunif(seq(-50,50, by = 0.1), 0, 100), 
       type = "l", col = "red")
abline(v = laply(dispersal.stan@inits, function(list_elem){list_elem$k}), 
       col = "blue")

# compare model estimates to capture data
par(mfrow = c(1,1))
plot_data = data.frame(model_data$Y, model_data$r, model_data$t, 
                       total_data$station_name)
colnames(plot_data) = c("Y", "r", "t", "name")

# Save a figure of model-data comparison to pdf for the supplement
pdf("./output/submodels/dispersal_data_comparison.pdf")
ggplot(plot_data, aes(x = r, y = Y, label = name)) + 
  geom_ribbon(aes(ymax = qnbinom(0.975, mu = t*exp(alpha)*exp(-r*beta), 
                                 size = k),
                  ymin = qnbinom(0.025, mu = t*exp(alpha)*exp(-r*beta), 
                                 size = k)),
              fill = "grey70", alpha = 0.6) +
  geom_line(aes(x = r, y = t*exp(alpha)*exp(-r*beta)), color = "red") +
  geom_text() + 
  theme_bw() + 
  ylab("Number Caught") + 
  xlab("Distance from Release") %>% print
dev.off()

# png("./output/submodels/dispersal_mod.png")
# ggplot(plot_data, aes(x = r, y = Y, label = name)) + 
#   geom_line(aes(x = r, y = beta/(2*pi)*exp(-r*beta)), color = "red") +
#   theme_bw() + 
#   ylab("Percent Mosquitoes") + 
#   xlab("Distance from Release")
# dev.off()

# Because we can, plot 95% confidence of curves
distances = 0:1000
dispersal_kerns = aaply(dispersal.dist$beta, 1, function(beta){
  beta/(2*pi)*exp(-distances*beta)
}) %>% t
matplot(distances, dispersal_kerns, type = 'l', col = "grey", lty = 1)
points(distances, beta/(2*pi)*exp(-distances*beta), type = "l", col = "black")