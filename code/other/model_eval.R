# Setup -----
library(plyr)
library(dplyr)
library(raster)
library(magrittr)
library(pracma)
library(ggplot2)
library(stringr)
library(rgdal)
library(pdp)
library(rgeos)
library(gbm)
library(dismo)
library(grid)

# Read in model and case data ----
model_case = readRDS("./data/cleaned/BRT_data/brt_data.rds")  
# Evaluating model ability to predict spillover with AIC and AUC ---- 
n_hypoth = 17
fpr_tpr_calc <- function(labels, predictions){
  df <- tibble(l = unlist(labels), p = unlist(predictions)) %>%
    arrange(desc(p)) %>% 
    mutate(N = length(unlist(labels)),#n(),
           N_pos = sum(l),
           N_neg = N - N_pos,
           pred = 1,
           pos = cumsum(pred),
           true = (pred == l),
           tp = cumsum(true),
           fp = pos - tp,
           tpr = tp/N_pos,
           fpr = fp/N_neg)
  return(dplyr::select(df, tpr,fpr))  
}

# Define the predictor variables
predictors <- c("mean_EnvRisk", "max_EnvRisk", 
                "mean_PhenomEnvRisk", "max_PhenomEnvRisk",
                "mean_ImmEnvRisk", "max_ImmEnvRisk", 
                "mean_PopRisk", "max_PopRisk"
                )
# For each predictor, run a GLM of spillover predicted by the predictor 
# and also calculate AUC from the predictor. Return the AIC from the GLM,
# the coefficient and p-value from the GLM, the FWER corrected pvalue, and AUC
model_evals <- adply(predictors, 1, function(pred){
  mod <- glm(substitute(spillover ~ i, list(i = as.name(pred))), 
             family = "binomial", data = model_case)
  fpr_tpr <- fpr_tpr_calc((na.omit(model_case)$cases > 0), 
                          na.omit(model_case)[,pred])
  return(c(pred = pred, 
           AIC = summary(mod)$aic, 
           pred_coef = summary(mod)$coef[2,1],
           pred_coef_pval = summary(mod)$coef[2,4],
           auc = pracma::trapz(fpr_tpr$fpr, fpr_tpr$tpr),
           adj_pred_coef_pval = min(summary(mod)$coef[2,4]*n_hypoth, 1)))
}, .id = NULL)

# Add columns saying which risk metric, and which municipality summary for 
# use in a table in the appendix
model_evals$summary = strsplit(as.character(model_evals$pred), 
                               split = "_") %>% 
  sapply("[", 1)
model_evals$mechanism = strsplit(as.character(model_evals$pred), 
                                 split = "_") %>% 
  sapply("[", 2)

# Make the mechanisms a factor so they can be ordered by something other than alphabetical
model_evals$mechanism = factor(model_evals$mechanism, 
                               levels = c("EnvRisk", 
                                          "PhenomEnvRisk", 
                                          "ImmEnvRisk", 
                                          "PopRisk"),
                               ordered = F)
write.csv(model_evals, "./output/results/model_evals.csv")

# Evaluating model + vaccine coverage ability to estimate cases given spillover  ----
# subset to just municipality-months where spillover occurred
spillover_sub = subset(na.omit(model_case), spillover == 1)
# Make a list of predictors on interest
predictors <- c("mean_EnvRisk", "max_EnvRisk", 
                "mean_PhenomEnvRisk", "max_PhenomEnvRisk",
                "mean_ImmEnvRisk", "max_ImmEnvRisk", 
                "mean_PopRisk", "max_PopRisk", "vax"
)
# For each predictor, regress # of cases on the predictor and extract the 
# R^2, adjusted R^2, predictor coefficient, predictor coefficient p-value.
# also calculate adjusted FWER pvalue and spearman correlation coefficient
case_preds <- adply(predictors, 1, function(pred){
  mod <- lm(substitute(cases ~ i, list(i = as.name(pred))), 
            data = spillover_sub)
  return(c(pred = pred, 
           r.sq = summary(mod)$r.squared,
           adj.r.sq = summary(mod)$adj.r.squared,
           coef = summary(mod)$coefficients[2,1],
           pval = summary(mod)$coefficients[2,4],
           adj_pval = min(summary(mod)$coefficients[2,4]*n_hypoth,1),
           spearman_cor = cor(model_case[[pred]], 
                              model_case$cases, 
                              method = "spearman",
                              use = "complete.obs")))
}, .id = NULL)
write.csv(case_preds, "./output/results/cases_preds.csv")

# Plot Figure 3: comparing mechanism evaluations ----
# Save figure 3 to pdf. Figure 3 compares AUC for the different mechanistic 
# models and municipality summary statistics. 
pdf("./output/results/Figure3.pdf", width = 6, height = 6)
par(mfrow = c(1,1), mar = c(6, 4, 4, 2) + 0.1)
plot(x = as.factor(subset(model_evals, summary  == "max")$mechanism), 
     y = subset(model_evals, summary  == "max")$auc, 
     xaxt = "n", las = 1,
     col = "red", pch = 16,
     xlab = "", ylab = "", bty = "l",
     type = "b", lty = "dashed", 
     ylim = c(0.45, 0.8),
     xlim = c(0.9, 4.1))
title(main = "Model performance by mechanism",
      font= 2, cex = 1.5)
mtext(side=1, line=4, "Mechanism",font=2,cex=1.3)
mtext(side=2, line=3, "AUC",font=2,cex=1.3)
points(x = as.factor(subset(model_evals, summary  == "mean")$mechanism), 
       y = subset(model_evals, summary  == "mean")$auc, 
       col = "blue", pch = 17,
       type = "b", lty = "dotted")
axis(1, at=1:4,line = 1,
     labels = c("Environmental\nRisk", "Periodic\nRisk", 
                "Immunological\nRisk", "Population-\nscaled\nRisk"),
     tick = F)
text(x=3.7, y=0.66, "Maximum risk",font=2, cex=0.9)
text(x=3.7, y=0.51, "Mean risk",font=2, cex=0.0)
dev.off()

# Figure 4a: seasonal pattern of model and case data ----
# In places with known municipality codes, calculate average, total, and max 
# of the periodic and environmental risk metrics and spillovers for each 
# municipality and month. Essentially, calculate the average seasonal trend 
# in risk metrics and number of spillovers for each municipality
seasonal_model_case <- ddply(model_case[model_case$CD_GEOCMU != "0", ] 
                              %>% dplyr::select(c(CD_GEOCMU, max_EnvRisk, 
                                               max_PhenomEnvRisk,
                                               month, cases, spillover,
                                               region)) %>% 
                               dplyr::rename(actual_cases = cases), 
                             .(CD_GEOCMU, month, region), summarize_if, 
                             .predicate = is.numeric, .funs = c(mean = mean, 
                                                                total = sum,
                                                                max = max), 
                             na.rm = T)
seasonal_model_case %<>% mutate(region = as.character(region))
seasonal_model_case_all = seasonal_model_case
# Taking the complete cases, calculate the total number of spillovers and average 
# environmental risk over each region 
seasonal_model_case = seasonal_model_case[complete.cases(seasonal_model_case),]
seasonal_regional_avg = ddply(seasonal_model_case, .(region, month), 
                                           summarize_at, 
                              .vars = c("max_EnvRisk_mean","spillover_total"), 
                              .funs = c(total = sum,
                                        mean = mean)) 
# Then calculate the correlation between region average seasonal risk and region 
# total seasonal spillovers
seasonal_regional_cor = daply(seasonal_regional_avg, .(region), 
                              function(reg_sub){
                                rho = cor(reg_sub$spillover_total_total, 
                                          reg_sub$max_EnvRisk_mean_mean)
                                out = paste0(reg_sub$region %>% 
                                               as.character %>% 
                                               unique, " (", round(rho, 2), ")")
                                if(is.na(rho)){
                                  out = reg_sub$region %>% 
                                    as.character %>% unique
                                }
                                out
                              })
# Plot Figure 4a ----
# Turn the above information into a figure with a panel for each region
pal = c("#FF0000", "#00A08A", "#F2AD00", "#F98400", "#5BBCD6") 
Fig4A <- ggplot(seasonal_model_case, aes(y = max_EnvRisk_mean, x = month, 
                                         group = CD_GEOCMU, color = region)) +
  geom_line(alpha = 0.2) + # Add a line for each municipality's seasonal average of max environmental risk
  geom_line(aes(y = max_EnvRisk_mean_mean, # Add a line for the region's average environmental risk(average is over the municipalities here)
                x = month, group = region), color = "white",
            data = seasonal_regional_avg) +
  geom_point(aes(y = spillover_total_total/110, # Add points for seasonal total number of spillover in each region 
                  x = month, group = region), color = "black",
             data = seasonal_regional_avg,
             alpha = 1) + 
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), 
        plot.title = element_text(face="bold"))+
  scale_x_discrete(breaks = c("1", "4", "7", "10")) + 
  scale_y_continuous(trans = "sqrt", 
    sec.axis = # Add a second y-axis for the number of spillovers
      sec_axis(~.*110, 
               name = "Number of municipality-months with spillover")) + 
  xlab("Month") + ylab("Average environmental risk") +
  ggtitle("a. Seasonal variation") + 
  guides(color = F, 
         size = guide_legend(title="Number of
                             \nmunicipality-months
                             \nwith spillovers")) +
  scale_colour_manual(values=pal) + 
  facet_wrap(~region, ncol = 1, scales = "fixed", 
             labeller=labeller(region = seasonal_regional_cor))

# Fibure 4c: Plot year patterns ----
# Repeat as above, but this time do calculations for each year so we get an
# interannual risk pattern and spillovers pattern
year_model_case <- ddply(model_case[model_case$CD_GEOCMU != "0", ] %>% 
                           mutate(spillover = (cases > 0)*1) %>% 
                               dplyr::select(c(CD_GEOCMU, max_EnvRisk, max_PhenomEnvRisk,
                                               region, year, cases, spillover)) %>% 
                               dplyr::rename(actual_cases = cases), 
                             .(CD_GEOCMU, year, region), summarize_if, 
                             .predicate = is.numeric, .funs = c(mean = mean, 
                                                                total = sum,
                                                                max = max), 
                         na.rm = T)

year_model_case %<>% mutate(region = as.character(region))
year_model_case_all = year_model_case
year_model_case = year_model_case[complete.cases(year_model_case),]
year_regional_avg = ddply(year_model_case, .(region, year), 
                              summarize_at, 
                              .vars = c("max_PhenomEnvRisk_mean",
                                        "spillover_total"), 
                              .funs = c(total = sum,
                                        mean = mean))
year_regional_cor = daply(year_regional_avg, .(region), function(reg_sub){
  rho = cor(reg_sub$spillover_total_total, reg_sub$max_PhenomEnvRisk_mean_mean)
  out = paste0(reg_sub$region %>% as.character %>% unique, 
               " (", round(rho, 2), ")")
  if(is.na(rho)){
    out = reg_sub$region %>% as.character %>% unique
  }
  out
})
# Plot Figure 4c ----
Fig4C <- ggplot(year_model_case, aes(y = max_PhenomEnvRisk_mean, x = year, 
                                         group = CD_GEOCMU, color = region)) +
  geom_line(alpha = 0.2) +
  geom_line(aes(y = max_PhenomEnvRisk_mean_mean,
                x = year, group = region), color = "white",
            data = year_regional_avg) +
  geom_point(aes(y = spillover_total_total/110, 
                 x = year, group = region), color = "black",
             data = year_regional_avg,
             alpha = 1) + 
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), 
        plot.title = element_text(face="bold"))+
  scale_x_continuous(breaks = c(2005, 2010, 2015)) + 
  scale_y_continuous(trans = "sqrt",
    sec.axis = 
      sec_axis(~.*110, 
               name = "Number of municipality-months with spillover")) + 
  xlab("Year") + ylab("Average periodic risk") +
  ggtitle("c. Interannual variation")+# in risk and spillover") + 
  guides(color = F, 
         size = guide_legend(title="Number of
                             \nmunicipality-months
                             \nwith spillovers")) +
  scale_colour_manual(values=pal) + 
  facet_wrap(~region, ncol = 1, scales = "fixed", 
             labeller=labeller(region = year_regional_cor))


# Figure 4b: Brazil's regions and municipalities with spillover ----
# Now lets make a panel with brazil's regions colored to match the lines from
# the previous two figures
# Read in the shapefile and add a state column to the dataframe
shape = readOGR(dsn = './data/raw/brazil_border_shapefiles/br_municipios')
shape@data$state = as.character(shape@data$CD_GEOCMU) %>% 
  substr(1,2) %>% mapvalues(from = c( 11,   12,   13,   14,   15,   
                                      16,   17,   21,   22,   23,   
                                      24,   25,   26,   27,   28,   
                                      29,   31,   32,   33,   35,   
                                      41,   42,   43,   50,   51,   
                                      52,   53) %>% as.character,
                            to = c('RO', 'AC', 'AM', 'RR', 'PA', 
                                   'AP', 'TO', 'MA', 'PI', 'CE', 
                                   'RN', 'PB', 'PE', 'AL', 'SE', 
                                   'BA', 'MG', 'ES', 'RJ', 'SP', 
                                   'PR', 'SC', 'RS', 'MS', 'MT', 
                                   'GO', 'DF'))
# Use the state name to make a region name column
shape@data$region = mapvalues(shape@data$state,
  from = c(c("RR", "AM", "AC", "RO", "PA", "AP", "TO"),
                     c("MA", "PI", "CE", "RN", "PB", "PE", "AL", "SE", "BA"),
                     c("MT", "GO", "DF", "MS"),
                     c("MG", "ES", "RJ", "SP"),
                     c("PR", "SC", "RS")),
            to = c(rep("North", 7), 
                   rep("Northeast", 9), 
                   rep("Central-West", 4), 
                   rep("Southeast", 4), 
                   rep("South", 3)))
# Merge by region
shape_region = gUnaryUnion(shape, id = shape@data$region)
shape_region = gSimplify(shape_region, tol = 0.1)

# Get the centroids of municipalities with spillover for plotting
spillover_sub = subset(na.omit(model_case), spillover == 1)
spillover_munis = shape[shape$CD_GEOCMU %in% spillover_sub$CD_GEOCMU, ]
spillover_cents = gCentroid(spillover_munis, byid = T)@coords

# Combine Figure 4 components and output to pdf ----
pdf("./output/results/Figure4.pdf", width = 10.5, height = 7)
par(mfrow=c(1,3))
plot.new() # Skip the first spot so we can put the brazil map in the middle
plot(shape_region, col = pal, adj = 0,
     xlim = shape_region@bbox[1,] + c(1, 0), 
     ylim = shape_region@bbox[2,])
title("b. Brazilian regions", cex.main = 1.7, line = -10)
points(spillover_cents[,1], spillover_cents[,2], pch = 16, cex = 0.6)

# create viewports for adding ggplots
vp.left <- viewport(height=unit(1, "npc"), width=unit(1/3, "npc"), 
                    just=c("left","top"), 
                    y=1, x=0)
vp.right <- viewport(height=unit(1, "npc"), width=unit(1/3, "npc"), 
                      just=c("left","top"), 
                      y=1, x=2/3)

# plot the ggplot using the print command
print(Fig4A, vp=vp.left)
print(Fig4C, vp=vp.right)
dev.off()

# Figure 5: boosted regression tree results  ----
# Read in the BRT fit with the lowest residual deviance (previously calculated) 
# and the test data
brt_fit <- readRDS("./output/BRT_fits/spillover_brt_tc10_lr0.001.rds")
test_data <- readRDS("./data/cleaned/BRT_data/brt_test.rds")

# calculate training AUC, just print it to the screen
train_roc = fpr_tpr_calc(brt_fit$gbm.call$dataframe[,brt_fit$gbm.call$gbm.y], 
                         brt_fit$fit)
pracma::trapz(train_roc$fpr, train_roc$tpr)
# calculate test AUC, again print to screen
test_preds <- predict.gbm(brt_fit, test_data, n.trees = brt_fit$n.trees)
test_roc = fpr_tpr_calc(test_data$spillover, 
                         test_preds)
pracma::trapz(test_roc$fpr, test_roc$tpr)

# Set up a dictionary to go from var.name to clean name for plotting
var.names.plotting = data.frame(
  var.names = c("max_EnvRisk", "lagged_max_EnvRisk", "pop_density", 
                "primates_max", "primates_mean", "phenom", "Tair_mean", 
                "precipitation", "month", "region", "vax", "fireArea",
                "fire_pct", "lagged_fireArea"),
  plotting.names = c("Maximum environmental risk", 
                     "Lagged maximum\nenvironmental risk",
                     "Log(Population density)",
                     "Maximum primate richness", 
                     "Average primate richness",
                     "Phenomenological primate\ndynamics", 
                     "Temperature (C)", "Precipitation (mm/hr)",
                     "Month", "Region", "Vaccine coverage (%)", 
                     "Fire Area", "Fire Percent", "Lagged Fire Area"),
  stringsAsFactors = F)

# get the data frame from the gbm call
brt_subset = brt_fit$gbm.call$dataframe

# j tracks the jth most important variable that we are looking at
# This can be run for j in 1:6 and mfrow = (2,3) for the main text figure
# and run with j in 1:14, mfrow = c(4,4) for the appendix figure, making 
# sure to change the pdf name. to switch between the two, just comment and 
# uncomment the appropriate lines in the following
pdf("./output/results/Figure5.pdf", width = 9.5, height = 5)
# pdf("./output/results/Figure5_appendix.pdf", width = 12, height = 9.5)
par(mfrow = c(2,3))
# par(mfrow = c(4,4))
par(mar=c(5, 5, 1.5, 2.5) + 0.1)
par(oma = c(0, 0, 0, 4))
par(cex.axis = 1.2, cex.lab = 1.2)
for (j in 1:6){
# for (j in 1:14){
  # Depending on where the plot is placed, set whether its left and right y 
  # axes will be labeled
  if( j %% par('mfrow')[2] == 1){left_ax = T
  } else{left_ax = F}
  if(j %% par('mfrow')[2] == 0){right_ax = T
  } else{right_ax = F}
  
  # Identify the name of the variable
  var.name = brt_fit$contributions$var[j] %>% as.character
  print(var.name)
  # Figure out which of the predictor names it was and get the data for it
  k <- match(brt_fit$contributions$var[j], brt_fit$gbm.call$predictor.names)
  var.data = brt_fit$gbm.call$dataframe[[var.name]]
  # also get the relative importance of the variable
  var.import = brt_fit$contributions[brt_fit$contributions$var == var.name, 2]
  # calculate the marginal effect of the variable using partial
  test_resp <- brt_fit %>%
    pdp::partial(pred.var = var.name, 
                 quantiles = T, probs = c(1:99/100),
                 plot=FALSE, train = brt_subset, n.trees = brt_fit$n.trees) 
  # if the variable is a factor, make a barplot and add it
  if(is.factor(var.data)) {
    test_sums = table(var.data)
    # make sure to save the barplot so you can access its x axis
    test_hist <- graphics::barplot(test_sums, col = "grey", 
                                   border = "white", las = 1, 
                                   ylab = ifelse(left_ax, "Density", ""),
                                   xlab = "")
    title(xlab = paste0(letters[j], ". ", 
                        mapvalues(var.name, 
                                  from = var.names.plotting$var.names, 
                                  to = var.names.plotting$plotting.names,
                                  warn_missing = F),
                        " (", round(var.import,1), "%)"),
          line = 3.5)
    par(new = TRUE)
    # Use the barplot's xaxis with the marginal effects calculated before
    plot(test_hist, 
         test_resp$yhat,
         xlab = "", ylab = "", axes = FALSE, bty = "n",
         type = "l")
    # Add left axis ticks, and a title if its on the far left
    axis(side=4, at = pretty(range(test_resp$yhat)), las = 1)
    mtext(ifelse(right_ax, "Marginal Effect", ""), side=4, line=4, cex = 0.75)
    box()
  } else { # If its not a factor, we'll make a histogram instead
    resp_x = test_resp[,1]
    if(var.name == "Tair_mean"){ # For air temperature, lets convert to the more understandable celsius
      var.data = var.data - 273.15
      resp_x = resp_x - 273.15}
    if(var.name == "pop_density"){ # For population density visibility purposes, use logs
      var.data = log(var.data)
      resp_x = log(resp_x)}
    # Save the histogram for use later, then plot it
    test_hist <- graphics::hist(var.data,
                                plot = F)
    plot(test_hist, freq = F, col = "grey", border = "white", 
         las = 1, 
         ylab = ifelse(left_ax, "Density", ""), xlab = "",
         main = NULL)
    title(xlab = paste0(letters[j], ". ", 
                        mapvalues(var.name, 
                                  from = var.names.plotting$var.names, 
                                  to = var.names.plotting$plotting.names,
                                  warn_missing = F), " (", 
                        round(var.import,1), "%)"),
          line = 3.5)
    par(new = TRUE)
    # Add the marginal effects
    plot(resp_x,
         test_resp$yhat,
         xlab = "", ylab = "", axes = FALSE, bty = "n",
         type = "l")
    # Add axis ticks, and an axis title if its on the far left
    axis(side=4, at = pretty(range(test_resp$yhat)), las = 1)
    mtext(ifelse(right_ax, "Marginal Effect", ""), side=4, line=4, cex = 0.75)
    box()
  }
}
dev.off()


# Figure 6: Lets look at the extended environmental risk in the municipalities ----
# where spillover happened in 2016-2018 and plot case data where available 
sub = readRDS("./data/cleaned/all_data.rds") %>% 
  base::subset(region == "Southeast")
# Add a state column and make it a ordered factor
sub$state = mapvalues(sub$state, 
                      from = c("RJ", "MG", "ES", "SP"),
                      to = c("Rio de Janeiro",
                             "Minas Gerais", 
                             "Espirito Santo",
                             "Sao Paulo"))
sub$state = factor(sub$state, levels = c("Minas Gerais", 
                                         "Espirito Santo",
                                         "Sao Paulo",
                                         "Rio de Janeiro"))
# For each state, calculate average risk over the municipalities for each time
sub_calc = ddply(sub, .(seq_ord, state), 
                 summarise_at, .vars = vars("max_EnvRisk"), 
                 .funs = c(mean), na.rm = T)
# Isolate the municipality-months with spillover to use in plotting 
sub_spillovers = subset(sub, spillover > 0)

# Now lets plot the risk metrics, faceting by region
Fig6 <- ggplot(sub_calc) + 
  geom_line(aes(y = max_EnvRisk, x = (seq_ord-1)/12 + 2000, group = CD_GEOCMU), 
            color = "grey", alpha = 0.5,
            data = sub) + # For every municipality, plot the risk metric over time
  geom_line(aes(x = (seq_ord-1)/12 + 2000, y = max_EnvRisk),
            data = sub_calc) + # For every region, plot the average of municipality risks
  geom_rect(data=data.frame(xmin=2016 + 11/12, xmax=2017 + 6/12, 
                            ymin=0, ymax=Inf),
            aes(xmin=xmin, xmax =xmax , ymin = ymin, ymax = ymax),
            color = NA, fill = "red",  alpha = 0.3) + # Add a box over the 2016-2017 outbreak
  geom_rect(data=data.frame(xmin=2017 + 10/12, xmax=2018 + 4/12, 
                            ymin=0, ymax=Inf),
            aes(xmin=xmin, xmax =xmax , ymin = ymin, ymax = ymax),
            color = NA, fill = "red", alpha = 0.3) + # And another one over the 2017 - 2018 outbreak
  geom_segment(x=2016+11/12,xend = 2016+11/12,
               y = 0, yend = 0.4, color = "blue",
               linetype = "dashed") + # Add a dashed line to show when our case data stop
  geom_point(aes(x = (seq_ord-1)/12 + 2000, y = max_EnvRisk), color = "red",
             data = sub_spillovers) + # Add points for all muni-months with spillover, showing the estimated risk when the spillover occured
  facet_wrap(~state) +
  coord_trans(y='sqrt', limx = c(2001, 2019), limy = c(0, 0.35)) + # Lets square-root transform the y-axis for visibility
  ylab("Environmental Risk") + 
  xlab("Year") + theme_bw()

# Now save the figure to a pdf
pdf("./output/results/Figure6.pdf", width = 9, height = 5)
print(Fig6)
dev.off()

# the end ----



