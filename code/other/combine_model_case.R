# This code combines the many sources of data to collect the model estimates, 
# aggregated to the municipality-month, municipality-month covariates, and case
# data for use in model-data comparisons and BRT

# Setup -----
library(plyr)
library(dplyr)
library(raster)
library(magrittr)
library(pracma)
library(ggplot2)
library(stringr)
library(rgeos)
library(rgdal)

# Read in case data and limit to municipalities in Brazil, 
# so ones with nonzero municipality code
source("./code/other/clean_case_data.R")
case.data = case.data[case.data$code != "0", ]

# read in the various model estimates, extract their type from the file name, 
# reshape them to be long form and save to a list ----
model_ests <- alply(list.files("./output/municipality_model_estimates", 
                               pattern = ".csv", full.names = T), 1, 
                    function(path){
                      data <- read.csv(path, stringsAsFactors = F) %>% 
                        dplyr::select(-c(system.index, .geo))
                      file_name = substring(path, 39, nchar(path)-18)
                      summary_type = strsplit(file_name, split = "_")[[1]][1]
                      est_type = strsplit(file_name, split = "_")[[1]][2]
                      str_lengths = nchar(setdiff(colnames(data),
                                                  c("CD_GEOCMU", "NM_MUNICIP")))
                      data %<>% reshape(idvar = "CD_GEOCMU", 
                                        varying = setdiff(colnames(data),
                                                          c("CD_GEOCMU", "NM_MUNICIP")), 
                                        v.names = paste(summary_type, 
                                                        est_type, sep = "_"), 
                                        times = substring(setdiff(colnames(data),
                                                                  c("CD_GEOCMU", "NM_MUNICIP")), 
                                                          first = str_lengths - 6), 
                                        direction = "long")
                      return(data)
                    })

# join all of the model estimates together, then add columns for year and month
all_model_ests = join(model_ests[[1]], model_ests[[2]], type = "full") %>%
  join(model_ests[[3]], type = "full") %>% 
  join(model_ests[[4]], type = "full") %>% 
  join(model_ests[[5]], type = "full") %>% 
  join(model_ests[[6]], type = "full") %>% 
  join(model_ests[[7]], type = "full") %>% 
  join(model_ests[[8]], type = "full") 
rm(model_ests)

# Split the times variable into year and month
time_splits = strsplit(all_model_ests$time, split = "_")
nsplits = sapply(time_splits, length)
all_model_ests$time = gsub("^\\_","", all_model_ests$time)
all_model_ests$month = strsplit(all_model_ests$time, split = "_") %>% 
  sapply(`[`,2) %>% as.numeric
all_model_ests$year = strsplit(all_model_ests$time, split = "_") %>% 
  sapply(`[`,1) %>% as.numeric

# merge case data and model estimates, first making datasets comparable
case.data %<>% dplyr::select(-c(municipality.infected, state.code, 
                                state.abbr, disease)) %>%
  mutate(month = as.numeric(factor(month, 
                                   levels = c("Jan", "Fev", "Mar", "Abr", 
                                              "Mai", "Jun", "Jul", "Ago", 
                                              "Set", "Out", "Nov", "Dez"), 
                                   ordered = T)),
         year = as.numeric(levels(year))[year],
         code = as.numeric(levels(code))[code])

model_case <- full_join(all_model_ests, case.data, by = c("month" = "month",
                                                          "year" = "year",
                                                          "CD_GEOCMU" = "code"))
# make the NA cases zeros
model_case$cases[which(is.na(model_case$cases))] = 0

# Read in population data ----
Sys.setlocale(category = "LC_ALL", locale = "C")
pop_data = read.csv("./data/cleaned/BRT_data/population_2000_2016.csv", 
                    sep = ";", header = F, quote = "", skip = 4, 
                    stringsAsFactors = F)
# Exclude lines after the totals
end <- which(pop_data[,1] == '\"Total\"') - 1
pop_data = pop_data[1:end,]
# Name the columns of the population data frame
colnames(pop_data) = c("name", paste("population", 2000:2016, sep = "."))

# Split the current name variable into code and name
pop_data$code = pop_data$name %>% strsplit(split = " ") %>% 
  sapply("[", 1) %>% substring(first = 2) 
pop_data$name = str_split(pop_data$name, pattern = " ", n = 2) %>% 
  sapply("[", 2)
pop_data %<>% dplyr::select(code, name, everything())


# Estimate 2017 population as 2016 population plus growth from 2015 to 2016
pop_data$population.2017 = (as.numeric(pop_data$population.2016) + 
                              (as.numeric(pop_data$population.2016) - 
                                 as.numeric(pop_data$population.2015))) %>%
  as.character
# Estimate 2018 population as 2016 population plus 2*growth from 2015 to 2016
pop_data$population.2018 = (as.numeric(pop_data$population.2016) + 
                              2*(as.numeric(pop_data$population.2016) - 
                                   as.numeric(pop_data$population.2015))) %>%
  as.character

# convert population data to long
pop_data_long = reshape(pop_data, direction = "long", 
                        varying = colnames(pop_data)[-(1:2)], 
                        idvar = colnames(pop_data)[1:2],
                        timevar = "year")
pop_data_long$population = as.numeric(pop_data_long$population)

# read in fire data ----
fire = read.csv("./data/cleaned/BRT_data/municipality_fireArea_200111_201809.csv", 
                stringsAsFactors = F) %>% 
  dplyr::select(-c("system.index", ".geo"))

# Name the columns then reshape to long
colnames(fire) = c("code", "name", 
                   paste("fireArea", substring(colnames(fire)[-(1:2)], 6, 12), 
                         sep = "."))
fire %<>% reshape(direction = "long", varying = colnames(fire)[-(1:2)],
                  idvar = colnames(fire)[1:2], timevar = "date")

# Use the date column to make month and year columns
fire$month = substr(fire$date, 6, 7) %>% as.numeric
fire$year = substr(fire$date, 1, 4) %>% as.numeric

# Read in municipality area data
muniArea = read.csv("./data/cleaned/BRT_data/municipality_Area.csv", 
                    stringsAsFactors = F) %>%
  dplyr::select(-c("system.index", ".geo"))
colnames(muniArea) = c("code", "name", "total_area")

# join fire area with municipality area data and calculate percent of 
# municipality with fires
fire %<>% left_join(muniArea) %>%
  mutate(fire_pct = fireArea/total_area) %>% 
  mutate(code = as.character(code))

# Read in primate data extracted from GEE----
primates_mean = read.csv("data/cleaned/BRT_data/primate_species_mean.csv", 
                         stringsAsFactors = F) %>% 
  dplyr::select(-c("system.index", ".geo")) 
colnames(primates_mean) = c("code", "name", "primates_mean")
primates_max = read.csv("data/cleaned/BRT_data/primate_species_max.csv", 
                        stringsAsFactors = F) %>% 
  dplyr::select(-c("system.index", ".geo")) 
colnames(primates_max) = c("code", "name", "primates_max")

# Read in temperature data and reshape to long format ----
temperature = read.csv("./data/cleaned/BRT_data/mean_temperature_200002_201612.csv",
                       stringsAsFactors = F) %>% 
  dplyr::select(-c("system.index", ".geo"))
temperature %<>% reshape(idvar = "CD_GEOCMU", 
                         varying = setdiff(colnames(temperature),
                                           c("CD_GEOCMU", "NM_MUNICIP")), 
                         v.names = "Tair_mean",
                         times = substring(
                           setdiff(colnames(temperature),
                                   c("CD_GEOCMU", "NM_MUNICIP")), 
                           first = 11),
                         direction = "long")

# Make year and month columns
temperature$year = strsplit(temperature$time, split = "_") %>% sapply("[[", 1) %>% as.numeric
temperature$month = strsplit(temperature$time, split = "_") %>% sapply("[[", 2) %>% as.numeric
temperature %<>% dplyr::select(-c(time, NM_MUNICIP))

# Read in precipitation data and reshape to long----
precip = read.csv("./data/cleaned/BRT_data/mean_precipitation_200001_201712.csv",
                  stringsAsFactors = F) %>% 
  dplyr::select(-c("system.index", ".geo"))
precip %<>% reshape(idvar = "CD_GEOCMU", 
                    varying = setdiff(colnames(precip), 
                                      c("CD_GEOCMU", "NM_MUNICIP")), 
                    v.names = "precipitation",
                    times = substring(setdiff(colnames(precip),
                                              c("CD_GEOCMU", "NM_MUNICIP")), 
                                      7, 12),
                    direction = "long")
# Add year and month columns
precip$year = substring(precip$time, 1, 4) %>% as.numeric
precip$month = substring(precip$time, 5, 6) %>% as.numeric
precip %<>% dplyr::select(-c(time, NM_MUNICIP))

# Read in vaccine coverage estimates, which is taken as modal coverage rate in the municipality from the raster from Freya Shearer----
vaccine_est = read.csv("./data/cleaned/BRT_data/municipality_vax_estimates_2001_2016.csv", 
                       stringsAsFactors = F) %>% dplyr::select(-c("system.index", ".geo"))
# Reshape to long
vaccine_est %<>% reshape(idvar = "CD_GEOCMU", 
                    varying = setdiff(colnames(vaccine_est),
                                      c("CD_GEOCMU", "NM_MUNICIP")), 
                    v.names = "vax",
                    times = substring(setdiff(colnames(vaccine_est),
                                              c("CD_GEOCMU", "NM_MUNICIP")), 
                                      27, 30),
                    direction = "long") %>% 
  dplyr::select(-NM_MUNICIP) %>% 
  mutate(time = as.numeric(time))

# join fire data, then primates data, then population data, then temperature, then precip ----
df_brt = model_case %>% mutate(code6 = substr(CD_GEOCMU, 1, 6),
                               spillover = (cases > 0)*1)

df_brt %<>% mutate(CD_GEOCMU = as.character(CD_GEOCMU))
df_brt %<>% full_join(dplyr::select(fire, -c(name)), 
                      by = c("CD_GEOCMU" = "code", 
                             "year" = "year",
                             "month" = "month")) %>% 
  full_join(dplyr::select(primates_mean, -c(name)) %>% 
              mutate(code = as.character(code)), 
            by = c("CD_GEOCMU" = "code")) %>% 
  full_join(dplyr::select(primates_max, -c(name)) %>% 
              mutate(code = as.character(code)), 
            by = c("CD_GEOCMU" = "code")) %>% 
  full_join(dplyr::select(pop_data_long, -c(name)), 
            by = c("code6" = "code", 
                   "year" = "year")) %>%   
  full_join(vaccine_est %>% 
              mutate(CD_GEOCMU = as.character(CD_GEOCMU)),
            by = c("CD_GEOCMU" = "CD_GEOCMU", 
                   "year" = "time")) %>% 
  full_join(temperature %>% 
              mutate(CD_GEOCMU = as.character(CD_GEOCMU)),
            by = c("CD_GEOCMU" = "CD_GEOCMU", 
                   "year" = "year",
                   "month" = "month")) %>% 
  full_join(precip %>% 
              mutate(CD_GEOCMU = as.character(CD_GEOCMU)),
            by = c("CD_GEOCMU" = "CD_GEOCMU", 
                   "year" = "year",
                   "month" = "month")) %>%
  mutate(CD_GEOCMU = as.character(CD_GEOCMU)) %>%
  mutate(pop_density = population/total_area)  

# Make a state variable by taking the first two characters of the municipality code
df_brt$state = substr(df_brt$CD_GEOCMU, 1, 2)
df_brt$state = mapvalues(df_brt$state, 
                         from = c( 11,   12,   13,   14,   15,   
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

# Make a region variable based on states in regions
df_brt$region = mapvalues(df_brt$state, 
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
# Remove the data with a zero municipality code
df_brt = df_brt[-which(df_brt$CD_GEOCMU == "0"),]

# to add lagged environmental risk and fire area, first add a sequence number,
# then iterate over each sequence to add lagged data
df_brt %<>% mutate(seq_ord = month + 12*(year-2000))
df_brt %<>% ddply(.(CD_GEOCMU), function(muni_sub){
  n = nrow(muni_sub)
  muni_sub = muni_sub[order(muni_sub$seq_ord),]
  muni_sub$lagged_max_EnvRisk = c(NA, muni_sub$max_EnvRisk[-n])
  muni_sub$lagged_fireArea = c(NA, muni_sub$fireArea[-n])
  return(muni_sub)
})

df_brt %<>% mutate(month = as.factor(month))
df_brt %<>% mutate(region = as.factor(region))

# Calculate the phenomenological primate curve from the parameter estimates
df_brt$phenom = -0.190384045495927*
  sin((as.numeric(as.character(df_brt$month))/12 + df_brt$year)*2*pi/7) +
  0.462335197559909*
  cos((as.numeric(as.character(df_brt$month))/12 + df_brt$year)*2*pi/7) + 1/2

# Save full data
saveRDS(df_brt, "./data/cleaned/all_data.rds")

# Save BRT data as all data between 2001 and 2016
brt_subset = df_brt[df_brt$year >= 2001 & df_brt$year <= 2016, ]
saveRDS(brt_subset, "./data/cleaned/BRT_data/brt_data.rds")
