# This function clean the yellow fever case data. Most useful when called within 
# another R script

# Set up libraries ----
library(sp)
library(plyr)
library(dplyr)
library(magrittr)
library(rgdal)
old_wd = getwd()

# Set the system locale so R can read the letters in the files
Sys.setlocale(category = "LC_ALL", locale = "C")

# Define a function to load the data based on the disease and year, this works for dengue as well -----
# returns the data in long format with a columns for year, disease, month, municipality, and number of cases
load_data <- function(disease, year){
  abbr = ifelse(disease == "dengue", "DENV", "YFV")
  filename <- paste0('data/raw/disease_data_by_month/', disease, 
                     '/download_by_year_of_first_symptoms/',abbr, 
                     '_Municipality_infected_', year, '.csv')
  data <- read.csv(filename, sep = ';', header = T, quote = "", skip = 4)
  end1 <- which(data[,1] == '\"Total\"') - 1
  end2 <- ncol(data) - 1
  data <- data[1:end1, 1:end2]
  colnames(data)[1] <- "municipality.infected"
  data %<>% reshape(direction = 'long', varying = colnames(data)[-1], 
                    idvar = colnames(data)[1], timevar = "month") %>% 
    cbind(year = year) %>% cbind(disease = disease)
}

# Run the load_data function for all the years where we have yellow fever cases -----
# and put it together in a dataframe
YFVyears <- cbind(disease = 'yellow_fever', c(2001:2010, 2013:2016))
case.data <- mdply(YFVyears, load_data, .id = NULL)[,-1]
colnames(case.data)[3] <- "cases"

# Clean up the case data data frame and order the months
case.data %<>% mutate(cases = as.integer(as.character(cases))) %>% 
  mutate(municipality.infected = as.character(municipality.infected)) %>% 
  mutate(code = NA) %>%
  mutate(month = factor(month, 
                        levels = c('Jan', "Fev", "Mar", "Abr", "Mai", "Jun", 
                                   "Jul", "Ago", "Set", "Out", "Nov", "Dez"), 
                        ordered = T))

# NAs in the case.data are from "-"s which are zeros
case.data[which(is.na(case.data), arr.ind = T)] <- 0

# read in municipality shapefile and dataset with municipality names paired with municipality codes
shape = readOGR(dsn = 'data/raw/brazil_border_shapefiles/br_municipios')
levels(shape$NM_MUNICIP) %<>% c("EXTERIOR")
levels(shape$CD_GEOCMU) %<>% c('0000000')
temp = rbind(shape@data, c("EXTERIOR", '0000000'))

# define state codes
state.dictionary = data.frame(state.abbr = c('RO', 'AC', 'AM', 'RR', 'PA', 'AP', 
                                             'TO', 'MA', 'PI', 'CE', 'RN', 'PB', 
                                             'PE', 'AL', 'SE', 'BA', 'MG', 'ES', 
                                             'RJ', 'SP', 'PR', 'SC', 'RS', 'MS',
                                             'MT', 'GO', 'DF'),
                              state.code = c(11, 12, 13, 14, 15, 16, 17, 21,
                                             22, 23, 24, 25, 26, 27, 28, 29, 
                                             31, 32, 33, 35, 41, 42, 43, 50,
                                             51, 52, 53) %>% as.character)

# Split municipality names from codes in the data and match the 
# municipality codes with their full version in the municipalities shapefile
whole.names <- case.data$municipality.infected %>% strsplit(split = ' ')
codes <- sapply(whole.names, `[`, 1)
codes %<>% substr(2,nchar(codes))
#names <- sapply(whole.names, `[`, 2)
#case.data$municipality.infected <- substr(names, 1, nchar(names)-1)
case.data$code <- temp$CD_GEOCMU[pmatch(codes, temp$CD_GEOCMU, duplicates.ok = T)]
case.data$state.code <- substr(codes, 1,2)
case.data %<>% left_join(state.dictionary)

# Remove everything unnecessary and return to the original working directory
rm(temp, whole.names, codes, state.dictionary, load_data, YFVyears)
setwd(old_wd)
