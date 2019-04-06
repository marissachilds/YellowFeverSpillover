# Clean and combine species occurences from different sources

# load libraries
library(reshape2)
library(mefa)
library(magrittr)
library(plyr)
library(dplyr)

# Importing and cleaning colandr data----
# which are species occurences extracted from published papers
# Import colandr data and subset to finished papers, then to the relevant columns 
colandr.data.raw <- read.csv(
  './data/raw/SDM/species_occurrence_data/colandr_mosquito_data.csv', 
  stringsAsFactors = F)
colandr.data <-subset(colandr.data.raw, 
                      data_extraction_screening_status=='finished')
colandr.data <- dplyr::select(colandr.data, c(citation_title, 
                                citation_authors,
                                citation_pub_year, 
                                Hg..janthinomys, 
                                Hg..Leucocelaenous, 
                                S..chloropterus,
                                location.format.))

# Reshape data to have one row for each species for each paper
colandr.data<-reshape2::melt(colandr.data, 
                             id.vars=c("citation_title","citation_authors",
                                       "citation_pub_year","location.format."))
names(colandr.data)<-c("Title","Authors","Pub_Year",
                       "Location_Format","Species","GPS_location")

# Then reshape data to have one row for each species in each 
# location in each paper, removing NA and data thats already in GBIF locations 
# (which occured when a species was not found in a paper)
temp <- rep(subset(colandr.data, select = -c(GPS_location)), 
            times = laply(strsplit(colandr.data$GPS_location, split = ";"), 
                          length))
colandr.data <-cbind.data.frame(
  temp, 
  GPS_location = unlist(strsplit(colandr.data$GPS_location, split = ";"))) %>%
  subset(GPS_location !="N/A") %>%
  subset(GPS_location !="gbif") %>% 
  mutate_if(is.factor, as.character) %>%
  mutate(Location_Format = as.factor(Location_Format))
colandr.data$Location_Format <- mapvalues(colandr.data$Location_Format, 
                                          from = c("d; /; m; /; s", "d; /; m", "d; e; c; i; m; a; l"), 
                                          to = c("dms", "dm", "decimal"))
colandr.data$Location_Format = as.character(colandr.data$Location_Format)

## Convert GPS locations to all be in decimal degrees
# Split up the coordinates (latitude and longitude parts of the coordinates are even spit up into seperate columns)
colandr.data$lat <- laply(strsplit(colandr.data$GPS_location, split = ","), function(x){x[[1]]})
colandr.data$long <- laply(strsplit(colandr.data$GPS_location, split = ","), function(x){x[[2]]})

# Construct a function that converts a lat or longitude into decimal degrees with +/- instead of E/W/N/S
convert_to_dd <- function(location_string, format = "decimal", coordinate){
  split_string <- strsplit(trimws(location_string), split = " ")[[1]]
  
  # for decimal formatted locations, extract decimal degrees and direction 
  if(format == "decimal"){
    dd <- gsub(pattern = "°", replacement = "", x = split_string[1]) %>% as.numeric
    direction <- split_string[2]
  }
  # for dm and dms formatted locations, extract direction and degrees, minutes, seconds then convert to decimal
  if(format == "dm" | format == "dms"){
    # extract direction and remove from the vector
    direction <- split_string[length(split_string)]
    split_string <- split_string[-length(split_string)]
    if (length(split_string) < 3){ # if its dm, add 0 seconds to we can treat them all the same
      split_string <- append(split_string, "0''")
    }
    # remove degree, minute, and second marks
    split_string <- gsub(pattern = "°|'", replacement = "", x = split_string) %>% as.numeric 
    # convert to decimal degrees
    dd <- (((split_string[3]/ 60) + split_string[2])/60) + split_string[1]
  }
  # flip the sign if its south or west
  if((coordinate == "latitude" & direction == "S") | 
     (coordinate == "longitude" & direction == "W") ){
    dd = -dd
  }
  return(dd)
}

# Apply the above function to calculate decimal degrees latitude and longitude
colandr.data$decimallatitude <- dplyr::select(colandr.data, c(lat, Location_Format)) %>% 
  transmute(location_string = lat, format = Location_Format) %>%
  maply(.fun = convert_to_dd,
        coordinate = "latitude",
        .expand = F)

colandr.data$decimallongitude <- dplyr::select(colandr.data, c(long, Location_Format)) %>% 
  transmute(location_string = long, format = Location_Format) %>%
  maply(.fun = convert_to_dd,
        coordinate = "longitude",
        .expand = F)

# Convert the species names
colandr.data$Species %<>% mapvalues(from = c("Hg..janthinomys", 
                                             "Hg..Leucocelaenous", 
                                             "S..chloropterus"), 
            to = c("Haemagogus janthinomys", 
                   "Haemagogus leucocelaenus", 
                   "Sabethes chloropterus"))

# Create and identifier from authors, publication year and publication title
colandr.data$id = paste(colandr.data$Authors, colandr.data$Pub_Year, colandr.data$Title)
colandr.data$source = "colandr"
  
# Importing and cleaning GBIF data ----
gbif.raw<-rbind(
  read.table('./data/raw/SDM/species_occurrence_data/Hg.janthinomys.gbif.csv',
             sep= "\t", header = T, stringsAsFactors = F),
  read.table('./data/raw/SDM/species_occurrence_data/Hg.leucoceleanus.gbif.csv',
             sep="\t", header = T, stringsAsFactors = F),
  read.table('./data/raw/SDM/species_occurrence_data/Sa.chloropterus.gbif.csv', 
             sep="\t", header = T, stringsAsFactors = F))

# subset to needed columns and add column for source
gbif.data <- select(gbif.raw, c(species, decimallatitude, decimallongitude, year, gbifid)) %>%
  mutate(source = "gbif")

# Combine the two data sources ----
all.data <- select(colandr.data, c(Species, decimallatitude, decimallongitude, Pub_Year, id, source)) %>%
  dplyr::rename(species = Species, 
                year = Pub_Year) %>%
    rbind(dplyr::rename(gbif.data, id = gbifid))

# Write out a csv file
write.csv(all.data, file = "./data/cleaned/SDM/species_occurrences.csv", row.names = F)
rm(colandr.data.raw, colandr.data, temp, convert_to_dd, gbif.raw, gbif.data, all.data)



