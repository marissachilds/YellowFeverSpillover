library(rgdal)
library(magrittr)

# read in IUCN shapefile
primates <- readOGR(dsn = "./data/raw/Primates/TERRESTRIAL_MAMMALS", layer = "TERRESTRIAL_MAMMALS")
primates %<>% subset(order_ == "PRIMATES")
YF_primates = subset(primates, genus %in% c("Ateles", "Aotus", "Alouatta", "Saimiri", 
                                            "Cebus", "Callicebus", "Callithrix", 
                                            "Saguinus", "Lagothrix"))

writeOGR(YF_primates, "data/cleaned/primates", "primates", driver = "ESRI Shapefile")