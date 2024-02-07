##
## Date Created: 2024-02-07
##
## Copyright (c) Rachel Roday, 2024
##
##      This code is considered intellectual property and any reproduction of 
##      of figures or analysis from the following script is a violation of US Copyright laws
##      "Under copyright law, source code is a literary work (like a book). 
##      And, just like any other writing, it is immediately copyrighted 
##      regardless of the author registering it with the U.S. Copyright Office."
##
## Email: rroday@udel.edu
##
##        
##
## ---------------------------
##
## Notes: This is the Final Script for Roday et al. 2024
##
##
## ---------------------------

########################  Working Dir and Packages #############################

# Set WD
setwd("C:/Users/RER/Documents/Masters UD/Code")


# Load required libraries
library(ggplot2)
library(glatos)
library(data.table) #setnames
library(ggplot2)
library(dplyr)
library(tidyverse)
library(lubridate)
library(network)
library(igraph)
library(ggraph)
library(tidygraph)
library(viridis)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(ggrepel)
library(maps)
library(sf)
library(mapdata)
library(gganimate) 
library(gifski)
library(plotrix) 
library(ggforce)
library(lawstat)
library(epitools)
library(vcd)
library(FSA)
library(ggsignif)
library(ggpubr)
library(tidyr)
library(dataRetrieval)
library(ggtext)
library(scales)
library(rcompanion)
library(tm)
library(cowplot)
library(rstatix)

########################  Fxns and Fish data frame ###############################

# Create required functions and dataframes

#make new function to extract second element from a hyphen-delimited string
#apply get_rsn() to each record in Receiver column and create a new column 
#named receiver_sn
get_rsn <- function(x) strsplit(x, "-")[[1]][2]

#NOTE that lat and long are release lat and long from 
rcv <- data.frame(
  glatos_array = c("BR1", "BR2", "CR1", "RC2","WC1", "CR2","RC1","WC2"),
  station = c("1","2"), 
  deploy_lat = c(39.75756, 39.73994, 39.68424, 39.7339, 39.70129, 39.68424, 39.71467, 39.69822),
  deploy_long = c(-75.55351, -75.5371, -75.6342, -75.63419, -75.65033, -75.63474, -75.64003, -75.66525),
  ins_serial_no = c("486334", "486336", "490249", "490250", "490263", "490252", "490251", "486335"),
  stringsAsFactors = FALSE) 



## Extract transmitter code and ID - make a new function to extract id from Transmitter
parse_tid <- function(x) strsplit(x, "-")[[1]][3]

#make a new function to extract codespace from Transmitter
parse_tcs <- function(x) {
  tx <- strsplit(x, "-")[[1]][1:2]
  return(paste(tx[1:2], collapse = "-"))} 


# Make Fish data frame
# This will change from year to year, as more fish are tagged. 
fsh <- data.frame(animal_id = c(13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 
                1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), 
  tag_code_space = c("A69-1604","A69-1604","A69-1604","A69-1604","A69-1604",
                     "A69-1604","A69-1604","A69-1604","A69-1604","A69-1604",
                     "A69-1604","A69-1604","A69-1604","A69-1604",
                     "A69-9007","A69-9007","A69-1602","A69-1602","A69-1602","A69-1602","A69-1602",
                     "A69-1602","A69-1602","A69-1602","A69-1602","A69-1602"),
  tag_id_code = c("32457","32458","32459","32460","32461","32462","32463",
                  "32464","32465","32466","32467","32468","32469","32470",
                  "15557","15558","49176", "49177","49178","49179","49180",
                  "49181","49182","49183","49184","49185"), 
  common_name = "American Shad", 
  release_date_time = as.POSIXct(c("2022-05-10 10:00:00", "2022-05-10 10:00:00",
                                   "2022-05-10 10:00:00", "2022-05-10 10:00:00",
                                   "2022-05-10 10:00:00", "2022-05-10 10:00:00",
                                   "2022-05-10 10:00:00", "2022-05-10 10:00:00",
                                   "2022-05-10 10:00:00", "2022-05-10 10:00:00", 
                                   "2022-05-10 10:00:00", "2022-05-10 10:00:00", 
                                   "2022-05-10 10:00:00", "2022-05-10 10:00:00",
                                   "2021-05-11 10:00:00", "2021-05-11 10:00:00",
                                   "2021-05-11 10:00:00", "2021-05-11 10:00:00",
                                   "2021-05-11 10:00:00", "2021-05-11 10:00:00",
                                   "2021-05-11 10:00:00", "2021-05-11 10:00:00",
                                   "2021-05-11 10:00:00", "2021-05-11 10:00:00",
                                   "2021-05-11 10:00:00", "2021-05-11 10:00:00"), 
                                 tz = "UTC"),
  recapture_date_time = as.POSIXct(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                                     NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                                     NA, NA, NA, NA, NA ),
                                   tz = "UTC"), 
  Sex = c("Female","Female","Female","Female","Female", "Male", "Female", "Male",
          "Male","Female","Female","Female", "Male", "Male","Male","Male","Female", "Male", "Female", "Male","Male","Female",
          "Male","Female","Female","Female"),
  Insertion_type = c("Gastric", "Gastric", "Surgical", "Surgical","Surgical",
                     "Surgical","Surgical","Surgical","Surgical","Surgical",
                     "Surgical", "Surgical","Surgical","Surgical", "Surgical", "Surgical","Surgical",
                     "Surgical","Surgical","Surgical","Surgical","Surgical", "Surgical", "Surgical","Surgical","Surgical"),
  Length.m = c(0.495, 0.540, 0.54, 0.545, 0.495, 0.495, 0.55, 0.52, 0.465, 0.530, 0.515, 0.555, 0.490, 0.455,
               0.472,0.493,0.593,0.484,0.522,0.493,0.475,0.517,0.446,0.540,0.521,0.525), 
  Fork.mm = c(440, 470, 490, 510, 445, 450, 495, 460, 425, 480, 470, 490, 445, 420,
              420, 437, 450, 443, 468, 437, 420, 460, 398, 480, 458, 468),
  Deploy.Year = c("2022", "2022","2022", "2022","2022", "2022","2022", "2022",
                  "2022", "2022","2022", "2022","2022", "2022","2021","2021","2021","2021","2021","2021",
                  "2021","2021","2021","2021","2021","2021"),
  stringsAsFactors = FALSE)  


# Create color palette to use for whole paper
colpal <- c("Atlantic" = "#4477AA",  "Bay" = "#FD5C70", "Brandywine 1" = "#228833", 
            "Brandywine 2" = "#AA3377", "Brandywine 3" = "#66CCEE",  "River" = "#E69F00", "#000000",
            "Brandywine" = "#8B79C3", "Canal" = "#000000",
            "Atlantic Ocean" = "#4477AA",  "Delaware Bay" = "#FD5C70", "Delaware River" = "#E69F00",
            "C&D Canal" = "#000000", "Brandywine " = "#8B79C3", "Brandywine River" = "#8B79C3")

# Creating 5 day latenent mortality window deemed survival gap
survival_gap <- as.duration(hms("120:00:00"))

# Creating shape files for country/world
world <- ne_coastline(scale = "medium", returnclass = "sf")
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
canada <- map_data("worldHires", "Canada") 

########################  Loading  data #########################################
# Go! Get! Your! Data

# Loading Data from Station BR1 (Receiver VR2Tx_486334)
BR1_2022 <- read.csv("VR2Tx_486334_20230119_2.csv", as.is = TRUE, check.names = FALSE, fileEncoding = "UTF-8-BOM")
BR1_2022$`Station Name` <- as.character(BR1_2022$`Station Name`)

# Loading Data from Ian Park -- All Reciever data with Ed's fish from 2021 and 2022 (Minus BR1)
# AllDetects_2021_2022 <- read.csv("Shad_2021_and_2022_data.csv",as.is = TRUE, 
#                                  check.names = FALSE, fileEncoding = "UTF-8-BOM")
AllDetects_2021_2022 <- read.csv("Copy of EdShad_FROMIANPARK.csv",as.is = TRUE, 
                                 check.names = FALSE, fileEncoding = "UTF-8-BOM")


# Loading more data from Ian Park AND all MATOS data
#MATOS <- read.csv("Acoustic Data from Ian and Dec 2021 MATOS Data Combined.xlsx - EdsShad.csv",
#as.is = TRUE, check.names = FALSE, fileEncoding = "UTF-8-BOM")
MATOS <- read.csv("Acoustic Data from Ian and Dec 2021 MATOS Data Combined THIS ONE.csv",
                  as.is = T, check.names =  F, fileEncoding = "UTF-8-BOM", stringsAsFactors = F)

## This is a data frame  that Ian gave us for receiver names 
sampling_stations <- read.csv("Ian Receiver Names and Locations 2021.csv", 
                              header = TRUE)

## This is another data frame from Ian that incldues river KM for all stations.
# yes there are 5 excel sheets that need to be uploaded. I dont make the rules.
RiverKM <- read.csv("Receiver RKM.csv", header = TRUE) %>% rename ("StationName" = "Station.Name")



dewayne <- read.csv("VUE_Export_Dewayne_Nov2023.csv",as.is = TRUE, 
                    check.names = FALSE, fileEncoding = "UTF-8-BOM")
ian <- read.csv("IAN2023.csv", as.is = TRUE, 
                check.names = FALSE, fileEncoding = "UTF-8-BOM") %>%
  dplyr::select(detection_timestamp_utc:`Sensor Precision`)

dewayneBatches <- read.csv("DewayneFullBatches.csv", as.is = T, check.names = F, 
                           fileEncoding = "UTF-8-BOM") %>%
  dplyr::select(detection_timestamp_utc:`Sensor Precision`)

########################  Tidying  data #########################################

# Full station name df
sampling_stations2 <- full_join(sampling_stations, RiverKM, by=c("StationName", "RKM")) %>%
  filter(!ID %in% c(561, 67,552,565,184,524,527,43,529,531,532,538,539, 522))
# Removed duplicate stations (stations that appeared more than once. Defaulted to higher RKM value for conservative distance measurments)
sampling_stations2$Y[sampling_stations2$StationName == "Christina River"] <- round(39.73251,5)
sampling_stations2$X[sampling_stations2$StationName == "LL# 2950 Delaware River Lighted Buoy 25"] <- round(-75.54714,5)
sampling_stations2$Y[sampling_stations2$StationName == "LL# 2950 Delaware River Lighted Buoy 25"] <- round(39.65608,5)
sampling_stations2$X[sampling_stations2$StationName == "LL# 2970 Delaware River Lighted Buoy 28"] <- -75.52269
sampling_stations2$Y[sampling_stations2$StationName == "LL# 2970 Delaware River Lighted Buoy 28"] <- 39.67289
sampling_stations2$X[sampling_stations2$StationName == "LL# 2995 Delaware River Lighted Buoy 29"] <- -75.51460
sampling_stations2$Y[sampling_stations2$StationName == "LL# 2995 Delaware River Lighted Buoy 29"] <- 39.70297
sampling_stations2$X[sampling_stations2$StationName == "LL# 3205 Marcus Hook Anchorage Buoy C"] <- -75.38810
sampling_stations2$Y[sampling_stations2$StationName == "LL# 3205 Marcus Hook Anchorage Buoy C"] <- 39.81054
sampling_stations2$X[sampling_stations2$StationName == "LL# 3505 Delaware River Lighted Buoy 68"] <- -75.16659
sampling_stations2$Y[sampling_stations2$StationName == "LL# 3505 Delaware River Lighted Buoy 68"] <- 39.88011
sampling_stations2$Region[sampling_stations2$StationName == "LL# 3000 Delaware River Buoy 56"] <- "River"



# Change the column name from Date and Time (UTC) to detection_timestamp_utc
# setnames(x, old, new)
setnames(AllDetects_2021_2022, "Date and Time (UTC)", "detection_timestamp_utc") 
setnames(BR1_2022, "Date and Time (UTC)", "detection_timestamp_utc")
setnames(dewayne, "Date and Time (UTC)", "detection_timestamp_utc")
setnames(MATOS, "Date and Time (UTC)", "detection_timestamp_utc")

# Parse date column
AllDetects_2021_2022$detection_timestamp_utc <- as.character(mdy_hms(AllDetects_2021_2022$detection_timestamp_utc))
MATOS$detection_timestamp_utc <- as.character(mdy_hms(MATOS$detection_timestamp_utc))

# Merge all datasets
dtc <- full_join(BR1_2022, AllDetects_2021_2022)
dtc1 <- full_join(dtc, MATOS)
dtc1.2 <- full_join(dtc1, ian)
dtc1.3 <- full_join(dtc1.2, dewayneBatches)
dtc2 <- full_join(dtc1.3, dewayne)


# Parsing the serial number of receiver 
dtc2$receiver_sn <- sapply(dtc2$Receiver, get_rsn)

# Changing  the names in the dtc2 dataframe so that its easier to work with
setnames(dtc2, c("Sensor Value", "Sensor Unit", "Station Name"),  c("sensor_value", "sensor_unit", "station_name"))

# Adding Lat/Long for stations 1 and 2 into the dtc2 dataframe so that Full Join is possible
dtc2$Longitude[dtc2$station_name == "2"] <- -75.5371
dtc2$Latitude[dtc2$station_name == "2"] <- 39.73994
dtc2$Longitude[dtc2$station_name == "1"] <- -75.55351
dtc2$Latitude[dtc2$station_name == "1"] <- 39.75756

#full join on receiver serial number to add receiver data to detections
# this isnt crucial, but its good to have it 
dtc3 <- full_join(dtc2, rcv, by = c("receiver_sn" = "ins_serial_no", "station_name" = "station",
                                    "Latitude" = "deploy_lat", "Longitude" = "deploy_long"))

#apply parse_tcs() to Transmitter and assign to transmitter_codespace
dtc3$transmitter_codespace <- sapply(dtc3$Transmitter, parse_tcs)

#apply parse_tid() to Transmitter and assign to transmitter_id
dtc3$transmitter_id <- sapply(dtc3$Transmitter, parse_tid)

## This join is CRUCIAl -> gets us the fish number by transmitter id space
#simple left join on codespace and id
dtc_final <- full_join(dtc3, fsh, by = c("transmitter_codespace" = "tag_code_space",
                                               "transmitter_id" = "tag_id_code"))

# Removing all rows were there is no Fish ID
dtc_final <- dtc_final[!is.na(dtc_final$transmitter_id),]
dtc_final <- dtc_final[!is.na(dtc_final$animal_id),]
dtc_final <- dtc_final[!is.na(dtc_final$detection_timestamp_utc),]

# Add a year column and filter out year 2020 and random receiver name 
dtc_final1 <- dtc_final %>% 
  mutate(Year = year(as.Date(detection_timestamp_utc))) %>%
  filter(!Receiver == "PROJ162",
         !Year == 2020)





## Cherry Island and New Castle were missing Lat/long and names
dtc_final1$station_name[dtc_final1$Receiver == "VR2AR-547881"] <- "Cherry Island B"
dtc_final1$Longitude[dtc_final1$station_name == "Cherry Island B"] <- -75.50951
dtc_final1$Latitude[dtc_final1$station_name == "Cherry Island B"] <- 39.70142 
dtc_final1$station_name[dtc_final1$Receiver == "VR2AR-547883"] <- "Cherry Island Range A"
dtc_final1$Longitude[dtc_final1$station_name == "Cherry Island Range A"] <- -75.51718
dtc_final1$Latitude[dtc_final1$station_name == "Cherry Island Range A"] <- 39.70366  
dtc_final1$station_name[dtc_final1$Receiver == "VR2AR-547882"] <- "New Castle Flats C"
dtc_final1$Longitude[dtc_final1$station_name == "New Castle Flats C"] <- -75.51934
dtc_final1$Latitude[dtc_final1$station_name == "New Castle Flats C"] <- 39.6712 
dtc_final1$Longitude[dtc_final1$station_name == "C&D Canal West"] <- -75.64542
dtc_final1$Latitude[dtc_final1$station_name == "C&D Canal West"] <- 39.55535
dtc_final1$Longitude[dtc_final1$station_name == "Christina River"] <- -75.53233
dtc_final1$Latitude[dtc_final1$station_name == "Christina River"] <- round(39.73251,5)
dtc_final1$Longitude[dtc_final1$station_name == "LL# 2950 Delaware River Lighted Buoy 25"] <- round(-75.54714,5)
dtc_final1$Latitude[dtc_final1$station_name == "LL# 2950 Delaware River Lighted Buoy 25"] <- round(39.65608,5)
dtc_final1$Longitude[dtc_final1$station_name == "DE River Gate 3A"] <- -75.53357
dtc_final1$Latitude[dtc_final1$station_name == "DE River Gate 3A"] <- 39.42175
dtc_final1$Longitude[dtc_final1$station_name == "DE River Gate 2A"] <- -75.51089
dtc_final1$Latitude[dtc_final1$station_name == "DE River Gate 2A"] <- 39.44077
dtc_final1$Longitude[dtc_final1$station_name == "DE River Gate 1A"] <- -75.51767
dtc_final1$Latitude[dtc_final1$station_name == "DE River Gate 1A"] <- 39.43398
dtc_final1$Longitude[dtc_final1$station_name == "LL# 3920 Upper DE River CLB 36"] <- -74.89752722
dtc_final1$Latitude[dtc_final1$station_name == "LL# 3920 Upper DE River CLB 36"] <- 40.07428


# Remove duplicate
sampling_stations2 <- sampling_stations2 %>%
  filter(!StationName == "LL# 3055 Cherry Is. Flats E. Channel 2")

# lets recode dtc_final1 so that the following join is successful
dtc_final1.5<-dtc_final1 %>% mutate(station_name = recode(station_name, "LL# 3050 Cherry Is Fl East Buoy 1" ="LL# 3050 Cherry Island Flats East Channel Buoy 1",
                                                          "LL# 3055 Cherry Is. Flats E. Channel 2" = "LL# 3055 Cherry Island Flats East Channel Buoy 2",
                                                          "LL# 3065 Cherry Is Flats E Chan. Bouy 4" = "LL# 3065 Cherry Island Flats East Channel Buoy 4",
                                                          "LL#2515 Delaware River Lighted Buoy #3" = "LL# 2515 Delaware River Lighted Buoy 3",
                                                          "LL#2460 Delaware Bay Lighted Buoy #46" = "LL# 2460 Delaware Bay Main Channel Lighted Buoy 46",
                                                          "LL#2472 Delaware Bay Lighter Buoy #49" = "LL# 2472 Delaware Bay Main Channel Lighted Buoy 49",
                                                          "2021 DE Bay Gate 017" = "DE Bay Gate 17",
                                                          "LL# 3055 Cherry Island Flats East Channel Buoy 2" = "LL# 3055 Cherry Island Flats East Channel Buoy 2",
                                                          "LL# 3000 Delaware River Buoy 56" = "LL# 3300 Delaware River Buoy 56",
                                                          "Mantua Creek Anchorage Buoy B" = "LL# 3395 Mantua Creek Anchorage Buoy B",
                                                          "LL# 3775 Channel Lighted Buoy 18" = "LL# 3775 Upper DE River CLB 18",
                                                          "LL# 4120 Channel Buoy 66" = "LL# 4120 Upper DE River CB 66"))
dtc_final1.5$Longitude[dtc_final1.5$station_name == "LL# 3055 Cherry Island Flats East Channel Buoy 2"] <- -75.48903
dtc_final1.5$Latitude[dtc_final1.5$station_name == "LL# 3055 Cherry Island Flats East Channel Buoy 2"] <- 39.72470


#unique(dtc_final1.5$station_name)

# Join to have all station names
full.stationz <- full_join(sampling_stations2, dtc_final1.5, 
                           by = c("StationName" ="station_name", "X" = "Longitude", 
                                                                    "Y" = "Latitude")) %>%
  distinct(X, Y, .keep_all = TRUE)      # distinct matches unique lats and longs 


# Except these regions were missing ugh
full.stationz$Region[full.stationz$StationName == "LL# 3685 Upper DE River CB 9"] <- "River"
full.stationz$Region[full.stationz$StationName == "C&D Canal West"] <- "Canal"
full.stationz$Region[full.stationz$StationName == "Christina River"] <- "Brandywine"
full.stationz$Region[full.stationz$StationName == "New Castle Flats B"] <- "River"
full.stationz$Region[full.stationz$StationName == "Cherry Island Range A"] <- "River"
full.stationz$Region[full.stationz$StationName == "LL# 2820 Bulkhead Shoal Channel B 8"] <- "River"
full.stationz$Region[full.stationz$StationName == "New Castle Flats C"] <- "River"
full.stationz$Region[full.stationz$StationName == "C&D Canal East"] <- "Canal"
full.stationz$Region[full.stationz$StationName == "LL# 2965 Delaware River Lighted Buoy 27"] <- "River"
full.stationz$Region[full.stationz$StationName == "Cherry Island B"] <- "River"
full.stationz$Region[full.stationz$StationName == "LL# 3055 Cherry Island Flats East Channel Buoy 2"] <- "River"
full.stationz$Region[full.stationz$StationName == "C&D Canal East"] <- "Canal"
full.stationz$RKM[full.stationz$StationName == "LL# 2965 Delaware River Lighted Buoy 27"] <- 110
full.stationz$RKM[full.stationz$StationName == "New Castle Flats B"] <- 110
full.stationz$RKM[full.stationz$StationName == "New Castle Flats C"] <- 110
full.stationz$RKM[full.stationz$StationName == "Cherry Island Range A"] <- 112
full.stationz$RKM[full.stationz$StationName == "Cherry Island B"] <- 112
full.stationz$RKM[full.stationz$StationName == "LL# 3055 Cherry Island Flats East Channel Buoy 2"] <- 115
full.stationz$RKM[full.stationz$StationName == "LL# 2820 Bulkhead Shoal Channel B 8"] <- 98
full.stationz$Region[full.stationz$StationName == "LL# 2630 Artificial Is Anchorage Buoy B"] <- "River"
full.stationz$Region[full.stationz$StationName == "LL# 2470 DE Bay Main Channel Buoy 48"] <- "River"
full.stationz$Region[full.stationz$StationName == "LL# 2760 Chesapeake and Delaware Canal Light 3"] <- "Canal"
full.stationz$RKM[full.stationz$StationName == "LL# 3075 Upper Junction Buoy"] <- 120
full.stationz$Region[full.stationz$StationName == "LL# 3075 Upper Junction Buoy"] <- "River"
full.stationz$RKM[full.stationz$StationName == "LL# 3075 Upper Junction Buoy CF"] <- 120
full.stationz$Region[full.stationz$StationName == "LL# 3075 Upper Junction Buoy CF"] <- "River"
full.stationz$RKM[full.stationz$StationName == "Delaware River Lighted Buoy 52"] <- 134
full.stationz$Region[full.stationz$StationName == "Delaware River Lighted Buoy 52"] <- "River"
full.stationz$RKM[full.stationz$StationName == "LL# 3685 Chanel Buoy 9"] <- 172
full.stationz$Region[full.stationz$StationName == "LL# 3685 Chanel Buoy 9"] <- "River"
full.stationz$RKM[full.stationz$StationName == "LL# 3920 Upper DE River CLB 36"] <- 187
full.stationz$Region[full.stationz$StationName == "LL# 3920 Upper DE River CLB 36"] <- "River"
full.stationz$RKM[full.stationz$StationName == "LL# 3920 Channel Lighted Buoy 36"] <- 187
full.stationz$Region[full.stationz$StationName == "LL# 3920 Channel Lighted Buoy 36"] <- "River"
full.stationz$RKM[full.stationz$StationName == "LL#4070 Channel Lighted Buoy 58"] <- 196
full.stationz$Region[full.stationz$StationName == "LL#4070 Channel Lighted Buoy 58"] <- "River"
full.stationz$RKM[full.stationz$StationName == "Delaware River Buoy 48"] <- 130
full.stationz$Region[full.stationz$StationName == "Delaware River Buoy 48"] <- "River"
full.stationz$RKM[full.stationz$StationName == "LL# 4170 Channel Lighted Buoy 76"] <- 202
full.stationz$Region[full.stationz$StationName == "LL# 4170 Channel Lighted Buoy 76"] <- "River"
full.stationz$RKM[full.stationz$StationName == "LL# 3815 Channel Lighted Buoy 26"] <- 181
full.stationz$Region[full.stationz$StationName == "LL# 3815 Channel Lighted Buoy 26"] <- "River"
full.stationz$RKM[full.stationz$StationName == "LL# 3990 Channel Lighted Buoy 48"] <- 192
full.stationz$Region[full.stationz$StationName == "LL# 3990 Channel Lighted Buoy 48"] <- "River"
full.stationz$RKM[full.stationz$StationName == "LL# 3130 Delaware River Lighted Buoy 42"] <- 123
full.stationz$Region[full.stationz$StationName == "LL# 3130 Delaware River Lighted Buoy 42"] <- "River"
full.stationz$RKM[full.stationz$StationName == "LL# 3645 Lighted Buoy 81"] <- 167
full.stationz$Region[full.stationz$StationName == "LL# 3645 Lighted Buoy 81"] <- "River"

# I only need the station specific guts for the next, and final, join :)
full.stations <- full.stationz %>% 
  dplyr::select(StationName:ID) 
full.stations$Region2 <- full.stations$Region
full.stations$Region2[full.stations$Region == "WED"] <- "Atlantic"
full.stations$Region2[full.stations$Region == "COX"] <- "Atlantic"


# Ok last join. Here We add back all the guts of each station (RKM, Region, coordinates, etc) 
# back to dtc_final2, which is our FINAL dataframe (bc its used in the next 1000 lines)
dtc_final1.6 <- left_join(dtc_final1.5, full.stations,
                          by = c("station_name" ="StationName", "Longitude" = "X", 
                                 "Latitude" = "Y")) 


# Some last edits in post
# dtc_final1.6$Region[dtc_final1.6$station_name == "Mantua Creek Anchorage Buoy B"] <- "River"
# dtc_final1.6$RKM[dtc_final1.6$station_name == "Mantua Creek Anchorage Buoy B"] <- 148
# dtc_final1.6$Region[dtc_final1.6$station_name == "LL# 3000 Delaware River Buoy 56"] <- "River"
# dtc_final1.6$RKM[dtc_final1.6$station_name == "LL# 3000 Delaware River Buoy 56"] <- 138
dtc_final1.6$Region[dtc_final1.6$station_name == "DCR-Cooling Water Channel-Outer"] <- "River"
dtc_final1.6$RKM[dtc_final1.6$station_name == "WEDGEPORT"] <- -1000
dtc_final1.6$RKM[dtc_final1.6$station_name == "202102_COX02"] <- -400
dtc_final1.6$Region[dtc_final1.6$station_name == "WEDGEPORT"] <- "Atlantic"
dtc_final1.6$Region[dtc_final1.6$station_name == "202102_COX02"] <- "Atlantic"
dtc_final1.6$Region[dtc_final1.6$station_name == "1"] <- "Brandywine 1"
dtc_final1.6$Region[dtc_final1.6$station_name == "2"] <- "Brandywine 2"
dtc_final1.6$Region[dtc_final1.6$station_name == "Christina River"] <- "Brandywine 3"
dtc_final1.6$RKM[dtc_final1.6$station_name == "Delaware River Lighted Buoy 52"] <- 134
dtc_final1.6$Region[dtc_final1.6$station_name == "Delaware River Lighted Buoy 52"] <- "River"

dtc_final1.6$RKM[dtc_final1.6$Region == "Brandywine 1"] <- 121
dtc_final1.6$RKM[dtc_final1.6$Region == "Brandywine 2"] <- 118.5
dtc_final1.6$RKM[dtc_final1.6$Region == "Brandywine 3"] <- 116.5

# False detections function
dtc_final1.7 <- false_detections(dtc_final1.6, tf=10800)
# make a Region 2 column (again) with Brandywines combined
dtc_final1.7 <- dtc_final1.7 %>%
  mutate(Region2 = str_replace_all(Region, "[[:space:]0-9]", "")) %>%
# Remove any instances of VR2W-137328 
  filter(!Receiver == "VR2W-137328")

# Removing Dead Fish 
dtc_final1.8 <- dtc_final1.7 %>%
  filter(!animal_id %in% c(11,23,17,14,19))

# Finding duplicate occurances
dtc_final2 <- dtc_final1.8 %>%
 distinct(detection_timestamp_utc, animal_id, .keep_all = T) 


########################  Figure 1 - Map of All Receivers  ######################

## Map all receivers
ggplot(data = world) +
  geom_sf(data = states, aes(fill = ID)) +
  geom_polygon(data = canada, aes(x=long, y = lat, group = group), 
               fill = "white", color = "black", linewidth = .05) +
  geom_point(data = full.stations, aes(y= Y, x= X, color= Region2 ), size = 3, alpha= .5) +
  geom_point(data = full.stations, aes(y= Y, x= X, color = Region2 ), size = 3, shape= 1) +
  annotation_scale(location = "bl", width_hint = 0.4) +
  coord_sf(xlim = c(-76.2, -65), ylim = c(38.5, 45))+
  scale_color_manual("Region",values=colpal, labels = c("Atlantic Ocean", "Delaware Bay", "Brandywine River", "C&D Canal", "Delaware River"))+
  theme_bw()+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", 
                                        size = 0.5), 
        panel.background = element_rect(fill = "#EAEAEB"),
        axis.text.x=element_text(angle = 60, hjust = 1)) +
  scale_y_continuous(breaks = seq(38.5, 45.0, 0.5)) +
  scale_x_continuous(breaks = seq(from = -76.5, to = -65, by=1)) +
  scale_fill_manual(values = c("delaware" = "white", "new jersey" = "white",
                               "pennsylvania" = "white", "maryland" = "white",
                               "new york" = "white", "massachusetts" = "white",
                               "connecticut" = "white", "rhode island"= "white",
                               "maine" = "white", "new hampshire" = "white",
                               "vermont" = "white"))+
  guides(fill="none") 

ggsave("AtlanticReceivers3.jpeg", plot = last_plot(), width = 7, height = 5,
       units = "in", dpi = 300, path = "C:/Users/RER/Documents/Masters UD/Code")



## Subset to see Brandywine Receivers
BR<-full.stations %>%  filter(Region %in% c("Brandywine", "Christina"))

# Map just the Delaware receivers
ggplot(data = world) +
  geom_sf(data = states, aes(fill = ID)) +
  geom_point(data = full.stations, aes(y= Y, x= X, color= Region), 
             size = 3, alpha= .5) +
  geom_point(data = full.stations, aes(y= Y, x= X, color = Region), 
             size = 3, shape= 1) +
  annotation_scale(location = "bl", width_hint = 0.4) +
  coord_sf(xlim = c(-76.2, -74.6), ylim = c(38.5, 40.4))+
  scale_color_manual("Region",values=colpal, labels = c("Delaware Bay", "Brandywine River", "C&D Canal", "Delaware River")) +
  theme_bw()+
  xlab("Longitude")+
  ylab("Latitude")+
  geom_text_repel(data = BR, aes(x = X, y = Y, label = glatos_array.x),
                  fontface = "bold",
                  nudge_x = c(.2, -.25),
                  nudge_y = c(0.1, -0.15)) +
  theme(panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", 
                                        size = 0.5), 
        panel.background = element_rect(fill = "#EAEAEB"),
        # legend.position = "none",
        axis.text.x=element_text(angle = 60, hjust = 1)) +
  scale_y_continuous(breaks = seq(38.5, 41.0, 0.5)) +
  scale_x_continuous(breaks = seq(from = -76.5, to = -74.5, by=0.5)) +
  scale_fill_manual(values = c("delaware" = "white", "new jersey" = "white",
                               "pennsylvania" = "white", "maryland" = "white"))+
  guides(fill="none")

ggsave("DelawareReceivers3.jpeg", plot = last_plot(), width = 5, height = 7,
       units = "in", dpi = 300, path = "C:/Users/RER/Documents/Masters UD/Code")

########################  Figure 2 - Length Stats #############################

# AVerage Lengths per year and sex
fsh %>%
  group_by(Deploy.Year) %>%
  summarise(TL.avg = mean(Fork.mm),
            TL.SD = sd(Fork.mm))


## movement across years
# Test fish  for normality          # We meet normality for both years!                     
shapiro.test(fsh$Fork.mm[fsh$Deploy.Year==2021])
shapiro.test(fsh$Fork.mm[fsh$Deploy.Year==2022])
# Test for differences in variance across years     # Variances are equal too. 
levene.test(fsh$Fork.mm, group = fsh$Deploy.Year)
# Time for a t test to test in length differences across years
# No difference in length between years :)
t.test(fsh$Fork.mm[fsh$Deploy.Year==2021], fsh$Fork.mm[fsh$Deploy.Year==2022])


# Box plot of fish length across years
fsh %>%
  ggplot(aes(x=Deploy.Year,y=Fork.mm))+
  geom_boxplot(notch = F, outlier.color = NA)+
  theme_classic() +
  labs(x = "Tagging Year", y="Fork Length (mm)") 
# geom_signif(comparisons = list(c("Female", "Male")), map_signif_level=F)

ggsave("ForkLengthbyYEAR.png", plot = last_plot(), width = 3, height = 3,
       units = "in", dpi = 300, path = "C:/Users/RER/Documents/Masters UD/Code")



########################  Figure 3 and 4  - All Detects ###########################

dtc_final2 %>%
  filter(Year == 2021, !Region == "Atlantic") %>%
  ggplot(aes(x=as.Date(detection_timestamp_utc), y=RKM, 
             color=Region, alpha=.2))+
  geom_point(size=3, alpha = .3)+
  geom_point(size =3, shape = 1)+
  labs(title = "", y= "River Kilometer", x="Date")+
  facet_wrap(~animal_id) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = 15)) +
  guides(alpha = "none", color=guide_legend(title="Region")) +
  geom_vline(xintercept = (as.Date("2021-05-11 10:00:00")), 
             colour = "red") +
  #geom_vline(xintercept = as.POSIXct(as.Date("2021-05-11 10:00:00")+survival_gap), colour = "black") +
  scale_color_manual("Region",values=colpal, labels = c("Delaware Bay", "Brandywine River 1", "Brandywine River 2", "Brandywine River 3", "Delaware River")) + 
  scale_x_date(date_breaks = "month", date_labels = "%B", limits = as.Date(c("2021-05-1", "2021-07-01")))



ggsave("2021Detectss2.2.jpeg", plot = last_plot(), width = 11, height = 7, units = "in",
       dpi = 300,  path = "C:/Users/RER/Documents/Masters UD/Code")

dtc_final2 %>%
  filter(Year == 2022, !Region == "Atlantic") %>%
  ggplot(aes(x=as.Date(detection_timestamp_utc), y=RKM, 
             color=Region, alpha=.2))+
  geom_point(size=3, alpha = .3)+
  geom_point(size =3, shape = 1)+
  labs(title = "", y= "River Kilometer", x="Date")+
  facet_wrap(~animal_id) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = 15)) +
  guides(alpha = "none", color=guide_legend(title="Region")) +
  geom_vline(xintercept = as.Date("2022-05-10 10:00:00"), colour = "red") +
  #geom_vline(xintercept = as.POSIXct(as.Date("2021-05-11 10:00:00")+survival_gap), colour = "black") +
  scale_color_manual("Region",values=colpal, labels = c("Delaware Bay", "Brandywine River 1", "Brandywine River 2", "Brandywine River 3", "C&D Canal", "Delaware River")) +
  scale_x_date(date_breaks = "month", date_labels = "%B", limits = as.Date(c("2022-05-1", "2022-08-01")))

ggsave("2022Detectss2.2.jpeg", plot = last_plot(), width = 11, height = 7, units = "in",
       dpi = 300,  path = "C:/Users/RER/Documents/Masters UD/Code")


#####  Fallback stats  ####

Fish7 <- dtc_final2 %>% filter(animal_id == 7) %>% arrange(detection_timestamp_utc) %>%
  select(detection_timestamp_utc, RKM, station_name, animal_id, Region, Region2)
Fish22 <- dtc_final2 %>% filter(animal_id == 22) %>% arrange(detection_timestamp_utc)%>%
  select(detection_timestamp_utc, RKM, station_name, animal_id, Region, Region2)
Fish24 <- dtc_final2 %>% filter(animal_id == 24) %>% arrange(detection_timestamp_utc)%>%
  select(detection_timestamp_utc, RKM, station_name, animal_id, Region, Region2)

Fish7 %>%
  filter(!Region == "Atlantic") %>%
  ggplot(aes(x=as.Date(detection_timestamp_utc), y=RKM, 
             color=Region, alpha=.2))+
  geom_point(size=3, alpha = .3)+
  geom_point(size =3, shape = 1)+
  labs(title = "", y= "River Kilometer", x="Date")+
  facet_wrap(~animal_id) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = 15)) +
  guides(alpha = "none", color=guide_legend(title="Region")) +
  geom_vline(xintercept = (as.Date("2021-05-11 10:00:00")), 
             colour = "red") +
  geom_vline(xintercept = (as.Date("2021-05-12 14:00:00")), colour = "black") +
  scale_color_manual("Region",values=colpal, labels = c("Delaware Bay", "Brandywine River 1", "Brandywine River 2", "Brandywine River 3", "Delaware River")) + 
  scale_x_date(date_breaks = "month", date_labels = "%B", limits = as.Date(c("2021-05-1", "2021-07-01")))

Fish22 %>%
  filter(!Region == "Atlantic") %>%
  ggplot(aes(x=as.Date(detection_timestamp_utc), y=RKM, 
             color=Region, alpha=.2))+
  geom_point(size=3, alpha = .3)+
  geom_point(size =3, shape = 1)+
  labs(title = "", y= "River Kilometer", x="Date")+
  facet_wrap(~animal_id) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = 15)) +
  guides(alpha = "none", color=guide_legend(title="Region")) +
  geom_vline(xintercept = (as.Date("2022-05-10 10:00:00")), 
             colour = "red") +
  geom_vline(xintercept = (as.Date("2022-05-11 14:00:00")), colour = "black") +
  scale_color_manual("Region",values=colpal, labels = c("Delaware Bay", "Brandywine River 1", "Brandywine River 2", "Brandywine River 3", "Delaware River")) + 
  scale_x_date(date_breaks = "month", date_labels = "%B", limits = as.Date(c("2022-05-1", "2022-07-01")))


########################  Figure 5 and 6  - Occupancy ##########################
# Set regions based on RKM
residency_detects <- dtc_final2%>%
  group_by(animal_id) %>%
  mutate(Region3 = removeNumbers(Region),
         segment = NA,
         detection_timestamp_utc = as.POSIXct(detection_timestamp_utc)) %>%
  rename("deploy_lat" ="Latitude", "deploy_long" = "Longitude") %>%
  dplyr::select(detection_timestamp_utc, station_name, Deploy.Year, Year, Region, Owner, RKM, Region3, segment, animal_id, deploy_lat, deploy_long) %>%
  mutate(segment = case_when((Region3 == "Brandywine " & RKM == 121) ~ "Brandywine 1",
                             (Region3 == "Brandywine " & RKM == 118.5) ~ "Brandywine 2",
                             (Region3 == "Brandywine " & RKM == 116.5) ~ "Brandywine 3",
                             (Region3 == "River" & RKM > 120) ~ "Delaware River",
                             (Region3 == "River" & RKM <= 120 & RKM >= 109) ~ "Delaware River",
                             (Region3 == "River" & RKM > 76 & RKM < 109) ~ "Delaware River",
                             (Region3 == "Bay") ~ "Delaware Bay",
                             (Region3 == "Atlantic") ~ "Atlantic Ocean",
                             (Region3 == "Canal") ~ "Delaware River",
                             .default = "NA")) %>%
  filter(segment != "NA")

# Run the glatos function to calculate unique detection events
FishDetects <- detection_events(residency_detects, location_col = "segment", time_sep = Inf, condense = TRUE) 

# Figure 5 - bar plot showing relative habitat occupancy
FishDetects %>%
  filter(!animal_id %in% c(11,14,17,19,23)) %>%
  ggplot(aes(x = factor(animal_id), y = res_time_sec, fill = location))+
  geom_bar(stat = 'identity', position = 'fill') +
  theme_bw() +
  labs(x = "Fish ID", y = "Proportional Occupancy", fill = "Region") + 
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(values = colpal) +
  theme(panel.grid = element_blank())

ggsave("OccupancyAll2.jpeg", plot = last_plot(), width = 7,  height =4,
       units = "in", dpi = 500,  path = "C:/Users/RER/Documents/Masters UD/Code")


# Reset locations in unique detection event dataframe so that Brandywine is combined
FishDetects2 <- FishDetects %>%
  mutate(location = recode(location, "Brandywine 1" = 'Brandywine River',
                           "Brandywine 2" = 'Brandywine River',
                           "Brandywine 3" = 'Brandywine River',
                           "Canal" = "Delaware River"))

# And plot! Bar plot of proportional occpancy 
FishDetects2 %>%
  filter(!animal_id %in% c(11,14,17,19,23)) %>%
  ggplot(aes(x = factor(animal_id), y = res_time_sec, fill = location))+
  geom_bar(stat = 'identity', position = 'fill') +
  theme_bw() +
  labs(x = "Fish ID", y = "Proportional Occupancy", fill = "Region") +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  theme(panel.grid = element_blank()) +
  scale_fill_manual(values = colpal)+
  theme(panel.grid = element_blank())



ggsave("OccupancyBRComb2.jpeg", plot = last_plot(), width = 7,  height =4,
       units = "in", dpi = 500,  path = "C:/Users/RER/Documents/Masters UD/Code")

####

# Remove Atlantic detections to avoid that shit
DEdetects <- dtc_final2 %>%
  filter(!Region2 == "Atlantic") %>%
  mutate(Region3 = recode(Region2, "Canal" = "River")) %>%
  rename("deploy_lat" ="Latitude", "deploy_long" = "Longitude") %>%
  mutate(detection_timestamp_utc = as.POSIXct(detection_timestamp_utc)) 
  

occupancy <- detection_events(DEdetects, location_col = "Region3", time_sep = Inf, condense = TRUE)
occupancy2 <- occupancy %>% mutate (res_time_min = res_time_sec/60, 
                                    res_time_hr = round(res_time_min/60,4),
                                    res_time_hr_log = log(res_time_hr),
                                    Year = year(first_detection)) %>%
  filter(is.finite(res_time_hr_log)) 

occupancy.avg <- occupancy2 %>%
  group_by(animal_id, location, Year) %>%
  summarise(total.time = sum(res_time_sec)) %>%
  ungroup() %>%
  group_by(location, Year) %>%
  summarise(avg.time = mean(total.time))


occupancy.avg %>%
ggplot(aes(x = factor(Year), y = avg.time, fill = location))+
  geom_bar(stat = 'identity', position = 'fill') +
  theme_bw() +
  labs(x = "Detection Year", y = "Average Proportional Occupancy", fill = "Region") + 
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(values = colpal) +
  theme(panel.grid = element_blank())


occupancy.avg %>%
  ggplot(aes(x = location, y = avg.time/3600, fill = factor(Year)))+
  geom_bar(stat = 'identity') +
  theme_bw() +
  labs(x = "Location", y = "Average Occupancy Time (hr)", fill = "Detection Year") + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(expand = expansion(mult = c(0.1, 0.1))) +
  scale_fill_manual(values = c("lightgray", "darkgrey")) +
  theme(panel.grid = element_blank()) +
  facet_wrap(~Year)


occupancy.avg2 <- occupancy2 %>%
  group_by(animal_id, location, Year) %>%
  summarise(total.time = sum(res_time_sec))


occupancy.avg2 %>%
  ggplot(aes(x = location, y = total.time/3600, fill = factor(Year)))+
  geom_boxplot(coef = 3)+
  theme_bw() +
  labs(x = "Location", y = "Total Occupancy Time (hr)", fill = "Detection Year") + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(expand = expansion(mult = c(0.1, 0.1)), labels = c("Delaware \nBay", 
                   "Brandywine \nRiver", "Delaware \nRiver")) +
  scale_fill_manual(values = c("lightgray", "darkgrey")) +
  theme(panel.grid = element_blank()) +
  facet_wrap(~Year)

ggsave("OccpuancyBoxplot.jpeg", plot = last_plot(), width = 7,  height = 4,
       units = "in", dpi = 300,  path = "C:/Users/RER/Documents/Masters UD/Code") 


occupancy.avg2021 <- occupancy.avg2 %>% filter(Year == 2021)
shapiro.test(occupancy.avg2021$total.time[occupancy.avg2021$location=="Brandywine"])
shapiro.test(occupancy.avg2021$total.time[occupancy.avg2021$location=="Bay"])
shapiro.test(occupancy.avg2021$total.time[occupancy.avg2021$location=="River"])


occupancy.avg2022 <- occupancy.avg2 %>% filter(Year == 2022)
shapiro.test(occupancy.avg2022$total.time[occupancy.avg2022$location=="Brandywine"])
shapiro.test(occupancy.avg2022$total.time[occupancy.avg2022$location=="Bay"])
shapiro.test(occupancy.avg2022$total.time[occupancy.avg2022$location=="River"])


occ2022 <- occupancy.avg2022 %>% ungroup %>% select(location, total.time)
mat1 <- data.matrix(occ2022)
friedman.test(mat1)

occ2021 <- occupancy.avg2021 %>% ungroup %>% select(location, total.time)
mat2 <- data.matrix(occ2021)
friedman.test(mat2)



########################  Figure 7 - Acceleration Data  ######################

# Filter for sensors that contain motion data
ac.df <- dtc_final2 %>% 
  #filter(transmitter_codespace %in% c("A69-9007", "A69-9002")) %>% 
  #filter(!sensor_unit == "ADC") %>%
  drop_na(sensor_value) %>%
  mutate(accel = case_when(
    sensor_unit == "ADC" ~ (sensor_value * 0.013588),
    startsWith(sensor_unit, "m") ~ sensor_value
  ))


# whats the maximum speed value? 3.5 m/s^2
max(ac.df$accel)


# Summary statistics 
accel.stats <- ac.df %>%
  group_by(animal_id) %>%
  summarize(mean = mean(accel),
            median = mean(accel),
            min = min(accel),
            max = max(accel),
            mode = mode(accel),
            sd = sd(accel))



Locations = c("Bay" = "Delaware Bay", "Brandywine" = "Brandywine River", "River" = "Delaware River")
# plot acceleration 
ac.df %>%
  drop_na(Region) %>%
  filter(!sensor_value == 0) %>%
  ggplot(aes(x = accel, fill = factor(animal_id)))+
  geom_histogram(binwidth = .4, color = "black", closed = "left") +
  facet_wrap(~Region2, nrow = 1, labeller = as_labeller(Locations)) +
  theme_bw() +
  labs(y = "Frequency", fill = "Fish ID") +
  scale_fill_manual(values = c("lightgray", "darkgray")) +
  xlab(expression(Acceleration~Rate~(m/s^2))) +
  theme(panel.grid = element_blank()) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                expand = expansion(mult = c(0, .1)))

# plot acceleration 
ac.df %>%
  ggplot(aes(x = sensor_value, fill = factor(animal_id)))+
  geom_histogram(binwidth = .2, color = "black") +
  facet_wrap(~Region, nrow = 1, labeller = as_labeller(strip_labels)) +
  theme_bw() +
  labs(y = "Frequency", fill = "Fish ID") +
  #scale_fill_manual(values = colpal) +
  xlab(expression(Acceleration~Rate~(m/s^2))) +
  theme(panel.grid = element_blank(),
        text = element_text(size = 14)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  scale_fill_manual(values = c("#96b78c", "#ad8cb7"))


ggsave("Acceler.Histo2.png", plot = last_plot(),width = 7, height = 3, 
       units = "in", dpi = 300, path = "C:/Users/RER/Documents/Masters UD/Code")


shapiro.test(ac.df$sensor_value[ac.df$Region == "River"])
# NOT normal
shapiro.test(ac.df$sensor_value[ac.df$Region == "Brandywine 1"])
# NORMAL
shapiro.test(ac.df$sensor_value[ac.df$Region == "Brandywine 2"])
# NOT normal
shapiro.test(ac.df$sensor_value[ac.df$Region == "Brandywine 3"])
# NOT normal
shapiro.test(ac.df$sensor_value[ac.df$Region == "Bay"])
# NORMAL

# We move onto non-parametric tests across groups for a Kruskal - Wallis test
kruskal.test(sensor_value ~ Region, data=ac.df)

dunnTest(sensor_value ~ Region, data=ac.df, method="bonferroni")

## Investigating animal ID
shapiro.test(ac.df$sensor_value[ac.df$animal_id == "1"])
# NOT normal
shapiro.test(ac.df$sensor_value[ac.df$animal_id == "2"])
# NOT normal

t.test(ac.df$sensor_value[ac.df$animal_id == "2"], ac.df$sensor_value[ac.df$animal_id == "1"])
# Different

########################  APPENDIX FIGURES  ####################################

####   Figure S1   ####
dtc_final1.6 %>%
  filter(!Deploy.Year %in% c(2023, 2022)) %>%
  filter(animal_id %in% c(11,23,17,14,19)) %>%
  filter(!Receiver == "PROJ162") %>%
  mutate(detection_timestamp_utc = as.POSIXct(detection_timestamp_utc)) %>%
  ggplot(aes(x=detection_timestamp_utc, y=animal_id, color=Region, alpha=.2))+
  geom_point(size=4)+
  labs(title = "",y= "Animal ID", x="Date", color = "Region")+
  scale_y_continuous(breaks = seq(1,26,2)) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank()) +
  guides(alpha = "none", color=guide_legend(title="Region")) +
  geom_vline(xintercept = as.POSIXct(as.Date("2021-05-11 10:00:00")), 
             colour = "red") +
  geom_vline(xintercept = as.POSIXct(as.Date("2021-05-11 10:00:00")+survival_gap), 
             colour = "black") +
  scale_color_manual(values=colpal)

ggsave("2021FISH.1.png", plot = last_plot(), width = 8, height = 5,
       units = "in", dpi = 300, path = "C:/Users/RER/Documents/Masters UD/Code")

####   Figure S2   ####
dtc_final1.6 %>%
  filter(!Deploy.Year %in% c(2023, 2021)) %>%
  filter(animal_id %in% c(11,23,17,14,19)) %>%
  mutate(detection_timestamp_utc = as.POSIXct(detection_timestamp_utc)) %>%
  ggplot(aes(x=detection_timestamp_utc, y=factor(animal_id), color=Region, alpha=.2))+
  geom_jitter(size=4, height = .1)+
  labs(title = "",y= "Animal ID", x="Date", color = "Region")+
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank()) +
  guides(alpha = "none", color=guide_legend(title="Region")) +
  geom_vline(xintercept = as.POSIXct(as.Date("2022-05-10 10:00:00")), 
             colour = "red") +
  geom_vline(xintercept = as.POSIXct(as.Date("2022-05-10 10:00:00")+survival_gap), 
             colour = "black") +
  scale_color_manual(values=colpal)


ggsave("2022FISH.1.png", plot = last_plot(), width = 8, height = 5,
       units = "in", dpi = 300, path = "C:/Users/RER/Documents/Masters UD/Code")

####   Figure S3   ####
dtc_final2 %>%
  filter(animal_id %in% c(1,7)) %>%
  mutate(Region3 = recode(Region, "Bay" = "Delaware Bay", "River" = "Delaware River",
                                    "Atlantic" = "Atlantic Ocean")) %>%
  ggplot(aes(x=as.Date(detection_timestamp_utc), y=RKM, color=Region3))+
  geom_point(size=3, alpha = .15)+
  geom_point(size =3, shape = 1) +
  labs(title = "", y= "River Kilometer", x="Date")+
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        text = element_text(size=15),
        axis.text.x = element_text(angle=45, hjust =1)) +
  guides(alpha = "none", color=guide_legend(title="Region")) +
  facet_wrap(~animal_id, scales = "free") +
  scale_color_manual(values=colpal)

ggsave("Fish1and7.2.jpeg", plot = last_plot(), width = 10,height = 4, units = "in",
       dpi = 300,  path = "C:/Users/RER/Documents/Masters UD/Code") 

####   Figure S4   #####

# Arrange by timestamp and animal to calculate movement, time between stations, and speed
movement <- dtc_final2 %>%
  group_by(animal_id) %>%
  arrange(detection_timestamp_utc, .by_group = TRUE) %>%
  mutate(Mvmt.m = (abs(RKM - lag(RKM)))*1000,
         Time.diff = as.numeric(ymd_hms(detection_timestamp_utc) - lag(ymd_hms(detection_timestamp_utc))),
         Move_rate = Mvmt.m/Time.diff,
         last.station = lag(station_name)) 

# Total distance Traveled 
movement2 <- movement %>%
  group_by(animal_id) %>%
  filter(!is.na(Mvmt.m)) %>%
  filter(!is.infinite(Move_rate)) %>%
  # Time at large
  summarise(Total.move.km = sum(Mvmt.m)/1000,
            Time.At.Large = as.numeric(as.duration(max(ymd_hms(detection_timestamp_utc)) - min(ymd_hms(detection_timestamp_utc)))/604800),
            avg.rate = mean(Move_rate, na.rm = T)) 

# Summary statistics
move.stats <- movement2 %>%
  filter(!animal_id %in% c(1,7)) %>%
  summarise(mean.time = mean(Time.At.Large),
            sd.time = sd(Time.At.Large),
            min = min(Time.At.Large),
            max = max(Time.At.Large),
            mean.move = mean(Total.move.km),
            sd.move = sd(Total.move.km),
            min.m = min(Total.move.km),
            max.m = max(Total.move.km))

# Cast the DF longer
movement2.1 <- movement2 %>%
  pivot_longer(!animal_id, names_to = "Variable", values_to = "Value") %>%
  filter(Variable == "Total.move.km") %>%
  mutate(Variable = recode(Variable, "Time.At.Large" = 'Time at large',
                           "Total.move.km" = "Total movement")) 
movement2.2 <- movement2 %>%
  filter(!animal_id %in% c(11,17,23,14,19)) %>%
  pivot_longer(!animal_id, names_to = "Variable", values_to = "Value") %>%
  filter(Variable == "Time.At.Large") %>%
  mutate(Variable = recode(Variable, "Time.At.Large" = 'Time at large',
                           "Total.move.km" = "Total movement"))

# Graph total distance traveled
ggplot(movement2.1, aes(x=factor(Variable), y=Value))+
  geom_bar(aes(x = factor(animal_id), y = Value), stat= "identity", fill ="lightgrey", color = "black", width = .8) +
  theme_bw()+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                expand = expansion(mult = c(0, .1)))+
  labs(x="Fish ID", y="Total distance traveled (km)") +
  theme(text=element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 


ggsave("Total.mvmt2.1.3.png", plot = last_plot(), width = 6,height = 5, units = "in", dpi = 300,
       path = "C:/Users/RER/Documents/Masters UD/Code")

# Graph time at large
ggplot(movement2.2, aes(x=factor(Variable), y=Value))+
  geom_bar(aes(x = factor(animal_id), y = Value), stat= "identity", fill ="lightgrey", color = "black", width = .8) +
  theme_bw()+
  labs(x="Fish ID", y = "Time at large (wks)") +
  theme(text=element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(expand = expansion(mult = c(0, .1)))

ggsave("Total.time.2.3.png", plot = last_plot(), width = 6,height = 5, units = "in", dpi = 300,
       path = "C:/Users/RER/Documents/Masters UD/Code")
####   Figure S6   ####

# Rename variables 
DelMemSite <- "01482100"
ChristinaConfluence <- "01480120"
Brandywine <- "01481500"
UpperBR <- "01481000"
discharge <- "00060"
gageH <- "00065"
temp <- "00010"
rain <- "00045"
ph <- "00400"
DO <- "00300"

# Join variables
siteNumbers <- c(ChristinaConfluence, Brandywine)
parameters <- c(discharge, gageH, temp, rain, ph, DO)

# Raw daily data:
DailyData <- readNWISuv(siteNumbers, parameters, "2021-03-01", "2022-8-31")
DailyData <- renameNWISColumns(DailyData)

# Rename Columns 
DailyData$site_name[DailyData$site_no == "01480120"] <- "USGS Station 01480120 \n (Christina Confluence)"
DailyData$site_name[DailyData$site_no == "01481500"] <- "USGS Station 01481500 \n(Brandywine River)"
DailyData$site_name[DailyData$site_no == "01481000"] <- "USGS Station 01481000 \n(Brandywine River at Chadds Ford)"

# Convert Flow to meters
#DailyData$Flow.m <- DailyData$Flow*0.0283168
DailyData$Flow.m <- DailyData$Flow_Inst*0.0283168

# Find hourly averages of all environmental data
USGS.hour<-DailyData %>%
  mutate(month = month(dateTime),
         year = year(dateTime),
         day = day(dateTime),
         hour2 = (hour(dateTime) - 4), # conversion to UTC
         hour = ifelse(hour2 < 0, hour2 + 24, hour2),
         minute = minute(dateTime),
         Date = paste(year, month, day, hour, minute, sep = "-"),
         Gage.m = .Downstream._GH_Inst*0.3048) %>%
  filter(month %in% c(4,5,6,7,8))%>%
  group_by(site_name, site_no, year, month, day, hour) %>%
  summarise(hourlyTemp = mean(Wtemp_Inst),
            hourlypH = mean(pH_Inst),
            hourlyFlow = mean(Flow.m),
            hourlyGage = mean(Gage.m),
            hourlyDO = mean(DO_Inst)) %>%
  mutate(Date = paste(year, month, day, sep = "-"),
         DateTime = ymd_h(paste(Date, hour, sep = " ")))


# Filter to just the Brandywine USGS location
USGS.BR <- USGS.hour %>%
  filter(site_no == "01481500")

# USGS.CR <- USGS.hour %>%
#   filter(site_no == "01480120")


# Create data frame for hourly resolution on frequency resolution
dtc.hour <- dtc_final2 %>%
  mutate(month = month(detection_timestamp_utc),
         year = year(detection_timestamp_utc),
         day = day(detection_timestamp_utc),
         hour = hour(detection_timestamp_utc), 
         minute = minute(detection_timestamp_utc),
         Date = paste(year, month, day, hour, minute, sep = "-")) %>%
  filter(month %in% c(4,5,6,7,8))%>%
  filter(Region2 == "Brandywine") %>%
  group_by(year, month, day, hour) %>%
  summarise(Detects = n())%>%
  mutate(Date = paste(year, month, day, sep = "-"),
         DateTime = ymd_h(paste(Date, hour, sep = " ")))

# Join detections and USGS water data
water.BR <- full_join(dtc.hour, USGS.BR, by = c("DateTime", "Date", "year", "month", "day", "hour")) 
# water.CR <- full_join(dtc.hour, USGS.CR, by = c("DateTime", "Date", "year", "month", "day", "hour")) 

# Find the lag time difference to make sure that intervals on time is the same
water.BR2 <- water.BR %>%
  mutate(Detects2 = replace(Detects, is.na(Detects), 0)) %>%
  ungroup() %>%
  arrange(DateTime) %>%
  group_by(year) %>%
  filter(!is.na(hourlyTemp),
         !is.na(hourlyDO),
         !is.na(hourlypH)) %>%
  mutate(lag = lag(DateTime, n = 1),
         timediff = DateTime - lag) %>%
  filter(!is.na(lag))


# Run cross correlation functions
ccf(as.ts(water.BR2$hourlyTemp), as.ts(water.BR2$Detects2)) # Significant on the negative, -12 - -34 hours
ccf(as.ts(water.BR2$hourlyDO), as.ts(water.BR2$Detects2))   # not signif
ccf(as.ts(water.BR2$hourlyFlow), as.ts(water.BR2$Detects2)) # Significant on the positive, 30 hours and  -12 hours
ccf(as.ts(water.BR2$hourlyGage), as.ts(water.BR2$Detects2)) # not signif
ccf(as.ts(water.BR2$hourlypH), as.ts(water.BR2$Detects2))   # Significant on the negative, -35 - -27 hours, -8 - +34 hours



# Make graph showing the different environmetnal data over time
WaterParams = c("hourlyTemp" = "Average Temperature (°C)", "hourlyDO" = "Average DO (mg L-1)", "hourlypH" = "Average pH", 
                "hourlyFlow" = "Average Discharge (m3 s-1) ", "hourlyGage" = "Average Gage Height (m)", "2021" = "2021", "2022" = "2022")

water.BR2 %>%
  pivot_longer(cols = c(hourlyTemp, hourlyDO, hourlypH, hourlyFlow, hourlyGage), names_to = "Water Parameter", values_to = "Values") %>%
 # filter(!month == 4) %>%
  ggplot(aes(x = DateTime, y = Values))+
  geom_line() +
  facet_grid(`Water Parameter` ~ year, scales = "free", labeller = as_labeller(WaterParams)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 13)) +
  labs(x = "Date", y = "")

ggsave("WaterParams.png", plot = last_plot(), width = 15,height = 11, units = "in", dpi = 300,
       path = "C:/Users/RER/Documents/Masters UD/Code")



water.BR2 %>%
  group_by(year, month) %>%
  summarise(avgTemp = mean(hourlyTemp),
            avgFlow = mean(hourlyFlow),
            avgpH = mean(hourlypH))



FLOW <- water.BR2 %>% 
  filter(year == 2021,
         month == 6)
# June 9
GONE <- dtc_final2 %>%
  mutate(month = month(detection_timestamp_utc)) %>%
  filter(animal_id %in% c(4,6,9),
         month == 6,
         station_name == 1) %>%
  dplyr::select(detection_timestamp_utc, station_name, animal_id) %>%
  group_by(animal_id) %>%
  summarise(max = max(detection_timestamp_utc))



FLOW %>%
  filter(month %in% c(5,6,7,8)) %>%
  ggplot(aes(x = DateTime, y = hourlyGage))+
  geom_point() +
  facet_wrap(~year, scales = "free") +
  labs(y = "Gage Height (m)")
