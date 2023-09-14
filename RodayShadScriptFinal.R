##
## Date Created: 2023-09-12
##
## Copyright (c) Rachel Roday, 2023
## Email: rroday@udel.edu
##
## ---------------------------
##
## Notes: This is the Final Script for Roday et al. 2024
##
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
  Deploy.Year = c("2022", "2022","2022", "2022","2022", "2022","2022", "2022",
                  "2022", "2022","2022", "2022","2022", "2022","2021","2021","2021","2021","2021","2021",
                  "2021","2021","2021","2021","2021","2021"),
  stringsAsFactors = FALSE)  


# Create color palette to use for whole paper
colpal <- c("Atlantic" = "#4477AA",  "Bay" = "#EE6677", "Brandywine 1" = "#228833", 
            "Brandywine 2" = "#AA3377", "Brandywine 3" = "#66CCEE",  "River" = "#E69F00", "#000000",
            "Brandywine" = "#66CCEE", "Canal" = "#228833")

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
setnames(MATOS, "Date and Time (UTC)", "detection_timestamp_utc")

# Parse date column
AllDetects_2021_2022$detection_timestamp_utc <- as.character(mdy_hms(AllDetects_2021_2022$detection_timestamp_utc))
MATOS$detection_timestamp_utc <- as.character(mdy_hms(MATOS$detection_timestamp_utc))

# Merge 3 datasets
dtc <- full_join(BR1_2022, AllDetects_2021_2022)
dtc2 <- full_join(dtc, MATOS)


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

# Add a year column
dtc_final1 <- dtc_final %>% 
  mutate(Year = lubridate::year(as.Date(detection_timestamp_utc)))


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
                                                          "LL# 3055 Cherry Island Flats East Channel Buoy 2" = "LL# 3055 Cherry Island Flats East Channel Buoy 2"))
dtc_final1.5$Longitude[dtc_final1.5$station_name == "LL# 3055 Cherry Island Flats East Channel Buoy 2"] <- -75.48903
dtc_final1.5$Latitude[dtc_final1.5$station_name == "LL# 3055 Cherry Island Flats East Channel Buoy 2"] <- 39.72470




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
dtc_final1.6 <- dtc_final1.6 %>% drop_na(Transmitter) 

# Some last edits in post
dtc_final1.6$Region[dtc_final1.6$station_name == "Mantua Creek Anchorage Buoy B"] <- "River"
dtc_final1.6$RKM[dtc_final1.6$station_name == "Mantua Creek Anchorage Buoy B"] <- 148
dtc_final1.6$Region[dtc_final1.6$station_name == "LL# 3000 Delaware River Buoy 56"] <- "River"
dtc_final1.6$RKM[dtc_final1.6$station_name == "LL# 3000 Delaware River Buoy 56"] <- 138
dtc_final1.6$Region[dtc_final1.6$station_name == "DCR-Cooling Water Channel-Outer"] <- "River"
dtc_final1.6$RKM[dtc_final1.6$station_name == "WEDGEPORT"] <- -1000
dtc_final1.6$RKM[dtc_final1.6$station_name == "202102_COX02"] <- -400
dtc_final1.6$Region[dtc_final1.6$station_name == "WEDGEPORT"] <- "Atlantic"
dtc_final1.6$Region[dtc_final1.6$station_name == "202102_COX02"] <- "Atlantic"
dtc_final1.6$Region[dtc_final1.6$station_name == "1"] <- "Brandywine 1"
dtc_final1.6$Region[dtc_final1.6$station_name == "2"] <- "Brandywine 2"
dtc_final1.6$Region[dtc_final1.6$station_name == "Christina River"] <- "Brandywine 3"

dtc_final1.6$RKM[dtc_final1.6$Region == "Brandywine 1"] <- 121
dtc_final1.6$RKM[dtc_final1.6$Region == "Brandywine 2"] <- 118.5
dtc_final1.6$RKM[dtc_final1.6$Region == "Brandywine 3"] <- 116.5

# False detections function
dtc_final1.7 <- false_detections(dtc_final1.6, tf=10800)


# Removing Dead Fish 
dtc_final1.8 <- dtc_final1.7 %>%
  filter(!animal_id %in% c(11,23,17,14,19))

# Finding duplicate occurances
DUPES <- dtc_final1.8 %>%
  distinct(detection_timestamp_utc, animal_id, station_name, .keep_all = T) %>%
  group_by(detection_timestamp_utc, animal_id) %>%
  mutate(FREQ = n()) %>%
  filter(FREQ > 1) %>%
  dplyr::select(detection_timestamp_utc, Year, animal_id, station_name, Latitude, Longitude, RKM, Region, Owner, FREQ) %>%
  mutate(avg.RKM = mean(RKM),
         number = 1,
         n = cumsum(number)) %>%
  ungroup() 

# Merge DUPES and dtc_final, remove all duplicates
dtc_final2 <- full_join(dtc_final1.8, DUPES) %>%
  filter(n %in% c(NA, 1))


########################  Figure 1 - Map of All Receivers  ######################

## Map all receivers
ggplot(data = world) +
  geom_sf(data = states, aes(fill = ID)) +
  geom_polygon(data = canada, aes(x=long, y = lat, group = group), 
               fill = "white", color = "black", linewidth = .05) +
  geom_point(data = full.stations, aes(y= Y, x= X, color= Region2 ), size = 3, alpha= .3) +
  geom_point(data = full.stations, aes(y= Y, x= X, color = Region2 ), size = 3, shape= 1) +
  annotation_scale(location = "bl", width_hint = 0.4) +
  coord_sf(xlim = c(-76.2, -65), ylim = c(38.5, 45))+
  scale_color_manual("Region",values=colpal)+
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

ggsave("AtlanticReceivers.jpeg", plot = last_plot(), width = 7, height = 5,
       units = "in", dpi = 300, path = "C:/Users/RER/Documents/Masters UD/Code")



## Subset to see Brandywine Receivers
BR<-full.stations %>%  filter(Region == "Brandywine")


# Map just the Delaware receivers
ggplot(data = world) +
  geom_sf(data = states, aes(fill = ID)) +
  geom_point(data = full.stations, aes(y= Y, x= X, color= Region), 
             size = 3, alpha= .3) +
  geom_point(data = full.stations, aes(y= Y, x= X, color = Region), 
             size = 3, shape= 1) +
  annotation_scale(location = "bl", width_hint = 0.4) +
  coord_sf(xlim = c(-76.2, -74.6), ylim = c(38.5, 40.4))+
  scale_color_manual("Region",values=colpal)+
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

ggsave("DelawareReceivers.jpeg", plot = last_plot(), width = 5, height = 7,
       units = "in", dpi = 300, path = "C:/Users/RER/Documents/Masters UD/Code")

########################  Figure 2 - Length Stats #############################

# AVerage Lengths per year and sex
fsh %>%
  group_by(Deploy.Year) %>%
  summarise(TL.avg = mean(Length.m),
            TL.SD = sd(Length.m))


## movement across years
# Test fish  for normality          # We meet normality for both years!                     
shapiro.test(fsh$Length.m[fsh$Deploy.Year==2021])
shapiro.test(fsh$Length.m[fsh$Deploy.Year==2022])
# Test for differences in variance across years     # Variances are equal too. 
levene.test(fsh$Length.m, group = fsh$Deploy.Year)
# Time for a t test to test in length differences across years
# No difference in length between years :)
t.test(fsh$Length.m[fsh$Deploy.Year==2021], fsh$Length.m[fsh$Deploy.Year==2022])


# Box plot of fish length across years
fsh %>%
  ggplot(aes(x=Deploy.Year,y=Length.m))+
  geom_boxplot(notch = F, outlier.color = NA)+
  theme_classic() +
  labs(x = "Tagging Year", y="Length (m)") 
# geom_signif(comparisons = list(c("Female", "Male")), map_signif_level=F)

ggsave("LengthbyYEAR.png", plot = last_plot(), width = 3, height = 3,
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
  scale_color_manual(values=colpal) +
  scale_x_date(date_breaks = "month", date_labels = "%B", limits = as.Date(c("2021-05-1", "2021-07-01")))



ggsave("2021Detectss2.jpeg", plot = last_plot(), width = 11, height = 7, units = "in",
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
  scale_color_manual(values=colpal)

ggsave("2022Detectss2.jpeg", plot = last_plot(), width = 11, height = 7, units = "in",
       dpi = 300,  path = "C:/Users/RER/Documents/Masters UD/Code")





########################  Figure 5 and 6  - Residency ##########################


# Redo detects dataframe with new Region simplification
detects2.1 <- dtc_final2 %>%
  mutate(Region2 = recode(Region, "Brandywine 1" = "Brandywine", "Brandywine 2" = "Brandywine",
                          "Brandywine 3" = "Brandywine")) %>%
  mutate(Region3 = recode(Region2, "Canal" = "River", "Atlantic" = "Other", "Bay" = "Other")) %>%
  rename("deploy_lat" ="Latitude", "deploy_long" = "Longitude") %>%
  mutate(detection_timestamp_utc = as.POSIXct(detection_timestamp_utc)) %>%
  filter(!is.na(Region3)) 


# Prepare data set for GLATOS syntax - Compress all detections into unique detection events (BY REGION)
residency <- detection_events(detects2.1, location_col = "Region3", time_sep = Inf, condense = TRUE)
residency2 <- residency %>% mutate (res_time_min = res_time_sec/60, 
                                     res_time_hr = round(res_time_min/60,4),
                                     res_time_hr_log = log(res_time_hr)) %>%
  filter(!res_time_hr == 0)


# movement Normality across Regions
shapiro.test(residency2$res_time_hr[residency2$location=="Brandywine"])
# NOT normal
shapiro.test(residency2$res_time_hr[residency2$location=="Other"])
# NOT normal
shapiro.test(residency2$res_time_hr[residency2$location=="River"])
# NOT normal

# movement Normality across Regions using log transformed data
shapiro.test(residency2$res_time_hr_log[residency2$location=="Brandywine"])
# NORMAL
shapiro.test(residency2$res_time_hr_log[residency2$location=="Other"])
# NORMAL
shapiro.test(residency2$res_time_hr_log[residency2$location=="River"])
# NORMAL. PROCEED

# Levene test for equal variance using log transformed data
levene.test(residency2$res_time_hr_log, group = residency2$location)
# UNEQUAL VARIANCE. Proceed with non-parametric test.
# Kruskal Wallis test
kruskal.test(res_time_hr_log ~ location, data=residency2)
# Significant
# Dunn post hoc test to determine between group significance 
dunnTest(res_time_hr_log ~ location, data=residency2, method="bonferroni")
# barely signifcant between Brandywine and 'Other' Region


# Plot BR vs Other VS River
residency2 %>%
  mutate(Year = lubridate::year(first_detection)) %>%
  ggplot(aes(x=factor(location, levels = c("River", "Other", "Brandywine")), y=res_time_hr_log))+
  geom_boxplot(color = "black", coef = 3) +
  geom_jitter(width = .1, alpha = .2, color = "red", size = 3)+
  theme_bw() +
  labs(title = "", x= "", y= "Residency (log hour)")+
  #facet_wrap(~Year) +
  theme(axis.text.x = element_text(angle = 45, hjust =1),
        text = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_signif(comparisons = list(c("River", "Brandywine")), 
              map_signif_level=TRUE) 

res.stats <- residency2 %>%
  group_by(location) %>%
  summarise(mean = mean(res_time_hr),
            sd = sd(res_time_hr),
            min = min(res_time_hr),
            max = max(res_time_hr))


ggsave("ResidencyFig5.1.jpeg", plot = last_plot(), width = 6,  height = 6,
       units = "in", dpi = 300,  path = "C:/Users/RER/Documents/Masters UD/Code") 




### Sub set to look specifically at Brandywine Locations
detectsBR <- detects2.1 
residencyBR <- detection_events(detectsBR, location_col = "Region", time_sep = Inf, condense = TRUE)
residencyBR2 <- residencyBR %>% mutate (res_time_min = res_time_sec/60, 
                                        res_time_hr = round(res_time_min/60,4),
                                        res_time_hr_log = log(res_time_hr)) %>%
  filter_all(all_vars(!is.infinite(.))) %>%
  filter(str_detect(location, "Brandywine"))
shapiro.test(residencyBR2$res_time_hr_log[residencyBR2$location == "Brandywine 1"]) # NOT normal
shapiro.test(residencyBR2$res_time_hr_log[residencyBR2$location == "Brandywine 2"]) # NORMAL
shapiro.test(residencyBR2$res_time_hr_log[residencyBR2$location == "Brandywine 3"]) # NOT normal

kruskal.test(res_time_hr_log ~ location, data=residencyBR2)
dunnTest(res_time_hr_log ~ location, data=residencyBR2, method="bonferroni")

# Residency BR Stats
BR.Res.Stats <- residencyBR2 %>%
  group_by(location) %>%
  summarise(mean = mean(res_time_hr),
            sd = sd(res_time_hr),
            min = min(res_time_hr),
            max = max(res_time_hr))

# Plot Brandywine Residency Graph
ggplot(residencyBR2, aes(x=location, y=res_time_hr_log))+
  geom_boxplot(color = "black", coef = 100) +
  geom_jitter(width = .1, alpha = .2, color = "red", size = 3)+
  theme_bw() +
  labs(title = "", x= "", y= "Residency (log hour)")+
  #facet_wrap(~Year) +
  theme(axis.text.x = element_text(angle = 45, hjust =1),
        text = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_signif(comparisons = list(c("Brandywine 2", "Brandywine 3")), 
              map_signif_level=TRUE) +
  geom_signif(comparisons = list(c("Brandywine 1", "Brandywine 3")), 
              map_signif_level=TRUE) +
  geom_signif(comparisons = list(c("Brandywine 1", "Brandywine 2")), 
              map_signif_level=TRUE)

ggsave("BRResidencyFigure6.1.jpeg", plot = last_plot(), width = 6, height = 6,
       units = "in", dpi = 300, path = "C:/Users/RER/Documents/Masters UD/Code") 


########################  Figure 7 - Movement and time at large #################

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
  #filter(!animal_id %in% c(1,7)) %>%
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


ggsave("Total.mvmt2.1.2.png", plot = last_plot(), width = 6,height = 5, units = "in", dpi = 300,
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

ggsave("Total.time.2.2.png", plot = last_plot(), width = 6,height = 5, units = "in", dpi = 300,
       path = "C:/Users/RER/Documents/Masters UD/Code")


########################  Figure 8 - Acceleration Data  ######################


# Filter for sensors that contain motion data
ac.df <- dtc_final2 %>% 
  filter(transmitter_codespace == "A69-9007") %>% 
  drop_na(sensor_value)

# whats the maximum speed value? 3.5 m/s^2
max(ac.df$sensor_value)


# Summary statistics 
accel.stats <- ac.df %>%
  group_by(animal_id) %>%
  summarize(mean = mean(sensor_value),
            median = mean(sensor_value),
            min = min(sensor_value),
            max = max(sensor_value),
            mode = mode(sensor_value),
            sd = sd(sensor_value))

# plot acceleration 
ac.df %>%
  ggplot(aes(x = sensor_value, fill = factor(animal_id)))+
  geom_histogram(binwidth = .2) +
  facet_wrap(~Region, nrow = 1) +
  theme_bw() +
  labs(y = "Frequency", fill = "Fish ID") +
  #scale_fill_manual(values = colpal) +
  xlab(expression(Acceleration~Rate~(m/s^2))) +
  theme(panel.grid = element_blank()) +
  scale_y_continuous(expand = expansion(mult = c(0, .1)))

ac.df %>%
  ggplot(aes(x = sensor_value, fill = Region))+
  geom_histogram(binwidth = .2) +
  facet_wrap(~animal_id, nrow = 1) +
  theme_bw() +
  labs(y = "Frequency") +
  scale_fill_manual(values = colpal) +
  xlab(expression(Acceleration~Rate~(m/s^2))) +
  theme(panel.grid = element_blank()) +
  scale_y_continuous(expand = expansion(mult = c(0, .1)))


ggsave("Acceler.Histo.png", plot = last_plot(),width = 11, height = 7, 
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
  ggplot(aes(x=as.Date(detection_timestamp_utc), y=RKM, color=Region))+
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

ggsave("Fish1and7.1.jpeg", plot = last_plot(), width = 12,height = 6, units = "in",
       dpi = 300,  path = "C:/Users/RER/Documents/Masters UD/Code") 

####   Figure S4   ####

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
siteNumbers <- c(ChristinaConfluence, Brandywine, UpperBR)
parameters <- c(discharge, gageH, temp, rain, ph, DO)

# Raw daily data:
DailyData <- readNWISdv(siteNumbers, parameters, "1980-01-01", "2022-12-31")
DailyData <- renameNWISColumns(DailyData)

# Rename Columns 
DailyData$site_name[DailyData$site_no == "01480120"] <- "USGS Station 01480120 \n (Christina Confluence)"
DailyData$site_name[DailyData$site_no == "01481500"] <- "USGS Station 01481500 \n(Brandywine River)"
DailyData$site_name[DailyData$site_no == "01481000"] <- "USGS Station 01481000 \n(Brandywine River at Chadds Ford)"

# Calculcate historic averages to compare to study conditions
monthly.means <- DailyData %>%
  mutate(month = month(Date),
         year = year(Date)) %>%
  group_by(site_name, month, year) %>%
  summarize(temp.mean = mean(Wtemp, na.rm = T),
            DO.mean = mean(DO, na.rm = T),
            Flow.mean = mean(Flow, na.rm = T))
historic.mean <- monthly.means %>%
  group_by(site_name, month) %>%
  summarize(historic.month.mean.temp = mean(temp.mean, na.rm = T),
            historic.month.sd.temp = sd(temp.mean, na.rm = T),
            historic.month.mean.DO = mean(DO.mean, na.rm = T),
            historic.month.sd.DO = sd(DO.mean, na.rm = T),
            historic.month.mean.Flow = mean(Flow.mean, na.rm = T),
            historic.month.sd.Flow = sd(Flow.mean, na.rm = T))


# plot historic averages
ggplot(historic.mean, aes(x = month, y = historic.month.mean.temp))+
  geom_line() +
  geom_line(aes(x = month, y = historic.month.mean.temp-historic.month.sd.temp)) +
  geom_line(aes(x = month, y = historic.month.mean.temp+historic.month.sd.temp)) +
  facet_wrap(~site_name) +
  theme_bw() +
  labs(x = "Month", y = "Historic Average Water Temperature") +
  scale_x_continuous(breaks = c(1,3,5,7,9,11))
ggplot(historic.mean, aes(x = month, y = historic.month.mean.DO))+
  geom_line() +
  geom_line(aes(x = month, y = historic.month.mean.DO-historic.month.sd.DO)) +
  geom_line(aes(x = month, y = historic.month.mean.DO+historic.month.sd.DO)) +
  facet_wrap(~site_name) +
  theme_bw() +
  labs(x = "Month", y = "Historic Average Water DO") +
  scale_x_continuous(breaks = c(1,3,5,7,9,11))
ggplot(historic.mean, aes(x = month, y = historic.month.mean.Flow))+
  geom_line() +
  geom_line(aes(x = month, y = historic.month.mean.Flow-historic.month.sd.Flow)) +
  geom_line(aes(x = month, y = historic.month.mean.Flow+historic.month.sd.Flow)) +
  facet_wrap(~site_name) +
  theme_bw() +
  labs(x = "Month", y = "Historic Average Water Flow") +
  scale_x_continuous(breaks = c(1,3,5,7,9,11))


# Plot USGS Water Temp
DailyData %>%
  mutate(month = month(Date),
         year = year(Date)) %>%
  filter(year %in% c(2021, 2022)) %>%
ggplot(aes(Date, Wtemp)) +
  geom_point() +
  geom_path() +
  facet_wrap(~site_name)+
  theme_bw() +
  labs(y = "Water Temperature (°C)")

ggsave("WaterTemp.jpeg", plot = last_plot(), width = 12,height = 6, units = "in",
       dpi = 300,  path = "C:/Users/RER/Documents/Masters UD/Code") 

####   Figure S5   ####

# Plot USGS Water Flow
DailyData %>%
  filter(!site_no == "01480120") %>%
  mutate(month = month(Date),
         year = year(Date)) %>%
  group_by(site_name, month, year) %>%
  mutate(mean.flow = mean(Flow)) %>%
  filter(year %in% c(2021, 2022)) %>%
ggplot(aes(Date, (mean.flow))) +
  geom_point() +
  geom_path() +
  facet_wrap(~site_name)+
  theme_bw() +
  labs(y = "Discharge (ft/m^3)")

ggsave("WaterTemp.jpeg", plot = last_plot(), width = 12,height = 6, units = "in",
       dpi = 300,  path = "C:/Users/RER/Documents/Masters UD/Code") 

# Plot USGS Water DO
DailyData %>%
  filter(!site_no == "01480120") %>%
  mutate(month = month(Date),
         year = year(Date)) %>%
  filter(year %in% c(2021, 2022)) %>%
  ggplot(aes(Date, (DO))) +
  geom_point() +
  geom_path() +
  facet_wrap(~site_name)+
  theme_bw() +
  labs(y = "DO")

##############################################################################
BRWater<-DailyData %>%
  mutate(month = month(Date),
         year = year(Date)) %>%
  filter(year %in% c(2021, 2022)) %>%
  filter(site_no == "01481500")

dtc_date <- dtc_final2 %>%
  mutate(Date = ymd(as.Date(detection_timestamp_utc)))

dtc_water <- left_join(dtc_date, BRWater)

dtc_water %>%
  filter(!Year == 2020) %>%
ggplot(aes(x = as.Date(detection_timestamp_utc))) +
  stat_bin(binwidth=1, position="identity") + 
  scale_x_date(breaks=date_breaks(width="1 month")) +
  facet_wrap(~Year, scales ="free") +
  labs(x = "Date", y = "Number of detections") +
  theme_bw() +
  geom_line(aes(x = Date, y = Wtemp), color = "red")+
  geom_line(aes(x = Date, y = Flow), color = "blue")


dtc_water2 <- dtc_water %>%
  filter(!Year == 2020) %>%
  group_by(Date) %>%
  summarise(frq = n(),
            Flow = Flow,
            Wtemp = Wtemp,
            DO = DO) %>%
  distinct()


ccf(dtc_water2$Flow, dtc_water2$frq)
ccf(dtc_water2$DO, dtc_water2$frq)
ccf(dtc_water2$Wtemp, dtc_water2$frq)

dtc_water2_pvt <- dtc_water2 %>%
  pivot_longer(cols = c(Flow, DO, Wtemp, frq), names_to = "Variable", values_to = "Value") %>%
  mutate(Year = year(Date))


ggplot(dtc_water2_pvt, aes(x = Date, y = Value, color = Variable)) +
  geom_line(size =2) +
  facet_wrap(~Year, scales = "free")+
  theme_bw()


