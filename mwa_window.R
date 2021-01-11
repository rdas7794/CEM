######################################################################
#------------------This file contains the code for ------------------#
#----------------CEM matching using momentum an trend----------------#
######################################################################

######################################################################
#------------------Installing packages and libraries------------------
######################################################################
options(java.parameters = "-Xmx8g")

library(mwa)
library(xtable)
library(cem)
library(GGally)
library(hrbrthemes)
library(viridis)
library(dplyr)
library(lubridate)
library(geosphere)
library(ggplot2)
library(cowplot)

require(cem)

######################################################################
#------------------Data  Loading and Pre-processing-------------------
######################################################################

setwd("C:/Users/Rohan Das/Documents/Rohan_Docs/Study/Semester 4/Vis Hiwi/jcr/rnr2/data")

data <- read.csv("indis2support.csv")

setwd("C:/Users/Rohan Das/Documents/Rohan_Docs/Study/Semester 4/Vis Hiwi/jcr/rnr2/results")
#Prep data
data <- data[,c(4:7,11:21)]
data <- na.omit(data)

d_timestep <- NULL

data$lon <- data$long
data$kind <- as.character(data$kind)
data$kind[data$kind == "Turn In" | data$kind == "Evidence Turn-in/Received" | data$kind == "ERW/Turn-in"] <- "Turn in"
data$side[data$kind == "Turn in"] <- "FRIEND"
data$timestamp <- paste(as.character(data$timestamp),"00:00:00")
data$event_type <- paste(data$kind, data$side, sep="-")

#Coding assistance to ISAF and types of violence
data$event_type[data$event_type == "Turn in-FRIEND"] <- "gb_assist"
data$event_type[data$event_type == "Indirect Fire-FRIEND" | data$event_type == "Close Air Support-FRIEND"] <- "br_indis"
data$event_type[data$event_type == "Direct Fire-FRIEND"] <- "br_sel"
data$event_type[data$event_type == "Indirect Fire-ENEMY" | data$event_type == "Mine Strike-ENEMY"] <- "rb_indis"
data$event_type[data$event_type == "Direct Fire-ENEMY"] <- "rb_sel"

#Generate subset of data to analyze 
data <- subset(data, event_type == "gb_assist" | event_type == "rb_indis" | event_type == "rb_sel" | event_type == "br_indis" | event_type == "br_sel")
data <- na.omit(data)

# Saving data timestamp and format to be added momentum & trend
d_timestep <- data$timestamp
data$timestamp <- format(as.Date(data$timestamp), "%d.%m.%Y")

#####################################################################
#-------Momentum & Trend calculation for Spatio-Temporal window------
#####################################################################

#Treatment events
treatment_data <- data %>%
  filter(event_type == "br_indis")

#Dependent_events
dependent_data <- data %>%
  filter(event_type == "gb_assist" )


treatment_data$momentum <- 0
treatment_data$trend <- 0

#Time & spatio Window initialization
delta_time <- 10
delta_s <- 10

for(i in 1:nrow(treatment_data)) {
  count <- 0
  treatment_row <- treatment_data[i,]
  previous_events <- NULL
  
  # Time Window Calculation
  time_window <- as.Date(as.character(treatment_row$timestamp), format = "%d.%m.%Y")
  
  time_window <- format(time_window - delta_time, "%d.%m.%Y")
  
  for(j in 1:nrow(dependent_data)) {
    dependent_row <- dependent_data[j,]
    
    # Calculating the spatial difference between treatment and dependent 
    geo_difference <- distm (c(treatment_row$long, treatment_row$lat),
                             c(dependent_row$long, dependent_row$lat),
                             fun = distHaversine)
    
    #Converting default meters to km
    geo_difference <- geo_difference/1000
    
    # The time and spatial window filtered
    if(dependent_row$timestamp < time_window && geo_difference < delta_s){
      # No. of dependent events before treatment event i.e. Momemntum
      count = count+1
      
      #All the dependent events before treatment events
      previous_events <- rbind(previous_events, dependent_row)
    }
    
  }
  
  #If there are only records in the previous events
  if(count>0){
    # Boundaries of the cylinder before the treatment events
    min_date <- as.Date(as.character(time_window), format = "%d.%m.%Y") 
    max_date <- as.Date(as.character(treatment_row$timestamp), format = "%d.%m.%Y")
    
    tmpTimes <- difftime(max_date, min_date)
    tmpTimes_mid <- tmpTimes/2
    
    mid_date <- as.Date(min_date) + tmpTimes_mid
    
    # Convert format of date for comparisons
    previous_events$timestamp <- as.Date(as.character(previous_events$timestamp), format = "%d.%m.%Y")
    
    # Trend calculation in terms of dependent events before conflicts 
    lower_cylinder <- 0
    upper_cylinder <- 0
    for(k in 1:nrow(previous_events)){
      mid_row <- previous_events[k,]
      
      if(mid_row$timestamp < mid_date){
        lower_cylinder = lower_cylinder +1
      }
      else{
        upper_cylinder = upper_cylinder +1
      }
    }
    
    # Trend Calculation
    trend <- upper_cylinder / lower_cylinder
    
    # Assign count of dependent events as momentum and trend 
    treatment_data[i, 18] <- count
    treatment_data[i, 19] <- trend
  }
  
}

data$momentum <- 0
data$trend <- 0

data$momentum[data$event_type == "br_indis"] <- treatment_data$momentum
data$trend[data$event_type == "br_indis"] <- treatment_data$trend

data$timestamp <- d_timestep

#########################################################################
#----------------------------Applying CEM--------------------------------
#########################################################################


# cem match: automatic bin choice
mat <- cem(treatment= "event_type", data=data)

#Adding matches to the data dataset
data$matches <- mat$matched

#########################################################################
#----------------------------Visualization--------------------------------
#########################################################################

exposure <- colnames(data[8:15])

hist_pak <- ggplot(data) +
  geom_histogram(aes(x=dist_pak, color = matches,  fill=matches),
                 breaks=mat$breaks$dist_pak,alpha=0.7, position="identity")+
  scale_x_continuous(breaks = mat$breaks$dist_pak)+ # Adds tick values to the plot
  theme(legend.position="top")

hist_pak <- hist_pak+scale_color_manual(values=c("#FFFFFF", "#FFFFFF"))+
  scale_fill_manual(values=c("#000000", "#FCF500"))

hist_pak
hist_kabul <- ggplot(data) +
  geom_histogram(aes(x=dist_kabul, color = matches,  fill=matches),
                 breaks=mat$breaks$dist_kabul,alpha=0.7, position="identity")+
  #scale_x_continuous(breaks = mat$breaks$dist_pak)+ # Adds tick values to the plot
  theme(legend.position="top")

hist_kabul <- hist_kabul+scale_color_manual(values=c("#FFFFFF", "#FFFFFF"))+
  scale_fill_manual(values=c("#000000", "#FCF500"))

hist_elevation <- ggplot(data) +
  geom_histogram(aes(x=hist_elevation, color = matches,  fill=matches),
                 breaks=mat$breaks$hist_elevation,alpha=0.7, position="identity")+
  #scale_x_continuous(breaks = mat$breaks$dist_pak)+ # Adds tick values to the plot
  theme(legend.position="top")

hist_elevation <- hist_elevation+scale_color_manual(values=c("#FFFFFF", "#FFFFFF"))+
  scale_fill_manual(values=c("#000000", "#FCF500"))

hist_elevation
hist_pop <- ggplot(data) +
  geom_histogram(aes(x=population, color = matches,  fill=matches),
                 breaks=mat$breaks$population,alpha=0.7, position="identity")+
  #scale_x_continuous(breaks = mat$breaks$population)+ # Adds tick values to the plot
  theme(legend.position="top")

hist_pop <- hist_pop+scale_color_manual(values=c("#FFFFFF", "#FFFFFF"))+
  scale_fill_manual(values=c("#000000", "#FCF500"))

hist_los <- ggplot(data) +
  geom_histogram(aes(x=lineofsight, color = matches,  fill=matches),
                 breaks=mat$breaks$lineofsight,alpha=0.7, position="identity")+
  #scale_x_continuous(breaks = mat$breaks$dist_pak)+ # Adds tick values to the plot
  theme(legend.position="top")

hist_los <- hist_los+scale_color_manual(values=c("#FFFFFF", "#FFFFFF"))+
  scale_fill_manual(values=c("#000000", "#FCF500"))

hist_gecon <- ggplot(data) +
  geom_histogram(aes(x=gecon, color = matches,  fill=matches),
                 breaks=mat$breaks$gecon,alpha=0.7, position="identity")+
  #scale_x_continuous(breaks = mat$breaks$dist_pak)+ # Adds tick values to the plot
  theme(legend.position="top")

hist_gecon <- hist_gecon+scale_color_manual(values=c("#FFFFFF", "#FFFFFF"))+
  scale_fill_manual(values=c("#000000", "#FCF500"))

# Merge different graphs in one
plot_grid(hist_pak, hist_kabul, hist_pop, hist_los, hist_gecon,
          labels = c('Distance to Pak', 'Dist to Kabul',
          'Population','LineofSight','Gecon'))

