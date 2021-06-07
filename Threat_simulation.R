library(rgeos) #to use the gCentroid and gTouches functions
library(devtools) #to use the source_url function
library(reshape2) #to use melt function
library(sp) #to use spDists function
library(dplyr) #to use left_join function
library(ggplot2) #To use the fortify function
library(rgdal) #To use the readOGR function
source_url("https://raw.githubusercontent.com/EnvironmentalDecisions/MarxanScrips/main/Downstream_Distances.R")

#Function to simulate distribution of threats--------------------------------------
Create_threat_distribution <- function(file, csv_longitudinal_distance,
                                       initial_units,
                                       type_of_prop = 1,
                                       expansion_speed = 0,
                                       periods = 10,
                                       time_step = 1,
                                       adjacency = TRUE){
  
  #Checking inputs
  if(is(file, "SpatialPolygonsDataFrame")){}
  else{
    stop("Incompatible file (SpatialPolygonsDataFrame required)", call. = FALSE)
  }
  if(is.data.frame(csv_longitudinal_distance)){}
  else{
    stop("Incompatible longitudinal distance file (DataFrame required)", call. = FALSE)
  }
  
  assertthat::assert_that(
    is.numeric(initial_units),
    all(type_of_prop %in% c(1,2,3,4)),
    is.numeric(expansion_speed),
    is.numeric(periods),
    is.numeric(time_step),
    is.logical(adjacency)
  )
  pus <- length(file$GridID)
  
  #Reading source function of longitudinal distances from github
  longitudinal_distance <- csv_longitudinal_distance
  colnames(longitudinal_distance) <- c("Var1", "Var2", "value")
  uniques_pu <- unique(longitudinal_distance$Var1)
  itself_distances <- data.frame(Var1 = uniques_pu,
                                 Var2 = uniques_pu,
                                 value = 0)
  longitudinal_distance <- rbind(longitudinal_distance, itself_distances)
  
  #Calculating radial distances
  centroids <- rgeos::gCentroid(shp_mitchell, byid = TRUE)
  radial_distance <- sp::spDists(centroids,longlat = FALSE)
  radial_distance <- reshape2::melt(radial_distance)
  
  #Calculating adjacency data
  adjacency_data <- rgeos::gTouches(shp_mitchell, byid = TRUE)
  rownames(adjacency_data) <- seq(1:pus)
  colnames(adjacency_data) <- seq(1:pus)
  adjacency_data <- reshape2::melt(adjacency_data)
  adjacency_data <- adjacency_data[adjacency_data$value == TRUE,]
  colnames(adjacency_data) <- c("Var1", "Var2", "adj")
  
  
  #Simulation-----------------------------------------------------
  threat_status <- data.frame(id=1:pus,
                              state = 0,
                              level = 0,
                              distance_to_focus = NA,
                              last_change = 0)
  
  #pu is the differents initial points where the threats appear
  for(pu in initial_units){
    threat_status$state[pu] <- 1
    threat_status$level[pu] <- 1
    threat_status$distance_to_focus <- 0
    threat_status$last_change <- 1
  }
  
  #Setting data frame of connections depending of type of propagation
  if(type_of_prop == 1){
    connections <- radial_distance
  }else if(type_of_prop == 2){
    connections <- longitudinal_distance
  }else if(type_of_prop == 3){
    connections <- longitudinal_distance
    colnames(connections) <- c("Var2","Var1","value")
  }else{
    connections <- longitudinal_distance
    connections_pu_aux <- longitudinal_distance[longitudinal_distance$Var2 == pu,]
    colnames(connections_pu_aux) <- c("Var2","Var1","value")
    connections <- rbind(connections, connections_pu_aux)
  }
  connections <- dplyr::left_join(connections, adjacency_data, by = c("Var1" = "Var1", "Var2" = "Var2"))
  
  
  #Running simulation
  
  for(period in 2:periods){
    units_threated <- which(threat_status$state > 0)
    
    for(pu in units_threated){
      connections_pu <- connections[connections$Var1 == pu,]
      connections_pu <- connections_pu[order(connections_pu$value),]
      
      i = 1
      while(connections_pu$value[i] <= expansion_speed || (isTRUE(connections_pu$adj[i] && adjacency == TRUE)))
      {
        id2 <- connections_pu$Var2[i]
        
        #The unit is only evaluated if it didn't change its level this period
        if(threat_status$last_change[id2] != period){
          if(threat_status$level[id2] != 0){
            if(period - threat_status$last_change[id2] == time_step){
              threat_status$level[id2] <- threat_status$level[id2] + 1
              threat_status$last_change[id2] <- period
            }
          }
          else{
            threat_status$state[id2] <- 1
            threat_status$level[id2] <- 1
            threat_status$last_change[id2] <- period
            
            #Calculating the minimun distance to focus units
            connections_pu_id2 <- connections[connections$Var2 == id2,]
            connections_pu_id2 <- connections_pu_id2[connections_pu_id2$Var1 %in% initial_units,]
            
            threat_status$distance_to_focus[id2] <- connections_pu_id2$value[which.min(connections_pu_id2$value)]
          }
        }
        i <- i + 1
        if(i > nrow(connections_pu))
        {
          break
        }
      }
    }
  }
  
  return(threat_status)
}

#Parameters---------------------------------------------------------------------

#type_of_prop:
# type of propagation (1,2,3, and 4) described below
# 1 = radial
# 2 = downstream
# 3 = upstream
# 4 = both

#expansion_speed:
# numeric parameter that indicates the velocity of propagation of the threat in a period 
# (mts/period)

#initial_units:
# list parameter of units where the threat appear

#time_step:
# Number of periods that the threat needs to increase its level by 1 on an unit 

#adjacency: 
# Boolean parameter (TRUE/FALSE) that indicates if the propagation is do it between 
# adjacent units. This could be use with complement the expansion_speed parameter

#How to use--------------------------------------------------------------------

path <- system.file("extdata/input_big/", package = "prioriactions")
shp_mitchell = raster::shapefile(paste0(path,"Fish_Mitchell.shp"))

file = read.csv("Longitudinal_distance.csv",)

threat1 <- Create_threat_distribution(shp_mitchell, csv_longitudinal_distance = file, initial_units = c(200))


threat1$distance_to_focus <- threat1$distance_to_focus/1000

threat1$distance_to_focus[1] <- NA
threat1$d <- 1/threat1$distance_to_focus


shp_mitchell$t1 <- threat1$d

tmap::tm_shape(shp_mitchell) +
  tmap::tm_fill("t1", pal = c("white", "dodgerblue4")) +
  tmap::tm_borders(col="black", lwd = 0.5)


