#' @include internal.R
NULL

#' Simulate invasive species
#'
#' Simulate the expansion of a threat from differents initial points and propagations. By default,
#' the output will contain values between zero and one.
#'
#' @param x [`SpatialPolygonsDataFrame-class`] object to use as
#'    a template.
#'
#' @param csv_longitudinal_distance `data.frame-class` data frame object (optional) that contain information of
#'    longitudinal distances between sites (for example inside the river).
#'
#' @param initial_units `list` list of sites where the threat appear.
#'
#' @param type_of_prop `integer` type of propagation. This is represent for a radial propagation (1),
#'  downstream propagation (2), upstream propagation (3) or downstream and upstream propagation (4). The
#'  id's of planning units of \code{\link{csv_longitudinal_distance}} must be the same than the pu's from
#'  x object.
#'
#' @param expansion_speed `numeric` velocity of the propagation the threat.
#'
#' @param periods `integer` number of periods of simulation
#'
#' @param time_step `integer` number of periods needs to increase one level of the threat
#'
#' @param adjacency `logical` TRUE indicates that the propagation is at least to adjacent units
#'
#' @name simulate_invasive_specie
#'
#' @return [`data.frame-class`] object.
#'
#' @examples
#' \dontrun{
#' # simulate threat
#' path <- system.file("extdata/input_big/", package = "prioriactions")
#' x = raster::shapefile(paste0(path,"Fish_Mitchell.shp"))
#'
#' st <- simulate_invasive_specie(x, initial_units = c(150))
#'
#' }
#'
#' @export
simulate_invasive_specie <- function(x, csv_longitudinal_distance = NULL,
                                     initial_units = c(1),
                                     type_of_prop = 1,
                                     expansion_speed = 0,
                                     periods = 10,
                                     time_step = 1,
                                     adjacency = TRUE){
  
  
  #Checking inputs
  if(type_of_prop != 1){
    if(is.data.frame(csv_longitudinal_distance)){
      #Reading source function of longitudinal distances from github
      longitudinal_distance <- csv_longitudinal_distance
      colnames(longitudinal_distance) <- c("Var1", "Var2", "value")
      uniques_pu <- unique(longitudinal_distance$Var1)
      itself_distances <- data.frame(Var1 = uniques_pu,
                                     Var2 = uniques_pu,
                                     value = 0)
      longitudinal_distance <- rbind(longitudinal_distance, itself_distances)
    }
    else{
      stop("Incompatible longitudinal distance file (DataFrame required)", call. = FALSE)
    }
  }
  
  assertthat::assert_that(
    is.numeric(initial_units),
    all(type_of_prop %in% c(1,2,3,4)),
    is.numeric(expansion_speed),
    is.numeric(periods),
    is.numeric(time_step),
    is.logical(adjacency),
    inherits(x, "SpatialPolygonsDataFrame")
  )
  pus <- length(x$GridID)
  
  #Calculating radial distances
  centroids <- rgeos::gCentroid(x, byid = TRUE)
  radial_distance <- sp::spDists(centroids,longlat = FALSE)
  radial_distance <- reshape2::melt(radial_distance)
  
  #Calculating adjacency data
  adjacency_data <- rgeos::gTouches(x, byid = TRUE)
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
    threat_status$distance_to_focus[pu] <- 0.0
    threat_status$last_change[pu] <- 1
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
        if(id2 == -1){
          break
        }
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
