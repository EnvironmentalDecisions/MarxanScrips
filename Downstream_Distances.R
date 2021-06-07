library(rgdal) #To use the readOGR function
library(ggplot2) #To use the fortify function

#Function to create longitudinal distances-----------------------------------------------------------------
Create_longitudinal_distances<-function(file, limit_dist = .Machine$integer.max, 
                                        limit_connect = .Machine$integer.max,
                                        export_csv = FALSE){
  #checking if the input is a shapefile
  if(is(file, "SpatialLinesDataFrame") || is(file, "SpatialPolygonsDataFrame")){
    df = fortify(file@data)
    
    #checking if the variables exists
    var = c("GridID","HydroID","NextDownID","Length")
    if(!all(var %in% colnames(df))){
      var_needed <- which(!var %in% colnames(df))
      stop(paste("The following variables are needed: ", var[var_needed]),call. = FALSE)
    }
    
    #changing hydros values by grids values
    if(df$GridID[1] != df$HydroID[1])
    {
      df<-df[order(df$GridID),]
      df$NextDownID<-match(df$NextDownID,df$HydroID)
      df$NextDownID[is.na(df$NextDownID)] <- -1 
      df<-df[,var]
    }
  }
  else if(is.data.frame(file)){
    df = file
  }
  else{
    stop("Incompatible file (SpatialLinesDataFrame or DataFrame required)", call. = FALSE)
  }
  
  #creating empty data frame to distances
  length_df <- length(df$GridID)
  data_distances<-data.frame(id=integer(0),
                             down=integer(0),
                             distance=double(0))  
  
  row_distances_data = 1
  for(i in 1:length_df)
  {
    #internal parameters
    row_df = i
    count_distance = 0
    connections = 0
    
    while(df$NextDownID[row_df] != -1)
    {
      #stop by number of connections
      if(connections == limit_connect || count_distance + df[row_df,3] >= limit_dist){
        break
      }
      
      #adding next down cell to distances data frame
      data_distances[row_distances_data,] <- list(df$GridID[i], df$NextDownID[row_df], count_distance + df[row_df,3])
      
      count_distance <- count_distance + df[row_df,3]
      row_df <- which(df$GridID==df$NextDownID[row_df])[1]
      
      
      row_distances_data = row_distances_data + 1
      connections = connections + 1
      
    }
    #adding last connection
    if(df$NextDownID[row_df] == -1 && connections != limit_connect){
      data_distances[row_distances_data,] <- list(df$GridID[i], df$NextDownID[row_df], count_distance + df[row_df,3])
      
      row_distances_data = row_distances_data + 1
    }
  }
  rownames(data_distances)<-c(1:nrow(data_distances))
  
  if(export_csv == FALSE){
    return(data_distances)
  }
  else{
    write.csv(data_distances, file="Longitudinal_distance.csv", row.names = FALSE)
  }
}

#How to use---------------------------------------------------------
#From .csv
#The .csv file must contain the following variables: GridID, NextDownID, Length

##Name_data_csv  = "Douro2.csv"
##file = read.csv(Name_data_csv)

#From shapefile data
#The shapefile must contain the following variables: GridID, HydroID, NextDownID, Length
##Name_data_shp = "River_500_Albers.shp"
##file = readOGR(Name_data_shp)

#Function
#limit_dist = maximum distance allowed between units
#limit_connect = maximum connections allowed between units
##Create_longitudinal_distances(file)



