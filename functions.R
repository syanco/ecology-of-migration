#Assign seasons to  data
assignSeasons <- function(x, start_breed, end_breed, start_wint, end_wint) {
  if (x >= start_breed & x <= end_breed) { 
    "summer"
  } else {
    if (x > end_breed & x < start_wint) { 
      "fall"
    } else {
      if (x >= start_wint | x <= end_wint) { 
        "winter"
      } else {
        if (x > end_wint & x < start_breed) {
          "spring" 
        } else {
          "unknown"
        }
      }
    }
  }
}

#function to turn scaling factors of a heterogeneous variance lme object (x) 
#into standard deviation or variance.
getHeteroSD <- function(x, variance = F) {
  if(variance == F) {
    out <- (c(1.0000000, 
              coef(x$modelStruct$varStruct, unconstrained=F))*x$sigma)
  } else {
    out <- (c(1.0000000, 
              coef(x$modelStruct$varStruct, unconstrained=F))*x$sigma)^2
  }
  return(out)
}

#vectorizable Var function
varCI <-function(x){
  ci <- quantile(replicate(10000, var(sample(na.omit(x), replace=TRUE))), 
                 probs=c(0.025, 0.975), na.rm = T) 
  return(ci)
}

bootMuCI <- function(x) {
  ci <- quantile(replicate(10000, mean(sample(na.omit(x), replace=TRUE))), 
                 probs=c(0.025, 0.975), na.rm = T) 
  return(ci)
}

#vectorizable SE function
se <- function(x) {
  sd(x)/sqrt(length(x))
}

#vectorizable mean function
samplemean <- function(x, d) {
  return(mean(x[d]))
}

#function to convert cell values to degrees celcius
nasaToCelcius <- function(temps) {
  temps*.02-273.15
}

#function to get last n characters from a string
getLastN <- function(string, n){
  substr(string, nchar(string)-n+1, nchar(string))
}

#extract spring moevements only
extractSpring <- function(birdID, move_data, doy_vec) {
  bird <- move_data[move_data$bird_year == birdID,]
  bird_spring <- bird[bird$event_doy >= as.numeric(min(doy_vec)) & 
                        bird$event_doy <= as.numeric(max(doy_vec)),]
  return(bird_spring)
}

#extract MODIS NDVI (and QA it) based on owl locations
getNDVI <- function(x, move_data, ndvi_brick, doy_index, QA_brick){
  #focal date
  date <- as.Date(move_data$date[x],
                  format = "%y-%m-%d")
  #find closest raster date
  difs <- index - date
  idx <- which(abs(difs) == min(abs(difs)))
  #pull raster for the first date
  r <- ndvi_brick[[idx[1]]]
  qar <- QA_brick[[idx[1]]]
  #pull first move location
  p <- as.matrix(move_data$geometry[[x]])
  #extract raster cell value at location
  
  qa <- raster::extract(qar, p)
  
  if(is.na(qa) == T) {
    return(NA)
    } else {
      if(qa <= 1) {
        cell <- raster::extract(r, p)
        return(cell)
        } else {
          return(NA)
        }
    }
}


#extract night temp data (and QA it) from MODIS based on tracking positions
getTemp <- function(x, move_data, temp_brick, doy_index, QA_brick){
  if(any(timesteps == as.Date(move_data$date[x],
                              format = "%y-%m-%d"))) {
    #pull raster for the first date
    r <- temp_brick[[which(doy_index == as.Date(move_data$date[x],
                                                format = "%y-%m-%d"))]]
    qar <- QA_brick[[which(doy_index == as.Date(move_data$date[x],
                                                format = "%y-%m-%d"))]]
    #pull first move location
    p <- as.matrix(move_data$geometry[[x]])
    
    qa <- raster::extract(qar, p)
    
    if(is.na(qa)) {
      cell <- raster::extract(r, p)
      return(cell)
    } else {
      if(qa == 0) {
        cell <- raster::extract(r, p)
        return(cell)
      } else {
        return(NA)
      }
    }
  } else {
    return(NA)
  }
  
}

#function to combine multiple sf objects into a single in a "rbind" type fashion
rbind_sf <- function(list) {
  if(length(list) > 0) {
    nogeo <- list
    for (i in 1:length(nogeo)) {
      st_geometry(nogeo[[i]]) <- NULL
    }  
    geos <- c()
    for (i in 1:length(list)) {
      geos[[i]] <- st_geometry(list[[i]])
    }
    b <- do.call(dplyr::bind_rows, nogeo)
    g<-do.call(c, geos)
    b_sf <- st_sf(b, geometry = g)
    return(b_sf)
  }
}

assignSeasonsInd <- function(df, time.df) {
  for(i in 1:nrow(df)) {
    start_breed <- time.df$first_on_terr_doy[time.df$bird == df$bird[i]]
    end_breed <- time.df$last_on_terr_doy[time.df$bird == df$bird[i]]
    start_wint <- time.df$wint_arr_doy[time.df$bird == df$bird[i]]
    end_wint <- time.df$wint_dep_doy[time.df$bird == df$bird[i]]
    if(is.na(start_breed)) {
      start_breed <- median(na.omit(time.df$first_on_terr_doy))}
    if(is.na(end_breed)) {
      end_breed <- median(na.omit(time.df$last_on_terr_doy))}
    if(is.na(start_wint)) {
      start_wint <- median(na.omit(time.df$wint_arr_doy))}
    if(is.na(end_wint)) {
      end_wint <- median(na.omit(time.df$wint_dep_doy))}
    x <- df$doy[i]
    if (x >= start_breed & x <= end_breed) { 
      df$season[i] <- "summer"
    } else {
      if (x > end_breed & x < start_wint) { 
        df$season[i] <- "fall"
      } else {
        if (x >= start_wint | x <= end_wint) { 
          df$season[i] <- "winter"
        } else {
          if (x > end_wint & x < start_breed) {
            df$season[i] <- "spring" 
          } else {
            df$season[i] <- "unknown"
          }
        }
      }
    }
  }
}

#make counterfactual simulated breeding locations based on winter location
makeHypoBreed <- function(bird, time_data, move_data, trackorigin, hyporigin,
                          thin = 1) {
  #get winter centroids
  arr <- time_data$wint_arr[time_data$id == bird]
  dep <- time_data$wint_dep[time_data$id == bird]
  if (is.na(arr)) { 
    arr <- as.Date(median(na.omit(time_data$wint_arr_doy)), origin = trackorigin)
    #arr <- median(as.Date(na.omit(time_data$wint_arr), format = "%m/%d/%Y"))
  }
  if (is.na(dep)) {
    dep <- as.Date(median(na.omit(time_data$wint_dep_doy)), origin = trackorigin)
  }
  
  wint_pts <- move_data[move_data$bird_year == bird &
                          move_data$date > as.Date(arr, format="%m/%d/%Y") &
                          move_data$date < as.Date(dep, format="%m/%d/%Y"),]
  
  if (nrow(wint_pts) == 0) {
    message("Move along - no points to interpolate.")
  } else {
    wint_pts <- st_combine(wint_pts)
    wint_cent <- suppressWarnings(st_centroid(wint_pts))
    #simulate data as though bird had stayed put after dep date
    
    ydayvec <- seq(yday(as.Date(dep, format = "%m/%d/%Y"))+1, 
                   yday(as.Date(arr, format = "%m/%d/%Y")), 
                   by = thin)
    hyp_dat <- st_sf("bird" = rep(bird, length(ydayvec)),
                     "event_doy" = ydayvec,
                     "date" = as.Date(ydayvec, origin = hyporigin),
                     "geometry" = st_geometry(wint_cent, length(ydayvec)))
    return(hyp_dat)}
}

#impute winter locations
makeWintVec <- function(bird, time_data, move_data, dist_thresh) {
  #define winter period
  start <- time_data$wint_arr[time_data$id == bird]
  stop <- time_data$wint_dep[time_data$id == bird]
  if (is.na(start)) { 
    start <- median(as.Date(na.omit(time_data$wint_arr), format = "%m/%d/%Y"))
  }
  if (is.na(stop)) {
    stop <- median(as.Date(na.omit(time_data$wint_dep), format = "%m/%d/%Y"))
  }
  wint_pts <- move_data[move_data$bird_year == bird &
                          move_data$date > as.Date(start, format="%m/%d/%Y") &
                          move_data$date < as.Date(stop, format="%m/%d/%Y"),]
  wint_pts <- wint_pts[order(as.Date(wint_pts$date)),]
  if (nrow(wint_pts) <= 1) {
    message("Move along - no points to interpolate.")
  } else {
    newdat <- list()
    for (i in 1:(nrow(wint_pts)-1)) {
      if (as.numeric(st_distance(x=wint_pts[i,], 
                                 y=wint_pts[i+1,])) <= dist_thresh) {
        if (length(seq(from=wint_pts$date[i]+1, 
                       to=(wint_pts$date[i+1]), 
                       by="day")) > 1) {
          dates <- seq(from=wint_pts$date[i]+1, to=(wint_pts$date[i+1]-1), by="day")
          newdat[[i]] <- st_sf("bird_year" = rep(bird, length(dates)),
                               "date" = dates,
                               "event_doy" = yday(dates),
                               "geometry" = st_geometry(st_centroid(
                                 x=wint_pts[i,],
                                 y=wint_pts[i+1,])),
                               length(dates))
        }
      }
    }
    combdat <- rbind_sf(plyr::compact(newdat))
    return(combdat)
  }
}

#impute breeding season locations
makeBreedVec <- function(bird, breed_loc, time_data, move_data, origin1, 
                         origin2 = NULL, thin = 1) {
  start <- time_data$first_on_terr_doy[time_data$id == bird]
  stop <- time_data$last_on_terr_doy[time_data$id == bird]
  if (is.na(start)) { 
    start <- median(na.omit(time_data$first_on_terr_doy))
  }
  if (is.na(stop)) {
    stop <- median(na.omit(time_data$last_on_terr_doy))
  }
  breed_dates <- seq(from = start, to = stop, by = thin)
  breed_dates <- breed_dates[!breed_dates %in% move_data$event_doy[
    move_data$bird_year == bird]]
  if (is.null(origin2)) {
    new_data <- data.frame("bird_year" = rep(bird, length(breed_dates)),
                           "event_doy" = breed_dates,
                           "date" = as.Date(breed_dates, origin = origin1),
                           "x" = rep(breed_loc$x[
                             breed_loc$bird_year == bird], 
                             length(breed_dates)),
                           "y" = rep(breed_loc$y[
                             breed_loc$bird_year == bird], 
                             length(breed_dates)))
  } else {
    new_data <- data.frame("bird_year" = rep(bird, 2*length(breed_dates)),
                           "event_doy" = c(breed_dates, breed_dates),
                           "date" = c(as.Date(breed_dates, origin = origin1),
                                      as.Date(breed_dates, origin = origin2)),
                           "x" = rep(breed_loc$x[
                             breed_loc$bird_year == bird], 
                             2*length(breed_dates)),
                           "y" = rep(breed_loc$y[
                             breed_loc$bird_year == bird], 
                             2*length(breed_dates)))
  }
  
  new_data_sf <- st_as_sf(new_data, 
                          coords = c("x", "y"),
                          crs = 26913) %>%
    st_transform(crs = 4326)
  return(new_data_sf)
}

#simulate counteractual winter locations
makeHypoWint <- function(bird, breed_loc, time_data, origin, thin = 1) {
  #if we have breeding arrival date...
  if(!is.na(time_data$last_mig[time_data$id == bird])){
    #...then set it + 1 year as the seq stop date
    stop <- as.Date(time_data$last_mig[time_data$id == bird],
                    format = "%m/%d/%Y")
  } else {
    #else take the median arrival DOY and set as th stop date 1 year after 
    #supplied origin
    stop <- as.Date(median(na.omit(time_data$first_on_terr_doy)),
                    origin = origin) + years(1)
  }
  
  #if we have a terr departure date, set as seq start
  if (!is.na(start <- time_data$first_mig[time_data$id == bird])) { 
    start <- as.Date(time_data$first_mig[time_data$id == bird],
                     format = "%m/%d/%Y")
  } else {
    #else acquire from sample median
    start <- as.Date(median(na.omit(time_data$first_mig_doy)), 
                     origin = origin)
  }
  #create a vector of days from start to stop
  breed_dates <- seq(from = start, to = stop, by = "day")
  breed_dates <- breed_dates[seq(1, length(breed_dates), by = thin)]
  
  #create new data frame with date, bird and location data
  new_data <- data.frame("bird_year" = rep(bird, length(breed_dates)),
                         "event_doy" = yday(breed_dates),
                         "date" = as.Date(breed_dates, origin = origin),
                         "x" = rep(breed_loc$x[
                           breed_loc$bird_year == bird], 
                           length(breed_dates)),
                         "y" = rep(breed_loc$y[
                           breed_loc$bird_year == bird], 
                           length(breed_dates)))
  
  #convert to spatial (sf) object and transform
  new_data_sf <- st_as_sf(new_data, 
                          coords = c("x", "y"),
                          crs = 26913) %>%
    st_transform(crs = 4326)
  
  return(new_data_sf)
}


#calculate variance on a rolling window
movingVar <- function(data, window = 10) {
  v <- c()
  for (i in (round(window)/2):(365-(round(window)/2))) {
    vec <- data$temp[data$doy > (i-(round(window)/2)) &
                       data$doy < (i + (round(window)/2))]
    if(is.null(vec)) {
      v[i] <- NA
    } else {
      v[i] <- var(na.omit(vec)) 
    }
  }
  return(v)
}


medfunc <- function(d, i) {
  d2 <- d[i]
  return(median(d2))
}

#Wald method to calc CIs for varianc estimates
varWaldCI <- function(x, level, df) {
  alph2 <- (1-level)/2
  qntls <- cbind(lower = qchisq(alph2, df, lower = F),
                 upper = qchisq(alph2, df))
  CI <- (x*parm/qntls)
  attr(CI, 'level') <- level
  return(CI)
}

#standalone legen function
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}