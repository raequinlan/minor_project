#' Extract vertical profiles of specified variables  from a specified latitude and longitude over a specified time-period from an eReefs or other EMS netcdf file.
#'
#' See also plot_ereefs_profile(), which relies on output from this function.
#'
#' Note that this function assumes consistent frequency of model output, even if the time-series extends across multiple output files (e.g.
#' multiple months of eReefs output). If the data in eta_stem is output on a different interval from the data in input_file, the function
#' will do its best, but surface elevation estimates may not exactly match the time-stamps in the main input file.
#'
#' @return a list containing a vector of dates, an array of surface elevations (eta), the vertical grid (z_grid) and a data frame of values.
#' @param var_name A vector of EMS variable names. Defailts to c('Chl_a_sum', 'TN'))
#' @param location_latlon Latitude and longitude of location to extract.  Defaults to c(-23.39189, 150.88852)
#' @param start_date Date from which to start extraction. Can be a date, or text formatted for as.Date(), or a (year, month, day) vector.
#'                   Defaults to c(2016, 02, 04). If start_date is a vector, 0.499 is added to the calculated date to bring the start 
#'                   as close to midday as possible.
#' @param end_date Date on which to end extraction, specified as for start_date. Defaults to c(2016, 03, 02).
#' @param input_file is the URI or file location of any of the EMS output files, 
#'        Defaults to a menu selection. Set to "choices" to see some other pre-defined options that
#'        can be used (codenames as used in https://research.csiro.au/ereefs/models/model-outputs/access-to-raw-model-output/ )
#' @param input_grid Either a list containing the coordinates of the cell corners (x_grid, y_grid and z_grid) or the name of the                                                        
#'      locally-stored or opendap-served netcdf file that contains these. If not specified, the function will first look for                                                            
#'      z_grid can be found in the first INPUT_STEM file, and if not found, will check whether the size                                                                                 
#'      of the variables in the input file corresponds to the size expected for GBR4 or GBR1, and load an appropriate                                                                   
#'      z grid from data files stored in this package. Alternatively, you can provide the location of a full                                                                            
#'      (not simple-format) ereefs netcdf output file such as                                                                                                                           
#'      "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_hydro_all/gbr4_all_2016-09.nc". 
#' @param eta_stem The URI or file location of the model output files that contains the surface elevation (eta), or the stem of that
#'       filename minus the
#'       date components of the filename in the case of GBR1 or GBR4 files, and ommitting the file extension, ".nc". Needed
#'       only if eta is not in the files indicated by input_stem (e.g. some GBR1 bgc files).
#' @param squeeze Whether to reduce the number of dimensions in the output profiles array if there is only one variable and/or 
#'       only one time-step. Default TRUE.
#' @export
get_ereefs_profile <- function(var_names=c('Chl_a_sum', 'TN'),
                               location_latlon=c(-23.39189, 150.88852),
                               start_date = c(2016, 02, 04),
                               end_date = c(2016, 03, 02),
                               input_file = "menu",
                               input_grid = NA,
                               eta_stem = NA,
                               squeeze = TRUE,
                               override_positive=FALSE)
{
  input_file <- substitute_filename(input_file)
  ereefs_case <- get_ereefs_case(input_file)
  input_stem <- get_file_stem(input_file)
  check_platform_ok(input_stem)
  if (!is.na(eta_stem)) {
    if (stringi::stri_endswith(eta_stem, fixed='.nc')) eta_stem <- get_file_stem(eta_stem) 
  }
  z_grid <- get_ereefs_grids(input_file, input_grid)[['z_grid']]
  nc <- ncdf4::nc_open(input_file)
  if (!is.null(nc$var[['latitude']])) {
    latitude <- ncdf4::ncvar_get(nc, 'latitude')
    longitude <- ncdf4::ncvar_get(nc, 'longitude')
  } else {
    latitude <- ncdf4::ncvar_get(nc, 'x_centre')
    longitude <- ncdf4::ncvar_get(nc, 'y_centre')
  }
  ncdf4::nc_close(nc)
  
  # Dates to plot
  if (is.vector(start_date)) {
    start_date <- as.Date(paste(start_date[1], start_date[2], start_date[3], sep='-')) + 0.499
  } else if (is.character(start_date)) {
    start_date <- as.Date(start_date)
  }
  start_day <- as.integer(format(start_date, '%d'))
  start_month <- as.integer(format(start_date, '%m'))
  start_year <- as.integer(format(start_date, '%Y'))
  
  if (is.vector(end_date)) {
    end_date <- as.Date(paste(end_date[1], end_date[2], end_date[3], sep='-')) + 0.499
  } else if (is.character(end_date)) {
    end_date <- as.Date(end_date)
  }
  end_day <- as.integer(format(end_date, '%d'))
  end_month <- as.integer(format(end_date, '%m'))
  end_year <- as.integer(format(end_date, '%Y'))
  
  if (start_date > end_date) {
    stop('start_date must preceed end_date')
  }
  
  if (start_year==end_year) {
    mths <- start_month:end_month
    years <- rep(start_year, length(mths))
  } else if ((start_year + 1) == end_year) {
    mths <- c(start_month:12, 1:end_month)
    years <- c(rep(start_year, 12 - start_month + 1), rep(end_year, end_month))
  } else {
    mths <- c(start_month:12, rep(1:12, end_year - start_year - 1), 1:end_month)
    years <- c(rep(start_year, 12 - start_month + 1), 
               rep((start_year + 1) : (end_year - 1), each=12),
               rep(end_year, end_month))
  }
  
  #var_list <- paste(var_names, collapse=",")
  
  if (ereefs_case == 4) {
    input_file <- paste0(input_stem, format(as.Date(paste(start_year, start_month, 1, sep='-')), '%Y-%m'), 
                         '.nc')
    if (!is.na(eta_stem)) etafile  <- paste0(eta_stem, format(as.Date(paste(start_year, start_month, 1, sep='-')), '%Y-%m'), 
                                             '.nc')
  } else if (ereefs_case == 1) {
    input_file <- paste0(input_stem, format(as.Date(paste(start_year, start_month, start_day, sep='-')), '%Y-%m-%d'), 
                         '.nc')
    if (!is.na(eta_stem)) etafile <- paste0(eta_stem, format(as.Date(paste(start_year, start_month, start_day, sep='-')), '%Y-%m-%d'), 
                                            '.nc')
  } else {
    input_file <- input_file
    if (!is.na(eta_stem)) etafile <- paste0(eta_stem, '.nc')
  }
  nc <- ncdf4::nc_open(input_file)
  if (!is.null(nc$var[['t']])) { 
    ds <- as.Date(ncdf4::ncvar_get(nc, "t"), origin = as.Date("1990-01-01"))
  } else {
    ds <- as.Date(ncdf4::ncvar_get(nc, "time"), origin = as.Date("1990-01-01"))
  }
  if (!is.na(eta_stem)) {
    nc3 <- ncdf4::nc_open(etafile)
    if (!is.null(nc3$var[['t']])) { 
      eta_ds <- as.Date(ncdf4::ncvar_get(nc3, "t"), origin = as.Date("1990-01-01"))
    } else {
      eta_ds <- as.Date(ncdf4::ncvar_get(nc3, "time"), origin = as.Date("1990-01-01"))
    }
  } else {
    eta_ds <- ds
  }
  if (length(ds)==1) {
    blank_length <- as.numeric(end_date - start_date + 1)
  } else {
    blank_length <- as.numeric(end_date - start_date + 1) / as.numeric(ds[2] - ds[1])
  }
  # Initialise the data frame with the right number of NAs
  blanks <- rep(NA, blank_length)
  dates <- as.Date(blanks)
  values <- array(blanks, dim=c(length(z_grid)-1, length(var_names), length(blanks)))
  colnames(values) <- var_names
  eta_record <- blanks
  
  if (is.integer(location_latlon)) {
    location_grid <- location_latlon
  } else { 
    # Find the nearest grid-points to the sampling location
    grid_ind <- (latitude - location_latlon[1])^2 + (longitude - location_latlon[2])^2 
    grid_ind <- which.min(grid_ind) 
    location_grid <- c(floor((grid_ind+dim(latitude)[1]-1)/dim(latitude)[1]),
                       (grid_ind+dim(latitude)[1]-1)%%dim(latitude)[1] + 1)
  }
  
  zat <- ncdf4::ncatt_get(nc, "botz")
  if (!is.null(zat$positive)) {
    if (zat$positive=="down") zsign <- -1 else zsign <- 1
    if (override_positive) zsign <- -zsign
  } else {
    zsign <-1
  }
  botz <- zsign * as.numeric(ncdf4::ncvar_get(nc, "botz", start=c(location_grid[2], location_grid[1]), count=c(1,1)))
  ncdf4::nc_close(nc)
  
  # Loop through monthly eReefs files to extract the data
  i <- 0
  mcount <- 0
  pb <- txtProgressBar(min = 0, max = length(mths), style = 3)
  for (month in mths) {
    mcount <- mcount + 1
    year <- years[mcount]
    if (mcount == 1) {
      from_record <- which.min(abs(ds - (start_date)))
      eta_from_record <- which.min(abs(eta_ds - (start_date)))
    } else {
      from_record <- 1
      eta_from_record <- 1
    }
    if ((start_year==end_year)&&(start_month==end_month)) {
      day_count <- end_day - start_day + 1
    } else if (mcount == 1) {
      day_count <- daysIn(as.Date(paste(year, month, 1, sep='-'))) - start_day + 1
    } else if (mcount == (length(mths))) {
      day_count <- end_day
    } else {
      day_count <- daysIn(as.Date(paste(year, month, 1, sep='-')))
    }
    if (ereefs_case == 4) { 
      fileslist <- 1
      input_file <- paste0(input_stem, format(as.Date(paste(year, month, 1, sep="-")), '%Y-%m'), '.nc')
      if (!is.na(eta_stem)) etafile <- paste0(eta_stem, format(as.Date(paste(year, month, 1, sep="-")), '%Y-%m'), '.nc')
    } else if (ereefs_case == 1) {
      fileslist <- 1:day_count
      day_count <- 1
    } else {
      fileslist <- 1
    }
    
    if ((start_date == end_date)|(length(ds)==1)) { # Assume we only want a single profile
      day_count <- 1
      eta_day_count <- 1
    } else if (length(ds)>1) {
      day_count <- day_count / as.numeric(ds[2]-ds[1])
      eta_day_count <- day_count / as.numeric(eta_ds[2]-eta_ds[1])
    }
    
    for (dcount in fileslist) {
      if (ereefs_case == 1) {
        input_file <- paste0(input_stem, format(start_date+dcount-1, sep="-", '%Y-%m-%d'), '.nc')
        if (!is.na(eta_stem)) etafile <- paste0(eta_stem, format(start_date+dcount-1, sep="-", '%Y-%m-%d'), '.nc')
      }
      #input_file <- paste0(input_file, '?', var_list, ',time,eta')
      nc <- ncdf4::nc_open(input_file)
      if (!is.na(eta_stem)) nc3 <- ncdf4::nc_open(etafile)
      if (ereefs_case == 0) {
        d <- ncdf4::ncvar_get(nc, "time", start=from_record, count=day_count)
        d <- as.Date(d, origin = as.Date("1990-01-01"))
      } else {
        d <- ds[from_record:(from_record + day_count - 1)]
      }
      if (!is.null(nc$var[['eta']])) { 
        eta <- ncdf4::ncvar_get(nc, 'eta', start=c(location_grid[2], location_grid[1],from_record), count=c(1,1,day_count))
      } else {
        if (is.na(eta_stem)) stop('eta not found in netcdf file. Please specify eta_stem.')
        eta <- ncdf4::ncvar_get(nc3, 'eta', start=c(location_grid[2], location_grid[1],eta_from_record), count=c(1,1,eta_day_count))
      }
      im1 = i+1
      i <- i + length(d)
      if (length(eta_ds)==length(ds)) {
        eta_record[im1:i] <- eta
      } else {
        if (length(eta_ds)<length(ds)) {
          if (dcount==1) warning(paste('Surface elevation (eta) in', etafile, 'is output less frequently than', var_names[1], 'in', input_file,
                                       '. Assuming eta always==0, though this is unlikely'))
          eta_record[im1:i] <- 0*c(im:i)
        } else {
          if (length(ds)==1) {
            if (start_date==end_date) {
              eta_record[im1:i] <- eta
              ind <- 1
            } else { 
              ind <- which.min(abs(eta_ds-ds[1])) 
            }
            if ((dcount==1)&((eta_ds[ind]-ds[1])>(1/48))) warning(paste('Surface elevation (eta) in', etafile, 'is output more frequently than', var_names[1], 'in', input_file,
                                                                        '. Using eta from closest times to output of', var_names[1], ', which is more than 30 mins away from desired time.'))
          } else {
            interval <- as.numeric(ds[2] - ds[1])/as.numeric(eta_ds[2] - eta_ds[1])
            ind <- seq(from=which.min(abs(eta_ds-ds[1])), to=length(eta), by=interval)
            if (dcount==1) {
              tgap <- abs(eta_ds[ind] - ds)
              if (max(tgap)>(1/48)) warning(paste('Surface elevation (eta) in', etafile, 'is output more frequently than', var_names[1], 'in', input_file,
                                                  '. Using eta from closest times to output of', var_names[1], ', which is sometimes more than 30 mins away from desired time.'))
            }
            eta_record[im1:i] <- eta[ind]
          }
        }
      }
      dates[im1:i] <- d
      z <- array(z_grid[2:length(z_grid)], dim=c(length(z_grid)-1, length(d)))
      zm1 <- array(z_grid[1:(length(z_grid)-1)], dim=c(length(z_grid)-1, length(d)))
      eta2 <- t(array(eta, dim=c(length(d), length(z_grid)-1)))
      wet <- (eta2 > zm1) & (z > botz)           # There is water in this layer
      dry <- !wet                                # There is no water in this layer
      
      for (j in 1:length(var_names)) {
        wc <- ncdf4::ncvar_get(nc, var_names[j], start=c(location_grid[2],location_grid[1],1,from_record), count=c(1,1,-1,day_count))
        wc[dry] <- NA
        if (dim(z)[2] == 1) wc <- array(wc, dim=dim(z))
        values[1:(length(z_grid)-1), j, im1:i] <- wc 
      }
      ncdf4::nc_close(nc)
      if (length(eta_stem)>1) ncdf4::nc_close(nc3)
      setTxtProgressBar(pb,mcount)
    }
  }
  close(pb)
  
  if (squeeze&(dim(values)[3] == 1)) {                          # Only one time-step
    values <- array(values, dim=dim(values)[c(1,2)])
    colnames(values) <- var_names
  }
  if (squeeze&(dim(values)[2] == 1)&(length(dim(values))==3)) { # Only one variable, but multiple time-steps
    values <- array(values, dim=dim(values)[c(1,3)])
    if (all(is.na(values))) warning('No wet cells in this profile. Either the location is a land cell or the positive attribute of botz is incorrect (use override_positive=TRUE) if this is the case)')
    return_list <- list(dates=dates, eta=eta_record, z_grid=z_grid, botz=botz, profiles=values)
    names(return_list) <- c(names(return_list)[1:4], var_names)
  } else { 
    if (all(is.na(values))) warning('No wet cells in this profile. Either the location is a land cell or the positive attribute of botz is incorrect (use override_positive=TRUE) if this is the case)')
    return_list <- list(dates=dates, eta=eta_record, z_grid=z_grid, botz=botz, profiles=values)
  }
  return(return_list)
}
