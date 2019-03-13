get_ereefs_slice <- function(var_names=c('Chl_a_sum', 'TN'),
                             location_latlon=data.frame(latitude=c(-20, -20), longitude=c(148.5, 149)),
                             target_date = c(2016, 02, 04),
                             input_file = "menu",
                             input_grid = NA,
                             eta_stem = NA,
                             robust = FALSE,
                             override_positive = FALSE)
{
  input_file <- substitute_filename(input_file)
  ereefs_case <- get_ereefs_case(input_file)
  input_stem <- get_file_stem(input_file)
  check_platform_ok(input_stem)
  
  # Remember which direction our slice should be facing
  if ((location_latlon[1,1]<location_latlon[2,1]) |
      ((location_latlon[1,1]==location_latlon[2,1])&(location_latlon[1,2]<location_latlon[2,2]))) {
    latlon_dir <- 1
  } else {
    latlon_dir <- -1
  }
  
  if (!is.na(eta_stem)) {
    if (stringi::stri_endswith(eta_stem, fixed='.nc')) eta_stem <- get_file_stem(eta_stem) 
  }
  grids <- get_ereefs_grids(input_file, input_grid)
  x_grid <- grids[['x_grid']]
  y_grid <- grids[['y_grid']]
  z_grid <- grids[['z_grid']] 
  nc <- ncdf4::nc_open(input_file)
  if (!is.null(nc$var[['latitude']])) {
    latitude <- ncdf4::ncvar_get(nc, 'latitude')
    longitude <- ncdf4::ncvar_get(nc, 'longitude')
  } else {
    latitude <- ncdf4::ncvar_get(nc, 'x_centre')
    longitude <- ncdf4::ncvar_get(nc, 'y_centre')
  }
  ncdf4::nc_close(nc)
  
  location_ll <- data.frame(latitude=NULL, longitude=NULL)
  location_edges <- data.frame(latitude=NULL, longitude=NULL)
  intersected <- NULL
  llind <- NULL
  for (i in 1:(dim(location_latlon)[1]-1)) {
    print(paste('transect section', i, 'of', dim(location_latlon)[1]-1))
    li <- find_intersections(location_latlon[i:(i+1),], x_grid, y_grid, latitude, longitude)
    if (dim(li[[1]])[1] > 0) {
      location_ll <- rbind(location_ll, li[[1]])
      location_edges <- rbind(location_edges, li[[2]])
      intersected <- rbind(intersected, li[[3]])
      llind <- c(llind, li[[4]])
    }
  }
  location_latlon <- location_ll
  
  if (dim(location_latlon)[1]==0) stop("The line segment given does not intersect any model cells on this grid.")
  i <- dim(location_latlon)[1]
  
  eta <- rep(NA, length(which(intersected)))
  botz <- rep(NA, length(which(intersected)))
  values <- array(NA, c(length(z_grid)-1, length(eta), length(var_names)))
  
  if (robust) {
    mydata <- get_ereefs_profile(var_names=var_names, location_latlon=as.numeric(location_latlon[1,c('latitude','longitude')]),
                                 start_date = target_date, end_date = target_date, 
                                 input_file = input_file, input_grid = input_grid, eta_stem = eta_stem, override_positive=override_positive)
    values[,1,] <- mydata$profiles
    eta[1] <- as.numeric(mydata$eta)
    botz[1] <- as.numeric(mydata$botz)
    grid_list <- mydata$grid_list
    
    for (i in 2:length(eta)) {
      print(paste('Extracting profile', i, 'of', length(eta)))
      mydata <- get_ereefs_profile(var_names=var_names, location_latlon=as.numeric(location_latlon[i,c('latitude','longitude')]),
                                   start_date = target_date, end_date = target_date, 
                                   input_file = input_file, input_grid = input_grid, eta_stem = eta_stem, override_positive=override_positive)
      values[,i,] <- mydata$profiles
      eta[i] <- as.numeric(mydata$eta)
      botz[i] <- as.numeric(mydata$botz)
    }
  } else {
    if (is.vector(target_date)) {
      target_date <- as.Date(paste(target_date[1], target_date[2], target_date[3], sep='-')) + 0.499
    } else if (is.character(target_date)) {
      target_date <- as.Date(target_date)
    }
    target_day <- as.integer(format(target_date, '%d'))
    target_month <- as.integer(format(target_date, '%m'))
    target_year <- as.integer(format(target_date, '%Y'))
    
    #var_list <- paste(var_names, collapse=",")
    
    location_grid <- cbind(floor((llind+dim(latitude)[1]-1)/dim(latitude)[1]),
                           (llind+dim(latitude)[1]-1)%%dim(latitude)[1] + 1)
    
    if (ereefs_case == 4) {
      input_file <- paste0(input_stem, format(as.Date(paste(target_year, target_month, 1, sep='-')), '%Y-%m'), 
                           '.nc')
      if (!is.na(eta_stem)) etafile  <- paste0(eta_stem, format(as.Date(paste(target_year, target_month, 1, sep='-')), '%Y-%m'), 
                                               '.nc')
    } else if (ereefs_case == 1) {
      input_file <- paste0(input_stem, format(as.Date(paste(target_year, target_month, target_day, sep='-')), '%Y-%m-%d'), 
                           '.nc')
      if (!is.na(eta_stem)) etafile <- paste0(eta_stem, format(as.Date(paste(target_year, target_month, target_day, sep='-')), '%Y-%m-%d'), 
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
    from_record <- which.min(abs(ds - (target_date)))
    eta_from_record <- which.min(abs(eta_ds - (target_date)))
    ########
    