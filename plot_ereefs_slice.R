#' Produces a coloured tile plot of a vertical slice already fetched from an eReefs or other EMS netcdf file.
#'
#' Relies on output from get_ereefs_slice().
#'
#' @param slice A list object as output by get_ereefs_slice(), containing dates, eta, z_grid, botz,
#'              a data frame of values and a data frame of latitudes and longitudes
#' @param scale_col Colours to use for low and high values. Default c("ivory", "hotpink").
#' @param scale_lim values for low and high limits of colourscale. Defaults to full range.
#' @return p handle for the generated figure
#' @export
plot_ereefs_slice <- function(slice, var_name='Chl_a_sum', scale_col=c("ivory", "navy"), scale_lim=NA) {
  numprofiles <- dim(slice$values)[2]
  layers <- length(slice$z_grid) - 1
  zmin <- array(slice$z_grid[1:layers], c(layers, numprofiles))
  zmax <- array(slice$z_grid[2:(layers+1)], c(layers, numprofiles))
  d <- rep(NA, numprofiles)
  for (i in 1:numprofiles) {
    zmin[zmin[,i]<slice$botz[i],i] <- slice$botz[i]
    zmin[zmin[,i]>slice$eta[i],i] <- slice$eta[i]
    zmax[zmax[,i]<slice$botz[i],i] <- slice$botz[i]
    zmax[zmax[,i]>slice$eta[i],i] <- slice$eta[i]
    d[i] <- earth.dist(slice$cell_intersections[i,'longitude'],slice$cell_intersections[i,'latitude'], slice$cell_intersections[i+1,'longitude'], slice$cell_intersections[i+1,'latitude']) 
  }
  
  d <- cumsum(d)
  dmin <- c(0, d[1:(numprofiles-1)])
  dmax <- d[1:numprofiles]
  dmin <- t(array(dmin, c(numprofiles, layers)))
  dmax <- t(array(dmax, c(numprofiles, layers)))
  
  ind <- which(!is.na(slice$values[, , var_name]))
  values <- slice$values[,, var_name]
  if (length(scale_lim)==1) {
    scale_lim[1] <- min(c(values[ind]))
    scale_lim[2] <- max(c(values[ind]))
  }
  
  mydata <- data.frame(xmin=dmin[ind], xmax=dmax[ind], ymin=zmin[ind], ymax=zmax[ind], z=values[ind])
  p <- ggplot2::ggplot(data=mydata,ggplot2:: aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, fill=z)) + 
    ggplot2::geom_rect() +
    ggplot2::scale_fill_gradient(name=var_name, low=scale_col[1], high=scale_col[2], limits=scale_lim, oob=scales::squish) +
    ggplot2::ylab('metres above msl') +
    ggplot2::xlab('kilometres from start of transect')
  plot(p)
  return(p)
}

#' Calculate rough distance in kilometers between two points
#'
#' Not exported. This is very approximate - a package is available if a more accurate distance is needed.
earth.dist <- function (long1, lat1, long2, lat2)
{
  rad <- pi/180
  a1 <- lat1 * rad
  a2 <- long1 * rad
  b1 <- lat2 * rad
  b2 <- long2 * rad
  dlon <- b2 - a2
  dlat <- b1 - a1
  a <- (sin(dlat/2))^2 + cos(a1) * cos(b1) * (sin(dlon/2))^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  R <- 6378.145
  d <- R * c
  return(d)
}

#' Produces a coloured rect plot of a vertical profile over time
#'
#' Relies on output from get_ereefs_profile().
#'
#' @param profileObj A list object as output by get_ereefs_profile(), containing dates, eta, z_grid, botz,
#'              and a data frame of values.
#' @param scale_col Colours to use for low and high values. Default c("ivory", "hotpink").
#' @param scale_lim values for low and high limits of colourscale. Defaults to full range.
#' @return p handle for the generated figure
#' @export
plot_ereefs_zvt <- function(slice, var_name='Chl_a_sum', scale_col=c("ivory", "hotpink"), scale_lim=NA) {
  numprofiles <- dim(slice$profiles)[3]
  layers <- length(slice$z_grid) - 1
  zmin <- array(slice$z_grid[1:layers], c(layers, numprofiles))
  zmax <- array(slice$z_grid[2:(layers+1)], c(layers, numprofiles))
  for (i in 1:numprofiles) {
    zmin[zmin[,i]<slice$botz,i] <- slice$botz
    zmin[zmin[,i]>slice$eta[i],i] <- slice$eta[i]
    zmax[zmax[,i]<slice$botz,i] <- slice$botz
    zmax[zmax[,i]>slice$eta[i],i] <- slice$eta[i]
  }
  d <- slice$dates
  dmin <- c(d[1]-(d[2]-d[1])/2, d[1:(length(d)-1)])
  dmin <- t(array(dmin, c(numprofiles, layers)))
  dmax <- c(d[2:length(d)], d[length(d)-1] + (d[length(d)] - d[length(d)-1])/2)
  dmax <- t(array(dmax, c(numprofiles, layers)))
  
  ind <- which(!is.na(slice$profiles[, var_name, ]))
  profiles <- slice$profiles[, var_name, ]
  if (length(scale_lim)==1) {
    scale_lim[1] <- min(c(profiles[ind]))
    scale_lim[2] <- max(c(profiles[ind]))
  }
  
  mydata <- data.frame(xmin=as.Date(dmin[ind], origin='1990-01-01'), xmax=as.Date(dmax[ind], origin='1990-01-01'), 
                       ymin=zmin[ind], ymax=zmax[ind], 
                       z=profiles[ind])
  p <- ggplot2::ggplot(data=mydata, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, fill=z)) + 
    ggplot2::geom_rect() +
    ggplot2::scale_x_date() +
    ggplot2::ylab('metres above msl') +
    ggplot2::scale_fill_gradient(name=var_name, low=scale_col[1], high=scale_col[2], limits=scale_lim, oob=scales::squish)
  plot(p)
  return(p)
}

#' Calculate rough distance in kilometers between two points
#'
#' Not exported. This is very approximate - a package is available if a more accurate distance is needed.
earth.dist <- function (long1, lat1, long2, lat2)
{
  rad <- pi/180
  a1 <- lat1 * rad
  a2 <- long1 * rad
  b1 <- lat2 * rad
  b2 <- long2 * rad
  dlon <- b2 - a2
  dlat <- b1 - a1
  a <- (sin(dlat/2))^2 + cos(a1) * cos(b1) * (sin(dlon/2))^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  R <- 6378.145
  d <- R * c
  return(d)
}

find_intersections <- function(location_latlon, x_grid, y_grid, latitude, longitude) {
  a <- (dim(x_grid) - 1)[1]
  b <- (dim(x_grid) - 1)[2]
  intersected <- array(FALSE, dim=c(a,b))
  gx <- c(x_grid[1:a, 1:b], x_grid[2:(a+1), 1:b], x_grid[2:(a+1), 2:(b+1)], x_grid[1:a, 2:(b+1)])
  gy <- c(y_grid[1:a, 1:b], y_grid[2:(a+1), 1:b], y_grid[2:(a+1), 2:(b+1)], y_grid[1:a, 2:(b+1)])
  gx <- array(gx, dim=c(a*b,4))
  gy <- array(gy, dim=c(a*b,4))
  # Find the grid cells intersected by the line between the two points.
  
  # First, define the line:
  lon1 <- location_latlon[1,'longitude']
  lon2 <- location_latlon[2,'longitude']
  lat1 <- location_latlon[1,'latitude']
  lat2 <- location_latlon[2,'latitude']
  
  # Define a line through the two points
  #Y = A*x+b ; A*x +b - Y = 0
  A <- (lat2 - lat1) / (lon2 - lon1)
  b = lat1 - (A * lon1)
  # For each edge, the line intersects the edge if the value of Ax+b-Y is +ve for one vertex and -ve for the other
  # Ignore exact hits on vertices.
  
  # Find grid-cells along this line
  c1 <- (A * gx[,1] + b - gy[,1]) > 0
  c2 <- (A * gx[,2] + b - gy[,2]) > 0
  c3 <- (A * gx[,3] + b - gy[,3]) > 0
  c4 <- (A * gx[,4] + b - gy[,4]) > 0
  
  intersected[c1!=c2] <- TRUE
  intersected[c2!=c3] <- TRUE
  intersected[c3!=c4] <- TRUE
  intersected[c4!=c1] <- TRUE
  
  # Exclude pesky NA values
  intersected[is.na(latitude)] <- FALSE
  intersected[is.na(c1+c2+c3+c4)] <- FALSE
  intersected[is.na(rowSums(gx)+rowSums(gy))] <- FALSE
  
  # Exclude grid-cells beyond the ends of the line segment
  minlat <- min(lat1, lat2)
  maxlat <- max(lat1, lat2)
  minlon <- min(lon1, lon2)
  maxlon <- max(lon1, lon2)
  intersected[((gx[,1]<minlon)&(gx[,2]<minlon)&(gx[,3]<minlon)&(gx[,4]<minlon))] <- FALSE
  intersected[((gx[,1]>maxlon)&(gx[,2]>maxlon)&(gx[,3]>maxlon)&(gx[,4]>maxlon))] <- FALSE
  intersected[((gy[,1]<minlat)&(gy[,2]<minlat)&(gy[,3]<minlat)&(gy[,4]<minlat))] <- FALSE
  intersected[((gy[,1]>maxlat)&(gy[,2]>maxlat)&(gy[,3]>maxlat)&(gy[,4]>maxlat))] <- FALSE
  
  llind <- which(intersected)
  gx <- gx[intersected, ]
  gy <- gy[intersected, ]
  
  # Find the locations where the edge intersections occur (could have done this all in one go, but this is easier to follow)
  #Ax + b = a + cx -kc
  c1 <- (gy[,2] - gy[,1]) / (gx[,2] - gx[,1])
  xi <- (gy[,1] - b -gx[,1]*c1) / (A -c1)
  d <- abs((xi<gx[,1]) + (xi<gx[,2]))
  xi[d>1] <- NA
  yi <- gy[,1] + (xi - gx[,1]) / (gx[,2] - gx[,1]) * (gy[,2] - gy[,1])
  
  c1 <- (gy[,3] - gy[,2]) / (gx[,3] - gx[,2])
  x <- (gy[,2] - b -gx[,2]*c1) / (A -c1)
  d <- abs((x<gx[,2]) + (x<gx[,3]))
  x[d>1] <- NA
  y <- gy[,2] + (gy[,3] - gy[,2]) * (x - gx[,2]) / (gx[,3] - gx[,2])
  xi <- c(xi, x)
  yi <- c(yi, y)
  
  c1 <- (gy[,4] - gy[,3]) / (gx[,4] - gx[,3])
  x <- (gy[,3] - b -gx[,3]*c1) / (A -c1)
  d <- abs((x<gx[,3]) + (x<gx[,4]))
  x[d>1] <- NA
  y <- gy[,3] + (gy[,4] - gy[,3]) * (x - gx[,3]) / (gx[,4] - gx[,3])
  xi <- c(xi, x)
  yi <- c(yi, y)
  
  c1 <- (gy[,1] - gy[,4]) / (gx[,1] - gx[,4])
  x <- (gy[,4] - b -gx[,4]*c1) / (A -c1)
  d <- abs((x<gx[,4]) + (x<gx[,1]))
  x[d>1] <- NA
  y <- gy[,4] + (gy[,1] - gy[,4]) * (x - gx[,4]) / (gx[,1] - gx[,4])
  xi <- c(xi, x)
  yi <- c(yi, y)
  yi <- yi[!is.na(xi)]
  xi <- xi[!is.na(xi)]
  
  if (maxlon>=minlon) {
    d <- sort(xi, index.return=TRUE)
    yi <- c(yi[d$ix][seq(1,length(xi),3)], yi[d$ix[length(xi)]])
    xi <- c(d$x[seq(1,length(xi),3)], d$x[length(xi)])
    if (longitude[llind[length(llind)]] < longitude[llind[1]] ) {
      yi <- yi[seq(length(xi), 1, -1)]
      xi <- xi[seq(length(xi), 1, -1)]
    }
  } else {
    d <- sort(yi, index.return=TRUE)
    xi <- c(xi[d$ix][seq(1,length(xi),3)], xi[d$ix[length(xi)]])
    yi <- c(d$x[seq(1,length(xi),3)], d$x[length(xi)])
    if (latitude[llind[length(llind)]] < latitude[llind[1]] ) {
      yi <- yi[seq(length(xi), 1, -1)]
      xi <- xi[seq(length(xi), 1, -1)]
    }
  }
  
  location_latlon <- data.frame(latitude=latitude[llind], longitude=longitude[llind])
  location_edges <- data.frame(latitude=yi, longitude=xi)
  return(list(location_latlon, location_edges, intersected, llind))
}