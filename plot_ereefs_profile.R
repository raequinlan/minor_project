#' Plots a single vertical profile using output from get_ereefs_profile()
#'
#' Relies on output from get_ereefs_profile().
#'
#' @param profileObj A list object as output by get_ereefs_profiles(), containing dates, eta, z_grid, botz and profiles
#' @param var_name The name of the variable to plot (must be a colname in profile$profiles). Default 'Chl_a_sum'.
#'        If profileObj contains only one variable, var_name is ignored and the content of the profile is shown.
#' @param target_date The target date (plot the profile closest in time to this).
#' @param p The handle of an existing figure, if you don't want to create a new figure
#' @return the handle of a figure containing the vertical profile plot
#' @examples 
#' \dontrun{
#' plot_ereefs_profile(get_ereefs_profile('TN'))
#' }
#' @export
plot_ereefs_profile <- function(profileObj, var_name='Chl_a_sum', target_date=c(2016,01,01), p=NA, colour='blue') {
  # Date to plot
  if (is.vector(target_date)) {
    target_date <- as.Date(paste(target_date[1], target_date[2], target_date[3], sep='-'))
  } else if (is.character(target_date)) {
    target_date <- as.Date(target_date)
  }
  day <- which.min(abs(target_date-profileObj$dates))
  if (names(profileObj)[5]=="profiles") { 
    colnum <- which(colnames(profileObj$profiles)==var_name)
    if (length(dim(profileObj$profiles))>2) {
      dind <- which.min(abs(profileObj$dates-target_date)) 
      values <- array(profileObj$profiles[, colnum, dind], length(profileObj$z_grid)-1)
      eta <- profileObj$eta[dind]
    } else {
      values <- array(profileObj$profiles[, colnum])
      eta <- profileObj$eta
    }
  } else { 
    # Only one variable
    var_name <- names(profileObj)[5]
    names(profileObj)[5] <- "profiles"
    if (is.null(dim(profileObj$profiles))) {
      # Only one time-step
      values <- profileObj$profiles
      eta <- profileObj$eta
    } else {
      dind <- which.min(abs(profileObj$dates-target_date)) 
      values <- array(profileObj$profiles[, dind], length(profileObj$z_grid)-1)
      eta <- profileObj$eta[dind]
    }
  }
  wet <- which(!is.na(values))
  values <- c(values[wet], values[max(wet)])
  z <- c(profileObj$botz, profileObj$z_grid[wet[1:length(wet)-1]+1], eta)
  mydata <- data.frame(z=z, values=values)
  if (length(p)==1) p <- ggplot2::ggplot(mydata)
  p <- p + ggplot2::geom_path(data=mydata, ggplot2::aes(x=values, y=z), colour=colour) + ggplot2::xlab(var_name) + ggplot2::ylab('metres above msl')
  print(p)
  return(p)
}