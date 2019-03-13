
library(devtools)
library(ncdf4)
library(ggplot2)
library(ereefs)




map_ereefs(var_name = "temp")

#manual ncdf package

map_ereefs(var_name = "temp", layer= -10)


janP <- get_ereefs_profile(var_names = "temp", location_latlon = c(-15, 153), start_date = c(2019, 01, 01), end_date = c(2019, 01, 31) )


jan10 <- get_ereefs_depth_specified_ts(var_names = "temp", depth = 10, start_date = c(2019, 01, 01), end_date = c(2019, 01, 31))

jan20 <- get_ereefs_depth_specified_ts(location_latlon = c(-15, 153), depth = 20, start_date = c(2019, 01, 01), end_date = c(2019, 01, 31))
janP2 <- get_ereefs_profile(var_names = "temp")
# profile 1 = bottom, increasing numbers = decrease depth

plot_ereefs_profile(profileObj = janP, var_name = "temp", target_date = c(2019, 01, 21))
#how to plot profile??
#what time of day does the model predict for?


map <- map_ereefs('temp', target_date = c(2019, 01, 21), layer = -5, input_file = 3 )
map + theme_classic()
#what exactly are the intergrated depths/how to plot map of different depths?
map_ereefs(jan10)

plot(janP$profiles)

depthI <- c(1:47)
p1 <- data.frame(janP$profiles, depthIs)

janP$profiles


#z_top or z_bot
# 5M LAYER 43

map_ereefs_movie('temp', start_date = c(2019, 01, 01), end_date = c(2019, 01, 31))

#taken from downloads
glider <- nc_open("IMOS_ANFOG_BCEOPSTUV_20160324T003515Z_SL248_FV01_timeseries_END-20160419T002111Z.nc")

glider_temp<-ncvar_get(glider, 'TEMP');
glider_dep<-ncvar_get(glider, 'DEPTH');
glider_time<-ncvar_get(glider, 'TIME')
glider_lat<-ncvar_get(glider, 'LATITUDE')
glider_long<-ncvar_get(glider, 'LONGITUDE')
glider_v <- data.frame(glider_time, glider_lat, glider_long, glider_dep, glider_temp)
glider_v1 <- cbind(glider_time, glider_lat, glider_long, glider_dep, glider_temp)

names(glider_v) <- c("time", "lat", "long", "depth", "temp")

glider_v$temp <- glider_v$glider_temp

glider_2 <- subset(glider_v, time ==  )

glider_1 <- subset(glider_v,glider_time == )

install.packages("dplyr")
library(dplyr)
glider_p3 <- select(filter(glider_TDT, glider_time == 24189.03),c(glider_dep, glider_temp))

qplot(x=glider_temp, y=glider_dep, data=glider_p1, geom="line", xlim = c(29, 29.8), ylim= rev(c(0, 30)))

qplot(x=glider_lat, y=glider_long, data=glider_v, geom="line")

ggplot(glider_v, aes(glider_long, glider_lat)) +  geom_path() + theme_bw()

subset.data.frame()

glider_p1 <- glider_v[1:100,]
plot(glider_p1$glider_temp, glider_p1$glider_dep, xlim = c(29, 29.8), ylim=rev(range(1,30)))



e2016_3<- get_ereefs_profile(var_names = c("temp", "time"), start_date = c(2016,03,24), end_date = c(2016, 03, 24), location_latlon = c(-17.57, 146.38) )
