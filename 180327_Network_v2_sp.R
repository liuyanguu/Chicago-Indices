options(digits = 4)
library(ggplot2)
library(ggmap)
library(plyr)
library(sp)
library(readxl)


radius <-  6371 # radius of earth: 6371 meters
radians <- function(deg) {(deg * pi) / 180.0}

# geo distance
geodist <- function(lat1, long1, lat2, long2, data){
  c0 = sin(radians(data[[lat1]]))*sin(radians(data[[lat2]])) +
    cos(radians(data[[lat1]]))*cos(radians(data[[lat2]]))*cos(radians(data[[long1]]-data[[long2]]))
  c = radius * acos(pmin(c0, 1)) # avoid NaN when c0>1 due to rounding
  # pmin = parallel min
  c
}
# get count
# weight is input as a vector every time
getC_Jobs <- function(lat, long, lat_dest, long_dest, dis,  data, weight)
{
  c = sin(radians(data[[lat_dest]]))*sin(radians(lat)) +
       cos(radians(data[[lat_dest]]))*cos(radians(lat))*
       cos(radians(data[[long_dest]] - long))
  c = radius * acos(pmin(c, 1)) # avoid NaN when c0>1 due to rounding
  # pmin = parallel min
  sum((c < dis) * data[[weight]])
}

source("D:/NYU Marron/R/Get Weighted SE and CI.R")

# Chicago Travel Record ---------------------------------------------------
mydata <- read.csv(file = "C:/Chicago/Chicago_OD_withlatlong.csv")
mydata$DisTravel <- geodist(data=mydata, lat1 = "h_cen_lat", long1 = "h_cen_long", 
                             lat2 = "w_cen_lat", long2 = "w_cen_long")
unlist(weighted_Output(mydata$DisTravel, mydata$S000))
hist(rep(mydata$DisTravel,mydata$S000), main='Distribution of Travel Distance')

# mydata[order(mydata$S000, decreasing = T),][c(1:13),]
# table(as.factor(mydata$S000))


# 3.23 sample from 3.5 M observations ----------------------------------------------------
# sampling
set.seed(123)
test3 <- mydata[sample(mydata$id, nrow(mydata)/100),]
test3 <- mydata
test3$DisTravel <- geodist(data=test3, lat1 = "h_cen_lat", long1 = "h_cen_long", 
                           lat2 = "w_cen_lat", long2 = "w_cen_long")
# Create the grid 
long_min <- min(mydata$h_cen_long)
long_max <- max(mydata$h_cen_long)
lat_min <- min(mydata$h_cen_lat)
lat_max <- max(mydata$h_cen_lat)
size_of_block <- 2000
library('geosphere')
long_n <- round(
                distGeo(c(long_min, 0.5*(lat_max+lat_min)), 
                        c(long_max, 0.5*(lat_max+lat_min)))/size_of_block)
lat_n <- round(
                distGeo(c(0.5*(long_min+long_max), lat_min),
                        c(0.5*(long_min+long_max), lat_max))/size_of_block)
long_grid <- seq(from = long_min, to = long_max, length.out = long_n)
lat_grid <- seq(from = lat_min, to = lat_max, length.out = lat_n)


# spatial join  ----------------------------------------------------------------
library('sp')
grid <- expand.grid(long = long_grid, lat = lat_grid) # 4550,3
grid$id <- 1:nrow(grid)
grid_df <- SpatialPixelsDataFrame(points=grid[c("long", "lat")], data=grid)
grid_poly <-  as(grid_df, "SpatialPolygonsDataFrame") 

# A second option is:
# grid2 <- grid
# coordinates(grid2) <- ~long + lat
# gridded(grid2) = TRUE
# grid_poly <- as(grid2, "SpatialPolygons")

# ggmap(chicago_map) + 
#   geom_point(data = grid2, aes(x = long, y = lat), color = "red", size = 2, alpha = 0.8) +
#   geom_point(data = grid, aes(x = long, y = lat), color = "blue", size = 0.5, alpha = 0.4) +
#   theme(legend.position="none",
#         axis.title.y = element_blank(), axis.title.x = element_blank()) 
# 

# 1. Dissimilarity index --------------------------------------------------
# spatial join by imposing test3 over grid by home coordinates
test_home <- test3

sys_time <- Sys.time()
coordinates(test_home) <- ~ h_cen_long + h_cen_lat
impose1 <- over(grid_poly, test_home[c("SE01", "SE02", "SE03")], fn = sum)
Sys.time() - sys_time

impose12 <- over(grid_poly, test_home[c("DisTravel")], fn = mean)

# combine test with grid
impose1[is.na(impose1)] <- 0
impose1$id <- impose12$id <-  1:nrow(impose1)
impose1$count <- impose1$SE01 + impose1$SE02 + impose1$SE03

homeBlock2 <- join_all(list(grid, impose1, impose12), by='id')
homeBlock2 <- homeBlock2[homeBlock2$count != 0, ]

homeBlock2$SE01_02 <- homeBlock2$SE01 + homeBlock2$SE02
homeBlock2$SE02_03 <- homeBlock2$SE02 + homeBlock2$SE03

p1 <- sum(homeBlock2$SE01)/sum(homeBlock2$count)
p2 <- sum(homeBlock2$SE01_02)/sum(homeBlock2$count)

# pi: percentage % of people working on lower-income jobs 
homeBlock2$pi_SE01 <- homeBlock2$SE01/homeBlock2$count
homeBlock2$pi_SE02 <- homeBlock2$SE01_02/homeBlock2$count

# the index for plotting is the abs difference from city average 
homeBlock2$Index_SE01 <- abs(homeBlock2$pi_SE01 - p1)
DissimilarityIndex_SE01 <- sum(homeBlock2$count*homeBlock2$Index_SE01)/
  (2*sum(homeBlock2$count)*p1*(1-p1))
DissimilarityIndex_SE01

homeBlock2$Index_SE02 <- abs(homeBlock2$pi_SE02 - p2)
DissimilarityIndex_SE02 <- sum(homeBlock2$count*homeBlock2$Index_SE02)/
  (2*sum(homeBlock2$count)*p2*(1-p2))
DissimilarityIndex_SE02

# 2. Job Segregation ------------------------------------------------------
test_work <- test3
coordinates(test_work) <- ~ w_cen_long + w_cen_lat
impose2 <- over(grid_poly, test_work[c("SE01", "SE02", "SE03")], fn = sum)
impose22 <- over(grid_poly, test_work[c("DisTravel")], fn = mean)
# combine test with grid
impose2[is.na(impose2)] <- 0
impose2$id <- impose22$id <- 1:nrow(impose2)
impose2$count <- impose2$SE01 + impose2$SE02 + impose2$SE03

workBlock2 <- join_all(list(grid, impose2, impose22), by='id')
workBlock2 <- workBlock2[workBlock2$count!=0, ]

workBlock2$SE01_02 <- workBlock2$SE01 + workBlock2$SE02
workBlock2$SE02_03 <- workBlock2$SE02 + workBlock2$SE03

p1 <- sum(workBlock2$SE01)/sum(workBlock2$count)
p2 <- sum(workBlock2$SE01_02)/sum(workBlock2$count)
# pi: percentage % of people working on lower-income jobs 
workBlock2$pi_SE01 <- workBlock2$SE01/workBlock2$count
workBlock2$pi_SE02 <- workBlock2$SE01_02/workBlock2$count
# the index for plotting is the abs difference from city average 
workBlock2$Index_SE01 <- abs(workBlock2$pi_SE01 - p1)
JobSegregationIndex_SE01 <- sum(workBlock2$count*workBlock2$Index_SE01)/
  (2*sum(workBlock2$count)*p1*(1-p1))
JobSegregationIndex_SE01

workBlock2$Index_SE02 <- abs(workBlock2$pi_SE02 - p2)
JobSegregationIndex_SE02 <- sum(workBlock2$count*workBlock2$Index_SE02)/
  (2*sum(workBlock2$count)*p2*(1-p2))
JobSegregationIndex_SE02

# 3. Job / Residence in each block -----------------------------------------

w1 <- workBlock2[c('id', 'count')]; names(w1) <- c('id', 'workCount')
h1 <- homeBlock2[c('id', 'count')]; names(h1) <- c('id', 'homeCount')
library(plyr)
work_home <- join_all(list(grid, w1, h1) , by = "id")
# fill NA with 0
work_home[is.na(work_home)] <- 0
work_home$totalCount <- work_home$workCount + work_home$homeCount
work_home <- work_home[work_home$totalCount!=0,]
p_city <- sum(work_home$workCount)/sum(work_home$totalCount)
# pi:  % of job 
work_home$pi <- work_home$workCount/work_home$totalCount
# the index for plotting is the abs difference from city average 
work_home$Abs_pi_p <- abs(work_home$pi - p_city)
JobHousingIndex <- sum(work_home$totalCount*work_home$Abs_pi_p)/
  (2*sum(work_home$totalCount)*p_city*(1-p_city))

summary(work_home$pi); hist(work_home$pi)
summary(work_home$Abs_pi_p); hist(work_home$Abs_pi_p)

# 4. Life - Work: people who work in the same block ---------------------------------
# people work within certain radius (2km): 

sum(test3$S000 * (test3$DisTravel<2))/sum(test3$S000)

# people who work within 2 km radius near his block.
test_home$workNear <- test_home$S000 * (test_home$DisTravel<2) 
impose4 <- over(grid_poly, test_home[c("workNear")], fn = sum)
# combine test with grid
impose4[is.na(impose4)] <- 0
impose4$id <- 1:nrow(impose4)
#
nearBlock2 <- join_all(list(grid, impose1, impose4), by='id')
nearBlock2$count <- nearBlock2$SE01+nearBlock2$SE02+nearBlock2$SE03
nearBlock2 <- nearBlock2[nearBlock2$count!=0, ]
nearBlock2$pi <- nearBlock2$workNear/nearBlock2$count
# summary(nearBlock2$pi); hist(nearBlock2$pi)


# 5.  Foraging Index ------------------------------------------------------
# sampling
set.seed(123)
test4 <- mydata[sample(mydata$id, nrow(mydata)/100),]
test4$DisTravel <- geodist(data=test4, lat1 = "h_cen_lat", long1 = "h_cen_long", 
                           lat2 = "w_cen_lat", long2 = "w_cen_long")
  
getC_Jobs_2 <- function(...){getC_Jobs(lat_dest = "w_cen_lat", long_dest="w_cen_long", 
                                       data=test4, weight = "S000", ...)}


# test4 <- mydata
test4$JobPassed <- mapply(getC_Jobs_2, 
                          lat = test4$h_cen_lat, long = test4$h_cen_long,
                          dis = test4$DisTravel)

# 
unlist(weighted_Output(test4$DisTravel, test4$S000)); summary(test4$DisTravel);hist(test4$DisTravel)
unlist(weighted_Output(test4$JobPassed, test4$S000)); summary(test4$JobPassed);hist(test4$JobPassed)

# foraging index done by: (weighted mean)/total jobs
sum(test4$JobPassed*test4$S000)/sum(test4$S000)/sum(test4$S000)
# use median?
median(rep(test4$JobPassed,test4$S000))/sum(test4$S000)

# create block
test4$NumJobPassed <- test4$JobPassed*test4$S000
coordinates(test4) <- ~h_cen_long + h_cen_lat
impose5 <- over(grid_poly, test4[c('NumJobPassed')], fn = sum)
# combine test with grid
impose5[is.na(impose5)] <- 0
impose5$id <- 1:nrow(impose5)

homeBlock_passed <- join_all(list(grid, impose1, impose5), by='id')
homeBlock_passed <- homeBlock_passed[homeBlock_passed$count!=0, ]
homeBlock_passed$aveJobPassed <-homeBlock_passed$NumJobPassed/homeBlock_passed$count
  

# ## > Figures  --------------------------------------------------------------------

library(ggplot2)
library(ggmap)
# the map of Chicago: 
chicago_map <- get_map(location = 'chicago_map', zoom = 9)

# Map 0. Show the 3 social economical (SE) groups ####
test3_long <- reshape(test3, direction = 'long',
                     varying = list(
                       c('SE01','SE02','SE03')
                     ),
                     timevar = 'SE_level',
                     times=c("SE01", "SE02", "SE03"),
                     v.names = c('SE')
)
test3_long <- test3_long[test3_long$SE!=0,]

# size, color and alpha by group
ggmap(chicago_map) + 
  geom_point(data = test3_long, aes(x = h_cen_long, y = h_cen_lat, 
                                       color = SE_level, size = SE_level, alpha = SE_level)) +
  scale_colour_manual(values = c("red", "yellow", "lightblue")) +
  scale_size_manual(values = c(3,2,1)/2) +
  scale_alpha_discrete(range = c(1,0.5,0.1)) + 
  theme(axis.title.y = element_blank(), axis.title.x = element_blank()) 

# equal size
ggmap(chicago_map) + 
  geom_point(data = test3_long, aes(x = h_cen_long, y = h_cen_lat, 
                                    color = SE_level,alpha = SE_level)) +
  scale_colour_manual(values = c("red", "yellow", "lightblue")) +
  scale_alpha_discrete(range = c(1,0.5,0.1))+
  theme(
    axis.title.y = element_blank(), axis.title.x = element_blank()) 

# Map 1. Dissimilarity Index #### 
# cut for even interval 
# map pi 
homeBlock2$RangeSE01 <- cut_interval(homeBlock2$pi_SE01, n=6)
homeBlock2$RangeSE02 <- cut_interval(homeBlock2$pi_SE02, n=6)
ggmap(chicago_map) + 
  geom_tile(data = homeBlock2, aes(x = long, y = lat, fill = RangeSE02), alpha = 0.8)+
  scale_fill_brewer(palette = "OrRd")+
  theme(legend.position = 'bottom',
        axis.title.y = element_blank(), axis.title.x = element_blank())
# map |pi - P|
homeBlock2$AbsSE01 <- cut_interval(homeBlock2$Index_SE01, n=6)
homeBlock2$AbsSE02 <- cut_interval(homeBlock2$Index_SE02, n=6)
ggmap(chicago_map) + 
  geom_tile(data = homeBlock2, aes(x = long, y = lat, fill = AbsSE02), alpha = 0.8)+
  scale_fill_brewer(palette = "OrRd")+
  theme(legend.position = 'bottom',
        axis.title.y = element_blank(), axis.title.x = element_blank())


# Map 2. Job Dissimilarity (Segregation) ####
workBlock2$RangeSE01 <- cut_interval(workBlock2$pi_SE01, n=6)
workBlock2$RangeSE02 <- cut_interval(workBlock2$pi_SE02, n=6)
ggmap(chicago_map) + 
  geom_tile(data = workBlock2, aes(x = long, y = lat, fill = RangeSE02), alpha = 0.8)+
  scale_fill_brewer(palette = "OrRd")+
  theme(legend.position = 'bottom',
        axis.title.y = element_blank(), axis.title.x = element_blank())

workBlock2$AbsSE01 <- cut_interval(workBlock2$Index_SE01, n=6)
workBlock2$AbsSE02 <- cut_interval(workBlock2$Index_SE02, n=6)
ggmap(chicago_map) + 
  geom_tile(data = workBlock2, aes(x = long, y = lat, fill = AbsSE02), alpha = 0.8)+
  scale_fill_brewer(palette = "OrRd")+
  theme(legend.position = 'bottom',
        axis.title.y = element_blank(), axis.title.x = element_blank())

# Map 3. Work-Housing Balance ####
work_home$Share_Jobs <- cut_interval(work_home$pi, n=6)
work_home$Abs_Jobs <- cut_interval(work_home$Abs_pi_p, n=6)
ggmap(chicago_map) + 
  geom_tile(data = work_home, aes(x = long, y = lat, fill = Abs_Jobs), alpha = 0.8)+
  scale_fill_brewer(palette = "OrRd")+
  theme(legend.position = 'bottom',
        axis.title.y = element_blank(), axis.title.x = element_blank())


# Map 4. % of people work close to home  ####
nearBlock2$Range_pi <- cut_number(nearBlock2$pi, n=6)
ggmap(chicago_map) + 
  geom_tile(data = nearBlock2, aes(x = long, y = lat, fill = Range_pi), alpha = 0.8)+
  scale_fill_brewer(palette = "OrRd")+
  theme(legend.position = 'bottom',
        axis.title.y = element_blank(), axis.title.x = element_blank())

# Use point?
ggmap(chicago_map) + 
  geom_point(data = nearBlock2, aes(long, lat, 
                                   color = pi),size = 1) +
  scale_colour_gradient(low='yellow', high = 'red') +
  theme(legend.position="none",
        axis.title.y = element_blank(), axis.title.x = element_blank()) 

# Map 5. Ave # of Job passed ####
homeBlock_passed$Log_Passed <- cut_interval(log(homeBlock_passed$aveJobPassed), n=6)
homeBlock_passed$Ave_Passed <- cut_number((homeBlock_passed$aveJobPassed), n=6)
ggmap(chicago_map) + 
  geom_tile(data = homeBlock_passed, aes(x = long, y = lat, fill = Ave_Passed), alpha = 0.8)+
  scale_fill_brewer(palette = "OrRd")+
  theme(legend.position = 'bottom',
        axis.title.y = element_blank(), axis.title.x = element_blank())



# Map more. --- Live and Work Density ---------------------------------------------------
# work_home$Range_log <- cut_interval(log(work_home$workCount), n=6)
levels(work_home$Range_log)
breaks_of_interval <- seq(0, max(log(work_home$workCount)),
                          length.out=7)
work_home$Range_log_work <- cut(log(work_home$workCount), breaks = breaks_of_interval)
work_home$Range_log_home <- cut(log(work_home$homeCount), breaks = breaks_of_interval)

## find the manual color ##
library("RColorBrewer")
display.brewer.pal(n = 8, name = 'OrRd')
color_scale <- brewer.pal(n = 8, name = "OrRd")[3:8]

# where the people live
ggmap(chicago_map) + 
  geom_tile(data = work_home, aes(x = long, y = lat, fill = Range_log_home), alpha = 0.8)+
  scale_fill_manual(values = color_scale)+
  theme(legend.position = 'bottom',
        axis.title.y = element_blank(), axis.title.x = element_blank())

# where the people work
ggmap(chicago_map) + 
  geom_tile(data = work_home, aes(x = long, y = lat, fill = Range_log_work), alpha = 0.8)+
  scale_fill_manual(values = color_scale)+
  theme(legend.position = 'bottom',
        axis.title.y = element_blank(), axis.title.x = element_blank())



# (if use test3)
# where people live:
ggmap(chicago_map) + 
  stat_density2d(aes(x = h_cen_long, y =h_cen_lat, fill=..level..), 
                 data=test3 ,geom="polygon", alpha=0.8) + 
  scale_fill_gradient(low = "yellow", high = "red")

# work:
ggmap(chicago_map) + 
  stat_density2d(aes(x = w_cen_long, y =w_cen_lat, fill=..level..), 
                 data=test3 ,geom="polygon", alpha=0.8) + 
  scale_fill_gradient(low = "yellow", high = "red")


# Average jobs passed - use point ----------------------------------------

ggmap(chicago_map) + 
  geom_point(data = homeBlock_passed, aes(x = long, y = lat, 
                                  color = aveJobPassed),size = 1) +
  scale_colour_gradient(low='white', high = 'red') +
  theme(legend.position="none",
        axis.title.y = element_blank(), axis.title.x = element_blank()) 



# 6. Clustering index --------------------------------------------------------
# 4/2 program is by myself
# get sum count of xj*cij for each xi 
getC_Cluster <- function(lat, long, lat_dest, long_dest, data, weight){
   c = sin(radians(data[[lat_dest]]))*sin(radians(lat)) +
    cos(radians(data[[lat_dest]]))*cos(radians(lat))*
    cos(radians(data[[long_dest]] - long))
   d = radius * acos(pmin(c, 1)) # avoid NaN when c0>1 due to rounding
   sum(exp(-d) * data[[weight]])
}

group1 <- "SE01_02"; group2 <- "SE03" # -> gives 1.0653
# group1 <- "SE01"; group2 <- "SE02_03" # -> gives 1.0144
homeBlock2$x_j <- mapply(function(...){getC_Cluster(lat_dest = "lat", long_dest="long", 
                                               data=homeBlock2, weight = group1, ...)},
                    lat = homeBlock2$lat, long = homeBlock2$long)
homeBlock2$y_j <- mapply(function(...){getC_Cluster(lat_dest = "lat", long_dest="long", 
                                  data=homeBlock2, weight = group2, ...)},
                    lat = homeBlock2$lat, long = homeBlock2$long)
homeBlock2$t_j <- mapply(function(...){getC_Cluster(lat_dest = "lat", long_dest="long", 
                                               data=homeBlock2, weight = "count", ...)},
                    lat = homeBlock2$lat, long = homeBlock2$long)
iipenalty <- 1- exp(-sqrt(0.6*(size_of_block/1000)^2))
Pxx <- sum((homeBlock2$x_j - homeBlock2[[group1]]*iipenalty)*homeBlock2[[group1]])
Pyy <- sum((homeBlock2$y_j - homeBlock2[[group2]]*iipenalty)*homeBlock2[[group2]])
Ptt <- sum((homeBlock2$t_j - homeBlock2$count*iipenalty)*homeBlock2$count)
X <- sum(homeBlock2[[group1]])
Y <- sum(homeBlock2[[group2]])
T <- sum(homeBlock2$count)
(IndexClustering<- (Pxx/X + Pyy/Y)/(Ptt/T))


## workBlock2
group1 <- "SE01_02"; group2 <- "SE03" # -> 1.136
# group1 <- "SE01"; group2 <- "SE02_03" # -> 1.051
workBlock2$x_j <- mapply(function(...){getC_Cluster(lat_dest = "lat", long_dest="long", 
                                                data=workBlock2, weight = group1, ...)},
                     lat = workBlock2$lat, long = workBlock2$long)
workBlock2$y_j <- mapply(function(...){getC_Cluster(lat_dest = "lat", long_dest="long", 
                                                data=workBlock2, weight = group2, ...)},
                     lat = workBlock2$lat, long = workBlock2$long)
workBlock2$t_j <- mapply(function(...){getC_Cluster(lat_dest = "lat", long_dest="long", 
                                                data=workBlock2, weight = "count", ...)},
                     lat = workBlock2$lat, long = workBlock2$long)
iipenalty <- 1- exp(-sqrt(0.6*(size_of_block/1000)^2))
Pxx <- sum((workBlock2$x_j - workBlock2[[group1]]*iipenalty)*workBlock2[[group1]])
Pyy <- sum((workBlock2$y_j - workBlock2[[group2]]*iipenalty)*workBlock2[[group2]])
Ptt <- sum((workBlock2$t_j - workBlock2$count*iipenalty)*workBlock2$count)
X <- sum(workBlock2[[group1]])
Y <- sum(workBlock2[[group2]])
T <- sum(workBlock2$count)
(IndexClustering<- (Pxx/X + Pyy/Y)/(Ptt/T))

# Output ------------------------------------------------------------------
# write.csv(test3, 
#           file = "D:/NYU Marron/Data3_R/test3_with home workd id.csv")
# 
# write.csv(homeBlock2, file = "D:/NYU Marron/Data3_R/Chicago_homeBlock.csv")
# write.csv(workBlock2, file = "D:/NYU Marron/Data3_R/Chicago_workBlock.csv")
# write.csv(work_home, file = "D:/NYU Marron/Data3_R/Chicago_home and work count.csv")
# write.csv(nearBlock2, file = "D:/NYU Marron/Data3_R/Chicago_Near Block.csv")

# test3 <- read.csv("C:/Chicago/Chicago Number of Jobs Passed.csv")
# test3 <- read.csv("D:/NYU Marron/Data3_R/test3_with home work id_0328.csv")
# work_home <- read.csv("D:/NYU Marron/Data3_R/Chicago_home and work count.csv")
homeBlock2 <- read.csv("D:/NYU Marron/Data3_R/Chicago_homeBlock.csv")
workBlock2 <- read.csv("D:/NYU Marron/Data3_R/Chicago_workBlock.csv")


# 7. Extra: Distance Travelled  -----------------------------------------------------
# make 3D maps for blog
# load homeBlock2 and workBlock2
# 3D
library("lattice")
# A single figure example for testing: 
wireframe(DisTravel ~ long * lat, data = homeBlock2,
          main = "Average Distance Travelled by Home Block",
          zlab = "Distance (km)",
          drape = TRUE, # add color 
          colorkey = F, # add color key
          pretty=T,     # not so useful in this case
          # adjust how to view the screen 
          screen = list(z = -80, x = -60), alpha = 0.6
)

# Draw two figures in one frame:
par.set <- list(axis.line = list(col = "transparent"),
           clip = list(panel = "off"))
# Ave Distance travelled by home block
print(wireframe(DisTravel ~ long * lat, data = homeBlock2, 
                main = "Distance Travelled by Where People Live",
                zlab = "Distance (km)",
                drape = TRUE, colorkey = TRUE,   
                screen = list(z = 0, x = -60),
                par.settings = par.set
                ),
      split = c(1,1,2,1), more = TRUE)
print(wireframe(DisTravel ~ long * lat, data = homeBlock2, 
                zlab = "Distance (km)",
                drape = TRUE, colorkey = TRUE, 
                screen = list(z = -80, x = -60),
                par.settings = par.set
                ),
      split = c(2,1,2,1))

# 
# # Ave Distance travelled by job block
# print(wireframe(DisTravel ~ long * lat, data = workBlock2, 
#                 main = "Distance Travelled by Where People Work",
#                 zlab = "Distance (km)",
#                 drape = TRUE, colorkey = TRUE,   
#                 screen = list(z = 0, x = -60),
#                 par.settings = par.set
# ),
# split = c(1,1,2,1), more = TRUE)
# print(wireframe(DisTravel ~ long * lat, data = workBlock2, 
#                 zlab = "Distance (km)",
#                 drape = TRUE, colorkey = TRUE, 
#                 screen = list(z = -80, x = -60),
#                 par.settings = par.set
# ),
# split = c(2,1,2,1))




