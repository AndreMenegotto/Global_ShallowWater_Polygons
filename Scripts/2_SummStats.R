
# ------------------------------------------------------------
# Script:   2_SummStats.R
# Author:   Andre Menegotto
# Purpose:  Summarise the shelf polygons data
# Updated:  12-11-2026
# ------------------------------------------------------------


# > Summary ----
library(terra)

# Load the files
polyShelf <- vect(here::here('Data','Shelf_Polygons','polyShelf.shp'))
dfShelf <- read.table(here::here('Data','Shelf_Polygons','polyShelfDF.txt'), h=T, sep = '\t')
polyShelf <- cbind(polyShelf, dfShelf[,2:ncol(dfShelf)])

# Shelf types and other submerged shallow features
table(dfShelf$type)
#Continental: 3   Insular: 1662   Other_Shallow: 4023


# Number of oceanic islands (considering only large islands in the Global Island Database)
x <- dfShelf[dfShelf$type!="Other_Shallow",]
round(sum(x$n_isl[4:nrow(x)], na.rm=T)/sum(x$n_isl, na.rm=T), 3)*100 #22.3%

# Oceanic shallow area (including Other_Shallow)
round(sum(dfShelf$area_0_200[4:nrow(dfShelf)])/sum(dfShelf$area_0_200), 3)*100 #10.8%

# Insular shelf area (excluding Other_Shallow)
round(sum(x$area_0_200[4:nrow(x)])/sum(x$area_0_200), 3)*100 #9.7%


# Number of insular shelves...
sum(!is.na(dfShelf$ALL_Uniq[dfShelf$type=="Insular"])) #...considering large islands in the Global Islands Database = 999
sum(!is.na(dfShelf$n_isl[dfShelf$type=="Insular"])) #...considering any island in the Global Islands Database = 1301
sum(is.na(dfShelf$n_isl[dfShelf$type=="Insular"]))  #...considering those not in it = 361

# Elevation in those 'islands' not present in the Global Islands Database
x <- which(dfShelf$type=="Insular" & is.na(dfShelf$n_isl))
quantile(dfShelf$maxElev[x]) #1   14   38   85  492
length(x) # n=361

# Elevation in those shelves indicated as insular by Global Islands Database, but with max topography below zero
x <- which(dfShelf$type=="Insular" & dfShelf$maxElev < 0)
rev(quantile(dfShelf$maxElev[x])) #-1   -1   -10  -53 -164
length(x) # n=109

# These would be classified as submerged if GID was not used?
all(!is.na(dfShelf$n_isl[x])) #True

# Any islands with no clear shelf
# These are probably really steep shelves and the topography raster only gets the surface and already deeper waters
# Because the GID polygons were merged to the shallow polygons...
#... and there is a discrepancy between the limits of the islands between the polygons and rasters
#Some deep cells in the raster overlap with the island polygon, allowing them enter the depth estimation
x <- which(dfShelf$type=="Insular" & dfShelf$depthMin < -200)
length(x) # n=1
dfShelf$depthMin[x] #-315

# These would be classified as submerged if GID was not used?
all(!is.na(dfShelf$n_isl[x])) #True

# How many shelves have maximum depths below -200m?
x <- which(dfShelf$depthMax < -200)
length(x) # n=914
table(dfShelf$type[x])
#Continental: 3   Insular: 634   Other_Shallow: 277
# The may occur because of the reasons above, for islands
# or because there are deep hollows, for the other features

# Shelves with no area?
x <- which(dfShelf$area_0_200==0)
table(dfShelf$type[x]) #Insular=1
# Again, islands with very steep shelves for the raster resolution, getting only the surface and deeper waters


# Other_Shallow minimum depth
x <- dfShelf[dfShelf$type=="Other_Shallow",]
rev(quantile(x$maxElev, na.rm=T)) #0  -87 -149 -183 -200

# Other_Shallow' minimum depth at the limit
x <- which(dfShelf$type=="Other_Shallow" & dfShelf$depthMin %in% -200)
length(x) #n=18


## Checking how much deeper cells in the border affect the polygons
x <- which(dfShelf$type=="Insular" & dfShelf$depthMax < -200)
length(x) #634

xF <- grep("DeepMargin", dfShelf$flag[x])
length(xF) #173
round(length(xF)/length(x)*100) #27%

flags <- unlist(strsplit(dfShelf$flag[x[xF]], "; "))
flags <- flags[grep("DeepMargin", flags)]
propDeep <- as.numeric(gsub("DeepMargin=", "", flags))

round(sum(propDeep <= .1) / length(propDeep) *100) #75%



# > Comparing with Harris et al (2014) ----

#. Shelves ----

## Load and prepare the shelves file
# Available at: https://doi.org/10.1016/j.margeo.2014.01.011
geomorhShelf <- vect(here::here('Data','1_raw','global seafloor geomorphic features map','Shelf.shp'))

# Sort by area
xArea <- expanse(geomorhShelf, "km")
geomorhShelf <- geomorhShelf[order(xArea, decreasing = T),]

# Continents are the geoms 1 to 5
plot(geomorhShelf[6,])

# Add shelf types
geomorhShelf$Type <- c(rep("Continental",5), rep("Insular", length(geomorhShelf)-5))
length(geomorhShelf)-5 #2575 insular shelves


## Matching shelves between data sets
x <- polyShelf[dfShelf$type!="Other_Shallow",]
rel <- relate(geomorhShelf, x, "intersects", pairs=T)

length(unique(rel[,1])) #1499
length(unique(rel[,1]))/length(geomorhShelf) #58% of theirs matching with ours
length(unique(rel[,2])) #1517
length(unique(rel[,2]))/length(x) #91% of ours matching with theirs

## Why many of theirs don't match ours?
gebco <- rast(here::here('Data','1_raw','GEBCO_2026','GEBCO_2026.nc'))

# Get those that don't match
pos <- setdiff(1:length(geomorhShelf), unique(rel[,1]))

# Get the maximum elevation of these shelves
altMax <- zonal(gebco, geomorhShelf[pos,], max)

# Most mismatches occur because their shelves are Other_Shallow in our data set...
#...and even these Other_Shallow may not entered our data set because they are deeper
sum(altMax$elevation <= 0)/length(altMax$elevation) #99.7%
sum(altMax$elevation <= -200)/length(altMax$elevation) #46.9%


## Matching full data set to include our Other_Shallow
relAll <- relate(geomorhShelf, polyShelf, "intersects", pairs=T)

length(unique(relAll[,1])) #2051
length(unique(relAll[,1]))/length(geomorhShelf) #79%

# Get those that don't match and their maximum elevation
pos <- setdiff(1:length(geomorhShelf), unique(relAll[,1]))
altMax <- zonal(gebco, geomorhShelf[pos,], max)

# Why some of their shelves don't match even our Other_Shallow?
# Most of them because they are deeper features, which did not enter our data set
sum(altMax$elevation <= -200)/length(altMax$elevation) #94.7%

# How about the other 5%?
# Check their area
x <- which(altMax$elevation > -200)
testArea <- numeric()
for(i in 1:length(x))
{
  # Crop the topography raster by the shelf polygon and delete deeper cells
  aa <- crop(gebco, geomorhShelf[pos[x[i]],])
  aa <- mask(aa, geomorhShelf[pos[x[i]],])
  aa <- ifel(aa < -200, NA, aa)
  
  # Isolate all existing polygons and calculate their maximum area
  bb <- disagg(aggregate(as.polygons(aa)))
  testArea[i] <- max(expanse(bb, "km"))
}

# All the shelves not matching are lower than 1km2 (our size threshold)
any(testArea>1) #FALSE



#. Seamounts ----

## Load and prepare the seamount file
geomorhSeamount <- vect(here::here('Data','1_raw','global seafloor geomorphic features map','Seamounts.shp'))
length(geomorhSeamount) #9951

# Get depth according to GEBCO
max_Depth_GEBCO <- zonal(gebco, geomorhSeamount, max)
geomorhSeamount$max_Depth_GEBCO <- max_Depth_GEBCO$elevation

# Compare depth from GEBCO and that available in the original data set
cor(geomorhSeamount$max_Depth_GEBCO, geomorhSeamount$Peak_Depth*-1) #0.86
plot(geomorhSeamount$max_Depth_GEBCO, geomorhSeamount$Peak_Depth*-1)

# How many seamounts are shallow according to each source?
sum(geomorhSeamount$Peak_Depth <= 200) #441 in the original file
sum(geomorhSeamount$max_Depth_GEBCO >= -200) #239 according to GEBCO

sum(geomorhSeamount$max_Depth_GEBCO > 0) #39 shelves seem to be an island
sum(geomorhSeamount$max_Depth_GEBCO >= -200 & geomorhSeamount$max_Depth_GEBCO < 0) #199 fitting the limits
geomorhSeamount$Type <- "Seamount"


## Matching the seamounts

# Get polygons within our criteria
x <- which(geomorhSeamount$max_Depth_GEBCO >= -200 & geomorhSeamount$max_Depth_GEBCO < 0)
sm <- geomorhSeamount[x,]

# Do they overlap?
rel <- relate(sm, polyShelf, "intersects", pairs=T)

length(unique(rel[,1])) #160 match with ours
table(dfShelf$type[unique(rel[,2])])

# How about the others?
# Check their area
x <- setdiff(1:length(sm), rel[,1])

areaS <- numeric()
for(i in 1:length(x))
{
  a <- crop(gebco, sm[x[i],])
  a <- mask(a, sm[x[i],])
  a <- ifel(a < -200, NA, a)
  b <- disagg(aggregate(as.polygons(a)))
  areaS[i] <- max(expanse(b, "km"))
}

# All have an area below 1km2
any(areaS>1) #FALSE


## Why our data set have so much more fully submerged shallow features?
# Because we don't have only seamounts, but also many other geomorphological features

# Load ridge file
ridges <- vect(here::here('Data','1_raw','global seafloor geomorphic features map','Ridges.shp'))

# Test overlap
rel <- relate(ridges, polyShelf, "intersects", pairs=T)

length(unique(rel[,1])) #409 ridge polygons matching our polygons
table(dfShelf$type[unique(rel[,2])]) #...which include insular shelves (587) and other_shallow features (435)



#### END ####