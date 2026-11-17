
# ------------------------------------------------------------
# Script:   1_ShelvesPolygons.R
# Author:   Andre Menegotto
# Purpose:  Define continental and insular shelf polygons
# Updated:  12-11-2026
# ------------------------------------------------------------

# Save packages
renv::init()

# > Creating shelves polygon ----
library(terra)
library(igraph)


#. Convert shelf raster to polygons ----

# Available at: https://www.gebco.net/data_and_products/gridded_bathymetry_data/
gebco <- rast(here::here('Data','1_raw','GEBCO_2026','GEBCO_2026.nc'))
gebco <- ifel(gebco < -200, NA, gebco)

# Set progress bar
pb <- txtProgressBar(min = 0, max = 180, style = 3)

# Convert the raster to polygon by slices
cont <- 0
myLat <- 90
polyShelf <- vect()
while(myLat > -90)
{
  land <- crop(gebco, ext(-180,180,myLat-1.05,myLat))
  tmp <- as.polygons(land, dissolve=T, values=F, round=F)
  
  if(length(tmp) > 0)
  {
    tmp <- aggregate(tmp)
    polyShelf <- terra::union(polyShelf, tmp)
    polyShelf <- aggregate(polyShelf)
  }
  
  myLat <- myLat-1
  cont <- cont+1
  setTxtProgressBar(pb, cont)
  
  if(cont%%5==0)
  {
    ## Save point
    #writeVector(x = polyShelf, filename = here::here('Data','2_tmp','polyShelf_0.shp'), overwrite=T)
  }
}
#polyShelf <- vect(here::here('Data','2_tmp','polyShelf_0.shp'))



#. Read islands polygons ----

# Available at: https://doi.org/10.1080/1755876X.2018.1529714
x <- here::here('Data','1_raw','USGSEsriWCMC_GlobalIslands_v3','v108','globalislandsfix.gdb')
gIsl <- vect(x, layer='USGSEsriWCMC_GlobalIslandsv2_BigIslands')
gIsl <- project(gIsl, 'epsg:4326')
nrow(gIsl) #22471

# Combine separated polygons from the same island (dateline)
uniPolys <- rbind(c(3251,22181),c(3252,22312), c(22182,22311))
#plot(rotate(gIsl[as.numeric(uniPolys),]))

for(i in 1:nrow(uniPolys))
{
  # For each set of polygons (unique island)
  # Extract polygons, attributes and find the biggest one
  islTMP <- gIsl[uniPolys[i,],]
  dtTMP <- values(islTMP)
  pos <- which.max(dtTMP$Area_Geode)
  
  # Keep id attributes of the biggest polygon and sum up the others values
  dtTMP_Unique <- dtTMP[pos,]
  dtTMP_Unique[,16:ncol(dtTMP_Unique)] <- colSums(dtTMP[,16:ncol(dtTMP)])
  
  # Aggregate polygons and add attributes
  islNew <- aggregate(islTMP)
  values(islNew) <- dtTMP_Unique
  
  # Add aggregated polygon to the global data set
  gIsl <- rbind(gIsl, islNew)
}

# Exclude isolated polygons
nrow(gIsl) #22474
gIsl <- gIsl[-c(as.numeric(uniPolys)),]
nrow(gIsl) #22468

# Reorder polygons based on area size
tmp <- values(gIsl)
polyOrder <- order(tmp$Area_Geode)
gIsl <- gIsl[polyOrder,]



#. Define the shelves ----

## Identify holes (they will disappear after aggregating the polygons)
holes <- disagg(fillHoles(polyShelf, inverse=T))
holes <- holes[order(expanse(holes, 'km'), decreasing = T)[1:5]]


## Relate the global islands to the shelves (to merge them)
polyShelf2 <- disagg(polyShelf)

posIsl <- relate(gIsl, polyShelf2, relation='intersects', pairs=T)

x <- sort(unique(posIsl[,1]))
length(gIsl)-length(x) #31


## Merge shelves and islands polygons
gIsl2 <- gIsl[x,]
polyShelf3 <- rbind(gIsl2, polyShelf)
polyShelf3 <- aggregate(polyShelf3, dissolve=T)

polyShelf3 <- fillHoles(polyShelf3)
polyShelf3 <- disagg(aggregate(polyShelf3, dissolve=T))


## Check for adjacent polygons (touching)
adP <- adjacent(polyShelf3, type='queen', pairs=T, symmetrical=T)

# Group all adjacent polygons
g <- igraph::graph_from_data_frame(adP, directed = FALSE)
cc <- igraph::components(g)
posList <- lapply(split(names(cc$membership), cc$membership), as.integer)


# Aggregate adjacent polygons into a single polygon
isol <- setdiff(1:length(polyShelf3), unique(as.numeric(adP)))
tmp <- polyShelf3[isol,]

for(i in 1:length(posList))
{
  xP <- aggregate(polyShelf3[posList[[i]],])
  tmp <- rbind(tmp, xP)
  print(i)
}

polyShelf4 <- tmp
pos <- order(expanse(polyShelf4, 'km'), decreasing = T)
polyShelf4 <- polyShelf4[pos,]


## Return the holes
cont <- polyShelf4[1,]
tmp <- geom(terra::union(cont, holes))
tmp[tmp[,'geom']>1,'hole'] <- 1
cont <- aggregate(vect(tmp, "polygons"))
#plot(cont, col='grey')

polyShelf5 <- rbind(cont, polyShelf4[-1,])
crs(polyShelf5) <- 'EPSG:4326'


## Aggregate polygons separated by crossing the dateline
rot5 <- rotate(polyShelf5, split=T)
adP <- adjacent(rot5, type='queen', pairs=T, symmetrical=T)

g <- igraph::graph_from_data_frame(adP, directed = FALSE)
cc <- igraph::components(g)
uniPolys <- lapply(split(names(cc$membership), cc$membership), as.integer)

polyShelf6 <- polyShelf5[-c(unlist(uniPolys)),]
for(i in 1:length(uniPolys))
{
  tmp <- aggregate(polyShelf5[uniPolys[[i]],])
  polyShelf6 <- rbind(polyShelf6, tmp)
}

pos <- order(expanse(polyShelf6, 'km'), decreasing = T)
polyShelf6 <- polyShelf6[pos,]
polyShelf6


## Relate the shelves to the global islands (to classify the shelves)
posIsl <- relate(gIsl, polyShelf6, relation='intersects', pairs=T)

nrow(posIsl)==length(unique(posIsl[,1])) #T
sum(duplicated(posIsl[,1])) #0 Perfect!

length(unique(posIsl[,2])) #995 Shelves
length(unique(posIsl[,1])) #22429 Islands identified
length(gIsl)-length(unique(posIsl[,1])) #39 islands with no shelf identified


## Organize data set
polyShelf6 <- polyShelf6[,-1]
polyShelf6$ID <- NA_integer_
polyShelf6$type <- c(rep('Continental',3), rep('Insular',length(polyShelf6)-3))
polyShelf6$area <- expanse(polyShelf6, 'km')

LatLong <- crds(centroids(polyShelf6))
polyShelf6$long <- LatLong[,'x']
polyShelf6$lat <- LatLong[,'y']

# Add ID of the islands in each shelf
GI_ALL_Uniq <- rep(NA_character_, length(polyShelf6))

gIslDF <- values(gIsl)
tmpShelves <- sort(unique(posIsl[,2]))
for(i in 1:length(tmpShelves))
{
  pos <- which(posIsl[,2]==tmpShelves[i])
  if(length(pos)>0)
  {
    GI_ALL_Uniq[tmpShelves[i]] <- paste(sort(gIslDF$ALL_Uniq[posIsl[pos,1]]), collapse=',')
  }
  print(i)
}
x <- unlist(strsplit(GI_ALL_Uniq,","))
sum(duplicated(na.omit(x))) #0

polyShelf6$ALL_Uniq <- GI_ALL_Uniq
polyShelf6$ALL_Uniq[polyShelf6$ALL_Uniq=='NA'] <- NA_character_

# Identify potential Seamounts and other fully submerged features
polyShelf6$type[is.na(polyShelf6$ALL_Uniq)] <- "Other_Shallow"

# Filter by size
pos <- which(polyShelf6$area >= 1)
polyShelf6 <- polyShelf6[pos,]



#. Find shelves within holes (seas) and select them ----

# get shelves in 'holes' (seas)
polyShelf2 <- disagg(polyShelf) # From the beginning
tmp <- relate(polyShelf2, holes, relation='within', pairs=T)
holeShelves <- polyShelf2[tmp[,1],]

# merge with islands without shelf
notF <- setdiff(1:length(gIsl), unique(posIsl[,1]))

z <- relate(gIsl[notF,], holeShelves, relation='intersects', pairs=T)
sum(duplicated(z[,1])) #0

holeShelves <- rbind(holeShelves, gIsl[notF[z[,1]],])
holeShelves <- aggregate(holeShelves, dissolve=T)
holeShelves <- fillHoles(holeShelves)
holeShelves <- disagg(aggregate(holeShelves, dissolve=T))

# Check for adjacent polygons (touching)
adP <- adjacent(holeShelves, type='queen', pairs=T, symmetrical=T)

g <- igraph::graph_from_data_frame(adP, directed = FALSE)
cc <- igraph::components(g)
posList <- lapply(split(names(cc$membership), cc$membership), as.integer)

# Aggregate adjacent polygons into a single polygon
isol <- setdiff(1:length(holeShelves), unique(as.numeric(adP)))
tmp <- holeShelves[isol,]

for(i in 1:length(posList))
{
  xP <- aggregate(holeShelves[posList[[i]],])
  tmp <- rbind(tmp, xP)
  print(i)
}

holeShelves <- tmp
holeShelves <- holeShelves[expanse(holeShelves, 'km') >= 1,]
holeShelves <- holeShelves[order(expanse(holeShelves, 'km'), decreasing = T),]

# Add metadata
holeShelves$ID <- NA_integer_
holeShelves$type <- 'Insular'
holeShelves$area <- expanse(holeShelves, 'km')

LatLong <- crds(centroids(holeShelves))
holeShelves$long <- LatLong[,'x']
holeShelves$lat <- LatLong[,'y']
holeShelves$ALL_Uniq <- NA_character_

# Relate to the missing global islands
tmp <- relate(gIsl[notF,], holeShelves, relation='intersects', pairs=T)

gIslDF <- values(gIsl[notF,])
zz <- sort(unique(tmp[,2]))
for(i in 1:length(zz))
{
  pos <- which(tmp[,2]==zz[i])
  if(length(pos)>0)
  {
    holeShelves$ALL_Uniq[zz[i]] <- paste(sort(gIslDF$ALL_Uniq[tmp[pos,1]]), collapse=',')
  }
  print(i)
}
holeShelves$type[is.na(holeShelves$ALL_Uniq)] <- "Other_Shallow"
holeShelves

# Check for adjacent polygons (touching)
# Some shelves in the holes are connected to the continent
polyShelf7 <- rbind(polyShelf6, holeShelves)
polyShelf7

adP <- adjacent(polyShelf7, type='queen', pairs=T, symmetrical=T)
adP # All adjacent polygons connected to the same continent

# Merge adjacent polygons
pos <- unique(as.numeric(adP))
polyShelf7$type[pos] #1 continent + some 'Other_Shallow' (no island here) 

zz <- aggregate(polyShelf7[pos,])
values(zz) <- values(polyShelf7)[1,]

polyShelf7 <- rbind(zz, polyShelf7[-pos,])
polyShelf7

# confirm classification is ok
sum(!is.na(polyShelf7$ALL_Uniq))-3 # Three polygons are the continents 
sum(polyShelf7$type=="Insular")

all(!is.na(polyShelf7$ALL_Uniq[polyShelf7$type=="Insular"]))
all(is.na(polyShelf7$ALL_Uniq[polyShelf7$type=="Other_Shallow"]))

x <- as.numeric(unlist(strsplit(polyShelf7$ALL_Uniq, ',')))
sum(duplicated(na.omit(x))) #0 Perfect!
#No island is associated with two or more shelves


## Save point
# polyShelfDF <- values(polyShelf7)
# write.table(x = polyShelfDF, file = here::here('Data','2_tmp','polyShelfDF_1.txt'), sep = '\t', row.names = F, col.names = T)
# writeVector(x = polyShelf7[,-c(2:ncol(polyShelfDF))], filename = here::here('Data','2_tmp','polyShelf_1.shp'), overwrite=T)

# polyShelf7 <- vect(here::here('Data','2_tmp','polyShelf_1.shp'))
# polyShelfDF <- read.table(here::here('Data','2_tmp','polyShelfDF_1.txt'), h=T, sep = '\t')
# polyShelf7 <- cbind(polyShelf7, polyShelfDF[,2:ncol(polyShelfDF)])

# all(is.na(polyShelfDF$ALL_Uniq[polyShelfDF$type %in% "Other_Shallow"])) #T
# all(!is.na(polyShelfDF$ALL_Uniq[!polyShelfDF$type %in% "Other_Shallow"])) #T



#. Re-organize attributes ----

# Until here, only geomorphologic features associated with large islands from (GID) were considered islands
# Everything else was considered Other_Shallow... Let's check this!

# Recalculate shelves area
polyShelf7$area <- NA_integer_
polyShelf7$area <- expanse(polyShelf7, 'km')
polyShelf7 <- polyShelf7[order(polyShelf7$area, decreasing = T)]

# Include ID
polyShelf7$ID <- 1:length(polyShelf7)

# Check colnames order
colnames(values(polyShelf7))

# Islands without shelf
x <- as.numeric(unlist(strsplit(polyShelf7$ALL_Uniq, ',')))
length(setdiff(values(gIsl)$ALL_Uniq, x)) #31 islands (GID) with no shelf found

setdiff(values(gIsl)$ALL_Uniq, x)
c("345738 361350  90143 340726 340709 340729 340879 340734 340711 340733 340809
340955  88122 340826 340928 340880 340712 340731 359126 340929 340771 340758
340752 340738 359127 340728 340708 345742 345733 340705 344013")

# Get number of island (from GID - big islands) by shelf
x <- strsplit(polyShelf7$ALL_Uniq, ',')
length(x)==length(polyShelf7)
polyShelf7$n_isl <- unlist(lapply(x, function(d) length(na.omit(d))))
polyShelf7$n_isl[which(is.na(x))] <- NA_integer_


## Get depth quantile and exclusive shelf area (leave the first polygon for later)
gebco <- rast(here::here('Data','1_raw','GEBCO_2026','GEBCO_2026.nc'))

n <- rep(NA,length(polyShelf7))

maxElev <- n
area_0_200 <- n
depthDF <- data.frame(depthMax=n, depth_q75=n, depth_q50=n, depth_q25=n, depthMin=n)

for(i in 2:length(polyShelf7))
{
  # Select polygon
  tmpS <- polyShelf7[i,]
  
  # Crop raster and get elevation
  x <- crop(gebco,tmpS)
  x <- suppressWarnings(mask(x,tmpS))
  maxElev[i] <- global(x, fun=max, na.rm=TRUE)[1,1]
  
  # Get depth inside the shelf
  x <- ifel(x > 0, NA, x)
  z <- global(x, fun=quantile, na.rm=TRUE)
  depthDF[i,] <- as.numeric(z[1,])
  
  # Get area of the cells between 0 and 200m
  x <- ifel(x < -200, NA, x)
  area_0_200[i] <- expanse(x, 'km')[1,2]

  # Control
  print(i)
}


## Get the same information for the main continental shelf (too heavy to do in one step)

# Split the polygon
left_ext  <- ext(-180, -25, -90, 90)
right_ext <- ext(-25, 180, -90, 90)

left_part  <- crop(polyShelf7[1,], left_ext)
right_part <- crop(polyShelf7[1,], right_ext)

# Get the data for each part
tmpS <- rbind(left_part, right_part)

elevT <- numeric()
areaT <- numeric()
valuT <- numeric()
for(i in 1:length(tmpS))
{
  x <- crop(gebco,tmpS[i])
  x <- suppressWarnings(mask(x,tmpS))
  elevT[i] <- global(x, fun=max, na.rm=TRUE)[1,1]

  x <- ifel(x > 0, NA, x)
  tmp <- values(x, na.rm=T)
  valuT <- rbind(valuT, tmp)
  
  x <- ifel(x < -200, NA, x)
  areaT[i] <- expanse(x, 'km')[1,2]
}

maxElev[1] <- max(elevT)
area_0_200[1] <- sum(areaT)
depthDF[1,] <- quantile(valuT[,1])


## Organize final data set
polyShelfDF <- data.frame(values(polyShelf7), area_0_200, depthDF, maxElev)

colnames(polyShelfDF)
polyShelfDF <- polyShelfDF[,c('ID','long','lat','type','area','area_0_200','depthMax','depth_q75','depth_q50','depth_q25','depthMin','maxElev','n_isl','ALL_Uniq')]


## Merge new dataframe to spatial polygons
polyShelf8 <- polyShelf7[,-c(2:7)]
polyShelf8 <- cbind(polyShelf8, polyShelfDF[,2:ncol(polyShelfDF)])


## Save point
# write.table(x = polyShelfDF, file = here::here('Data','2_tmp','polyShelfDF_2.txt'), sep = '\t', row.names = F, col.names = T)
# writeVector(x = polyShelf8[,-c(2:ncol(polyShelfDF))], filename = here::here('Data','2_tmp','polyShelf_2.shp'), overwrite=T)

# polyShelf8 <- vect(here::here('Data','2_tmp','polyShelf_2.shp'))
# polyShelfDF <- read.table(here::here('Data','2_tmp','polyShelfDF_2.txt'), h=T, sep = '\t')
# polyShelf8 <- cbind(polyShelf8, polyShelfDF[,2:ncol(polyShelfDF)])



#. Check apparent submerged features that are actually insular shelves ----
polyShelf9 <- polyShelf8

## Any Other_Shallow above water...
sum(polyShelf9$maxElev[polyShelf9$type=="Other_Shallow"]>0, na.rm = T) #583
# These are probably small islands (< 1km), or maybe islands missing in the Global Island Database


## Correct classification

# If it is above water, it must be an island...
pos <- which(polyShelf9$type=="Other_Shallow" & polyShelf9$maxElev>0)
polyShelf9$type[pos] <- "Insular"

# Let's open the Small Islands Dataset to check if these new insular shelves belong to them
x <- here::here('Data','1_raw','USGSEsriWCMC_GlobalIslands_v3','v108','globalislandsfix.gdb')
sIsl <- vect(x, layer='USGSEsriWCMC_GlobalIslandsv2_SmallIslands')
sIsl <- project(sIsl, 'epsg:4326')
length(sIsl) #269,391

vsIsl <- vect(x, layer='USGSEsriWCMC_GlobalIslandsv2_VerySmallIslands')
vsIsl <- project(vsIsl, 'epsg:4326')
length(vsIsl) #77,534

# Check if the shelf overlap with small islands polygons
id_SIs <- relate(polyShelf9[pos,], sIsl, relation="intersects", pairs=T)
id_VSIs <- relate(polyShelf9[pos,], vsIsl, relation="intersects", pairs=T)
newIsl <- unique(c(id_SIs[,1], id_VSIs[,1]))

# If there are evidence of small islands according to the global islands database...
# Inform that there are zero large islands (but we know there are small islands)
# So we are sure this is an insular shelf
polyShelf9$n_isl[pos[newIsl]] <- 0

# Otherwise, leave the number of large islands as NA.
# This indicates it is probably a insular shelf according to the topography raster...
# but there is no evidence of islands, even small ones, according to the global islands database

# Most very small islands occur together with small islands (in this subset of possible insular shelves)
setdiff(unique(id_SIs[,1]), unique(id_VSIs[,1]))
setdiff(unique(id_VSIs[,1]), unique(id_SIs[,1]))


## Any submerged feature associated with known small islands...
pos <- which(polyShelf9$type=="Other_Shallow")

# Check if the shelf overlap with small islands polygons
id_SIs <- relate(polyShelf9[pos,], sIsl, relation="intersects", pairs=T)
id_VSIs <- relate(polyShelf9[pos,], vsIsl, relation="intersects", pairs=T)
newIsl <- unique(c(id_SIs[,1], id_VSIs[,1]))

# Update classification
polyShelf9$type[pos[newIsl]] <- "Insular"
polyShelf9$n_isl[pos[newIsl]] <- 0

# These are not insular shelf according to the topography raster
# but are associated with small islands according to the global islands database...



#. Add flags for the apparently incoherent data ----
polyShelf9$flag <- NA_character_

# Insular shelves with no positive topography in GEBCO
polyShelf9$flag[which(polyShelf9$type == "Insular" & polyShelf9$maxElev <= 0)] <- "NoGEBCOLand"

# Potential insular shelves (positive topography in GEBCO) with no described islands in the Global Islands Database
polyShelf9$flag[which(polyShelf9$type == "Insular" & is.na(polyShelf9$n_isl))] <- "NoGIDMatch"

# No shelf area (too steep shelf - from island to deep sea)
polyShelf9$flag[polyShelf9$area_0_200==0] <- "NoShelfArea"

# Shelf border deeper than -200m due to mismatch between GEBCO derived shelf and GID islands polygons
gebco <- rast(here::here('Data','1_raw','GEBCO_2026','GEBCO_2026.nc'))

pos <- which(polyShelf9$type=="Insular" & polyShelf9$depthMax < -200)
length(pos) #634

areaTab <- data.frame(ID=polyShelf9$ID[pos], areaTot=NA, areaClean=NA)
for(i in 1:nrow(areaTab))
{
  # Select polygon
  tmpS <- polyShelf9[pos[i],]
  
  # Crop raster
  x2 <- crop(gebco,tmpS)
  x2 <- suppressWarnings(mask(x2,tmpS))
  areaTab$areaTot[i] <- expanse(x2, "km")[1,2]
  
  # Get area only for the cells above -200m
  x3 <- ifel(x2 < -200, NA, x2)
  x3 <- as.polygons(x3)
  if(length(x3)>0)
  {
    x4 <- aggregate(x3)
    x5 <- fillHoles(x4)
    areaTab$areaClean[i] <- expanse(x5, "km")
  }
  else
  {
    areaTab$areaClean[i] <- 0
  }
  
  print(i)
}

propDeep <- round(1-(areaTab$areaClean/areaTab$areaTot), 2)
zz <- which(propDeep>0)

flagDeep <- paste("DeepMargin=", sprintf("%1.2f", propDeep[zz]), sep="")
polyShelf9$flag[pos[zz]] <- ifelse( is.na(polyShelf9$flag[pos[zz]]),
                                    flagDeep,
                                    paste(polyShelf9$flag[pos[zz]], flagDeep, sep = "; "))


## Save final shapefile
# polyShelfDF <- values(polyShelf9)
# write.table(x =  polyShelfDF, file = here::here('Data','Shelf_Polygons','polyShelfDF.txt'), sep = '\t', row.names = F, col.names = T)
# writeVector(x = polyShelf9[,-c(2:ncol(polyShelfDF))], filename = here::here('Data','Shelf_Polygons','polyShelf.shp'), overwrite=T)


#### END ####