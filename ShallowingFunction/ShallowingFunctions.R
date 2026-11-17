
# ------------------------------------------------------------
# Script:   ShallowingFunctions.R
# Author:   Andre Menegotto
# Purpose:  Change the limits of the shelf polygons to shallower depths
# Updated:  11-11-2026
# ------------------------------------------------------------

## Packages
require(terra)
require(dplyr)
require(igraph)


## Change the depth limits for one or many shelf polygons
ShelfShallowing <- function(gebcoRast, PolyShelves, newTHR, GID=NULL)
{
  # gebcoRast = GEBCO raster object
  # PolyShelves = The polygon for one or many shelves from the polyShelf database
  # newTHR = A new threshold shallower than the -200 defined initially
  # GID = List of spatial vector returned from the function 'readGID'

  ## Set up progress bar
  cat(paste("Limitng the shelves to",abs(newTHR)*-1,"m...\n", sep=" "))
  pb <- txtProgressBar(min = 0, max = nrow(PolyShelves), style = 3)
  
  ## Run the procedure for each polygon
  resu <- list()
  for(i in 1:nrow(PolyShelves))
  {
    resu[[i]] <- ShelfShallow(gebcoRast = gebcoRast, singlePolyShelf = PolyShelves[i,], newTHR = newTHR, GID = GID)

    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  ## Remove NULLs and aggregate all the polygons into a single object
  resu <- resu[!sapply(resu, is.null)]
  resu <- do.call('rbind', resu)
  
  ## Update column's name
  if(!is.null(resu))
  {
    names(resu)[6] <- paste("area_0_",abs(newTHR),sep="")
  }
  
  ## Return output
  return(resu)
}


## Main function: Change the depth limits for a single shelf polygon
ShelfShallow <- function(gebcoRast, singlePolyShelf, newTHR, GID=NULL)
{
  # gebcoRast = GEBCO raster object
  # singlePolyShelf = The polygon for one shelf from the polyShelf database
  # newTHR = A new threshold shallower than the -200 defined initially
  # GID = List of spatial vector returned from the function 'readGID'
  
  ## Save attributes
  att <- values(singlePolyShelf)
  is_continent <- att$type=="Continental"
  
  ## Adjust the newTHR
  if(newTHR > 0){newTHR <- newTHR*-1}
  
  ## Control for deeper places
  if(att$maxElev < newTHR){return(NULL)}
  
  
  ## Extract the associated raster values and exclude grids deeper than the new threshold
  newCrop <- crop(gebcoRast, singlePolyShelf)
  newCrop <- suppressWarnings(mask(newCrop, singlePolyShelf))
  newRast <- ifel(newCrop < newTHR, NA, newCrop)
  #plot(newRast)
  
  
  ## Extract the polygons
  if(att$area > 10^5)
  {
    # For larger rasters use this adapted function that extract the polygon by latitudinal slices
    newPoly <- as.LargePolygons(myRast = newRast)
  }
  else
  {
    newPoly <- aggregate(as.polygons(newRast, dissolve=T, values=F, round=F))
  }
  newPoly <- fillHoles(newPoly)
  newPoly <- disagg(aggregate(newPoly, dissolve=T))
  #plot(newPoly)
  
  
  ## If GID is provided, include GID island polygons
  if(!is.null(GID) & !is.na(att$ALL_Uniq))
  {
    # Get islands unique ID at GID
    isl <- unlist(strsplit(unique(att$ALL_Uniq), ","))
    isl <- as.numeric(isl)
    
    # Find the respective GID polygons
    p <- which(values(GID$gIsl)$ALL_Uniq %in% isl)
    tmpIsl <- GID$gIsl[p,]
    
    rel <- relate(tmpIsl, newPoly, relation='intersects', pairs=T)
    UniqRel <- unique(rel[,1])
    
    # Merge all polygons
    mergPolys <- rbind(newPoly, tmpIsl[UniqRel,])
    mergPolys <- aggregate(mergPolys, dissolve=T)
    mergPolys <- fillHoles(mergPolys)
    newPoly <- disagg(aggregate(mergPolys, dissolve=T))
  }
  
  
  ## Check for adjacent polygons (touching)
  newPoly <- mergTouch(myPoly = newPoly)
  
  
  ## Reorder the selected polygons (area > 1km)
  pArea <- expanse(newPoly, 'km')
  pArea[pArea<1] <- NA
  newPoly <- newPoly[order(pArea, na.last = NA, decreasing = T),]
  #plot(newPoly)
  
  if(length(newPoly)==0){return(NULL)}
  
  
  ## Recalculate attributes...
  att <- do.call("rbind", replicate(length(newPoly), att, simplify=F))
  att$ID <- 1:length(newPoly)
  
  # Update coordinates
  LatLong <- crds(centroids(newPoly))
  att$long <- LatLong[,'x']
  att$lat <- LatLong[,'y']
  
  # Update depth
  att$maxElev <- zonal(newCrop, newPoly, 'max', na.rm=TRUE)[,1]
  
  newCrop <- ifel(newCrop > 0, NA, newCrop)
  for(i in 1:length(newPoly))
  {
    x <- crop(newCrop, newPoly[i,])
    x <- suppressWarnings(mask(x, newPoly[i,]))
    att[i,c("depthMax","depth_q75","depth_q50","depth_q25","depthMin")] <- global(x, fun=quantile, na.rm=TRUE)
  }
  
  # Update area
  newCrop <- ifel(newCrop < newTHR, NA, newCrop)
  for(i in 1:length(newPoly))
  {
    x <- crop(newCrop, newPoly[i,])
    x <- suppressWarnings(mask(x, newPoly[i,]))
    att$area_0_200[i] <- expanse(x, 'km')[1,2]
  }
  att$area <- expanse(newPoly, 'km')
  
  # Update classification (first round)
  att$type <- ifelse(att$maxElev > 0, "Insular", "Other_Shallow")
  
  
  ## Update insular information based on Global Islands Database
  if(!is.null(GID))
  {
    ## Update ID and number of large islands
    
    # Restart the island count
    att$n_isl <- NA_integer_
    
    # Get islands unique ID at GID
    isl <- unlist(strsplit(unique(att$ALL_Uniq), ","))
    isl <- as.numeric(isl)
    
    # If a ID is present (i.e. the shelf polygon overlaps with an GID island polygon)
    if(any(!is.na(isl)))
    {
      # Find the respective GID polygons
      p <- which(values(GID$gIsl)$ALL_Uniq %in% isl)
      tmpIsl <- GID$gIsl[p,]
      #plot(tmpIsl, add=T, col="red")
      
      # Check the overlap between the GID and the shelf polygons
      rel <- relate(newPoly, tmpIsl, relation='intersects', pairs=T)
      UniqRel <- unique(rel[,1])
      
      # Update the association and the number of islands
      att$ALL_Uniq <- NA_character_
      for(i in 1:length(UniqRel))
      {
        pos <- which(rel[,1] == UniqRel[i])
        att$ALL_Uniq[UniqRel[i]] <- paste(sort(tmpIsl$ALL_Uniq[rel[pos,2]]), collapse=',')
        att$n_isl[UniqRel[i]] <- length(pos)
        att$type[UniqRel[i]] <- "Insular" # just in case the shelf is apparently submerged but overlaps with an island
      }
    }
    
    
    ## Update indication for the presence of small islands
    
    # If classified as insular but no large islands is present
    pos <- which(is.na(att$ALL_Uniq) & att$maxElev>=0)
    
    if(length(pos)>0)
    {
      # Check overlap with small and very small islands
      relS <- relate(newPoly[pos,], GID$sIsl, relation='intersects', pairs=T)
      relVS <- relate(newPoly[pos,], GID$vsIsl, relation='intersects', pairs=T)
      rel <- rbind(relS, relVS)
      UniqRel <- unique(rel[,1])
      att$n_isl[pos[UniqRel]] <- 0
    }
    
    # Even if classified as Other_Shallow, check the presence of small islands
    pos <- which(att$type=="Other_Shallow")
    
    if(length(pos)>0)
    {
      # Check overlap with small and very small islands
      relS <- relate(newPoly[pos,], GID$sIsl, relation='intersects', pairs=T)
      relVS <- relate(newPoly[pos,], GID$vsIsl, relation='intersects', pairs=T)
      rel <- rbind(relS, relVS)
      UniqRel <- unique(rel[,1])
      
      if(length(UniqRel)>0)
      {
        att$type[pos[UniqRel]] <- "Insular"
        att$n_isl[pos[UniqRel]] <- 0
      }
    }
    
    
    ## Update quality flags
    att$flag <- NA_character_
    
    # Insular shelves with no positive topography in GEBCO
    att$flag[which(att$type == "Insular" & att$maxElev <= 0)] <- "NoGEBCOLand"
    
    # Potential insular shelves (positive topography in GEBCO) with no described islands in the Global Islands Database
    att$flag[which(att$type == "Insular" & is.na(att$n_isl))] <- "NoGIDMatch"
    
    # No shelf area (too steep shelf - from island to deep sea)
    att$flag[att$area_0_200==0] <- "NoShelfArea"
    
    # Shelf border deeper than -200m due to mismatch between GEBCO derived shelf and GID islands polygons
    # Check for insular shelves with deep regions
    pos <- which(att$type=="Insular" & att$depthMax < newTHR)
    
    if(length(pos)>0)
    {
      propDeep <- deepArea(gebcoRast = gebcoRast, subPolys = newPoly[pos,], subAtt = att[pos,], newTHR = newTHR)
      zz <- which(propDeep>0)
      
      if(length(zz)>0)
      {
        flagDeep <- paste("DeepMargin=", sprintf("%1.2f", propDeep[zz]), sep="")
        
        att$flag[pos[zz]] <- ifelse( is.na(att$flag[pos[zz]]), flagDeep, paste(att$flag[pos[zz]], flagDeep, sep = "; "))
      }
    }
    
  }
  else
  {
    att$ALL_Uniq <- "Dismissed"
    att$n_isl <- "Dismissed"
    att$flag <- "Dismissed"
  }
  
  # Add attributes
  if(is_continent){att$type[1] <- "Continental"}
  newPoly <- cbind(newPoly, att)
  
  ## Return
  return(newPoly)
}


## Read and adjust the polygons from the Global Island Database
readGID <- function(GID_dir)
{
  #GID_dir = Path to the GID folder (including the GID folder)
  
  ## Get the general path 
  x <- paste(GID_dir, '/v108/globalislandsfix.gdb', sep='')
  
  ## Read and adjust large islands
  gIsl <- vect(x, layer='USGSEsriWCMC_GlobalIslandsv2_BigIslands')
  gIsl <- project(gIsl, 'epsg:4326')
  #nrow(gIsl) #22471
  
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
  #nrow(gIsl) #22474
  
  # Exclude isolated polygons
  gIsl <- gIsl[-c(as.numeric(uniPolys)),]
  #nrow(gIsl) #22468
  
  # Reorder polygons based on area size
  tmp <- values(gIsl)
  polyOrder <- order(tmp$Area_Geode)
  gIsl <- gIsl[polyOrder,]
  
  
  ## Read and adjust small islands
  sIsl <- vect(x, layer='USGSEsriWCMC_GlobalIslandsv2_SmallIslands')
  sIsl <- project(sIsl, 'epsg:4326')
  #length(sIsl) #269,391
  
  vsIsl <- vect(x, layer='USGSEsriWCMC_GlobalIslandsv2_VerySmallIslands')
  vsIsl <- project(vsIsl, 'epsg:4326')
  #length(vsIsl) #77,534
  
  
  ## Return list with all three sizes
  return(list(gIsl=gIsl, sIsl=sIsl, vsIsl=vsIsl))
}


## Polygonize large raster by latitudinal slices
as.LargePolygons <- function(myRast)
{
  ## Get the latitudinal limits
  yMax <- ceiling(ext(myRast)[4])
  yMin <- floor(ext(myRast)[3])
  
  ## Set up progress bar
  cat('\nConverting large raster to polygon:\n')
  pb <- txtProgressBar(min = 0, max = yMax-yMin, style = 3)
  
  # Convert the raster to polygon by slices
  cont <- 0
  myLat <- yMax
  myPoly <- vect()
  while(myLat > yMin)
  {
    land <- crop(myRast, ext(-180,180,myLat-1.05,myLat))
    tmp <- as.polygons(land, dissolve=T, values=F, round=F)
    
    if(length(tmp) > 0)
    {
      tmp <- aggregate(tmp)
      myPoly <- terra::union(myPoly, tmp)
      myPoly <- aggregate(myPoly)
    }
    
    myLat <- myLat-1
    cont <- cont+1
    setTxtProgressBar(pb, cont)
  }
  close(pb)
  
  ## Return output
  return(myPoly)
}


## Aggregate all adjacent polygons into a single object
mergTouch <- function(myPoly)
{
  # Check for adjacent polygons (touching)
  adP_tmp <- adjacent(myPoly, type='queen', pairs=T, symmetrical=T)
  
  if(nrow(adP_tmp)>0)
  {
    # Find all polygons adjacent
    g <- igraph::graph_from_data_frame(adP_tmp, directed = FALSE)
    cc <- igraph::components(g)
    posList <- lapply(split(names(cc$membership), cc$membership), as.integer)
    
    # Aggregate adjacent polygons into a single polygon
    isol <- setdiff(1:length(myPoly), unique(as.numeric(unlist(posList))))
    tmp <- myPoly[isol,]
    
    for(i in 1:length(posList))
    {
      xP <- aggregate(myPoly[posList[[i]],])
      tmp <- rbind(tmp, xP)
    }
    
    myPoly <- tmp
  }
  
  # Return output
  return(myPoly)
}


## Get the proportion of external area located in deep waters (< newTHR)
deepArea <- function(gebcoRast, subPolys, subAtt, newTHR)
{
  areaTab <- data.frame(ID=subAtt$ID, areaTot=NA, areaClean=NA)
  for(i in 1:nrow(areaTab))
  {
    # Select polygon
    tmpS <- subPolys[i,]
    
    # Crop raster
    x2 <- crop(gebcoRast,tmpS)
    x2 <- suppressWarnings(mask(x2,tmpS))
    areaTab$areaTot[i] <- expanse(x2, "km")[1,2]
    
    # Get area only for the cells above new threshold
    x3 <- ifel(x2 < newTHR, NA, x2)
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
  }
  
  # Get the proportion and return the output
  propDeep <- round(1-(areaTab$areaClean/areaTab$areaTot), 2)
  return(propDeep)
}


#### END ####