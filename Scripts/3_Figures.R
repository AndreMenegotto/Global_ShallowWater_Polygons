
# ------------------------------------------------------------
# Script:   3_Figures.R
# Author:   Andre Menegotto
# Purpose:  Summarise the shelf polygons data
# Updated:  12-11-2026
# ------------------------------------------------------------


# > World map ----
library(terra)

# Load the files
polyShelf <- vect(here::here('Data','Shelf_Polygons','polyShelf.shp'))
dfShelf <- read.table(here::here('Data','Shelf_Polygons','polyShelfDF.txt'), h=T, sep = '\t')
polyShelf <- cbind(polyShelf, dfShelf[,2:ncol(dfShelf)])

# Project the polygons
aa <- rotate(polyShelf[1:3,], longitude=-20, split=T, left=F)
aa <- crop(aa, ext(-20.00001, 339.99999, -90, 90))
contisEE <- project(aa, "+proj=eqearth +lon_0=160 +datum=WGS84 +units=m +no_defs")

pI <- which(dfShelf$type=="Insular")
aa <- rotate(polyShelf[pI,], longitude=-20, split=T, left=F)
aa <- crop(aa, ext(-20.00001, 339.99999, -90, 90))
pI_EE <- project(aa, "+proj=eqearth +lon_0=160 +datum=WGS84 +units=m +no_defs")

pO <- which(dfShelf$type=="Other_Shallow")
aa <- rotate(polyShelf[pO,], longitude=-20, split=T, left=F)
aa <- crop(aa, ext(-20.00001, 339.99999, -90, 90))
pO_EE <- project(aa, "+proj=eqearth +lon_0=160 +datum=WGS84 +units=m +no_defs")

# Plot the global map
g <- graticule(360, 180, crs='+proj=eqearth')
plot(g, labels=F, lwd=.5, col=NA)

plot(pO_EE, add=T, col="#F17720", border="#F17720")
plot(pI_EE, add=T, col="#00A7E1", border="#00A7E1")
plot(contisEE, add=T, col="#EBEBEB", border=NA)
plot(g, labels=F, lwd=.5, add=T)



# > Close ups ----

# Read GEBCO
gebco <- rast(here::here('Data','1_raw','GEBCO_2026','GEBCO_2026.nc'))

# Read the Global Islands Dataset
x <- here::here('Data','1_raw','USGSEsriWCMC_GlobalIslands_v3','v108','globalislandsfix.gdb')
gIsl <- vect(x, layer='USGSEsriWCMC_GlobalIslandsv2_BigIslands')
gIsl <- project(gIsl, 'epsg:4326')

sIsl <- vect(x, layer='USGSEsriWCMC_GlobalIslandsv2_SmallIslands')
sIsl <- project(sIsl, 'epsg:4326')

vsIsl <- vect(x, layer='USGSEsriWCMC_GlobalIslandsv2_VerySmallIslands')
vsIsl <- project(vsIsl, 'epsg:4326')

# Function to adjut the plot area
ajustar_ext <- function(e, largura = 4.5, altura = 3.5)
{
  # Original dimensons
  dx <- e$xmax - e$xmin
  dy <- e$ymax - e$ymin
  
  # Extension center
  cx <- (e$xmin + e$xmax) / 2
  cy <- (e$ymin + e$ymax) / 2
  
  # Correction of the horizontal dimension by the latitude
  dx_real <- dx * cos(cy * pi / 180)
  
  # Desired proportion
  proporcao <- largura / altura
  
  # Current proportion
  proporcao_atual <- dx_real / dy
  
  # Adjust the dimension
  if(proporcao_atual > proporcao)
  {
    # Increase the height
    dy_novo <- dx_real / proporcao
    dx_novo <- dx
  }
  else
  {
    # Increase the width
    dx_novo <- dy * proporcao / cos(cy * pi / 180)
    dy_novo <- dy
  }
  
  # New extension
  ext(
    cx - dx_novo / 2,
    cx + dx_novo / 2,
    cy - dy_novo / 2,
    cy + dy_novo / 2
  )
}


## Plot B
x <- polyShelf[dfShelf$type=="Insular" & dfShelf$n_isl > 0,]
i=74

# adjust extension and crop the shelf
ex <- ajustar_ext(ext(x[i,])+.05)
x2 <- crop(gebco, ex)
x2 <- mask(x2, x[i,])

# Plot the shelf and the islands
plot(x2, box=T, axes=F, legend=F, #las=1, pax=list(cex.axis=1.5, tcl=.3, mgp = c(3, .5, .5)),
     breaks=c(-10000, -200, 10000), col=c("darkblue","lightblue"))

zI <- relate(x[i,], gIsl, relation="intersects", pairs=T)
if(nrow(zI)>0){plot(gIsl[zI[,2],], add=T, col="#E2CA76", lwd=.5)}

zSI <- relate(x[i,], sIsl, relation="intersects", pairs=T)
if(nrow(zSI)>0){plot(sIsl[zSI[,2],], add=T, col="#E2CA76", lwd=.5)}

zVSI <- relate(x[i,], vsIsl, relation="intersects", pairs=T)
if(nrow(zVSI)>0){plot(vsIsl[zVSI[,2],], add=T, col="#E2CA76", lwd=.5)}

# Add the contour of other shelves and isobaths
plot(polyShelf, add=T)

x3 <- crop(gebco, ex, mask=F)
contour(x3, levels=seq(-10000, -1000, 1000), col="grey", add=T)

# Print the shelf centroid
round(crds(centroids(x[i,])), 2) #117.59 -5.54


## Plot C
x <- polyShelf[dfShelf$type=="Insular" & dfShelf$n_isl > 0,]
i=189

# adjust extension and crop the shelf
ex <- ajustar_ext(ext(x[i,])+.05)
x2 <- crop(gebco, ex)
x2 <- mask(x2, x[i,])

# Plot the shelf and the islands
plot(x2, box=T, axes=F, legend=F,
     breaks=c(-10000, -200, 10000), col=c("darkblue","lightblue"))

zI <- relate(x[i,], gIsl, relation="intersects", pairs=T)
if(nrow(zI)>0){plot(gIsl[zI[,2],], add=T, col="#E2CA76", lwd=.5)}

zSI <- relate(x[i,], sIsl, relation="intersects", pairs=T)
if(nrow(zSI)>0){plot(sIsl[zSI[,2],], add=T, col="#E2CA76", lwd=.5)}

zVSI <- relate(x[i,], vsIsl, relation="intersects", pairs=T)
if(nrow(zVSI)>0){plot(vsIsl[zVSI[,2],], add=T, col="#E2CA76", lwd=.5)}

# Add the contour of other shelves and isobaths
plot(polyShelf, add=T)

x3 <- crop(gebco, ex, mask=F)
contour(x3, levels=seq(-10000, -1000, 1000), col="grey", add=T)

# Print the shelf centroid
round(crds(centroids(x[i,])), 2) #-149.39 -17.68


## Plot D
x <- polyShelf[dfShelf$type=="Insular" & dfShelf$depthMax < -200,]
i=100

# adjust extension and crop the shelf
ex <- ajustar_ext(ext(x[i,])+.05)
x2 <- crop(gebco, ex)
x2 <- mask(x2, x[i,])

# Plot the shelf and the islands
plot(x2, box=T, axes=F, legend=F,
     breaks=c(-10000, -200, 10000), col=c("darkblue","lightblue"), border=NA)

zI <- relate(x[i,], gIsl, relation="intersects", pairs=T)
if(nrow(zI)>0){plot(gIsl[zI[,2],], add=T, col="#E2CA76", lwd=.5)}

zSI <- relate(x[i,], sIsl, relation="intersects", pairs=T)
if(nrow(zSI)>0){plot(sIsl[zSI[,2],], add=T, col="#E2CA76", lwd=.5)}

zVSI <- relate(x[i,], vsIsl, relation="intersects", pairs=T)
if(nrow(zVSI)>0){plot(vsIsl[zVSI[,2],], add=T, col="#E2CA76", lwd=.5)}

# Add deeper areas within the shelf that overlap with the island polygon
x2D <- ifel(x2 < -200, 1, NA)
x2D <- as.polygons(x2D)
x2D <- intersect(x2D, gIsl[zI[,2],])
plot(x2D, col="darkred", border=NA, add=T)

# Add the contour of other shelves and isobaths
plot(polyShelf, add=T)

x3 <- crop(gebco, ex, mask=F)
contour(x3, levels=seq(-10000, -1000, 1000), col="grey", add=T)

# Print the shelf centroid
round(crds(centroids(x[i,])), 2) #-57.69 -64.11



# > Max Elevation Histograms ----

# Elevation in those 'islands' not present in the Global Islands Database
x <- which(dfShelf$type=="Insular" & is.na(dfShelf$n_isl))
quantile(dfShelf$maxElev[x]) #1  14  38  85  492
h <- hist(dfShelf$maxElev[x], breaks=seq(0,500,25))

topoCol <- c("#ffffe5","#f7fcb9","#d9f0a3","#addd8e","#78c679","#41ab5d","#238443","#006837","#004529")
topoCol <- colorRampPalette(topoCol)(20)

barplot(h$counts, space=0, ylim=c(0,150), axes=F, col=topoCol)
axis(side=2, at=seq(0,150,30), labels=seq(0,150,30), las=1, cex.axis=1.25)
axis(side=1, at=seq(0,20,4), labels=seq(0,500,100), cex.axis=1.25)


# Elevation in those shelves indicated as insular by Global Islands Database, but with max topography below zero
x <- which(dfShelf$type=="Insular" & dfShelf$maxElev < 0)
rev(quantile(dfShelf$maxElev[x])) #-1   -1   -10  -53  -164
h <- hist(dfShelf$maxElev[x], breaks=seq(-200,0,10))

brewerBlue <- rev(c("#f7fbff","#deebf7","#c6dbef","#9ecae1","#6baed6","#4292c6","#2171b5","#08519c","#08306b"))
brewerBlue <- colorRampPalette(brewerBlue)(20)

barplot(h$counts, space=0, ylim=c(0,60), axes=F, col=brewerBlue)
axis(side=2, at=seq(0,60,15), labels=seq(0,60,15), las=1, cex.axis=1.25)
axis(side=1, at=seq(0,20,4), labels=seq(-200,0,40), cex.axis=1.25)


# Other_Shallow minimum depth
x <- dfShelf[dfShelf$type=="Other_Shallow",]
rev(quantile(x$maxElev, na.rm=T)) #0  -87 -149 -183 -200
h <- hist(x$maxElev, breaks=seq(-200,0,10))

barplot(h$counts, space=0, ylim=c(0,800), axes=F, col=brewerBlue)
axis(side=2, at=seq(0,800,200), labels=seq(0,800,200), las=1, cex.axis=1.25)
axis(side=1, at=seq(0,20,4), labels=seq(-200,0,40), cex.axis=1.25)



#### END ####