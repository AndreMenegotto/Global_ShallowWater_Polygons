# Global Dataset of Insular Shelves and Other Submerged Shallow-Water Features
André Menegotto, Márcia Neunschwander Kurtz, Thomas J. Webb, Ana
Margarida Coelho Santos

## Overview

This repository provides a global spatial dataset of insular shelves and
other submerged shallow-water features, together with an R function for
modifying shelf boundaries and the R code used to generate the dataset.

The dataset consists of a shapefile containing 5,688 objects spatial
features, including three continental shelf polygons, 1662 insular shelf
polygons, and 4023 polygons representing other submerged shallow-water
features. Continental and insular shelf boundaries are defined by the
-200 m isobath, following the 1958 Convention on the Continental Shelf.

For each polygon, the dataset provides the following attributes:

- A unique identifier (ID);
- Centroid coordinates (longitude and latitude);
- Feature classification, distinguishing continental shelves, insular
  shelves, and other submerged shallow-water features;
- Total polygon area (km²);
- Polygon area within the 0 to -200 m depth range (km²);
- Depth quantiles (maximum, 75th, 50th, 25th, and minimum; m);
- Maximum elevation (m);
- Number of large islands (area \> 1 km²);
- A list of the identifiers of large islands, as recorded in the Global
  Islands Database.

A quality flag is also provided for insular shelf polygons to identify
four potential issues:

- Insular shelves with no corresponding island in the Global Islands
  Database;
- Insular shelves lacking associated positive topography;
- Insular shelves lacking areas within 0 to -200 m;
- The proportion of each polygon’s area affected by spatial mismatches
  between the topography-derived shelf boundary and the corresponding
  island boundary. These mismatches may result in portions of shelf
  polygons extending into areas deeper than -200 m.

Because the -200 m shelf boundary may extend beyond the depth range
relevant to many shallow-water organisms and research applications, we
also provide an R function that allows users to redefine the boundaries
of the original polygons according to shallower depth limits of their
choice.

## Repository structure

The repository is organized into three main directories:

**Data/**

Contains the spatial dataset and its associated documentation.

- Shelf_Polygons/: Contains the shapefile `polyShelf.shp`, the attribute
  table `polyShelfDF.txt`, and `data_dictionary.csv`, which describes
  the variables, including their definitions, units, and data formats.

**ShallowingFunction/**

Contains `ShallowingFunctions.R`, an R script defining the function for
modifying the depth limits of the original shelf polygons.

A PDF user guide is also provided in this directory, describing the
function’s inputs, outputs, and usage.

**Scripts/**

Contains the R scripts used to generate the dataset, calculate summary
statistics, and produce figures. These scripts are provided for
transparency and reproducibility but are not required to use the final
dataset.

**Other files and directories**

- `.Rproj`: RStudio project file for opening the repository as an
  RStudio project.
- `renv.lock`: Records the R package dependencies and their versions to
  support reproducibility.
- `LICENSE`: Specifies the terms under which the R scripts and
  shallowing function may be used, modified, and redistributed
  (GPL-3.0).
- `LICENSE-DATA`: Specifies the terms under which the spatial dataset
  and associated data may be used and redistributed (CC BY 4.0).

## Usage

**Accessing the dataset**

The spatial dataset can be accessed in R using either the `terra` or
`sf` package, as well as in any GIS software that supports the shapefile
format.

To load the spatial data in R:

``` r
library(terra)

polyShelf <- vect(here::here('Data','Shelf_Polygons','polyShelf.shp'))
```

The attribute table containing additional information for each polygon
is provided separately in `polyShelfDF.txt`. It can be loaded and
combined with the spatial dataset using cbind(), as the records in both
files are in the same order. Both files also contain a unique polygon
identifier (ID), which can be used to match records when combining the
datasets using other software or packages.

For example:

``` r
dfShelf <- read.table(here::here('Data','Shelf_Polygons','polyShelfDF.txt'), h=T, sep = '\t')
polyShelf <- cbind(polyShelf, dfShelf[,2:ncol(dfShelf)])

# Shelf types and other submerged shallow features
table(polyShelf$type)
```

      Continental       Insular Other_Shallow 
                3          1662          4023 

The resulting spatial object can be used for spatial analyses,
visualization, and data filtering. The available attributes allow users
to select polygons according to their characteristics and research
requirements.

``` r
# Selecting islands
insShelves <- polyShelf[polyShelf$type=="Insular",]
length(insShelves)
```

    [1] 1662

``` r
# Selecting only large islands with unique identifiers in the Global Island Database
insShelves <- polyShelf[polyShelf$type=="Insular" & !is.na(polyShelf$ALL_Uniq),]
length(insShelves)
```

    [1] 999

**Modifying shelf boundaries**

The `ShallowingFunctions.R` script provides a function to redefine the
boundaries of the original shelf polygons using shallower depth limits.

To make the function available in R:

``` r
source("ShallowingFunction/ShallowingFunctions.R")
```

For detailed instructions on the function’s arguments, required inputs,
and outputs, please refer to the PDF user guide provided in the
ShallowingFunction/ directory.

## Requirements

The spatial dataset can be accessed using standard GIS software and does
not require R.

For users working in R, the following packages may be used:

- `'terra'`: Reading, handling, and processing spatial vector data.
- `'tidyterra'`: Provides tidyverse-style methods for manipulating and
  visualizing spatial data objects from `'terra'`.
- `'sf'`: An alternative package for reading, handling, and processing
  spatial vector data.
- `'dplyr'`: An optional package for manipulating and combining
  attribute tables.

Additional packages may be required to run the scripts used to generate
the dataset or to apply the shallowing function (e.g. `'igraph'`).

The package dependencies and their versions are recorded in `renv.lock`.
To restore the computational environment used in the project, install
`'renv'` and run:

``` r
install.packages("renv")
renv::restore()
```

This will install the package versions recorded in the lockfile.

## Data Availability

The complete dataset and associated documentation are available in this
repository and archived on Zenodo.

Zenodo DOI: \[DOI\]

## Citation

If you use this dataset or the associated R function in your research,
please cite the archived dataset and the corresponding publication:
\[Citation\]

## License

The spatial dataset and its associated documentation are made available
under the Creative Commons Attribution 4.0 International (CC BY 4.0)
license. The dataset incorporates data from GEBCO and the Global Islands
Database. When using or redistributing the dataset, please provide
appropriate attribution to the dataset authors and acknowledge the
original data sources. Full license information and source data
references are provided in LICENSE-DATA.

The R scripts and the shallowing function are distributed under the GNU
General Public License v3.0 (GPL-3.0). The full terms of the software
license are provided in LICENSE.
