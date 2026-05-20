################################# Azimuth for Atlas ref ######################################
BE CAREFULL if cell in the query dataset that are not represented in the reference,
they will project to the MOST SIMILAR cell in the reference.  SO IT CAN MASK THE DISCOVERY of NEW CELLTYPE

https://satijalab.github.io/azimuth/articles/run_azimuth_tutorial.html


# install Azimuth directly from main branch
update.packages(oldPkgs = c("withr", "rlang"))
#remotes::install_github('satijalab/azimuth', ref = 'master')

devtools::install_github("satijalab/seurat-data", "seurat5")

# #library(Azimuth) DOES NOT WORK !!!
# # returns a Seurat object named fetusref
# #fetusref <- LoadData("fetusref")
# #data("fetusref")
# # Get SeuratData installation path
# 
# # The RunAzimuth function can take a Seurat object as input
# fetusref <- RunAzimuth(fetusref, reference = "fetusref")

# remove.packages(c("SeuratData", "Azimuth", "rappdirs"))
# .libPaths()  # shows your library paths
# 
# install.packages("rappdirs")        # fix the corrupt dependency first
# install.packages("Seurat")          # core Seurat package
# install.packages("SeuratData")      # Seurat reference datasets

library(Seurat)
library(SeuratData)

# Install dataset you want
Avdata <- AvailableData()
options(timeout = 1800)  # increase timeout to 30 minutes

InstallData("fetusref") # Installing human feotus ref for having neurblast + schwan cells
OR download directly through Zenodo --> There is the .RDS files !!
  
  
  # check all installed datasets
  available <- InstalledData()
available[grep("fetus", available$Dataset), ]


# FINDING the folder of installed terminal : 
system()::"locate" fetusref.DeuratData
setwd("/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/fetusref.SeuratData")  # Set working directory up to the file that contains

# your .rds and .rdb files. 
fetus <- readRDS("azimuth/ref.rds")

DimPlot(fetus,group.by = "annotation.l1", label = TRUE)  + NoLegend()
DimPlot(fetus,group.by = "annotation.l2", label = TRUE)  + NoLegend()

DimPlot(reduction= "refDR" ,fetus,group.by = "annotation.l1", label = TRUE)  + NoLegend() # It is PCA redcution.



