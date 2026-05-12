 ################################################################################
################################# MAPPING ######################################
################################################################################


################################# RESSOURCES ######################################
https://satijalab.org/seurat/articles/integration_mapping
https://satijalab.org/seurat/articles/multimodal_reference_mapping.html

################################# General idea ######################################
Seurat also supports the projection of reference data (or meta data) onto a query object.
While many of the methods are conserved (both procedures begin by identifying anchors), there are two important distinctions between data transfer and integration:
- In data transfer, Seurat does not correct or modify the query expression data.
- In data transfer, Seurat has an option (set by default) to project the PCA structure of a reference onto the query,
instead of learning a joint structure with CCA. 
We generally suggest using this option when projecting data between scRNA-seq datasets.


--> MapQuery() is a wrapper around three functions: TransferData(), IntegrateEmbeddings(), and ProjectUMAP(). 
TransferData() is used to transfer cell type labels and impute the ADT values;
IntegrateEmbeddings() is used to integrate reference with query by correcting the query’s projected low-dimensional embeddings; 
ProjectUMAP() is used to project the query data onto the UMAP structure of the reference. 



###############################  Coding using MAPQUERY #################################
# Import both dataset in a single session
reference # In Vivo
query  # In Vitro

##### Re-running the basic 
reference <- reference %>% 
  NormalizeData %>% # by default takes it all
  FindVariableFeatures %>% # 2000 by default
  ScaleData(vars.to.regress = c('nCount_RNA', 'percent.mt',"CC.diff")) %>% # by default use the FindVariableFeatures
  RunPCA %>% 
  RunUMAP( dims = 1:50, return.model = T) # return.model = MANDATORY !


query <- query %>% 
  NormalizeData %>% # by default takes it all
  FindVariableFeatures %>% # 2000 by default
  ScaleData(vars.to.regress = c('nCount_RNA', 'percent.mt',"CC.diff")) %>% # by default use the FindVariableFeatures
  RunPCA %>% 
  RunUMAP( dims = 1:50, return.model = T) # return.model = MANDATORY !

anchors <- FindTransferAnchors(
  reference = reference,
  query = query,
  reference.reduction = "pca"
)


# by orig.ident
# query <- MapQuery(anchorset = anchors, reference = reference, query = query,
#     refdata = reference$orig.ident, 
#     reference.reduction = "pca", 
#     reduction.model = "umap")

# by BAM labeling
query <- MapQuery(anchorset = anchors, reference = reference, query = query,
                  refdata = reference$BAM, 
                  reference.reduction = "pca", 
                  reduction.model = "umap")

################################# IF MapQuery do not work DO THAT instead : ######################################

# IF MapQuery do not work :  
pancreas.query <- TransferData(anchorset = pancreas.anchors, reference = pancreas.ref, query = pancreas.query,
                               refdata = list(celltype = "celltype"))
pancreas.query <- IntegrateEmbeddings(anchorset = pancreas.anchors, reference = pancreas.ref, query = pancreas.query,
                                      new.reduction.name = "ref.pca")
pancreas.query <- ProjectUMAP(query = pancreas.query, query.reduction = "ref.pca", reference = pancreas.ref,
                              reference.reduction = "pca", reduction.model = "umap")


# TransferData() returns a matrix with predicted IDs and prediction scores, which we can add to the query metadata. 
# it classify the query cells based on reference data (a vector of reference cell type labels)
predictions <- TransferData(anchorset = anchors, 
                            refdata = reference$orig.ident)

# You can save it in metadata
query <- AddMetaData(query, metadata = predictions[,c(1,ncol(predictions))])

query <- IntegrateEmbeddings(
  anchorset = anchors,
  reference = reference,
  query = query, 
  new.reduction.name = "projected.pca"
)

query <- ProjectUMAP(
  query = query, 
  query.reduction = "projected.pca", 
  reference = reference, 
  reference.reduction = "pca", 
  reduction.model = "umap",
  reduction.name = 'projected.umap'
)

# Ploting 

p1 <- DimPlot(reference, reduction = "umap", group.by = "orig.ident", label = TRUE, label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(query, reduction = "ref.umap", group.by = "BAM", label = TRUE,
              label.size = 3, repel = TRUE) + ggtitle("Query transferred labels")

png(filename = paste0("report/3_NBEI1_integration_clustering_DOX/Mapping_In_Vivo_on_In_VITRO_4.png"), width = 2400, height = 900, res=100 )
print(p2)
dev.off()




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













