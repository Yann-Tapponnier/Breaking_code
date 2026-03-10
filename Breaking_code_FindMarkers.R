## Finding markers
Presentation of the function stored in Emile's' function

###### AddMarkers to seurat with findAllMarkers. ######
# remove the NA cluster temporarly
# Store the info for cerebro
AddMarkerGenesToSeurat # Presentation 
                          function(seuratObj, 
                                   group_by,  # the meta.data you want to use
                                   filterOutNAcells= T, # it will remove the NA from clustering
                                   pct_in_min=0, pct_out_min=0, # 
                                   logFCposMin=0, logFCnegMax=-0)

                            
################################### Finding Markers ##############################
# Visualisation of meta.data names 
  names(join_obj_4@meta.data)
# Load the clustering name to activate it transiently as "active ident" to be able to calculate cluster from it
meta_data_names_to_use_as_idents <- c(  tail(names(int_obj1@meta.data),3), "orig.ident")


# Loop to automatically calculate with different clustering / meta.data

for (name in meta_data_names_to_use_as_idents){
      join_obj_4 <- AddMarkerGenesToSeurat(seuratObj = join_obj_4,
                                         group_by = name,
                                         filterOutNAcells = T ) # IT REMOVE TEMPORALLY THE NA clusters
                                  
}


################################### Manually Finding Markers ##############################
FindAllMarkers(object,
               assay = NULL, # default RNA
               features = NULL, # geneset to test / Default all genes
               group.by = NULL, # wich clusters to test against each other / Default = acive Idents.
               logfc.threshold = 0.1,
               test.use = "wilcox",
               slot = "data",
               min.pct = 0.01,
               random.seed = 1,
               min.cells.feature = 3,
               min.cells.group = 3)

D0_VS_D2_FindMarkers(object = obj, ident.1 = "D0", ident.2 = "D2", group.by = "orig.ident")
