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

############################### Dependencies 
library(future)
options(future.globals.maxSize = 10000 * 1024^2)  # 10 GiB

install.packages('BiocManager')
BiocManager::install('glmGamPoi')

###############################  Coding using MAPQUERY #################################
# Import both dataset in a single session
reference # In Vivo
query  # In Vitro


##### Re-running the basic 
reference <- reference %>% 
  NormalizeData %>% # by default "logNormalize" takes it all 
  FindVariableFeatures %>% # 2000 by default
  ScaleData(vars.to.regress = c('nCount_RNA', 'percent.mt',"CC.diff")) %>% # by default use the FindVariableFeatures
  RunPCA %>% 
  RunUMAP( dims = 1:50, return.model = T) # return.model = MANDATORY only for REFERENCE (it is the space where query will be projected) !

                      ############# Not necessary ###############
                      # query <- query %>% 
                      #   NormalizeData %>% # by default takes it all
                      #   FindVariableFeatures %>% # 2000 by default
                      #   ScaleData(vars.to.regress = c('nCount_RNA', 'percent.mt',"CC.diff")) %>% # by default use the FindVariableFeatures
                      #   RunPCA %>% 
                      #   RunUMAP( dims = 1:50, return.model = T) 

anchors <- FindTransferAnchors(
  reference = reference,
  query = query,
  normalization.method = "LogNormalize",  #Default but can be SCT
  reference.reduction = "pca",
  dims = 1:30  # 1:50 c'est beaucoup, 30 est standard
)


# by orig.ident
query <- MapQuery(anchorset = anchors,
                  reference = reference, 
                  query = query,
                  refdata = reference$orig.ident,
                  reference.reduction = "pca",
                  reduction.model = "umap") # Can be NOT added to have the LABEL TRANSFERT

# by BAM labeling
# query <- MapQuery(anchorset = anchors, 
#                   reference = reference,
#                   query = query,
#                   refdata = reference$BAM, 
#                   reference.reduction = "pca", 
#                   reduction.model = "umap")




################### Different Normalization : SCT VS LogNormalized ############################
Ex : reference (from a paper) is SCT transformed a
     query (my data) is LogNormalize
   
# If reference have both : reference@assay :  "RNA", "SCT"  
# You can do Choice 1 or 2 !!!   )))
     
## Chemin 2 — SCTransformer ta query (plus rigoureux car on ne modifie pas la ref qui vient d'un papier)
query <- SCTransform(query, vst.flavor = "v2", variable.features.n = 2000) # remplace findVariableFeature NormalizeData and scaleData / Default 2000
DefaultAssay(query) <- "SCT"
query <- RunPCA(query)
# Pas besoin de RunUMAP sur la query !

anchors <- FindTransferAnchors(
  reference = reference,
  query = query,
  normalization.method = "SCT",
  reference.reduction = "pca",
  dims = 1:30
)
            
query <- MapQuery(
  anchorset = anchors,
  reference = reference,
  query = query,
  referencedata = list(celltype = "celltype"),
  reference.reduction = "pca",
  reduction.model = "umap"
)


###################. Trouble Shooting ############################

                    ################### Different Normalization : SCT VS LogNormalized ############################
                    Ex : reference (from a paper) is SCT transformed a
                    query (my data) is LogNormalize
                    
                    # If reference have both : reference@assay :  "RNA", "SCT"  
                    # You can do Choice 1 or 2 !!!   )))
                    
                    ## Chemin 2 — SCTransformer ta query (plus rigoureux car on ne modifie pas la ref qui vient d'un papier)
                    variable.features.n = 2000) # remplace findVariableFeature NormalizeData and scaledata
                     DefaultAssay(query) <- "SCT"
                    query <- RunPCA(query)
                    # Pas besoin de RunUMAP sur la query !

                              # ERREUR --> Checking
                              # 1. Features communes
                              common_features <- intersect(rownames(reference@assays$SCT), 
                                                           rownames(query@assays$SCT))
                              length(common_features)  # combien ?
                              
                              # 2. Vérifier le PCA de la ref
                              dim(reference@reductions$pca@cell.embeddings)
                              
                              # 3. Vérifier les variable features de la ref
                              length(VariableFeatures(reference))
                              
                              # 4. Vérifier les variable features de la query
                              sum(VariableFeatures(reference) %in% rownames(query@assays$SCT))
                              
                              # Le SCT de la ref a bien une scale.data ?
                              dim(reference@assays$SCT@scale.data)
                              
                              # Le PCA de la ref vient bien de SCT ?
                              reference@reductions$pca@assay.used
                              
                              # La query a bien un SCT avec scale.data ?
                              dim(query@assays$SCT@scale.data)
                              
                              # DefaultAssay des deux ?
                              DefaultAssay(reference)
                              DefaultAssay(query)
                              
                              # Les 2000 features du SCT ref sont bien dans le scale.data de la query ? 531 in commun = too few
                              sum(rownames(reference@assays$SCT@scale.data) %in% 
                                    rownames(query@assays$SCT@scale.data))
                              
                  # 531 in commun = too few
                              
                  #  Increasing variable feature of query to 3000
                  query <- SCTransform(query, vst.flavor = "v2", variable.features.n = 3000) # remplace findVariableFeature NormalizeData and scaleData
                  DefaultAssay(query) <- "SCT"
                  query <- RunPCA(query) 
                  
                  # Still failing with 3000 variable features in the query (but more in commun with the ref : 531)              
                  
                    anchors <- FindTransferAnchors(
                                reference = reference,
                                query = query,
                                normalization.method = "SCT",
                                reference.reduction = "pca",
                                dims = 1:30,
                                features = VariableFeatures(reference), # ici on force l'utilisation des variable features de la ref (et pas de la query)
                                verbose = T
                              )
                              
                              query <- MapQuery(
                                anchorset = anchors,
                                reference = reference,
                                query = query,
                                referencedata = list(celltype = "celltype"),
                                reference.reduction = "pca",
                                reduction.model = "umap"
                              )
                        


## Chemin 2 — Forcer RNA (rapide, à tester en premier)
                        
## Chemin 2 — Forcer RNA (rapide, à tester en premier)
# Ref — recalculer sur RNA
DefaultAssay(reference) <- "RNA"
reference <- reference %>%
  NormalizeData %>%
  FindVariableFeatures %>%
  ScaleData() %>% #ScaleData(vars.to.regress = c('nCount_RNA', 'percent.mt', "CC.diff")) %>%
  RunPCA %>%
  RunUMAP(dims = 1:30, return.model = TRUE)

# Query — pipeline LogNormalize
DefaultAssay(query) <- "RNA"
query <- query %>%
  NormalizeData %>%
  FindVariableFeatures %>%
  ScaleData(vars.to.regress = c('nCount_RNA', 'percent.mt', "CC.diff")) %>%
  RunPCA

# Anchors
anchors <- FindTransferAnchors(
  reference = reference,
  query = query,
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  dims = 1:30
)

# MapQuery
query <- MapQuery(
  anchorset = anchors,
  reference = reference,
  query = query,
  refdata = reference$seurat_clusters,
  reference.reduction = "pca",
  reduction.model = "umap"
)



###################### NO UMAP recalculation  : Transfert CellTypes and Clusters without Projection !!!
tu veux juste transférer les labels (celltype, clusters...)

# Retire simplement reduction.model de MapQuery()
query <- MapQuery(
  anchorset = anchors,
  reference = ref,
  query = query,
  refdata = list(celltype = "celltype"),
  reference.reduction = "pca"
  # pas de reduction.model ici
)


Pas besoin de RunUMAP du tout dans ce cas. ✅

ObjectifRunUMAP nécessaire ?Transférer labels seulement❌ NonProjeter query sur UMAP ref✅ Oui (avec return.model = TRUE)
Tu veux faire quoi exactement — juste le transfert de labels, ou aussi la projection UMAP ?
  



  
  
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









