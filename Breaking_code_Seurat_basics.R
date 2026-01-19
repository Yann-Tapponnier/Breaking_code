colnames(seurat_obj) # = UMI
rownmaes(seurat_obj) # = Gene names

colnames(seurat_obj@meta.data) # All entries / categories of your metadata eg: orig.ident, nFeature_RNA, etc..
rownmaes(seurat_obj@meta.data) # UMI !! it is the contrary of SPARSE DATA

GetAssayData(seurat_obj,assay = "RNA", slot = "counts" ) # COUNT matrix
GetAssayData(seurat_obj,assay = "RNA", slot = "data" ) # normalized counts (log2 transformed) (divided by number of read for that cell + multiply by the mean numebr read per cell)
seurat_obj[['meta.data_name']]





################################################################################################################################################################################################
######################################################################### SUBSET() function on Seurat Object ###################################################################################
################################################################################################################################################################################################

| Argument       | Type                         | Default | Description                                                                                                             |
  | -------------- | ---------------------------- | ------- | ----------------------------------------------------------------------------------------------------------------------- |
  | **x**          | Seurat object                | —       | The Seurat object to subset.                                                                                            |
  | **subset**     | expression                   | `NULL`  | Logical expression to choose cells or features (e.g., `nFeature_RNA > 200`). Evaluated using metadata or feature stats. |
  | **cells**      | character vector             | `NULL`  | Vector of cell names to keep.                                                                                           |
  | **features**   | character vector             | `NULL`  | Vector of gene/feature names to keep.                                                                                   |
  | **idents**     | character / factor / numeric | `NULL`  | Identity classes to keep (uses `Idents(object)`).                                                                       |
  | **invert**     | logical                      | `FALSE` | If TRUE, inverts the selection (keeps everything *except* selected cells/features).                                     |
  | **downsample** | integer                      | `NULL`  | Randomly subsamples each identity class to this number of cells.                                                        |
  | **seed**       | integer                      | `NULL`  | Random seed for downsampling (defaults to 1 internally).                                                                |
  | **...**        | passed to internal functions | —       | Not typically used by end-users.                                                                                        |
  

  
  
# Subseting 
Idents(merged_obj) <- "orig.ident"

ds_merged <- subset(merged_obj, idents = c("D0", "D2", "D5"))
##### Downsampling by the active ident
ds_merged <- subset(merged_obj, downsample = 3555, seed = 1 ) # By default there is an automatic seed = 1

# checking the number of cells per experiments 
table(ds_merged@meta.data$orig.ident)




#############################################################################################
################################### Manipulating metadata ###################################
############################################################################################# 

table(int_obj1@meta.data$orig.ident)
# creating the new metadata column 
int_obj1@meta.data$IdentForCellChat <- NA
# Doing it manually
int_obj1@meta.data$IdentForCellChat[int_obj1@meta.data$orig.ident=="-Dox"] <- "-DOX"
int_obj1@meta.data$IdentForCellChat[int_obj1@meta.data$orig.ident=="+Dox"] <- "+DOX"

#### Doing it automatically
identsIwant <- names(table(int_obj1@meta.data$SingleR_pruned.labels))
for (name in identsIwant){
  int_obj1@meta.data$IdentForCellChat[int_obj1@meta.data$SingleR_pruned.labels== name] <- name
}

# For checking the output 
table(int_obj1@meta.data$IdentForCellChat, useNA = "ifany")


