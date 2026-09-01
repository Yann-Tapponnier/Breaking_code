## SingleR
SingleR is used to label cells
The point here is to check whether a particular identity is attributed to the outlayers cells before removing them

SingleR to label cells
https://bioconductor.org/packages/release/bioc/vignettes/SingleR/inst/doc/SingleR.html

https://bioconductor.org/packages/release/data/experiment/vignettes/celldex/inst/doc/userguide.html

############################################################################################################
BiocManager::install('SingleR')
library(SingleR)

BiocManager::install('celldex')
library(celldex)

browseVignettes("celldex")  # to visualise all datasets
surveyReferences() # To visualize the dataset presents in celldex
temp <- surveyReferences()

temp$title

#### 
hpca  <- HumanPrimaryCellAtlasData()
hpca  <- fetchReference("hpca", "2024-02-26") # the same a above


####### Human bulk RNA-seq
Blueprint (Martens and Stunnenberg 2013) and ENCODE projects (The ENCODE Project Consortium 2012).
ref <- fetchReference("blueprint_encode", "2024-02-26")

### Immunological Genome Project (ImmGen)
immgen <- fetchReference("immgen", "2024-02-26") # Mouse

#### Novershtern hematopoietic data
ref <- fetchReference("novershtern_hematopoietic", "2024-02-26") # Human

### Monaco immune data (Monaco et al. 2019)
ref <- fetchReference("monaco_immune", "2024-02-26") # Human

### Bulk RNA-seq data of sorted mouse cell populations
ref <- fetchReference("mouse_rnaseq", "2024-02-26") # Mouse



#######################################################################################################################
BiocManager::install('scRNAseq') 
library("scRNAseq")

data from (La Manno et al. 2016) especially hESCs
hESCs <- LaMannoBrainData('human-es')
sceM <- MuraroPancreasData()

*Basically SingleR* works like this : Given a reference dataset of samples (single-cell or bulk) with known labels, it labels new cells from a test dataset based on similarity to the reference. 
*Advanced SingleR :* the package also provides more advanced functionality that includes the use of multiple references simultaneously, manipulating the cell ontology and improving performance on big datasets. 

http://127.0.0.1:29244/library/celldex/doc/userguide.html









#############################################################################################################
# RESETING THE METADATA TO AVOID BUGS
# names(pre_obj@meta.data)
# pre_obj@meta.data <- pre_obj@meta.data[1:86]

# loading a reference from celldex package
### Immunological Genome Project (ImmGen)
immgen <- fetchReference("immgen", "2024-02-26") # Mouse
mouse_rnaseq <- fetchReference("mouse_rnaseq", "2024-02-26") # Mouse
hpca  <- HumanPrimaryCellAtlasData()


####### CHANGE THE REFERENCE BASED ON THE SPECIES OF YOUR DATASET BY " SEARCH / REMPLACE IN THE CODE IN R" ####### 


# Annotation based on Spearman correlation for each cell
temp_single_R<- SingleR(test = GetAssayData(pre_obj), # the data to annotate (join_obj[,1:20] here 20 first cell to test)
                        ref = mouse_rnaseq,
                        labels = mouse_rnaseq$label.main) # choose the labeling output : main / fine / ontology


table(temp_single_R$labels) # table counts the unique entry per category
table(temp_single_R$pruned.labels) # pruned. is above certain threshold !! generates NA

pre_obj@meta.data <- cbind(pre_obj@meta.data, SingleR_mouse_rnaseq = temp_single_R$pruned.labels) # assigning to meta data with the name : SingleR_pruned.labels
pre_obj$SingleR_mouse_rnaseq[is.na(pre_obj$SingleR_mouse_rnaseq)] <- "NA" # Assigning manually STRINGS as a category called NA

table(pre_obj$SingleR_mouse_rnaseq,useNA = "always" )

png(filename = paste0("report/2_RT_EV1_pre_integration_cells_to_keep/",Object_name,"/pre_obj_UMAP_SingleR_mouse_rnaseq.png"), width = 1200, height = 900, res = 100)
  print(DimPlot(pre_obj, reduction = "umap", group.by ='SingleR_mouse_rnaseq', label = T, repel = T))
dev.off()

DimPlot(pre_obj, reduction = "umap", group.by ='SingleR_mouse_rnaseq', label = T, repel = T) 


###### Removing tiny clusters (<10Cells)
# cleaning the groups with less than 20 cells 
tiny_cluster <- names(which(table(pre_obj$SingleR_mouse_rnaseq)<=20))
# Getting the UMI
tiny_cluster_UMI <-  pre_obj@meta.data %>% filter(SingleR_mouse_rnaseq %in%  tiny_cluster) %>% rownames()
# Assinging NA to them
pre_obj$SingleR_mouse_rnaseq[names(pre_obj$SingleR_mouse_rnaseq) %in% tiny_cluster_UMI] <- "NA" # assigning false NA because SEURAT BUGS in DECEMBER 2025

table(pre_obj$SingleR_mouse_rnaseq,useNA = "always" )


## Ploting again
png(filename = paste0("report/2_RT_EV1_pre_integration_cells_to_keep/",Object_name,"/pre_obj_UMAP_SingleR_mouse_rnaseq_curated.png"), width = 1200, height = 900, res = 100)
  DimPlot(pre_obj, reduction = "umap", group.by ='SingleR_mouse_rnaseq', label = T, repel = T) 
dev.off()

DimPlot(pre_obj, reduction = "umap", group.by ='SingleR_mouse_rnaseq', label = T, repel = T)

qsave(pre_obj, paste0( "output_seurat/2_",Object_name,"_join_vanilla_Reg_nCount_MT_CCdiff.qs"), nthreads = 12)



