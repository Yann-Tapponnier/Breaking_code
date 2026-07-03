## SingleR
SingleR is used to label cells
The point here is to check whether a particular identity is attributed to the outlayers cells before removing them
```{r, eval =F}
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

#### 
hpca  <- HumanPrimaryCellAtlasData()
hpca  <- fetchReference("hpca", "2024-02-26") # the same a above


#######
Blueprint (Martens and Stunnenberg 2013) and ENCODE projects (The ENCODE Project Consortium 2012).
ref <- fetchReference("blueprint_encode", "2024-02-26")

### Immunological Genome Project (ImmGen)
immgen <- fetchReference("immgen", "2024-02-26") # Mouse

#### Novershtern hematopoietic data
ref <- fetchReference("novershtern_hematopoietic", "2024-02-26") # Human

### Monaco immune data (Monaco et al. 2019)
ref <- fetchReference("monaco_immune", "2024-02-26") # Human




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
# loading a reference from celldex package
hpca <- HumanPrimaryCellAtlasData()

# Annotation based on Spearman correlation for each cell
pred_join <- SingleR(test = GetAssayData(join_obj), # the data to annotate (join_obj[,1:20] here 20 first cell to test)
                     ref = hpca,
                     labels = hpca$label.main) # choose the labeling output : main / fine / ontology

table(pred_join$labels) # table counts the unique entry per category
table(pred_join$pruned.labels) # pruned. is above certain threshold !!
join_obj@meta.data <- cbind(join_obj@meta.data,SingleR_pruned.labels = pred_join$pruned.labels) # assigning to meta data with the name : SingleR_pruned.labels


DimPlot(join_obj, reduction = "umap", group.by ='SingleR_pruned.labels') 
DimPlot(join_obj, reduction = "umap", group.by ='SingleR_pruned.labels', split.by ='orig.ident' )
DimPlot(join_obj, reduction = "umap", group.by ='SingleR_pruned.labels', split.by ='SingleR_pruned.labels' )


names(join_obj@meta.data)







