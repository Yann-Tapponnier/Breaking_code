#####################################################################################
  EMILE_S function use WilcoxAUC that calculates the geometric_mean
              
              FINDALLMARKERS = log ratio of Arithmetic means
              WILCOXAUC = Geometric mean

              
              Find ALl markers is more SPECIFIC
              WILCOXAUC is more SENSITIVE
#####################################################################################

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
names(int_obj1@meta.data)
# Load the clustering name to activate it transiently as "active ident" to be able to calculate cluster from it
meta_data_names_to_use_as_idents <- c( grep("dbscan", colnames(int_obj1@meta.data), value = TRUE), "orig.ident")

# Loop to automatically calculate with different clustering / meta.data

for (name in meta_data_names_to_use_as_idents){
  int_obj1 <- AddMarkerGenesToSeurat(seuratObj = int_obj1,
                                     group_by = name,
                                     filterOutNAcells = T, # IT REMOVE TEMPORALLY THE NA clusters
                                     logFCposMin = 0, # 
                                     logFCnegMax = -0) # 
  
}


view(int_obj1@misc$marker_genes$cerebro_seurat$optics_dbscan100_minptseps_cl0.6)
markers <- FindAllMarkers(int_obj1, group.by = "optics_dbscan100_minptseps_cl0.6", logfc.threshold = 0.1,min.pct = 0.01)
hist(markers$avg_log2FC)                   


DimPlot(int_obj1, group.by = "optics_dbscan100_minptseps_cl0.6")
FeaturePlot(int_obj1, features = "MYT1L")


#### the best for cluster 3
# According to Wilcoxauc the best for cluster 3
FeaturePlot(int_obj1, features = c("SCG2","LINGO2","TUBA1A", "TMEM158","THY1","PCDH7","TNC"))

# Accordign to FindAllMarkers 
FeaturePlot(int_obj1, features = c("LINC01931", 'SYNDIG1',	'RIPOR2','MYT1L',"POU3F2","STMN2" ))







####### TESTING WILCOXAUC VS FINDALLMARKERS
int_obj1_temp <- NormalizeData(int_obj1, scale.factor = 1)
POU3F2 <- data.frame(expr=GetAssayData(int_obj1_temp,assay ="RNA", layer = "data")["POU3F2",],
                     groups=int_obj1_temp$optics_dbscan100_minptseps_cl0.6)

POU3F2_1 <- ((POU3F2 %>% filter(groups==1) %>% select(expr) %>% exp) - 1) %>% unlist %>% mean
POU3F2_not_1 <- ((POU3F2 %>% filter(groups!=1) %>% select(expr) %>% exp) - 1) %>% unlist %>% mean
POU3F2_2 <- ((POU3F2 %>% filter(groups==2) %>% select(expr) %>% exp) - 1) %>% unlist %>% mean
POU3F2_not_2 <- ((POU3F2 %>% filter(groups!=2) %>% select(expr) %>% exp) - 1) %>% unlist %>% mean

log(POU3F2_1/POU3F2_not_1)
log(POU3F2_2/POU3F2_not_2)

MYT1L <- AverageExpression(int_obj1,assays ="RNA", features = "MYT1L", group.by = "optics_dbscan100_minptseps_cl0.6")
MYT1L$RNA
POU3F2 <- AverageExpression(int_obj1,assays ="RNA", features = "POU3F2", group.by = "optics_dbscan100_minptseps_cl0.6")
?NormalizeData
POU3F2 <- AverageExpression(int_obj1_temp,assays ="RNA", layer = "counts", features = "POU3F2", group.by = "optics_dbscan100_minptseps_cl0.6")
POU3F2$RNA

test <-  FindAllMarkers(int_obj1, logfc.threshold = 0.5849, min.pct = 0, group.by ="optics_dbscan100_minptseps_cl0.6", verbose = T )
View(test)

int_obj1@misc$enriched_pathways$cerebro_seurat2_enrichr$optics_dbscan100_minptseps_cl0.6
View(int_obj1@misc$marker_genes$cerebro_seurat$optics_dbscan100_minptseps_cl0.6)
View(int_obj1@misc$marker_genes$cerebro_seurat$orig.ident)
)

table( int_obj1@misc$marker_genes$cerebro_seurat$optics_dbscan100_minptseps_cl0.6$optics_dbscan100_minptseps_cl0.6)
names( int_obj1@misc$marker_genes$cerebro_seurat$optics_dbscan100_minptseps_cl0.6)
names( int_obj1@misc$marker_genes$cerebro_seurat)


################################################################################################################################
################################ New test for pairwise comparison ################################################################
###############################         FINDMARKERS               ################################################################
################################ wilcoxauc with groups_use.       ################################################################
################################ emile function on subset.        ################################################################
################################################################################################################################
ds_objf1 <- qread ('/Users/administrateur/Desktop/Bio_info/Single_Cell/RT_in_vivo/output_seurat/7_RT_in_vivo_integration_join_vanilla_big_cluster_only_scaled_nCount_mt_WO_ANSES_DS_f1_IN_PROGRESS.qs', nthreads = 12)
############### 

png(filename = "report/7_RT_in_vivo_WO_ANSES_ds_Signatures/ds_objf1_umap_DBSCAN_500_cl0.8.png", width = 1500, height = 900, res=100)
DimPlot(ds_objf1, reduction = "umap", group.by = c( "optics_dbscan500_minptseps_cl0.8"), label = T)
dev.off()

############# FINDMARKERS ###############      
#### Changing active idents
Idents(ds_objf1) <- "optics_dbscan500_minptseps_cl0.8"


Repro_signature_markers_5_vs_7_FM <- FindMarkers(
  object = ds_objf1,
  ident.1 = 5, #Target
  ident.2 = 7 # Control
)

# Descending order = Cl5_up
Cl5_up_FM_top1000_FC <- Repro_signature_markers_5_vs_7_FM %>%
  arrange(desc(avg_log2FC)) %>% 
  slice_head(n = 1000) %>% 
  rownames() %>% 
  unlist()
# Ascending order = Cl7_up
Cl7_up_FM_top1000_FC <- Repro_signature_markers_5_vs_7_FM %>%
  arrange(avg_log2FC) %>% 
  slice_head(n = 1000) %>% 
  rownames() %>% 
  unlist()

############# wilcoxauc with groups_use ############### 
Repro_signature_markers_5_vs_7_Wilcox <- wilcoxauc(ds_objf1,
                                                   group_by= "optics_dbscan500_minptseps_cl0.8",
                                                   groups_use = c(5,7) )


Cl5_up_wc_top1000_FC <- Repro_signature_markers_5_vs_7_Wilcox %>% 
  filter(group == 5) %>% 
  arrange(desc(logFC)) %>% 
  slice_head(n = 1000) %>% 
  select(feature) %>% 
  unlist()
Cl7_up_wc_top1000_FC <- Repro_signature_markers_5_vs_7_Wilcox %>%
  filter(group == 7) %>% 
  arrange(desc(logFC)) %>% 
  slice_head(n = 1000) %>% 
  select(feature) %>% 
  unlist()

Cl5_up_wc_top1000_AUC <- Repro_signature_markers_5_vs_7_Wilcox %>% 
  filter(group == 5) %>% 
  arrange(desc(auc)) %>% 
  slice_head(n = 1000) %>% 
  select(feature) %>% 
  unlist()
Cl7_up_wc_top1000_AUC <- Repro_signature_markers_5_vs_7_Wilcox %>%
  filter(group == 7) %>% 
  arrange(desc(auc)) %>% 
  slice_head(n = 1000) %>% 
  select(feature) %>% 
  unlist()

####### Emile function after subseting 2 cluster of interest  ######## 
temp <- subset(ds_objf1, idents = c(5,7))
temp@misc$marker_genes$cerebro_seurat <- NULL
temp <- AddMarkerGenesToSeurat(temp, group_by = "optics_dbscan500_minptseps_cl0.8")


Cl5_up_em <-  temp@misc$marker_genes$cerebro_seurat$optics_dbscan500_minptseps_cl0.8_upreg_only %>% filter (optics_dbscan500_minptseps_cl0.8 == 5)
Cl7_up_em <-  temp@misc$marker_genes$cerebro_seurat$optics_dbscan500_minptseps_cl0.8_upreg_only %>% filter (optics_dbscan500_minptseps_cl0.8 == 7)

Cl5_up_em_top1000_FC <- Cl5_up_em %>% arrange(desc(logFC)) %>% 
  slice_head(n = 1000) %>% 
  select(gene) %>% 
  unlist()

Cl7_up_em_top1000_FC <- Cl7_up_em %>% arrange(desc(logFC)) %>% 
  slice_head(n = 1000) %>% 
  select(gene) %>% 
  unlist()

########### Comparing the methods ################# 
# Emile_on_Subset vs Wilcox_goups_use
length(intersect(Cl5_up_em_top1000_FC, Cl5_up_wc_top1000_FC)) # Intersect are equal !!! 1000
length(intersect(Cl7_up_em_top1000_FC, Cl7_up_wc_top1000_FC)) # Intersect are almost good 997

# FindMakers VS Willcox on FC
length(intersect(Cl5_up_FM_top1000_FC, Cl5_up_wc_top1000_FC)) # 189 only 
length(intersect(Cl7_up_FM_top1000_FC, Cl7_up_wc_top1000_FC)) # 282
# FindMakers VS Willcox on AUC
length(intersect(Cl5_up_FM_top1000_FC, Cl5_up_wc_top1000_AUC)) # 241
length(intersect(Cl7_up_FM_top1000_FC, Cl7_up_wc_top1000_AUC)) # 321


########## Library 
library(VennDiagram)
,
Emile_FC = Cl5_up_em_top1000_FC
####### 
venn.diagram(
  x = list(
    Seurat = Cl5_up_FM_top1000_FC,
    Wilcox_FC = Cl5_up_wc_top1000_FC,
    Wilcox_AUC = Cl5_up_wc_top1000_AUC
  ),
  filename = "report/7_RT_in_vivo_WO_ANSES_ds_Signatures/venn_signature_Seurat_Wilcox2.png",
  fill = c("skyblue", "salmon","yellow"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  resolution = 300
)

# PLOTING THE TOP 9 to show the diffenrece in expression pattern
topToPlot <- c("Cl5_up_wc_top1000_FC", "Cl7_up_wc_top1000_FC", "Cl5_up_wc_top1000_AUC", "Cl7_up_wc_top1000_AUC", "Cl5_up_FM_top1000_FC", "Cl7_up_FM_top1000_FC")

for (names in topToPlot){
  genes = get(names)
  png(filename = paste0("report/7_RT_in_vivo_WO_ANSES_ds_Signatures/ds_objf1_WC_vs_FM_",names,".png"), width = 2400, height = 1600, res=100)
    print(FeaturePlot(ds_objf1, reduction= "umap", features = genes[1:9]))
  dev.off()
}
