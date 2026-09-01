################################## SLIGSHOT ###########################################
BiocManager::install("slingshot")

# You need clusters, you can use ORIG.IDENT but it is bisased
library(slingshot)
source('/Users/administrateur/Desktop/Bio_info/My_codes/Emile_s_functions_single_cell.R')


varfeats <- VariableFeatures(seurat_obj)

# REQUIRE the convertion from SEURAT_obj to SingleCellExperiment_Obj
seurat_obj_filt <- seurat_obj[varfeats,]
sce_filt <- as.SingleCellExperiment(seurat_obj_filt)

sce_filt <- slingshot(sce_filt, 
                      clusterLabels = 'dbscan_100_minpts',  # Choose the clustering you want, dbscan / Orig.ident etc
                      reducedDim = 'UMAP', # Choose the embeding you want (UMAP, DiffMap, PCA)
                      start.clus='1') # NEED to choose 
                      #end.clus= 'Transdiff_D14)

# Slingshot_obj
Slingshot_obj <- SlingshotDataSet(sce_filt)
Slingshot_obj@lineages # show the progressing  through clusters
Slingshot_obj@adjacency # 
Slingshot_obj@curves # select some curve to plot or not
Slingshot_obj@slingParams # 



# To visualize the trajectories (without passing to Cerebro):
library(RColorBrewer)
n <- length(unique(sce_filt$optics_dbscan10_minptseps_cl0.4))
my_colors <- colorRampPalette(brewer.pal(n,'Set3'))
# my_colors(n)  --> pick n colors in the palette
plotcol <- my_colors(n)[as.numeric(sce_filt$colors2)] # assign the 4 colors to the different entry ! 



png(filename = "Report/7_side_project_ReproTransfo_D2T16/Trajectory_UMAP_DO.png", width = 1600, height = 1200 )
  plot(reducedDims(sce_filt)$UMAP, col = plotcol, pch=16, asp = 1)
  lines(SlingshotDataSet(sce_filt), lwd=2, col='black')
  title("Trajectory on UMAP - DO")
dev.off()  

# OR
library(RColorBrewer)
colors <- colorRampPalette(brewer.pal(11,'Spectral'))
plotcol <- colors(length(unique(sce_filt$dbscan_100_minpts))[as.numeric(sce_filt$dbscan_100_minpts)])

plot(reducedDims(sce_filt)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce_filt), lwd=2, col='black')








### the name "monocle2" only is recognized by CEREBRO It is a trick to save the trajectory and visualised it.
seurat_obj@misc$trajectories$monocle2[['trajectory_name']] <- MassageSlingshotResult_3d_and_2d(sce_filt)





#IF you want you can assign colors to metadata entry : 
    ##### Creating a slot in metadata (colData) populated with numerical value for the category/group we want to plot.
    sce_filt@colData # is real metadata
  sce_filt@colData$colors <- 1
  n <- 0
  for (i in unique(sce_filt$orig.ident)) {
    n <- n+1
    sce_filt@colData$colors[sce_filt$orig.ident == i] <- n
  } # is real metadata
  
  sce_filt@colData # is real metadata
  sce_filt@colData$colors2 <- 1
  n <- 0
  for (i in unique(sce_filt$optics_dbscan10_minptseps_cl0.4))) {
    n <- n+1
    sce_filt@colData$colors2[sce_filt$optics_dbscan10_minptseps_cl0.4 == i] <- n
  } # is real metadata
  
  table(sce_filt@colData$colors)  # = 4 for orgi.idents   
  table(sce_filt@colData$colors2 ) # = 7 (6+NA) for optics_dbscan10_minptseps_cl0.4
  
  
  
  
  ##### Downstream Analysis 
  # 4.1 Identifying temporally dynamic genes
  # BiocManager::install("tradeSeq")
  library(tradeSeq)
  
  # fit negative binomial GAM
  sce <- fitGAM(sce_filt)
  
  # test for dynamic expression
  ATres <- associationTest(sce)
  
  
  
  
  
  