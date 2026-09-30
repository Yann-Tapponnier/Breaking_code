################################## SLIGSHOT ###########################################
BiocManager::install("slingshot")

# You need clusters, you can use ORIG.IDENT but it is bisased
library(slingshot)
source('/Users/administrateur/Desktop/Bio_info/My_codes/Emile_s_functions_single_cell.R')

varfeats <- VariableFeatures(OSKM_KRASsub0)

# REQUIRE the convertion from OSKM_KRASsub0 to SingleCellExperiment_Obj
OSKM_KRASsub0_filt <- OSKM_KRASsub0[varfeats,]
sce_filt <- as.SingleCellExperiment(OSKM_KRASsub0_filt)

ClusterLabels <- 'orig.ident'  # Choose the clustering you want, dbscan / Orig.ident etc
sce_filt <- slingshot(sce_filt, 
                      clusterLabels = ClusterLabels,  # Choose the clustering you want, dbscan / Orig.ident etc
                      reducedDim = 'UMAP', # #CAPITAL FOR scexperiement / Choose the embeding you want (UMAP, DiffMap, PCA)
                      start.clus='KRAS_D2') # NEED to choose 
#end.clus= 'Transdiff_D14)

# Slingshot_obj
Slingshot_obj <- SlingshotDataSet(sce_filt)
        Slingshot_obj@lineages # show the progressing  through clusters
        Slingshot_obj@adjacency # 
        Slingshot_obj@curves # select some curve to plot or not
        Slingshot_obj@slingParams # 


################################## Assigning the color ###########################################

        # To visualize the trajectories (without passing to Cerebro):
# easy version
library(RColorBrewer)
n <- length(unique(sce_filt[[ClusterLabels]]))
my_colors <- colorRampPalette(brewer.pal(n,'Accent'))
# my_colors(n)  --> pick n colors in the palette
plotcol <- my_colors(n)[as.numeric(sce_filt[[ClusterLabels]])] # assign the n colors to the different entry, like remplacing the cluster name by the hexa ! 

                     # VERSION COLORISATION COMPLICATED by INDICES : here I worked on a subseted object so the number of cluster remaining does no correspond to the INDICES of clusters
                     # indices : 2 3 4 5 7 8 9 10 (it miss 1 and 6) --> so the color assigned to 1 2 3 4 5 6 7 8 and 9 and 10 will be NA !!
                     library(RColorBrewer)

                     clusters <- factor(sce_filt[[ClusterLabels]])
                     n <- nlevels(clusters)
                     my_colors <- colorRampPalette(brewer.pal(8, "Set3"))
                     plotcol <- my_colors(n)[as.numeric(clusters)]
                     table(plotcol)

png(filename = "Report/7_side_project_ReproTransfo_D2T16/Trajectory_UMAP_DO.png", width = 1600, height = 1200 )
  plot(reducedDims(sce_filt)$UMAP, col = plotcol, pch=16, asp = 1)
  lines(SlingshotDataSet(sce_filt), lwd=2, col='black')
  title("Trajectory on UMAP - D2")
dev.off()  

                  ## MANUAL COLORATION :
                  # IF you want you can assign colors to metadata entry : 
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






# OR colorRampPalette
library(RColorBrewer)
# Easy version
colors <- colorRampPalette(brewer.pal(11,'Spectral'))
plotcol <- colors(length(unique(sce_filt$dbscan_100_minpts))[as.numeric(sce_filt$dbscan_100_minpts)])

          # # VERSION COLORISATION COMPLICATED by INDICES
          clusters <- factor(sce_filt[[ClusterLabels]])
          n <- nlevels(clusters)
          colors <- colorRampPalette(brewer.pal(n,'Spectral'))
          plotcol <- colors(n)[as.numeric(clusters)]
          
png(filename = "Report/7_side_project_ReproTransfo_D2T16/Trajectory_UMAP_DO.png", width = 1600, height = 1200 )         
  plot(reducedDims(sce_filt)$UMAP, col = plotcol, pch=16, asp = 1)
  lines(SlingshotDataSet(sce_filt), lwd=2, col='black')
dev.off()

                  

################################## Ploting the lines ###########################################
# Ploting only Line 1 and 3
png(filename = paste0("report/3_RA_integration_clustering/Transdiff/Slingshot_Transdiff_2.png"), width = 1200, height = 900, res = 100 )
plot(reducedDims(sce_filt)$UMAP, col = plotcol, pch=16, asp = 1)
  lines(Slingshot_obj@curves[[1]], lwd = 2)
  lines(Slingshot_obj@curves[[3]], lwd = 2)
dev.off()

# Ploting only Line 1 and 3
png(filename = paste0("report/3_RA_integration_clustering/Transdiff/Slingshot_Transdiff_3.png"), width = 1200, height = 900, res = 100 )
plot(reducedDims(sce_filt)$UMAP, col = plotcol, pch=16, asp = 1)
  lines(Slingshot_obj@curves[[1]], lwd = 2)
  lines(Slingshot_obj@curves[[2]], lwd = 2)
dev.off()


# Other way to plot directly without extracting the Slingshot obj before
plot(reducedDims(sce_filt)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce_filt), lwd=2, col=c('black', "blue"), lineage = c(1,3))


################################### ADDING ARROWS TO THE END OF EACH LINEAGE ###########################################
        # Créating manually arrow at the tip of the lines
        sds <- SlingshotDataSet(sce_filt)
        lineages <- c(1, 2)  # Specify the lineages you want to plot
        cols <- c('red', 'blue')
        names(cols) <- lineages
        Lwd = 4
        
        # Plot de base
        png(filename = "Report/4_RT_EV1_further_analysis/OSKM_KRASsub0/OSKM_KRASsub0_UMAP_D2_3.png", width = 1600, height = 1200 )     
        plot(reducedDims(sce_filt)$UMAP, col = plotcol, pch = 16, asp = 1)
        lines(sds, lwd = Lwd , col = cols, lineage = lineages)
        
        # Récupération des courbes ajustées
        crvs <- slingCurves(sds)
        
        # Ajout d'une pointe de flèche à l'extrémité de chaque lignée
        for (i in lineages) {
          crv <- crvs[[i]]
          pts <- crv$s[crv$ord, ]          # points ordonnés selon le pseudotemps
          n <- nrow(pts)
          
          arrows(x0 = pts[n - 1, 1], y0 = pts[n - 1, 2],
                 x1 = pts[n, 1],     y1 = pts[n, 2],
                 length = 0.15, angle = 20,
                 col = cols[as.character(i)], lwd = Lwd )
        }
        dev.off()



################################## Exporting the line to Seurat Obj ###########################################
### the name "monocle2" only is recognized by CEREBRO It is a trick to save the trajectory and visualised it.
seurat_obj@misc$trajectories$monocle2[['trajectory_name']] <- MassageSlingshotResult_3d_and_2d(sce_filt)


        
        
################################## Downstream Analysis  ###########################################
  # 4.1 Identifying temporally dynamic genes
  # BiocManager::install("tradeSeq")
  library(tradeSeq)
  
  # fit negative binomial GAM
  sce <- fitGAM(sce_filt)
  
  # test for dynamic expression
  ATres <- associationTest(sce)
  
  
  # topgenes <- rownames(ATres[order(ATres$pvalue), ])[1:500]
  # Order on pVal is not good enough 
  expr_frac <- rowMeans(assays(sce)$counts > 0)
  ATres$expr_frac <- expr_frac[rownames(ATres)]
  
  ATres_filt <- subset(
    ATres,
    pvalue < 0.05 &
      expr_frac > 0.1
  )
  
  topgenes <- rownames(
    ATres_filt[order(ATres_filt$pvalue), ]
  )[1:250]
  
  
  pst.ord <- order(sce$slingPseudotime_3, na.last = NA)
  #heatdata <- assays(sce)$counts[topgenes, pst.ord]
  heatdata <- as.matrix(assays(sce)$counts[topgenes, pst.ord, drop = FALSE])
  #heatclus <- sce$GMM[pst.ord]
  heatclus <- sce$orig.ident[pst.ord]
  n <- length(levels(heatclus))
  
  heatmap(log1p(heatdata), Colv = NA,
          ColSideColors = brewer.pal(n,"Set2")[heatclus])
  
  
  
  
  