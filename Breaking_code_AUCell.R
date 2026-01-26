#################################################################################################
##################### Highlighting the differentiated cell to label #####################
library(AUCell)

### DEFINTING A LIST OF MARKERS
          #### Manual definitinon of list
          Genesets <- list(AT2=c("Sftpc", "Lyz2"),
                           Ciliated_cell = c("Foxj1",'Cd24a'),
                           Club_cell = c("Scgb1a1", "Scgb3a2"),
                           AT1 = c("Aqp5" , "Pdpn")
          )

#OR
          #### Extracting list from dataframe using split
          NB_Markers <- read.xlsx(xlsxFile =  "data/List_genes_on_interest_NB.xlsx", colNames = T)
          # Making it a list using "split"
          List_NB_markers <- split (NB_Markers$symbol,NB_Markers$cell_type)
                          
                      # Individual ploting
                          # for (i in 1:length(List_NB_markers)){
                          #   name <- names(List_NB_markers)[i] # retreiving the name
                          #   Dim <- ceiling(sqrt(length(List_NB_markers[[i]]))) # Taking the square root rounded upper (ceiling) for miminum number of slot for plots
                          #   
                          #   png(filename = paste0("report/3_NK_integration_clustering/int_obj1_Manual_individual_Markers_UMAP",name,".png"), width = Dim*400, height = Dim*300 )
                          #   print(FeaturePlot(int_obj1, features = List_NB_markers[[i]], pt.size = 1, ncol= Dim)+
                          #           plot_annotation(title =paste0("UMAP Feature Plot -", name)) &       # Adding the plot name of the category
                          #           theme( plot.title = element_text(hjust = 0.5, face = "bold"))
                          #   )
                          #   dev.off()
                          # }
          
          
#OR 
          
          # Extracting column as independant list
          List_Makers <- setNames(
            as.list(omalley),     # turn each column into a list element
            colnames(omalley)     # assign column names to the list thanks to setNames
          )

          
          

          
Counts <- GetAssayData(RT_3_ds4, assay = "RNA", layer= "counts") #AUCell is not afected if using COUNTS or DATA !!
cells_rankings <-AUCell_buildRankings(Counts, nCores = 12, plotStats = TRUE)

# Calculate enrichment scores
cells_AUC <- AUCell_run(exprMat = Counts,
                        geneSets = Genesets, 
                        aucMaxRank=nrow(cells_rankings)*0.05)

# View AUC scores for cells
head(getAUC(cells_AUC))

    # Plot histogram
    png(filename = paste0("report/3_NK_integration_clustering/int_obj1_Markers_AUCscore_histogramm.png"), width = 1600, height = 1200 )
        par(mfrow = c(2, 2)) # to save it 2x2
        for (i in seq(1:length(Genesets))){
          name <- names(Genesets)[i]
          
          hist(getAUC(cells_AUC)[name, ], main = paste0("AUC for ",name))
          
        }
par(mfrow = c(1, 1)) # to revert toward normal plotting  

### Check
auc_matrix <- as.data.frame(t(getAUC(cells_AUC)))  # Transpose to match Seurat: cells = rows
answer <- all(rownames(auc_matrix) %in% colnames(RT_3_ds4))  # Should be TRUE
print( "testing is all rownames of AUC_MATRIX are the colnames of ", obj_name,' and the answer is :', answer)

# Intergrating the result in Metadata
RT_3_ds4 <- AddMetaData(RT_3_ds4, metadata = auc_matrix)

# ploting
png(filename = paste0("report_seurat/2_inference_project_downsampling_D4/RT_3_ds4_UMAP_Signatures.png"), width = 1800, height = 1350 )
FeaturePlot(RT_3_ds4, features = c("AT2", "Ciliated_cell", "Club_cell", "AT1"), pt.size = 1)
dev.off()





###################### Puting NA below a certain threshold
##### semi-Automation ######## 
threshold <-  c (0.5)
sig_names <- c("AT2", "Ciliated_cell", "Club_cell", "AT1")


for (sig in sig_names) {
  message("Plotting ", sig)
  
  # Create a new metadata column with NA for low-score cells
  highlight_col <- paste0(sig, "_highlight")
  RT_3_ds4[[highlight_col]] <- "NA"
  
  # Safely extract the metadata column as a numeric VECTOR !!!
  sig_vector <- FetchData(RT_3_ds4, vars = sig)[, 1]
  
  # Create a highlight column: keep values > threshold, else NA
  RT_3_ds4[[highlight_col]] <- ifelse(sig_vector > threshold, sig_vector, 0)# Returns NA as real number
  
  # Plot, grey cells with NA values
  p <- FeaturePlot(RT_3_ds4, features = highlight_col) + 
    ggtitle(paste0(sig, " (AUC > ", threshold, ")"))
  
  png(filename = paste0("report_seurat/2_inference_project_downsampling_D4/RT_3_ds4_UMAP_Signature",sig,".png"), width = 1800, height = 1350 )
  print(p)  # show the plot
}





##############  Automation ############## 
##### Ploting with threshold
 
# Define the quantiles and keep probs as names
probs <- c(0, 0.25, 0.5, 0.75, 0.9, 0.99)
thresholds <- quantile(int_obj1$NK_Cl3_up_top100_FC, probs = probs)
names(thresholds) <- probs  # assign the probs as names

sig_names <- c("NK_Cl3_up_top100_FC")

library(patchwork)

for (sig in sig_names) {
  
  message("Plotting all thresholds for ", sig)
  
  plots_list <- list()  # store plots for this signature
  
  # Safely extract the metadata column as numeric
  sig_vector <- FetchData(int_obj1, vars = sig)[, 1]
  
  for (threshold_name in names(thresholds)) {
    
    threshold_value <- thresholds[threshold_name]  # the numeric value
    
    # Create a temporary highlight column
    highlight_col <- paste0(sig, "_highlight_", threshold_name)
    int_obj1[[highlight_col]] <- ifelse(sig_vector > threshold_value, sig_vector, 0)
    
    # Create the FeaturePlot
    p <- FeaturePlot(int_obj1, features = highlight_col) +
      ggtitle(paste0(sig, " (quantile ", threshold_name, ")")) +
      NoLegend()
    
    plots_list[[threshold_name]] <- p
  }
  
  # Combine all plots with patchwork
  ncol_layout <- ceiling(sqrt(length(thresholds)))
  combined_plot <- wrap_plots(plots_list, ncol = ncol_layout)
  
  # Save as a single PNG
  png(filename = paste0("report/3_NBEI1_integration_clustering_DOX/int_obj1_UMAP_Signature_", sig, ".png"),
      width = ncol_layout*600, height = ncol_layout*400, res = 100)
  print(combined_plot)
  dev.off()
}

