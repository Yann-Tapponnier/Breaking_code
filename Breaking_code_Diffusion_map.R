
# DIFFUSION MAP 
######################################################## PCA = FAST VERSION ###########################################################################################
ElbowPlot(int_obj1)


###### Calculating diffusion map on PCA axis (FAST VERSION)
# npcs = number of PCA componant to take --> YOU HAVE to  play around manually
#number_of_PCA = 10
number_of_PCA = seq(from = 2, to =10, by = 1) #  # Try low number based on Elbow plot - Over a certain point the increase just colapse the plot


for (npcs in number_of_PCA){
  data4diffmap <- int_obj1[['pca']]@cell.embeddings[,1:npcs]
  
  options(Matrix.warnDeprecatedCoerce = 0) # transform and error back in a warning to not crash the calculation !
  diffmappo <- destiny::DiffusionMap(data = data4diffmap, verbose = T)
  
  # Do not forget to assign a unique name to the diffusion reduction
  int_obj1[[paste('diffmap',npcs)]] <- CreateDimReducObject(embeddings = 10000*diffmappo@eigenvectors, key ='diffmap', assay = 'RNA') #BE CAREFULL to have a unique name form
  # 10000 is just to adapt the scale to cerebro
}


#### Manual Print
DimPlot(int_obj1, reduction = 'diffmap_10', group.by = c('orig.ident',"Phase"))
DimPlot(int_obj1, reduction = paste('diffmap',npcs,sep="_"), group.by = c(names(tail(int_obj1@meta.data,3)),'orig.ident'))

### Printing out each plot from the loop in 1 file
embeds <- c(names(tail(int_obj1@reductions, length(number_of_PCA))))
p <- NULL
for (i in embeds){
  p[[i]] <- DimPlot(int_obj1, reduction = i, group.by = 'orig.ident')
}

png(filename = "report/BLABALBLABLA/int_obj1_diffmap_PCA.png", width = 2400, height = 2000 )
CombinePlots(plots = p, ncol = 3 )
dev.off() 


####################################################### High Variable Gene = LONG VERSION ##########################################################################################
###### Calculating diffusion map on the most High Variable Features (LONG BUT BETTER VERSION)

################ calculating and taking the 500 variable genes #############################################
int_obj1 <- FindVariableFeatures(int_obj1, nfeatures = 500) # 500genes 3 min for 13 000 Cells / 1000 genes 9 mins / 2000 genes 36min 
top500hvg <- VariableFeatures(int_obj1)
data4diffmap <- t(GetAssayData(int_obj1,)[top500hvg,]) # transpose / GetAssayData gather the matrix normalised data (Gene X Cell)
# and retreive only 500 row 

# Caclulating the Diffusion map
options(Matrix.warnDeprecatedCoerce = 0)# transform and error back in a warning to not crash the calculation !
diffmappo <- destiny::DiffusionMap(data = as.matrix(data4diffmap), verbose = T)
int_obj1[['diffmap500hvg']] <- CreateDimReducObject(embeddings = 10000*diffmappo@eigenvectors, key = "diffmap500hvg", assay = 'RNA')


#DimPlot(int_obj1, reduction = 'diffmap500hvg', group.by = c(tail(names(int_obj1@meta.data),2),'orig.ident'))

png(filename = "report/BLABALBLABLA/int_obj1_diffusion_map_hvg_500.png", width = 1600, height = 1200 )
DimPlot(int_obj1, reduction = 'diffmap500hvg', group.by ='orig.ident')
dev.off()


################# calculating and taking the 1000 variable genes ################
int_obj1 <- FindVariableFeatures(int_obj1, nfeatures = 1000) # 1000genes 3 min for 13 000 Cells / 1000 genes 9 mins / 2000 genes 36min 
top1000hvg <- VariableFeatures(int_obj1)
data4diffmap <- t(GetAssayData(int_obj1, )[top1000hvg,]) # transpose / GetAssayData gather the matrix normalised data (Gene X Cell)
# and retreive only 1000 row 

### Caclulating the Diffusion map  
options(Matrix.warnDeprecatedCoerce = 0)# transform and error back in a warning to not crash the calculation !
diffmappo <- destiny::DiffusionMap(data = as.matrix(data4diffmap), verbose = T)
int_obj1[['diffmap1000hvg']] <- CreateDimReducObject(embeddings = 10000*diffmappo@eigenvectors, key = "diffmap1000hvg", assay = 'RNA')

### Plotings
png(filename = "report/BLABALBLABLA/int_obj1_diffusion_map_hvg_1000.png", width = 2400, height = 2000 )
DimPlot(int_obj1, reduction = 'diffmap1000hvg', group.by = c(tail(names(int_obj1@meta.data),3),'orig.ident'))
dev.off()

png(filename = "report/BLABALBLABLA/int_obj1_diffusion_map_hvg_1000_2.png", width = 2000, height = 1000 )
DimPlot(int_obj1, reduction = 'diffmap1000hvg', group.by = c('orig.ident','Phase'))
dev.off()

png(filename = "report/BLABALBLABLA/diffusion_map_hvgs.png", width = 2000, height = 1000 )
DimPlot(int_obj1, reduction = 'diffmap500hvg', group.by = 'orig.ident') +
  DimPlot(int_obj1, reduction = 'diffmap1000hvg', group.by = 'orig.ident')
dev.off() 


DimPlot(int_obj1, reduction = 'diffmap1000hvg', group.by = c(tail(names(int_obj1@meta.data),3),'orig.ident')
        
       
        
        
        ####################   Calculation on MNN reconstructed ASSAY ########################## 
        As MNN reconstructed contains Data slot, we can play around
        FindVariable feature can search RNA assay, 
        Then intersect the one present in MNN.reconstructed (data slot).
        --> there we can do Diff map with approximatly the TOP 500 or 1000 !!!
          
        ######################################################## PCA = FAST VERSION ###########################################################################################
        # DIFFUSION MAP on MNN reconstructed
        ElbowPlot(int_obj1)
        
        ###### Calculating diffusion map on PCA axis (FAST VERSION)
        # npcs = number of PCA componant to take --> YOU HAVE to  play around manually
        #number_of_PCA = 10
        number_of_PCA = seq(from = 2, to =10, by = 1) #  # Try low number based on Elbow plot - Over a certain point the increase just colapse the plot
        
        
        for (npcs in number_of_PCA){
          data4diffmap <- int_obj1[['integrated.mnn']]@cell.embeddings[,1:npcs]
          
          options(Matrix.warnDeprecatedCoerce = 0) # transform and error back in a warning to not crash the calculation !
          diffmappo <- destiny::DiffusionMap(data = data4diffmap, verbose = T)
          
          # Do not forget to assign a unique name to the diffusion reduction
          int_obj1[[paste0('diffmap',npcs)]] <- CreateDimReducObject(embeddings = 10000*diffmappo@eigenvectors, key ='diffmap', assay = 'integrated.mnn') #BE CAREFULL to have a unique name form
          # 10000 is just to adapt the scale to cerebro
        }
        
        DimPlot(int_obj1, reduction = 'diffmap10', group.by = c('orig.ident',"Phase"))
        DimPlot(int_obj1, reduction = paste('diffmap',npcs), group.by = c(names(tail(int_obj1@meta.data,3)),'orig.ident'))
        
        # Printing out each plot from the loop in 1 file
        embeds <- c(names(tail(int_obj1@reductions, length(number_of_PCA))))
        p <- NULL
        for (i in embeds){
          p[[i]] <- DimPlot(int_obj1, reduction = i, group.by = 'orig.ident')
        }
        
        png(filename = "report/BLABALBLABLA/int_obj1_diffmap_PCA.png", width = 2400, height = 2000 )
        CombinePlots(plots = p, ncol = 3 )
        dev.off() 
        
    
        ###################################################### High Variable Gene = LONG VERSION ##########################################################################################
        ###### Calculating diffusion map on the most High Variable Features (LONG BUT BETTER VERSION)
                                Calculation on MNN corrected assays !!!!)

        
        ###### Calculating diffusion map on the 500 most variable features (LONG BUT BETTER VERSION)
        # calculating and taking the 500 variable genes 
        int_obj1 <- FindVariableFeatures(int_obj1, nfeatures = 500) # 500genes 3 min for 13 000 Cells / 1000 genes 9 mins / 2000 genes 36min 
        top500hvg <- VariableFeatures(int_obj1)
        # Selecting the intersection of top500 in mnn.reconstructed data slot
        mnn.reconstructed.hvg <- rownames(int_obj1@assays$mnn.reconstructed)
        top500intersect <- intersect(mnn.reconstructed.hvg, top500hvg)
        
        data4diffmap <- t(GetAssayData(int_obj1,assay = "mnn.reconstructed")[top500intersect,]) # transpose / GetAssayData gather the matrix normalised data (Gene X Cell)
        # and retreive only 500 row 
        
        # Caclulating the Diffusion map
        options(Matrix.warnDeprecatedCoerce = 0)# transform and error back in a warning to not crash the calculation !
        diffmappo <- destiny::DiffusionMap(data = as.matrix(data4diffmap), verbose = T)
        int_obj1[['diffmap500hvg_mnn']] <- CreateDimReducObject(embeddings = 10000*diffmappo@eigenvectors, key = "diffmap500hvg_mnn", assay = 'integrated.mnn')
        
        #DimPlot(int_obj1, reduction = 'diffmap500hvg', group.by = c(tail(parameters,2),'orig.ident'))
        #DimPlot(int_obj1, reduction = 'diffmap500hvg', group.by = c(tail(names(int_obj1@meta.data),2),'orig.ident'))
        
        png(filename = "report/BLABALBLABLA/int_obj1_diffusion_map_hvg_500_mnn.png", width = 1600, height = 1200 )
        DimPlot(int_obj1, reduction = 'diffmap500hvg_mnn', group.by ='orig.ident')
        dev.off()
        
    
        ################# calculating and taking the 1000 variable genes ################
        int_obj1 <- FindVariableFeatures(int_obj1, nfeatures = 1000) # 500genes 3 min for 13 000 Cells / 1000 genes 9 mins / 2000 genes 36min 
        top1000hvg <- VariableFeatures(int_obj1)
        # Selecting the intersection of top1000 in mnn.reconstructed data slot
        mnn.reconstructed.hvg <- rownames(int_obj1@assays$mnn.reconstructed)
        top1000intersect <- intersect(mnn.reconstructed.hvg, top1000hvg)
        
        data4diffmap <- t(GetAssayData(int_obj1,assay = "mnn.reconstructed")[top1000intersect,]) # transpose / GetAssayData gather the matrix normalised data (Gene X Cell)
        # and retreive only 1000 row 
        # Caclulating the Diffusion map
        options(Matrix.warnDeprecatedCoerce = 0)# transform and error back in a warning to not crash the calculation !
        diffmappo <- destiny::DiffusionMap(data = as.matrix(data4diffmap), verbose = T)
        int_obj1[['diffmap1000hvg_mnn']] <- CreateDimReducObject(embeddings = 10000*diffmappo@eigenvectors, key = "diffmap1000hvgmnn", assay = 'integrated.mnn')
        
        
        
        ### Plotings
        png(filename = "report/BLABALBLABLA/int_obj1_diffusion_map_hvg_1000_mnn.png", width = 2400, height = 2000 )
        DimPlot(int_obj1, reduction = 'diffmap1000hvg_mnn', group.by = c(tail(names(int_obj1@meta.data),3),'orig.ident'))
        dev.off()
        
        png(filename = "report/BLABALBLABLA/int_obj1_diffusion_map_hvg_1000_mnn_2.png", width = 2000, height = 1000 )
        DimPlot(int_obj1, reduction = 'diffmap1000hvg_mnn', group.by = c('orig.ident','Phase'))
        dev.off()
        
        png(filename = "report/BLABALBLABLA/diffusion_map_hvgs_mnn.png", width = 2000, height = 1000 )
        DimPlot(int_obj1, reduction = 'diffmap500hvg_mnn', group.by = 'orig.ident') +
          DimPlot(int_obj1, reduction = 'diffmap1000hvg_mnn', group.by = 'orig.ident')
        dev.off() 