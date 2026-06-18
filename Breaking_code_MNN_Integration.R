###### Importing datasets
### rt_obj comming from RT_in_vivo_continuous/RT_in_vivo_continuous_CC_diff_reg_2.Rmd
rt_obj1 <- readRDS("/Users/administrateur/Desktop/Bio_info/Single_Cell/RT_in_vivo_continuous/data/Emile_data/UpdatedReproTransfoInVivoQueryVanillaNoDoubletsRegressednGenesAndMTpct.RDS")
int_obj1 <- qread('output_seurat/5_RT_IVC2_Reg_nCount_MT_CCdiff_global_integration_MNN_IN_PROGRESS.qs', nthreads = 12)


###### Merging object #####
# T16 = 13912¢ and D2 = 2082¢ !! WARING NOT UNIQUE NAME !
merged_obj <-  merge (x= T16, y=D2)
merged_obj@meta.data$batch <- "RT1"
merged_obj@meta.data$batch[merged_obj$orig.ident %in% c("D0","D2")] <- "RT2"

pre_obj <- merged_obj
pre_obj[['RNA']] <- JoinLayers(merged_obj[['RNA']])


# Calculating Scale and PCA 
# becaus MNN require PCA. SOOOOOOO Scaling will affect MNN integration.
pre_obj <- pre_obj %>% 
  NormalizeData %>% # by default takes it all
  FindVariableFeatures %>% # 2000 by default
  ScaleData(vars.to.regress = c('nCount_RNA', 'percent.mt')) %>% # vars.to.regress c('nCount_RNA','percent.mt','CC.diff') # by default use the FindVariableFeatures
  RunPCA %>% 
  RunUMAP( dims = 1:50)



###### Too disjoined, need integration ###### 
pre_obj_MNN <- pre_obj
pre_obj_MNN[['RNA']] <- split(pre_obj_MNN[["RNA"]], f = pre_obj_MNN$batch) 

pre_obj_MNN <- IntegrateLayers(object = pre_obj_MNN  , 
                               method =  FastMNNIntegration, 
                               new.reduction = "integrated.mnn", 
                               verbose = T)
                                # d=50 = dimensions to take in account
                                # k=30 = number of Neighbors to take in KNN, larger = bigger papatoid BUT splits more appart large group from each other... 



# ReRunning UMAP on integrated data...
pre_obj_MNN <- RunUMAP(pre_obj_MNN, dims = 1:50, reduction = "integrated.mnn" )



