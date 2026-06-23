Downloading the total DB in TSV from :
  https://panglaodb.se/markers.html?cell_type=%27choose%27

PanglaoDB_markers_27_Mar_2020 <- read_delim("~/Desktop/Bio_info/Single_Cell/Ressources/PanglaoDB/PanglaoDB_markers_27_Mar_2020.tsv", delim = "\t", escape_double = FALSE,  trim_ws = TRUE)
table(PanglaoDB_markers_27_Mar_2020$`organ`)

# Selecting the organ we want 
panglaoDB_brain_Hs <- PanglaoDB_markers_27_Mar_2020 %>% filter(organ == "Brain") %>% filter (grepl('Hs', species))
panglaoDB_brain_Mm <- PanglaoDB_markers_27_Mar_2020 %>% filter(organ == "Brain") %>% filter (grepl('Mm', species))

panglaoDB_lung_Hs <- PanglaoDB_markers_27_Mar_2020 %>% filter(organ == "Lungs") %>% filter(grepl('Hs',species))
panglaoDB_lung_Mm <- PanglaoDB_markers_27_Mar_2020 %>% filter(organ == "Lungs") %>% filter(grepl('Mm',species))

panglaoDB_immune_Hs <- PanglaoDB_markers_27_Mar_2020 %>% filter(organ == "Immune system") %>% filter(grepl('Hs',species))
panglaoDB_immune_Mm <- PanglaoDB_markers_27_Mar_2020 %>% filter(organ == "Immune system") %>% filter(grepl('Mm',species))

### Needs to be a named list
brain_Hs <- split (panglaoDB_brain_Hs$`official gene symbol`,panglaoDB_brain_Hs$`cell type`)
names(brain_Hs) <- paste0("PanglaoDB_",names(brain_Hs))
  qsave (brain_Hs, "~/Desktop/Bio_info/Single_Cell/Ressources/PanglaoDB/Signatures_list_PanglaoDB_Brain_Hs.qs")

brain_Mm <- split (panglaoDB_brain_Mm$`official gene symbol`,panglaoDB_brain_Mm$`cell type`)
names(brain_Mm) <- paste0("PanglaoDB_",names(brain_Mm))
brain_Mm <- lapply(brain_Mm, stringr::str_to_title)# I put the gene name in title case to match with my data
  qsave (brain_Mm, "~/Desktop/Bio_info/Single_Cell/Ressources/PanglaoDB/Signatures_list_PanglaoDB_Brain_Mm.qs")
  

lung_Hs <- split (panglaoDB_lung_Hs$`official gene symbol`,panglaoDB_lung_Hs$`cell type`)
names(lung_Hs) <- paste0("PanglaoDB_",names(lung_Hs))
  qsave (lung_Hs, "~/Desktop/Bio_info/Single_Cell/Ressources/PanglaoDB/Signatures_list_PanglaoDB_Lungs_Hs.qs")

lung_Mm <- split (panglaoDB_lung_Mm$`official gene symbol`,panglaoDB_lung_Mm$`cell type`)
names(lung_Mm) <- paste0("PanglaoDB_",names(lung_Mm))
lung_Mm <- lapply(lung_Mm, stringr::str_to_title)# I put the gene name in title case to match with my data
  qsave (lung_Mm, "~/Desktop/Bio_info/Single_Cell/Ressources/PanglaoDB/Signatures_list_PanglaoDB_Lungs_Mm.qs")
  
immune_Hs <- split (panglaoDB_immune_Hs$`official gene symbol`,panglaoDB_immune_Hs$`cell type`)
names(immune_Hs) <- paste0("PanglaoDB_",names(immune_Hs))
  qsave (immune_Hs, "~/Desktop/Bio_info/Single_Cell/Ressources/PanglaoDB/Signatures_list_PanglaoDB_immune_Hs.qs")
  
immune_Mm <- split (panglaoDB_immune_Mm$`official gene symbol`,panglaoDB_immune_Mm$`cell type`)
names(immune_Mm) <- paste0("PanglaoDB_",names(immune_Mm))
immune_Mm <- lapply(immune_Mm, stringr::str_to_title)# I put the gene name in title case to match with my data
  qsave (immune_Mm, "~/Desktop/Bio_info/Single_Cell/Ressources/PanglaoDB/Signatures_list_PanglaoDB_immune_Mm.qs")

# Reimporting it
panglaoDB_lungs <- readRDS("~/Desktop/Bio_info/Single_Cell/Ressources/PanglaoDB/Signatures_list_PanglaoDB_lung.RDS")



# Using it
seu_sc_f <- AddSignatures2CerebroAlt(obj = seu_sc_f,
                                     signatures = panglaoDB_brain,
                                     sourceAssay = 'RNA',
                                     mouse = T , 
                                     calculateUCell = F