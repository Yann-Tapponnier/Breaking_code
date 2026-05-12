####################################################################################
# Multiple strategies to get the final format --> List of genesets 
"BLABLA" = c("gene1", "gene2", "gene3"))
################################################################################################
############  all strategies ############

# From Column to List
# Making it a list from each colomn using "setNames"
List_Makers <- setNames(
  as.list(omalley),     # turn each column into a list element
  colnames(omalley)     # assign column names to the list
)


# From Sheets to list
# Getting sheetsNames
sheets <- getSheetNames("/Users/administrateur/Desktop/Bio_info/Data and supp data from paper/Belmonte_2025/ScienceDirect_files_30Jan2026_10-10-32.888/1-s2.0-S0092867425008530-mmc3.xlsx")

# Importing all sheets as a list 
Belmonte_all <- lapply(
  sheets,
  function(s) read.xlsx("/Users/administrateur/Desktop/Bio_info/Data and supp data from paper/Belmonte_2025/ScienceDirect_files_30Jan2026_10-10-32.888/1-s2.0-S0092867425008530-mmc3.xlsx", sheet = s)
)


# from groups to list
### Creating a list of gene sets from groups of a LONG TABLE for the AUC score calculation
result <- Polo %>%
  group_by(module) %>%
  summarise(genes = list(geneName)) %>%
  deframe()



########################## EXEMPLES ##########################

########################################################################################################################
############## From O'Malley 2013 sup table4  ##########

omalley <- read.xlsx("/Users/administrateur/Desktop/Bio_info/Data and supp data from paper/o_maley_2013/O_Malley_2013_sup table4.xlsx")[,2:6] # eliminating the first col which countains nothing.

# Making it a list from each colomn using "setNames"
List_Makers <- setNames(
  as.list(omalley),     # turn each column into a list element
  colnames(omalley)     # assign column names to the list
)

# Convert mouse symbols to uppercase first
List_Makers_upper <- lapply(List_Makers, toupper)

int_obj1 <- AddSignatureAUCscoreSeuratObject(int_obj1, 
                                             Probs = c(0,0.66),
                                             Genesets = List_Makers_upper,
                                             Output_path = "report/4_hiPAR_further_analysis/",
                                             Plot = T ,
                                             verbose = T
)



########################################################################################################################
############### PLoting mesenchymal drift siganture BELMONTE ##############
Belmonte <- read.xlsx("/Users/administrateur/Desktop/Bio_info/Data and supp data from paper/Belmonte_2025/ScienceDirect_files_30Jan2026_10-10-32.888/1-s2.0-S0092867425008530-mmc3.xlsx")

# Getting sheetsNames
sheets <- getSheetNames("/Users/administrateur/Desktop/Bio_info/Data and supp data from paper/Belmonte_2025/ScienceDirect_files_30Jan2026_10-10-32.888/1-s2.0-S0092867425008530-mmc3.xlsx")

# Importing all sheets as a list 
Belmonte_all <- lapply(
  sheets,
  function(s) read.xlsx("/Users/administrateur/Desktop/Bio_info/Data and supp data from paper/Belmonte_2025/ScienceDirect_files_30Jan2026_10-10-32.888/1-s2.0-S0092867425008530-mmc3.xlsx", sheet = s)
)

# Assigning names to the list elements based on sheet names
names(Belmonte_all) <- sheets
sheets
print(Belmonte_all[[sheets[1]]]
      Aging_signatures = Belmonte_all[[sheets[2]]],
      Fibroblast_subtype_signatures = Belmonte_all[[sheets[3]]],
      Reprog_cell_state_signatures = Belmonte_all[[sheets[4]]]
      
      
      int_obj1 <- AddSignatureAUCscoreSeuratObject(int_obj1, 
                                                   Genesets = list(Belmonte_MD_signatures = Belmonte_all$MD_signatures$MD.score,
                                                                   Belmonte_Aging_signatures_up = Belmonte_all$Aging_signatures$Age.up,
                                                                   Belmonte_Aging_signatures_down = Belmonte_all$Aging_signatures$Age.down,
                                                                   Belmonte_Reprog_cell_state_signatures_fibroblast = Belmonte_all$Reprog_cell_state_signatures$Fibroblast,
                                                                   Belmonte_Reprog_cell_state_signatures_PartialReprog = Belmonte_all$Reprog_cell_state_signatures$PartialReprog,
                                                                   Belmonte_Reprog_cell_state_signatures_EarlyPluripotency= Belmonte_all$Reprog_cell_state_signatures$EarlyPluripotency,
                                                                   Belmonte_Reprog_cell_state_signatures_Pluripotency= Belmonte_all$Reprog_cell_state_signatures$Pluripotency,
                                                                   Belmonte_Reprog_cell_state_signatures_NonReprog= Belmonte_all$Reprog_cell_state_signatures$NonReprog
                                                   ),
                                                   Output_path = "report/4_hiPAR_further_analysis/",
                                                   Plot = T ,
                                                   verbose = T)
      
      
      
      
########################################################################################################################
############# Importing signature from Liu_Polo 2020 ############# 

Polo <- read.xlsx("/Users/administrateur/Desktop/Bio_info/Data and supp data from paper/Liu_Polo_2020/MOESM3_ESM/Supplementary Table 3.xlsx")

### Creating a list of gene sets from groups of a LONG TABLE for the AUC score calculation
result <- Polo %>%
  group_by(module) %>%
  summarise(genes = list(geneName)) %>%
  deframe()

int_obj11 <- AddSignatureAUCscoreSeuratObject(int_obj11, 
                                              Genesets = result,
                                              Output_path = "report/11_hiPAR_further_analysis/Liu_Polo_signatures/",
                                              Plot = T ,
                                              verbose = T)
