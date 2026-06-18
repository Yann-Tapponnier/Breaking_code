######################################################################################
##################################### MSiG DB ########################################
######################################################################################
library(tidyverse)

# Package vignette
https://igordot.github.io/msigdbr/
# Broad Institute page
https://www.gsea-msigdb.org/gsea/msigdb


# install.packages("msigdbr")
library(msigdbr)

# Checking what the is inside
msigdbr_collections()

print(n = 25, msigdbr_collections())
gs_collection gs_subcollection  gs_collection_name                     num_genesets
1 C1            ""                "Positional"                                    302
2 C2            "CGP"             "Chemical and Genetic Perturbations"           3538
3 C2            "CP"              "Canonical Pathways"                             19
4 C2            "CP:BIOCARTA"     "BioCarta Pathways"                             292
5 C2            "CP:KEGG_LEGACY"  "KEGG Legacy Pathways"                          186
6 C2            "CP:KEGG_MEDICUS" "KEGG Medicus Pathways"                         658
7 C2            "CP:PID"          "PID Pathways"                                  196
8 C2            "CP:REACTOME"     "Reactome Pathways"                            1787
9 C2            "CP:WIKIPATHWAYS" "WikiPathways"                                  885
10 C3            "MIR:MIRDB"       "miRDB"                                        2377
11 C3            "MIR:MIR_LEGACY"  "MIR_Legacy"                                    221
12 C3            "TFT:GTRD"        "GTRD"                                          505
13 C3            "TFT:TFT_LEGACY"  "TFT_Legacy"                                    610
14 C4            "3CA"             "Curated Cancer Cell Atlas gene sets "          148
15 C4            "CGN"             "Cancer Gene Neighborhoods"                     427
16 C4            "CM"              "Cancer Modules"                                431
17 C5            "GO:BP"           "GO Biological Process"                        7583
18 C5            "GO:CC"           "GO Cellular Component"                        1042
19 C5            "GO:MF"           "GO Molecular Function"                        1855
20 C5            "HPO"             "Human Phenotype Ontology"                     5748
21 C6            ""                "Oncogenic Signature"                           189
22 C7            "IMMUNESIGDB"     "ImmuneSigDB"                                  4872
23 C7            "VAX"             "HIPC Vaccine Response"                         347
24 C8            ""                "Cell Type Signature"                           866
25 H             ""                "Hallmark"                                       50



[1]"HALLMARK_ADIPOGENESIS"                      "HALLMARK_ALLOGRAFT_REJECTION"               "HALLMARK_ANDROGEN_RESPONSE"                
[4] "HALLMARK_ANGIOGENESIS"                      "HALLMARK_APICAL_JUNCTION"                   "HALLMARK_APICAL_SURFACE"                   
[7] "HALLMARK_APOPTOSIS"                         "HALLMARK_BILE_ACID_METABOLISM"              "HALLMARK_CHOLESTEROL_HOMEOSTASIS"          
[10] "HALLMARK_COAGULATION"                       "HALLMARK_COMPLEMENT"                        "HALLMARK_DNA_REPAIR"                       
[13] "HALLMARK_E2F_TARGETS"                       "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" "HALLMARK_ESTROGEN_RESPONSE_EARLY"          
[16] "HALLMARK_ESTROGEN_RESPONSE_LATE"            "HALLMARK_FATTY_ACID_METABOLISM"             "HALLMARK_G2M_CHECKPOINT"                   
[19] "HALLMARK_GLYCOLYSIS"                        "HALLMARK_HEDGEHOG_SIGNALING"                "HALLMARK_HEME_METABOLISM"                  
[22] "HALLMARK_HYPOXIA"                           "HALLMARK_IL2_STAT5_SIGNALING"               "HALLMARK_IL6_JAK_STAT3_SIGNALING"          
[25] "HALLMARK_INFLAMMATORY_RESPONSE"             "HALLMARK_INTERFERON_ALPHA_RESPONSE"         "HALLMARK_INTERFERON_GAMMA_RESPONSE"        
[28] "HALLMARK_KRAS_SIGNALING_DN"                 "HALLMARK_KRAS_SIGNALING_UP"                 "HALLMARK_MITOTIC_SPINDLE"                  
[31] "HALLMARK_MTORC1_SIGNALING"                  "HALLMARK_MYC_TARGETS_V1"                    "HALLMARK_MYC_TARGETS_V2"                   
[34] "HALLMARK_MYOGENESIS"                        "HALLMARK_NOTCH_SIGNALING"                   "HALLMARK_OXIDATIVE_PHOSPHORYLATION"        
[37] "HALLMARK_P53_PATHWAY"                       "HALLMARK_PANCREAS_BETA_CELLS"               "HALLMARK_PEROXISOME"                       
[40] "HALLMARK_PI3K_AKT_MTOR_SIGNALING"           "HALLMARK_PROTEIN_SECRETION"                 "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"  
[43] "HALLMARK_SPERMATOGENESIS"                   "HALLMARK_TGF_BETA_SIGNALING"                "HALLMARK_TNFA_SIGNALING_VIA_NFKB"          
[46] "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"         "HALLMARK_UV_RESPONSE_DN"                    "HALLMARK_UV_RESPONSE_UP"                   
[49] "HALLMARK_WNT_BETA_CATENIN_SIGNALING"        "HALLMARK_XENOBIOTIC_METABOLISM"


####### gathering database
HALLMARKS_HS <- msigdbr(species = "human", collection = "H")

Apoptosis_HS <- HALLMARKS_HS %>% filter(gs_name == "HALLMARK_APOPTOSIS") %>%  pull(gene_symbol)


# For Mouse
genesetmH <- msigdbr(species = "mouse", collection = "H") # (gather the human MSigDB with ortholog mapping to mouse )
genesetMM <- msigdbr(species = "mouse", collection = "MH", db_species = "MM")

### Testing the difference between MH and H (converted to mouse)
apoptosisMM <- genesetMM %>% filter(gs_name == "HALLMARK_APOPTOSIS") %>%  pull(gene_symbol)
apoptosismH <- genesetmH %>% filter(gs_name == "HALLMARK_APOPTOSIS") %>%  pull(gene_symbol)

commun <- intersect(apoptosisMM,apoptosismH) # 198 over 200 are commun

Best practice = go for  msigdbr(species = "mouse", collection = "MH", db_species = "MM")

######## Use the homemade function next
AddSignatureAUCscoreSeuratObject() to calculate and plot the signature score 





                      ex of use : 
                        HALLMARKS <- c(
                          "HALLMARK_APOPTOSIS",
                          "HALLMARK_G2M_CHECKPOINT",
                          "HALLMARK_IL2_STAT5_SIGNALING",
                          "HALLMARK_IL6_JAK_STAT3_SIGNALING",
                          "HALLMARK_INFLAMMATORY_RESPONSE",
                          "HALLMARK_INTERFERON_ALPHA_RESPONSE",
                          "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                          "HALLMARK_KRAS_SIGNALING_DN",
                          "HALLMARK_KRAS_SIGNALING_UP",
                          "HALLMARK_MITOTIC_SPINDLE",
                          "HALLMARK_MYC_TARGETS_V1",
                          "HALLMARK_MYC_TARGETS_V2",
                          "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
                          "HALLMARK_TGF_BETA_SIGNALING",
                          "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
                        )
                      
                      
                      HALLMARKS <- "HALLMARK_APOPTOSIS"
                      Genesets <- list()
                      ##### Checking that all entry where written corretly #########
                      if (!all(HALLMARKS %in% msigmm$gs_name)) {
                        stop(
                          "Some HALLMARK names are incorrect:\n",
                          paste(setdiff(HALLMARKS, msigmm$gs_name), collapse = ", ")
                        )
                      } else {
                        print("All HALLMARK names are written correctly.")
                      }
                      
                      ####### Fetching the lists ########## 
                      for (HM in HALLMARKS) {
                        Genesets[[HM]] <- msigmm %>% # Genesets[[HM]] adding HM here will name te list : list(HM = pull (HM))
                          filter(gs_name == HM) %>% 
                          pull(gene_symbol) # Pull allow to not unlist the column and return a vector diretly
                      }
                      
                      ########### RUNNING ALL HALLMARKS ################## 
                      ds_obj <- AddSignatureAUCscoreSeuratObject (Seurat_Obj = ds_obj,
                                                                    Genesets = Genesets,
                                                                    Output_path = "report/8_RT_in_vivo_WO_ANSES_ds_Signatures_clean/" ,
                                                                    Plot = TRUE,
                                                                    verbose = T)
                      
                      
