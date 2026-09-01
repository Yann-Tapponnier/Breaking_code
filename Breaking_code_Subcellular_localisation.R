########################################################
THE AIM HERE IS TO GET ACESS TO SUBCELLULAR LOCALISATION
########################################################


                          BROKEN
                          # Using HUMAN PROTEIN ATLAS
                          BiocManager::install("HPAanalyze")
                          library(HPAanalyze)
                          
                          
                          # downloadList :
                              'Normal tissue'
                              'Pathology'
                              'Subcellular location'
                              'RNA tissue'
                              'RNA cell line'
                              'RNA transcript tissue'
                              'RNA transcript cell line'
                              'all': download everything
                              'histology': same as c('Normal tissue', 'Pathology', 'Subcellular location')
                              'rna': same as c('RNA tissue', 'RNA cell line')
                              'isoform': same as c('RNA transcript tissue', 'RNA transcript cell line')
                          
                          # Charger les données de localisation subcellulaire
                          hpa_data <- hpaDownload(downloadList = "Subcellular location", version = "latest")
                          
                          
                          
                          
                          BROKEN
                          #BiocManager::install("AnnotationHub")
                          library(AnnotationHub)
                          
                          ah <- AnnotationHub()
                          
                          query(ah, "Human Protein Atlas")
                          
                          hpa_subcell <- ah[["AHXXXXX"]]  # remplace par l’ID exact affiché
                          
                          
                          
install.packages("UniprotR")
library(UniprotR)

# Obtenir la localisation subcellulaire pour un gène spécifique

GetSubcellular_location("POU5F1", directorypath = NULL)
                          



############
library(org.Hs.eg.db)
library(GO.db)


#### EVIDENCE LEVEL ####
#Code	Signification	Fiabilité
IDA	Expérience directe	⭐⭐⭐
IC	Inféré par curateur	⭐⭐
TAS	Déclaration d’auteur	⭐⭐
IEA	Inféré automatiquement	⭐
NAS	Déclaration non traçable	⭐
ISS	Similarité de séquence	⭐–⭐⭐



res <- select(
  org.Hs.eg.db,
  keys = c("POU5F1"),
  columns = c("GO", "ONTOLOGY"),
  keytype = "SYMBOL"
) %>% filter(ONTOLOGY == "CC") %>% 
  as.data.frame()

go2term <- select(
  GO.db,
  keys = res$GO,
  columns = "TERM",
  keytype = "GOID"
)

res$TERM <- go2term$TERM[match(res$GO, go2term$GOID)]





######## Final alternative ########
Downloading the Transmembrane list on Uniprot

Va sur UniProt → Advanced Search
Filtre :
  Organism: Homo sapiens (Human)
Subcellular location: contains “Cell membrane” ou contient “Transmembrane”
Exporte : "/Users/administrateur/Desktop/Bio_info/Single_Cell/Ressources/Uniprot/uniprotkb_taxonomy_id_9606_2026_02_06.tsv"
  
Membrane_Uniprot <- read.table(
  "/Users/administrateur/Desktop/Bio_info/Single_Cell/Ressources/Uniprot/uniprotkb_taxonomy_id_9606_2026_02_06.tsv",
  sep = "\t",
  header = TRUE,
  quote = "",
  fill = TRUE,
  comment.char = ""
)

library(dplyr)
Trans_membrane_only <- Membrane_Uniprot %>%
  filter(Transmembrane != "") %>%   # garde seulement les protéines transmembranaires
  distinct(Gene.Names, .keep_all = TRUE)  # garde une seule ligne par gène

Trans_membrane_only_curated <- Trans_membrane_only %>% dplyr::select(Gene.Names, Entry, Transmembrane)
write.table(Trans_membrane_only, "/Users/administrateur/Desktop/Bio_info/Single_Cell/Ressources/Uniprot/Transmembrane_protein_list.tsv", sep="\t", row.names = F) #9631
write.table(Trans_membrane_only_curated, "/Users/administrateur/Desktop/Bio_info/Single_Cell/Ressources/Uniprot/Transmembrane_protein_list_curated.tsv", sep="\t", row.names = F)

Cell_membrane_only <- Membrane_Uniprot %>% 
        filter(grepl("Cell membrane", Subcellular.location..CC.)) %>%
        distinct(Gene.Names, .keep_all = TRUE) # garde une seule ligne par gène --# 5821
write.table(Cell_membrane_only, "/Users/administrateur/Desktop/Bio_info/Single_Cell/Ressources/Uniprot/Cell_membrane_only.tsv", sep="\t", row.names = F)

