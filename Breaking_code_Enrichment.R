#### ENRICHMENT !!!
https://yulab-smu.top/biomedical-knowledge-mining-book/index.html
The enrichplot package implements several visualization methods to help interpreting enrichment results.

# Nice detailed tuto here:
  https://www.melbournebioinformatics.org.au/tutorials/tutorials/seurat-go/seurat-go/

#enrichplot implements several visualization methods
supports visualizing enrichment results obtained from DOSE (Yu et al. 2015), clusterProfiler (Yu et al. 2012; Wu et al. 2021), ReactomePA (Yu and He 2016) and meshes (Yu 2018). Both over representation analysis (ORA) and gene set enrichment analysis (GSEA) are supported.
BiocManager::install("enrichplot")

# Disease Ontology Semantic and Enrichment analysis
BiocManager::install("DOSE")

# Gene Ontology only and can be done off line
BiocManager::install("clusterProfiler")

library(clusterProfiler)
clusterProfiler::enrichGO()
--> can plot the dot plot of enrichment !!
  ##### EnrichR
  --> online picking all response from enricheR as a table1
# super comprehensive and need connection
install.packages("enrichR")


###### Converting to another format
Blabla  <- mapIds(org.Mm.eg.db,
                  keys = df$ENSEMBL, # What to transform in your data
                  keytype = "ENSEMBL", # Which input format
                  column = "ENTREZID", # Which output format (searching in Database)
                  multiVals = "first") %>%  # What to do if muliple names
                  discard(is.na)

ego <- clusterProfiler::enrichGO(gene          = df$gene, # the signature to test --> SHOULD BE ENTREZID format
                                 universe      = names(geneList), # Whole list of total potential genes detected in your experiment
                                 OrgDb         = org.Mm.eg.db,
                                 ont           = "CC", # CC Cell Compartement, MF molecular Function, BP biological Process
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 0.01,
                                 qvalueCutoff  = 0.05,
                                 readable      = TRUE)

#
ego <- clusterProfiler::enrichKEGG(gene          = gene, # the signature to test
                                   organism      = "hsa"
                                   keyType       = "kegg",
                                   universe      = names(geneList), # Whole list of total potential genes detected in your experiment
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.01,
                                   qvalueCutoff  = 0.05,
                                   readable      = TRUE)

# GENE LIST SHOULD BE SORTED FROM HIGH TO LOW !!! 
gene is a vector of logfold_change sorted + names !
  http://yulab-smu.top/biomedical-knowledge-mining-book/universal-api.html
TERM2GENE : correspondance of GO term correspond to which Genesets

ego <- clusterProfiler::GSEA(geneList         = gene, # the signature to test
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.01,
                             qvalueCutoff  = 0.05,
                             readable      = TRUE)


# MSIGDBR can queries DB of geneset of GSEA.
https://www.gsea-msigdb.org/gsea/msigdb
generally we use H (Hallmark cancer) / C2 (curated gene sets) / C5 (ontology gene sets)

install.packages('msigdbdf', repos = 'https://igordot.r-universe.dev')
library(msigdbr)
# Fetching DB
genesets <- msigdbr(species = "mouse", collection = "C2", subcollection = "CGP")
genesets




################# Basic Cerebro tricks #############
int_obj1 <- cerebroApp::getEnrichedPathways(object = int_obj1, 
                                            marker_genes_input = int_obj1@misc$marker_genes$cerebro_seurat, # Name of list of markers gene table . Default is taking the last markers calculated to the object ?? maybe..
                                            databases = c("GO_Biological_Process_2018", "GO_Cellular_Component_2018","GO_Molecular_Function_2018", "KEGG_2016", "WikiPathways_2016", "Reactome_2016","Panther_2016", "Human_Gene_Atlas", "Mouse_Gene_Atlas"),
                                            adj_p_cutoff = 0.05,
                                            max_terms = 100,
                                            URL_API = "http://maayanlab.cloud/Enrichr")



int_obj1 <- cerebroApp::getEnrichedPathways(object = int_obj1)

names(int_obj1@misc$marker_genes$cerebro_seurat)


