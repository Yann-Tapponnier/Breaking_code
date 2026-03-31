library(hgu133a.db)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

Blabla  <- mapIds(hgu133a.db,
                   keys = df$ENSEMBL, # What to transform in your data
                   keytype = "ENSEMBL", # Which input format
                   column = "PROBEID", # Which output format (searching in Database)
                   multiVals = "first") %>%  # What to do if muliple names
                      discard(is.na)

### Types of Annotatino to use
columns(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
"ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"     "GENETYPE"     "GO"           "GOALL"       
"IPI"          "MAP"          "OMIM"         "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"         "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"       "UNIPROT"     



