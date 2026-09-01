# Example: subsetting a Seurat object using built-in downsampling

set.seed(123)

pbmc.sub <- subset(
  x = pbmc,
  cells = cells_to_keep                                                        # names of the cells 
  subset    = nFeature_RNA > 300 & cluster %in% optics_dbscan500_cl0.8,        # LOGICAL EXPRESSION !! QC filtering (metadata) OR in the group you define
  idents    = c("D0", "D2", "D4")                                              # select clusters / identities from the active IDENTS 
  features  = c("POU5F1", "SOX2", "KLF4", "MYC"),                              # keep variable genes only
  downsample = 500                                                             # max cells per identity using active ident
  invert = FALSE                                                               # include or reject the selection
)