###################################################################################################
                                #### Genomic Ranges Objects ####
###################################################################################################
https://bioconductor.org/packages/devel/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.html
BiocManager::install("GenomicRanges")
library(GenomicRanges)

## Basic Accessors for GRanges Objects
| Function             | Description                      |
  | -------------------- | -------------------------------- |
  | seqnames(gr)       | Chromosome names                 |
  | ranges(gr)         | Start and end as IRanges       |
  | start(gr), end(gr) | Access or modify start/end       |
  | strand(gr)         | Strand (+, -, or *)        |
  | width(gr)          | Width of each range              |
  | mcols(gr)          | Metadata columns (extra columns) |
  | length(gr)         | Number of ranges                 |
  | gr[i]              | Subset like a vector             |
  
ex : 
  mcols(gr)$score <- 
# Subseting
  gr[gr$gene == "A"]
  subset(gr, gene =='A')




#### Common Operations on "ranges/peaks/intervals"
| Function             | Purpose                                         |
  | -------------------- | ----------------------------------------------- |
  | subset(gr, ...)    | Filter based on metadata or coordinates         |
  | resize(gr, width)  | Resize ranges to new width (from start/mid/end) |
  | flank(gr, width)   | Get flanking regions                            |
  | promoters(gr, ...) | Get promoter regions (upstream of TSS)          |
  | trim(gr)           | Trim to avoid negative positions                |
 
#### Set-Like & Genomic Operations
  | Function                     | Description                                  |
  | ---------------------------- | -------------------------------------------- |
  | findOverlaps(gr1, gr2)     | Find overlapping ranges                      |
  | countOverlaps(gr1, gr2)    | Count overlaps per range                     |
  | subsetByOverlaps(gr1, gr2) | Keep ranges in gr1 overlapping gr2       |
  | intersect(gr1, gr2)        | Intersection (shared regions)                |
  | union(gr1, gr2)            | Union of two GRanges                         |
  | setdiff(gr1, gr2)          | Ranges in gr1 not in gr2                 |
  | reduce(gr)                 | Merge overlapping/adjacent ranges            |
  | disjoin(gr)                | Split overlapping ranges into disjoint parts |

  
##### Exporting to BED files
  library(rtracklayer)
export(gr, "output.bed")
gr2 <- import("input.bed")




######### GRangesList: Groups of Genomic Ranges ##################
## 📦 **What is a GRangesList?**

A GRangesList is a **list of GRanges objects**, used to group ranges—for example:

* Exons grouped by transcript
* Peaks grouped by sample
* Regions grouped by condition or feature

---

## 🧱 **Creating a GRangesList**
library(GenomicRanges)

# From multiple GRanges
gr1 <- GRanges("chr1", IRanges(1, 100))
gr2 <- GRanges("chr2", IRanges(200, 300))
grl <- GRangesList("txA" = gr1, "txB" = gr2)

# From a data frame with grouping
df <- data.frame(chr = c("chr1", "chr1", "chr2"),
                 start = c(1, 100, 200),
                 end = c(50, 150, 250),
                 group = c("tx1", "tx1", "tx2"))
grl <- split(makeGRangesFromDataFrame(df), df$group)


## 🧰 **Key Accessors and Operations**

| Function            | Description                               |
| ------------------- | ----------------------------------------- |
| names(grl)        | Get or set names of each group            |
| length(grl)       | Number of groups                          |
| elementNROWS(grl) | Number of elements (ranges) in each group |
| grl[[1]]          | Access individual GRanges object        |
| unlist(grl)       | Flatten into a single GRanges           |
| grl[["tx1"]]      | Access by name                            |
| mcols(grl)        | Access or assign metadata for groups      |

---

## 🔁 **Apply-like Operations**

Use lapply, end, start, etc. across groups:


start(grl)               # List of start positions
width(grl)               # List of widths
lapply(grl, reduce)      # Reduce each group


---

## 🔁 **Functional Utilities**

| Function                   | Purpose                                    |
| -------------------------- | ------------------------------------------ |
| endoapply(grl, fun)      | Like lapply, but returns a GRangesList |
| unlist(grl)              | Flatten to one GRanges with group metadata |
| reduce(grl)              | Reduce within each group                   |
| disjoin(grl)             | Disjoin each group                         |
| subsetByOverlaps(grl, x) | Subset within each group                   |


# GRangesList
grl <- GRangesList("ctrl" = ctrl_gr5, "ntr" = ntr_gr5)
gr <- unlist(grl)

# Find overlaps between the gr and the grl 
findOverlaps(gr,grl)
countOverlaps(gr,grl)
subsetByOverlaps(gr,grl)

## 🧬 **Set-Like and Overlap Functions**
Work similarly to GRanges, often after flattening:

findOverlaps(unlist(grl), gr_other)
subsetByOverlaps(grl, query)


#To keep group structure:
endoapply(grl, function(gr) subsetByOverlaps(gr, query))


---

## 📤 **Export/Import**

For export, you’ll often unlist and include a group column:

 
gr <- unlist(grl)
gr$group <- rep(names(grl), elementNROWS(grl))
rtracklayer::export(gr, "grouped_regions.bed")
 

---

## 🧪 Example Use Case: Exons by Transcript

 
library(GenomicFeatures)
txdb <- makeTxDbFromGFF("genes.gtf")
exons_by_tx <- exonsBy(txdb, by = "tx", use.names = TRUE)  # returns GRangesList