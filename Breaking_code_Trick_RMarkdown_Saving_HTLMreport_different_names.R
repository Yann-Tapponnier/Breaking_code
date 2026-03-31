

##############################################################################################################
################################### TRICK RMarkdown to save the html output with DIFFERENT NAMES ###################################
############################################################################################################## 
# TRICK to save the html output with different names
--> RUNIT IN THE CONSOLE
rmarkdown::render("a2_NBEI_Integration_SC_ATAC.Rmd",
                  output_file = paste0("a2_NBEI_Integration_SC_ATAC_",sample_name,".html"))

sample_name <- 'E0_PBS_ATAC'
Min_nCount_peaks = 10000 
Max_nCount_peaks = 100000
Min_TSS.enrichment = 3.5
Max_blacklist_ratio = 1.03
Max_nucleosome_signal = 3
Min_pct_reads_in_peaks = 30

sample_name <- 'E6_DOX_ATAC'
Min_nCount_peaks = 800 
Max_nCount_peaks = 30000
Min_TSS.enrichment = 3
Max_blacklist_ratio = 0.01
Max_nucleosome_signal = 3
Min_pct_reads_in_peaks = 40