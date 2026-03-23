library(GSVA)            

# Loading the TPM Matrix
TPM_TN_TP_symbol <- readRDS("/Users/administrateur/CerebroData/ChIPSeq_Yap_ActiveMotiv/data/working_data/TPM_for_GSVA/TPM_TN_TP_symbol.RDS")

##### Gene names format needs to be in the same format (Symbol)
gsva_TN_TP <- GSVA::gsva(gsvaParam( TPM_TN_TP_symbol,
                                    list(r_Net_specifi_prom_genes))) # Named Vector list of you gene names 


plot <- reshape2::melt(gsva_TN_TP)

ggplot(plot, aes(x=Var2, y=value, fill=Var2)) + 
  geom_boxplot(width=0.4)+
  #geom_point(position = position_jitterdodge(jitter.width = 0.0, dodge.width = 0.8))+
  ggtitle('TN TP') +  
  xlab("") +
  ylab("GSVA value")+
  theme_classic() +
  #ylim(manual_range)+ 
  theme(legend.position = "none") + theme(axis.text = element_text(size = 12))+scale_x_discrete(guide = guide_axis(n.dodge = 2))

