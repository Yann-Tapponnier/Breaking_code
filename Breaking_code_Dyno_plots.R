
###### Nice plot Dyno pagkace
https://dynverse.org/
  https://dynverse.org/dyno/
  https://github.com/dynverse/dynmethods#list-of-included-methods
60 available tools (includes other people package)


# install.packages("devtools")
devtools::install_github("dynverse/dyno")

As they were a lot of fetshing from github I created a personal token
Personal token Github
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



library(dyno); library(tidyverse) 
#############Constructing the Dyno object requires 4 Elements + 1 Embeding of interest
        # 1) basic dataset
        # It requires Counts = the counts and the expression = Data = normalised counts
        dataset <- wrap_expression(counts=Matrix::t(GetAssayData(int_obj1, layer = "counts")),
                                   expression=Matrix::t(GetAssayData(int_obj1, layer = "data")))
        
        
        # 2) # Milestone_network is milestone of the trajectory of intereset.
        # HERE it is dummy table just to run plot_dimred()
        milestone_network <- tribble(
          ~from, ~to, ~length, ~directed,
          "A", "B", 1, FALSE,
          "B", "C", 2, FALSE,
          "B", "D", 1, FALSE,
        )
        #Show the table
        milestone_network
        
        # 3)
        # It is like progression accros the milestone, like 
        # Like milestone it is require and Dummy
        
        progressions <- milestone_network %>%
          sample_n(length(dataset$cell_ids), replace = TRUE, weight = length) %>%
          mutate(
            cell_id = dataset$cell_ids,
            percentage = runif(n())
          ) %>%
          select(cell_id, from, to, percentage)
        
        # Printing it out 
        progressions
        # Give you the "from where to where" and the % is the % of progression among this segment of road !
        
        #4)
        # Dummy again for the bifurcation
        divergence_regions <- tribble(
          ~divergence_id, ~milestone_id, ~is_start,
          "1", "A", TRUE,
          "1", "B", FALSE,
          "1", "C", FALSE
        )
        
        divergence_regions
        

################################################################################################
###### You need all the 4 elements to create the trajectory object and plot it correctly
trajectory <- add_trajectory(
  dataset,
  milestone_network = milestone_network,
  divergence_regions = divergence_regions,
  progressions = progressions
)

# Adding the groups/cluster for coloring
trajectory <- add_grouping(trajectory, 
                           grouping = data.frame(group_id=as.character(int_obj1$orig.ident),
                                                 cell_id=colnames(int_obj1)))

# LOAD the dr Dimension Reduction of Interest 
DR <- data.frame(cell_id=colnames(int_obj1), int_obj1[['umap']]@cell.embeddings[,1:2])
# dr[,2] <- dr[,2]*2




################################################################################################
Ploting
#png(filename = "~/Analyses/scRNAseq/Oiseau/Oiseau_small_basic.png", width = 5, height = 5, units = 'in', res=300)
plot_dimred(trajectory, 
                  dimred = DR, 
                  plot_trajectory = F, 
                  hex_cells = FALSE, 
                  color_cells = 'grouping', # feature or grouping
                  feature_oi = NULL,
                  size_cells = 3,
                  border_radius_percentage = 0.5,
                  alpha_cells = 1,
                  color_milestones = "given")+
  scale_color_manual(values = IBM_palette)+
  guides(colour = guide_legend(override.aes = list(size=10))) # larger texte for legend

#dev.off()

##### Adding hex_cells 
plot_dimred(trajectory, 
                  dimred = DR, 
                  plot_trajectory = F, 
                  hex_cells = ifelse(length(trajectory$cell_ids) > 10000, 50, FALSE), 
                  color_cells = 'grouping', # feature or grouping
                  feature_oi = NULL,
                  size_cells = 1,
                  border_radius_percentage = 0.9,
                  alpha_cells = 1,
                  color_milestones = "given")+
  scale_color_manual(values = IBM_palette)+
  guides(colour = guide_legend(override.aes = list(size=10))) # larger texte for legend




##################### Loop for feature plot #####################
Plots <- list()
genes <- c("Pou5f1", "Sox2", "Nanog", "Klf4", "Sftpc", "Lyz2" )
for (i in genes){
  Plots[[i]] <- plot_dimred(trajectory, 
                            dimred = DR, 
                            plot_trajectory = F, 
                            hex_cells = F, 
                            color_cells = 'feature', # feature or grouping
                            feature_oi = i ,
                            size_cells = 0.5,
                            border_radius_percentage = 0.6,
                            alpha_cells = 1)
  #scale_color_manual(values = IBM_palette)+
  #guides(colour = guide_legend(override.aes = list(size=10))) # larger texte for legend
}
# Printing out
png(filename = "balbala.png", width = 2400, height = 2000 )
    CombinePlots(plots = Plots, ncol = 3 )
dev.off() 



      ####### Changing the scale of the colors using squish !!!
      library(scales) # for squish

      Plots_scaled <- list()
      for (i in genes){
        Plots_scaled[[i]] <- plot_dimred(trajectory, 
                                         dimred = DR, 
                                         plot_trajectory = F, 
                                         hex_cells = F, 
                                         color_cells = 'feature', # feature or grouping
                                         feature_oi = i ,
                                         size_cells = 0.5,
                                         border_radius_percentage = 0.6,
                                         alpha_cells = 1)+
          #scale_color_manual(values = IBM_palette)+
          #guides(colour = guide_legend(override.aes = list(size=10))) # larger texte for legend
          scale_colour_distiller(
            palette = "RdYlBu",
            limits = c(0, 1.5), # fixing the limit to plots
            oob = squish, # Reasign the scale insid the "limits"
            direction = -1  # or -1 to reverse the palette
          )
      }
      
      CombinePlots(plots = Plots, ncol = 3 )
      CombinePlots(plots = Plots_scaled, ncol = 3 )









####################################################################################################
#################### Khroma package and other for the colorblind palettes #########################

library(khroma)

bright <- color("bright")
plot_scheme(bright(7), colours = TRUE, names = TRUE, size = 0.9)
bright_6 <- bright(6)
bright_6 <- as.character(bright_6)
plot_scheme(bright_6, colours = TRUE, names = TRUE, size = 0.9)

highcontrast <- color("high contrast")
plot_scheme(highcontrast(3), colours = TRUE, names = TRUE, size = 0.9)
highcontrast_3 <- highcontrast(3)

vibrant <- color("vibrant")
plot_scheme(vibrant(7), colours = TRUE, names = TRUE, size = 0.9)
vibrant_6 <- vibrant(6)

###### IBM pallette x5
IBM_palette <- c('#648FFF', '#785EF0', '#DC267F', '#FE6100', '#FFB000', '#808080' )
IBM_palette_3 <- IBM_palette[c(1,3,5)]



