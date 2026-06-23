############################################################################## 
################################### plotly ################################### 
##############################################################################
PlotlyForLassoSelection  <- function(Seurat_obj,
                                     Reduction = "umap",
                                     Color.by = "orig.ident",
                                     Threshold = 0,
                                     Colors = "viridis") {
  
          ## 🔍 QC — vérifications des arguments
          if (!inherits(Seurat_obj, "Seurat")) {
            stop("❌ Seurat_obj must be a Seurat object")
          }
          
          if (!is.character(Color.by)) {
            stop("❌ Color.by must be a single character string (metadata name)")
          }
          
          if (!Color.by %in% colnames(Seurat_obj@meta.data)) {
            stop(paste0(
              "❌ Color.by = '", Color.by,
              "' not found in Seurat object metadata"
            ))
          }
          
          if (!is.numeric(Threshold) || length(Threshold) != 1 ||
              Threshold < 0 || Threshold > 1) {
            stop("❌ Threshold must be a single numeric value between 0 and 1")
          }
  # #0 Libraris and QC
  require(plotly)
  require(Seurat)
  require(shiny)
  
  
  
  # #1 Conversion of Seurat Object
  umap_df <- as.data.frame(Embeddings(Seurat_obj, reduction = Reduction))
  umap_df$cell <- colnames(Seurat_obj)
  if (Threshold==0){
  umap_df$cluster <- (Seurat_obj[[Color.by]]) %>% unlist()
  }else{
    umap_df$cluster <- (Seurat_obj[[Color.by]]) %>% unlist() 
    thr <- quantile(umap_df$cluster, Threshold, na.rm = TRUE)
    umap_df$cluster[umap_df$cluster <= thr] <- 0
  }
  
  # #2 — Créer l’interface Shiny (UI)  
  ui <- fillPage(
    h3("UMAP interactive (lasso)"), 
    plotlyOutput("umap_plot", height = "85vh"), # Height adjusted for fillPage of the umap plot. 85% for umap 15% for selected UMIs
    actionButton("save", "Sauvegarder la sélection"),
    verbatimTextOutput("info")
  )
  
  
  #3 — Créer le serveur Shiny
  server <- function(input, output, session) {
    ## 🔁 Sélection réactive unique
    selected_cells <- reactive({
      sel <- event_data("plotly_selected", source = "umap_source")
      if (is.null(sel)) return(NULL)
      unique(sel$key)
    })
  
    ## 1️⃣ Plot UMAP
    output$umap_plot <- renderPlotly({
      plot_ly(
        data = umap_df,
        x = ~umap_1,
        y = ~umap_2,
        key = ~cell,
        type = "scatter",
        mode = "markers",
        color = ~cluster,
        colors = Colors,
        marker = list(size = 5),
        source = "umap_source"
      ) %>%
        layout(dragmode = "lasso")
    })
    
    
    ## 2️⃣ Print des cellules sélectionnées
    output$info <- renderPrint({
      cells <- selected_cells()
      if (is.null(cells)) return("Aucune sélection")
      cells
    })
    
    
    ## 3️⃣ Nombre total de cellules sélectionnées
    output$n_cells <- renderText({
      cells <- selected_cells()
      if (is.null(cells)) return("0 cellule")
      paste(length(cells), "cellules sélectionnées")
      message("UMI récupérés pour ", length(cells), " cellules")
    })
    
    observeEvent(input$save, {
      assign(
        "cells_selected",
        selected_cells(),
        envir = .GlobalEnv
      )
    })
  }

#4 — Lancer l’application Shiny
shinyApp(ui = ui, server = server)
}



                                      # #### Exemple of launch :
                                      # PlotlyForLassoSelection(
                                      #   Seurat_obj = int_obj1,
                                      #   Reduction = "umap",
                                      #   Color.by = "NK_Cl3_up_top100_FC",
                                      #   Threshold = 0.9)
                                      # 
