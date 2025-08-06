library(shiny)
library(ggplot2)
library(viridis)
library(dplyr)
library(proxy)
library(readxl)
library(shiny)
library(vegan)
library(plotly)
library(data.table)
library(tools)
library(bslib)
library(gridExtra)
if (!requireNamespace("Cardinal", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("Cardinal")
}
library(Cardinal)
library(SummarizedExperiment)
library(MALDIquant)
library(MALDIquantForeign)



options(shiny.maxRequestSize = 100*1024^2)

# -------- FUNZIONI --------

read_any_data <- function(file) {
  ext <- tolower(tools::file_ext(file$name))
  if (ext == "rds") return(readRDS(file$datapath))
  if (ext %in% c("csv", "txt", "tsv")) return(read.table(file$datapath, header = TRUE))
  if (ext %in% c("xls", "xlsx")) return(readxl::read_excel(file$datapath))
  stop("Unsupported file format:", ext)
}

fun_extractData <- function(zip_path) {
  unzip_dir <- tempfile()
  dir.create(unzip_dir)
  unzip(zip_path$datapath, exdir = unzip_dir)
  
  imzml_file <- list.files(unzip_dir, pattern = "\\.imzML$", full.names = TRUE)
  if (length(imzml_file) != 1) stop("❌ The zip archive must contain exactly one .imzML file")
  
  spectra <- importImzMl(imzml_file, centroided = TRUE)
  
  n_pixel <- length(spectra)
  coor_x <- numeric(n_pixel)
  coor_y <- numeric(n_pixel)
  intensity <- vector("list", n_pixel)
  
  for (i in seq_len(n_pixel)) {
    coor_x[i] <- spectra[[i]]@metaData$imaging$pos[[1]]
    coor_y[i] <- spectra[[i]]@metaData$imaging$pos[[2]]
    intensity[[i]] <- spectra[[i]]@intensity
  }
  
  name_mz <- round(spectra[[1]]@mass, 3)
  featureMatrix <- t(do.call(cbind, intensity))
  colnames(featureMatrix) <- name_mz
  coordinates <- data.frame(X = coor_x, Y = coor_y)
  
  list(featureMatrix = featureMatrix, coordinates = coordinates, mz = name_mz)
}

downsample_grid <- function(Data_IM, Data_XY, cell_size) {
  Data_IM <- as.data.frame(Data_IM)
  Data_IM[] <- lapply(Data_IM, as.numeric)
  
  cell_x <- floor(Data_XY$X / cell_size)
  cell_y <- floor(Data_XY$Y / cell_size)
  cell_id <- paste(cell_x, cell_y, sep = "_")
  
  df_full <- cbind(data.frame(cell_id = cell_id, X = Data_XY$X, Y = Data_XY$Y), Data_IM)
  
  summary_df <- df_full %>%
    group_by(cell_id) %>%
    summarise(
      X = mean(X), Y = mean(Y),
      across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop"
    )
  
  col_IM <- setdiff(colnames(summary_df), c("cell_id", "X", "Y"))
  IM_reduced <- as.matrix(summary_df[, col_IM])
  XY_reduced <- as.data.frame(summary_df[, c("X", "Y")])
  
  list(IM_reduced = IM_reduced, XY_reduced = XY_reduced)
}

# -------- UI --------

ui <- navbarPage(
  title        = "Clustering for Mass Spectrometry Data",
  id           = "main_nav",
  theme        = bs_theme(bootswatch = "lux"),
  windowTitle  = "MSI Clustering Suite",
  
  # ===== HOME TAB =====
  tabPanel(
    title = "Home",
    value = "home",
    fluidPage(
      h2("Welcome to the MSI Clustering Suite"),
      p("This Shiny application lets you explore and cluster imaging mass-spectrometry data, either as an IM matrix with XY coordinates or as imzML archives, in a reproducible and interactive way."),
      tags$hr(),
      h4("Quick start?"),
      tags$ol(
        tags$li("Download the example files or use your own data."),
        tags$li("Go to the Workspace tab."),
        tags$li("Choose a distance metric, a clustering algorithm and k."),
        tags$li("Press “Assign clusters” to view and explore the results.")
      ),
      tags$p(
        tags$a(href = "FinalMatrix_Peptides.RDS", download = NA, "⬇️ Download IM example (RDS)"),
        " · ",
        tags$a(href = "FinalMatrix_Coor.RDS", download = NA, "⬇️ Download XY example (RDS)")
      ),
      tags$hr(),
      h4("Developed by"),
      tags$blockquote("Giulia Capitoli & Alex Pietrarota"),
      actionLink("jump_to_analysis", "▶︎ Go to Workspace")
    )
  ),
  
  # ===== WORKSPACE TAB =====
  tabPanel(
    title = "Workspace",
    value = "analysis",
    tabsetPanel(id = "analysis",
                
      # === RDS Upload ===
      tabPanel("Upload IM + XY Files",
               fluidPage(
                 sidebarLayout(
                   sidebarPanel(
                     radioButtons("mode_rds", "Classification type:",
                                  choices = c("Show only IM clusters" = "solo_cluster",
                                              "Classify new tissue" = "tessuto")),
                     
                     actionLink("show_note",  "⬆️ Upload Information"),
                     
                     fileInput("original_IM", "Upload IM file", accept = c(".rds", ".csv", ".txt", ".xlsx")),
                     fileInput("original_XY", "Upload XY file", accept = c(".rds", ".csv", ".txt", ".xlsx")),
                     conditionalPanel("input.mode_rds == 'tessuto'",
                                      fileInput("file_new_im", "Matrice IM nuovo tessuto", accept = c(".rds", ".csv", ".txt", ".xlsx")),
                                      fileInput("file_new_xy", "Coordinate XY nuovo tessuto", accept = c(".rds", ".csv", ".txt", ".xlsx"))
                     ),
                     numericInput("cell_size", "Grid cell size (1-30)", value = 3, min = 1, max = 30),
                     
                     actionLink("show_gloss", "🌏 General Information"),
                     
                     selectInput("method", "Distance method",
                                 choices = c("angular", "bray", "canberra", "chord", "correlation",
                                             "cosine", "divergence", "eDice", "euclidean", "eJaccard", "geodesic",
                                             "hellinger", "kullback", "manhattan",
                                             "soergel", "squared_euclidean",
                                             "supremum", "wave", "whittaker"),
                                 selected = "euclidean"),
                     selectInput("algorithm", "Clustering algorithm",
                                 choices = c("Hierarchical Clustering", "K-means"),
                                 selected = "Hierarchical Clustering"),
                     conditionalPanel(
                       condition = "input.algorithm == 'Hierarchical Clustering'",
                       selectInput("hclust_method", "Hierarchical clustering method",
                                   choices = c("complete", "single", "ward.D", "ward.D2", "average"),
                                   selected = "complete")
                     ),
                     conditionalPanel(
                       condition = "input.algorithm == 'K-means' && !(input.method == 'euclidean' || input.method == 'squared_euclidean')",
                       numericInput("mds_var", "Desired explained variance (%)", min = 50, max = 99, value = 80)
                     ),
                     checkboxInput("weight_space", "Apply spatial weighting", value = FALSE),
                     conditionalPanel(
                       condition = "input.weight_space == true",
                       helpText("ℹ️ Spatial penalization adjusts molecular distances considering the physical proximity of pixels: nearby points are considered more similar, even with identical molecular profiles.")
                     ),
                     numericInput("k", "Number of clusters (2-20)", value = 5, min = 2, max = 20),
                     actionButton("run_rds", "Assign clusters")
                   ),
                   mainPanel(
                     plotlyOutput("plot_rds_orig"),
                     conditionalPanel(
                       condition = "input.mode_rds == 'tessuto'",
                       plotlyOutput("plot_rds_new"),
                       downloadButton("downloadClusters", "Download clustering report (.csv)")
                     ),
                     conditionalPanel(
                       condition = "input.algorithm == 'K-means' && !(input.method == 'euclidean' || input.method == 'squared_euclidean')",
                       plotOutput("scree_plot")
                     )
                   )
                 )
               )
      ),
      
      # === imzML Upload ===
      tabPanel("Upload imzML Files",
               fluidPage(
                 sidebarLayout(
                   sidebarPanel(
                     radioButtons("mode_zip", "Classification type:",
                                  choices = c("Show only IM clusters" = "solo_cluster", "Classify new tissue" = "tessuto")),
                     
                     actionLink("show_note_zip", "⬆️ Upload Information"), 
                     
                     fileInput("original_zip", "Original tissue (.zip)", accept = ".zip"),
                     conditionalPanel("input.mode_zip == 'tessuto'",
                                      fileInput("new_zip", "New tissue (.zip)", accept = ".zip")
                     ),
                     numericInput("cell_size", "Grid cell size (1-30)", value = 3, min = 1, max = 30),
                     
                     actionLink("show_gloss", "🌏 General Information"),
                     
                     selectInput("method", "Distance method",
                                 choices = c("angular", "braun-blanquet", "bray", "canberra", "chord", "correlation",
                                             "cosine", "divergence", "eDice", "euclidean", "eJaccard", "faith", "fager", "geodesic",
                                             "hamman", "hellinger", "jaccard", "kulczynski2", "kullback", "manhattan", "michael",
                                             "mozley", "ochiai", "russel", "simple matching", "simpson", "soergel", "squared_euclidean",
                                             "supremum", "tanimoto", "wave", "whittaker"),
                                 selected = "euclidean"),
                     selectInput("algorithm", "Clustering algorithm",
                                 choices = c("Hierarchical Clustering", "K-means"),
                                 selected = "Hierarchical Clustering"),
                     conditionalPanel(
                       condition = "input.algorithm == 'Hierarchical Clustering'",
                       selectInput("hclust_method", "Hierarchical clustering method",
                                   choices = c("complete", "single", "ward.D", "ward.D2", "average"),
                                   selected = "complete")
                     ),
                     conditionalPanel(
                       condition = "input.algorithm == 'K-means' && !(input.method == 'euclidean' || input.method == 'squared_euclidean')",
                       numericInput("mds_var", "Desired explained variance (%)", min = 50, max = 99, value = 80)
                     ),
                     checkboxInput("weight_space", "Apply spatial weighting", value = FALSE),
                     conditionalPanel(
                       condition = "input.weight_space == true",
                       helpText("ℹ️ Spatial penalization adjusts molecular distances considering the physical proximity of pixels: nearby points are considered more similar, even with identical molecular profiles.")
                     ),
                     numericInput("k", "Number of clusters (2-20)", value = 5, min = 2, max = 20),
                     actionButton("run_zip", "Assign clusters")
                   ),
                   mainPanel(
                     plotlyOutput("plot_zip_orig"),
                     conditionalPanel(
                       condition = "input.mode_zip == 'tessuto'",
                       plotlyOutput("plot_zip_new"),
                       downloadButton("downloadClusters", "Download clustering report (.csv)")
                     ),
                     conditionalPanel(
                       condition = "input.algorithm == 'K-means' && 
                 !(input.method == 'euclidean' || input.method == 'squared_euclidean')",
                       plotOutput("scree_plot_zip")
                     )
                   )
                 )
               )
      )
    )
  ),
  
  # ===== PRIVACY TAB =====
  tabPanel(
    title = "Privacy & Impact",
    value = "privacy",
    fluidPage(
      h3("Privacy and Impact"),
      p("This application is designed to help analyze and cluster mass spectrometry imaging (MSI) data. The results of this analysis are intended solely for scientific research and should not be used for other purposes without appropriate consent."),
      p("RStudio collects data in line with their Privacy Policy for the necessary functioning of their cloud products, including our hosting provider, shinyapps.io."),
      p("If you have any concerns about data privacy or would like more information, please contact us at: alex.pietrarota@unimib.it.")
    )
  )
)


# ------------   SERVER   --------------

server <- function(input, output, session) {
  
  observeEvent(input$jump_to_analysis, {
    updateNavbarPage(session, inputId = "main_nav", selected = "analysis")
  })
  
    # ---- Modale: File structure e requisiti (.rds) ----
    observeEvent(input$show_note, {
      req(input$main_nav, input$analysis)
      if (input$main_nav == "analysis" && input$analysis == "Upload IM + XY Files") {
        showModal(modalDialog(
          title = "Important Note",
          HTML(
            "<h4><b>Supported File Formats</b></h4>
          <p>The supported file formats are: <code>.rds</code>, <code>.RDS</code>, 
          <code>.csv</code>, <code>.txt</code>, <code>.xlsx</code></p>

          <h4><b>Files Structure</b></h4>
          <p>The IM file must have rows representing the pixels or cells, and columns 
          representing the m/z values (specific mass-to-charge ratios).</p>
          <p>Please ensure that the IM and XY files have the same number of rows.</p>

          <h4><b>Column Names</b></h4>
          <p>It is recommended that the IM file has column names for each m/z. 
          Additionally, the XY file must have columns named <code>'X'</code> for the 
          X-axis and <code>'Y'</code> for the Y-axis coordinates.</p>

          <h4><b>Spatial Downsampling</b></h4>
          <p>
          Downsampling can be applied 
          to reduce computation. The tissue is divided into square cells based on the selected <code>cell size</code>, and points within each cell are averaged both in spatial 
          position and intensity. The <code>cell size</code> controls the grid resolution: 
          smaller values preserve more detail, while larger values result in stronger 
          aggregation and fewer representative points.
          </p>
          <p>
          <b>Note:</b> If <code>cell size</code> is set to 1, no reduction is performed 
          and the original data are preserved, even for large datasets.
          </p>

          <h4><b>Example Datasets</b></h4>
          <p>For example datasets, refer to the Home tab.</p>"
          ),
          easyClose = TRUE,
          footer = modalButton("Close")
        ))
      }
    })
    
    # ---- Modale: Info generali sulle distanze ----
    observeEvent(input$show_gloss, {
      req(input$main_nav, input$analysis) 
      if (input$main_nav == "analysis") {
        showModal(modalDialog(
          title = "General Information",
          HTML(
            "<h4><b>Note</b></h4>
          <p><i>Consult <code>proxy::dist()</code> and <code>vegan::vegdist</code> documentation for the complete list of distance methods.</i></p><br>

          <h4><b>Non-Euclidean Distance for K-means</b></h4>
          <p>When a non-Euclidean distance method is selected for K-means, dimensionality
          reduction is applied. The user can specify the percentage of explained variance
          to be preserved during the reduction. Additionally, a screeplot will be
          displayed to visualize the cumulative explained variance.</p><br>"
          ),
          easyClose = TRUE,
          footer = modalButton("Close")
        ))
      }
    })
  
  observeEvent(input$show_note_zip, {
    if (input$main_nav == "analysis" && input$analysis == "Upload imzML Files") {
      showModal(modalDialog(
        title = "Upload Information (.zip)",
        HTML("
        <h4><b>SUPPORTED FILE FORMATS</b></h4>
        <p>Only compressed <code>.zip</code> files are accepted for uploading tissue data.</p>
        <p>Each <code>.zip</code> archive <b>must contain exactly two files</b>:</p>
        <ul>
          <li>One <code>.imzML</code> file</li>
          <li>One <code>.ibd</code> file</li>
        </ul>
        <p>The application will automatically extract molecular intensity (IM) and spatial coordinate (XY) data from the archive.</p>

        <h4><b>SPATIAL DOWNSAMPLING</b></h4>
        <p>
        Downsampling can be applied to reduce computation. The tissue is divided into square cells based on the selected <code>cell size</code>, and points within each cell are averaged both in spatial 
        position and intensity. The <code>cell size</code> controls the grid resolution: 
        smaller values preserve more detail, while larger values result in stronger 
        aggregation and fewer representative points.
        </p>
        <p>
        <b>Note:</b> If <code>cell size</code> is set to 1, no reduction is performed 
        and the original data are preserved, even for large datasets.
        </p>

        <h4><b>Example Datasets</b></h4>
        <p>For example datasets, refer to the Home tab.</p>
      "),
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
    }
  })
  
  
  
  observeEvent(input$run_rds, {
    req(input$original_IM, input$original_XY)
    
    withProgress(message = "Processing data...", value = 0, {
      
      # ---- Lettura e downsampling ----
      res1_IM <- read_any_data(input$original_IM)
      res1_XY <- read_any_data(input$original_XY)
      down1 <- downsample_grid(res1_IM, res1_XY, input$cell_size)
      original_IM <- down1$IM_reduced
      original_XY <- down1$XY_reduced
      
      if (nrow(original_IM) != nrow(original_XY)) {
        showNotification("❌ IM and XY must have the same number of rows",type = "error")
        return()
      }
      
      incProgress(0.3, detail = "Computing molecular distance...")
      
      # ---- Calcolo distanza molecolare ----
      dist_method <- input$method
      dist_matrix <- if (tolower(dist_method) == "squared_euclidean") {
        base_dist <- dist(original_IM, method = "euclidean")
        as.dist(as.matrix(base_dist)^2)
      } else if (tolower(dist_method) == "bray") {
        vegan::vegdist(original_IM, method = "bray")
      } else {
        proxy::dist(original_IM, method = dist_method)
      }
      
      # ---- Penalizzazione spaziale (opzionale) ----
      incProgress(0.5, detail = "Applying spatial weighting...")
      
      if (input$weight_space) {
        spatial_dist <- dist(original_XY, method = "euclidean")
        weight_matrix <- 1 - exp(-as.matrix(spatial_dist))
        
        dist_mat <- as.matrix(dist_matrix)
        if (!all(dim(dist_mat) == dim(weight_matrix))) {
          showNotification("❌ Matrix dimensions do not match", type = "error")
          return()
        }
        
        dist_mat <- dist_mat * weight_matrix
        dist_matrix <- as.dist(dist_mat)
      }
      
      # ---- Clustering ----
      incProgress(0.7, detail = "Clustering...")
      
      if (input$algorithm == "K-means") {
        if (input$method == "euclidean" || input$method == "squared_euclidean") {
          clusters <- kmeans(original_IM, centers = input$k, nstart = 10)$cluster
          
        } else {
          incProgress(0.75, detail = "Running MDS...")
          
          dist_mat <- proxy::dist(original_IM, method = input$method)
          mds_result <- cmdscale(dist_mat, eig = TRUE, k = nrow(original_IM) - 1)
          
          eig_vals <- mds_result$eig[mds_result$eig > 0]
          
          if (length(eig_vals) == 0) {
            showNotification("⚠️ No positive eigenvalues found in MDS. Check your data or method.", type = "error")
            return()
          }
          
          var_explained <- cumsum(eig_vals) / sum(eig_vals)
          k_mds <- which(var_explained >= input$mds_var / 100)[1]
          
          if (is.na(k_mds)) {
            showNotification("⚠️ No dimension reaches the desired explained variance.", type = "warning")
            return()
          }
          
          mds_coords <- mds_result$points[, 1:k_mds, drop = FALSE]
          clusters <- kmeans(mds_coords, centers = input$k)$cluster
          
          output$scree_plot <- renderPlot({
            plot(var_explained * 100, type = "b", pch = 19,
                 xlab = "Number of dimensions", ylab = "Cumulative explained variance (%)",
                 main = "Multidimensional Scaling", sub = "Screeplot",
                 ylim = c(0, 100), xlim = c(1, min(length(var_explained), 50)))
            abline(h = input$mds_var, col = "red", lty = 2)
          })
  
          mds_coords <- mds_result$points[, 1:k_mds]
          km <- kmeans(mds_coords, centers = input$k)
          clusters <- km$cluster
        }
      } else if (input$algorithm == "Hierarchical Clustering") {
        hc <- hclust(dist_matrix, method = input$hclust_method)
        clusters <- cutree(hc, k = input$k)
      }
      
      original_XY$cluster <- as.factor(clusters)
      
      output$plot_rds_orig <- renderPlotly({
        
        p <- ggplot(original_XY, aes(
          x = X, y = Y, color = cluster,
          text = paste0("X: ", X, "<br>Y: ", Y, "<br>Cluster: ", cluster)
        )) +
          geom_point(shape = 15, size = 1.5) +
          scale_color_viridis_d(option = "plasma") +
          labs(title = "Clusters in the original tissue", x = "X", y = "Y") +
          coord_fixed() +
          theme_minimal()
        
        ggplotly(p, tooltip = "text")
      })
      
      
      # ---- Classificazione nuovo tessuto ----
      if (input$mode_rds == "tessuto") {
        incProgress(0.8, detail = "Classifying new tissue...")
        req(input$file_new_im, input$file_new_xy)
        
        res2_IM <- read_any_data(input$file_new_im)
        res2_XY <- read_any_data(input$file_new_xy)
        down2 <- downsample_grid(res2_IM, res2_XY, input$cell_size)
        new_IM <- down2$IM_reduced
        new_XY <- down2$XY_reduced
        
        colnames(original_IM) <- round(as.numeric(colnames(original_IM)), 1)
        colnames(new_IM) <- round(as.numeric(colnames(new_IM)), 1)
        common_mz <- intersect(colnames(original_IM), colnames(new_IM))
        
        if (length(common_mz) == 0) {
          new_XY$cluster_plot <- "unclassified"
          output$plot_rds_new <- renderPlotly({
            new_XY$cluster_plot <- as.factor("unclassified")
            
            p <- ggplot(new_XY, aes(
              x = X, y = Y, color = cluster_plot,
              text = paste0("X: ", X, "<br>Y: ", Y, "<br>Status: unclassified")
            )) +
              geom_point(size = 1.5, shape = 15) +
              scale_color_manual(name = "cluster", values = c("unclassified" = "grey80")) +
              labs(
                title = "Clusters in the new tissue",
                subtitle = "No matching m/z values",
                x = "X", y = "Y"
              ) +
              coord_fixed() +
              theme_minimal()
            
            ggplotly(p, tooltip = "text")
          })
          showNotification("⚠️ No common m/z values between the two tissues. Classification not performed.")
          return()
        }
        
        valid_idx <- complete.cases(new_IM[, common_mz])
        classified_IM <- new_IM[valid_idx, common_mz, drop = FALSE]
        classified_XY <- new_XY[valid_idx, , drop = FALSE]
        unclassified_XY <- new_XY[!valid_idx, , drop = FALSE]
        
        method_for_proxy <- if (tolower(input$method) == "squared_euclidean") "euclidean" else input$method
        dist_to_new <- proxy::dist(classified_IM, original_IM[, common_mz], method = method_for_proxy)
        assigned <- apply(as.matrix(dist_to_new), 1, function(row) clusters[which.min(row)])
        classified_XY$cluster_plot <- as.factor(assigned)
        
        if (nrow(unclassified_XY) > 0) {
          unclassified_XY$cluster_plot <- "non_classificato"
        }
        
        combined_df <- rbind(classified_XY, unclassified_XY)
        
        output$plot_rds_new <- renderPlotly({
          combined_df$cluster_plot <- as.factor(combined_df$cluster_plot)
          p <- ggplot(combined_df, aes(
            x = X, y = Y,
            text = paste0("X: ", X, "<br>Y: ", Y, "<br>Cluster: ", cluster_plot)
          )) +
            geom_point(aes(color = cluster_plot), shape = 15, size = 1.5, show.legend = TRUE) +
            scale_color_manual(
              name = "cluster",
              values = viridis::viridis(input$k, option = "plasma"),
              drop = FALSE
            ) +
            labs(
              title = "Clusters in the new tissue",
              x = "X", y = "Y"
            ) +
            coord_fixed() +
            theme_minimal()
          
          # Conversione in Plotly
          ggplotly(p, tooltip = "text")
        })
      }
      incProgress(1, detail = "Completed!")
      showNotification("Classification completed ✅", type = "message", duration = 8)
    })
    output$downloadClusters <- downloadHandler(
      filename = function() {
        paste0("clusters_", input$mode_rds, "_", Sys.Date(), ".csv")
      },
      content = function(file) {
        # Check che i cluster esistano
        if (!exists("assigned") || is.null(assigned)) {
          showNotification("❌ No cluster assignment found.", type = "error")
          return(NULL)
        }
        
        # 1. Dati clustering
        cluster_df <- data.frame(Pixel = seq_along(assigned), Cluster = assigned)
        
        # Riga vuota per separazione
        empty_row <- data.frame(Pixel = "", Cluster = "")
        
        # 2. Report descrittivo
        report_df <- data.frame(
          Pixel = c("=== CLUSTERING REPORT ===",
                    "Input type:",
                    "Cell size",
                    "Clustering method:",
                    "Distance method:",
                    "Spatial penalization:",
                    "Number of clusters:",
                    if (input$algorithm == "K-means" && !(input$method %in% c("euclidean", "squared_euclidean"))) "MDS dimension:" else NULL,
                    "Date:"),
          Cluster = c("",
                      input$mode_rds,
                      input$cell_size,
                      input$algorithm,
                      input$method,
                      ifelse(input$weight_space, "Yes", "No"),
                      input$k,
                      if (input$algorithm == "K-means" && !(input$method %in% c("euclidean", "squared_euclidean"))) input$mds_var else NULL,
                      format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
          stringsAsFactors = FALSE
        )
        
        # 3. Salva su file
        final_df <- rbind(cluster_df, empty_row, report_df)
        write.csv(final_df, file, row.names = FALSE)
      },
      contentType = "text/csv"
    )
  })




  # ----------------- ZIP -----------------
  
  observeEvent(input$run_zip, {
    withProgress(message = "Processing data...", value = 0, {
      
      req(input$original_zip)
      
      res <- fun_extractData(input$original_zip)
      down <- downsample_grid(res$featureMatrix, res$coordinates, input$cell_size)
      original_IM <- down$IM_reduced
      original_XY <- down$XY_reduced
      
      if (nrow(original_IM) != nrow(original_XY)) {
        showNotification("❌ IM e XY devono avere lo stesso numero di righe", type = "error")
        return()
      }
      
      incProgress(0.3, detail = "Computing molecular distance...")
      
      dist_method <- input$method
      dist_matrix <- if (tolower(dist_method) == "squared_euclidean") {
        base_dist <- dist(original_IM, method = "euclidean")
        as.dist(as.matrix(base_dist)^2)
      } else if (tolower(dist_method) == "bray") {
        vegan::vegdist(original_IM, method = "bray")
      } else {
        proxy::dist(original_IM, method = dist_method)
      }
      
      incProgress(0.5, detail = "Applying spatial weighting...")
      
      if (input$weight_space) {
        spatial_dist <- dist(original_XY, method = "euclidean")
        weight_matrix <- 1 - exp(-as.matrix(spatial_dist))
        
        dist_mat <- as.matrix(dist_matrix)
        if (!all(dim(dist_mat) == dim(weight_matrix))) {
          showNotification("Errore: dimensioni incompatibili tra le matrici di distanza", type = "error")
          return()
        }
        
        dist_mat <- dist_mat * weight_matrix
        dist_matrix <- as.dist(dist_mat)
      }
      
      incProgress(0.7, detail = "Clustering...")
      
      set.seed(123)
      if (input$algorithm == "K-means") {
        
        if (tolower(input$method) %in% c("euclidean", "squared_euclidean")) {
          clusters <- kmeans(original_IM, centers = input$k, nstart = 10)$cluster
          
        } else {
          # MDS per distanze non euclidee (solo autovalori positivi)
          mds_result <- cmdscale(dist_matrix, eig = TRUE, k = nrow(original_IM) - 1)
          eigvals <- mds_result$eig[mds_result$eig > 0]
          
          if (length(eigvals) == 0) {
            showNotification("⚠️ Nessun autovalore positivo trovato in MDS. Controlla i dati.", type = "error")
            return()
          }
          
          var_explained <- cumsum(eigvals) / sum(eigvals)
          k_mds <- which(var_explained >= input$mds_var / 100)[1]
          
          if (is.na(k_mds)) {
            showNotification("⚠️ Nessuna dimensione raggiunge la varianza richiesta.", type = "warning")
            return()
          }
          
          mds_coords <- mds_result$points[, 1:k_mds, drop = FALSE]
          clusters <- kmeans(mds_coords, centers = input$k)$cluster
          
          # Screeplot aggiornato
          output$scree_plot_zip <- renderPlot({
            plot(var_explained * 100, type = "b", pch = 19,
                 xlab = "Number of dimensions", ylab = "Cumulative explained variance (%)",
                 main = "Multidimensional Scaling", sub = "Screeplot",
                 ylim = c(0, 100), xlim = c(1, min(length(var_explained), 50)))
            abline(h = input$mds_var, col = "red", lty = 2)
          })
        }
        
      } else {
        hc <- hclust(dist_matrix, method = input$hclust_method)
        clusters <- cutree(hc, k = input$k)
      }
      
      output$plot_zip_orig <- renderPlotly({
        original_XY$cluster <- as.factor(clusters)
        
        p <- ggplot(original_XY, aes(
          x = X, y = Y, color = cluster,
          text = paste0("X: ", X, "<br>Y: ", Y, "<br>Cluster: ", cluster)
        )) +
          geom_point(shape = 15, size = 1.5) +
          scale_color_viridis_d(option = "plasma") +
          labs(title = "Clusters in the original tissue", x = "X", y = "Y") +
          coord_fixed() +
          theme_minimal()
        
        ggplotly(p, tooltip = "text")
      })
      
      # ---- Classificazione nuovo tessuto (ZIP) ----
      if (input$mode_zip == "tessuto") {
        incProgress(0.8, detail = "Classifying new tissue...")
        
        req(input$new_zip)
        
        res2 <- fun_extractData(input$new_zip)
        down2 <- downsample_grid(res2$featureMatrix, res2$coordinates, input$cell_size)
        new_IM <- down2$IM_reduced
        new_XY <- down2$XY_reduced
        
        colnames(original_IM) <- round(as.numeric(colnames(original_IM)), 1)
        colnames(new_IM) <- round(as.numeric(colnames(new_IM)), 1)
        common_mz <- intersect(colnames(original_IM), colnames(new_IM))
        
        if (length(common_mz) == 0) {
          new_XY$cluster_plot <- as.factor("unclassified")
          
          output$plot_zip_new <- renderPlotly({
            p <- ggplot(new_XY, aes(
              x = X, y = Y,
              text = paste0("X: ", X, "<br>Y: ", Y, "<br>Status: unclassified")
            )) +
              geom_point(aes(color = cluster_plot), shape = 15, size = 1.5, show.legend = TRUE) +
              scale_color_manual(
                name = "cluster",
                values = c("unclassified" = "grey80"),
                drop = FALSE
              ) +
              labs(
                title = "Clusters in the new tissue",
                subtitle = "No compatible m/z values",
                x = "X", y = "Y"
              ) +
              coord_fixed() +
              theme_minimal()
            
            ggplotly(p, tooltip = "text")
          })
          
          showNotification("⚠️ No common m/z values between the two tissues. Classification not performed.", type = "warning")
          return()
        }
        
        valid_idx <- complete.cases(new_IM[, common_mz])
        classified_IM <- new_IM[valid_idx, common_mz, drop = FALSE]
        classified_XY <- new_XY[valid_idx, , drop = FALSE]
        unclassified_XY <- new_XY[!valid_idx, , drop = FALSE]
        
        dist_to_new <- proxy::dist(classified_IM, original_IM[, common_mz], method = input$method)
        assigned <- apply(as.matrix(dist_to_new), 1, function(row) clusters[which.min(row)])
        classified_XY$cluster_plot <- as.character(assigned)
        
        if (nrow(unclassified_XY) > 0) {
          unclassified_XY$cluster_plot <- "non_classificato"
        }
        
        combined_df <- rbind(classified_XY, unclassified_XY)
        
        livelli_cluster <- as.character(1:input$k)
        combined_df$cluster_plot <- factor(combined_df$cluster_plot, levels = c(livelli_cluster, "non_classificato"))
        
        cluster_colors <- viridis::viridis(input$k, option = "plasma")
        names(cluster_colors) <- livelli_cluster
        cluster_colors <- c(cluster_colors, "non_classificato" = "grey80")
        
        output$plot_zip_new <- renderPlotly({
          p <- ggplot(combined_df, aes(
            x = X, y = Y, color = cluster_plot,
            text = paste0("X: ", X, "<br>Y: ", Y, "<br>Cluster: ", cluster_plot)
          )) +
            geom_point(shape = 15, size = 1.5) +
            scale_color_manual(values = cluster_colors, drop = FALSE, name = "cluster") +
            labs(
              title = "Clusters in the new tissue",
              x = "X", y = "Y"
            ) +
            coord_fixed() +
            theme_minimal()
          
          ggplotly(p, tooltip = "text")
        })
      }
      # <-- Alla fine della procedura -->
      incProgress(1, detail = "Completed!")
      showNotification("Classification completed ✅", type = "message", duration = 8)
    }) # fine withProgress
    
    
    output$downloadClusters <- downloadHandler(
  filename = function() {
    paste0("clusters_", input$mode_zip, "_", Sys.Date(), ".csv")
  },
  content = function(file) {
    # Check che i cluster esistano
    if (!exists("assigned") || is.null(assigned)) {
      showNotification("❌ No cluster assignment found.", type = "error")
      return(NULL)
    }

    # 1. Dati clustering
    cluster_df <- data.frame(Pixel = seq_along(assigned), Cluster = assigned)

    # Riga vuota per separazione
    empty_row <- data.frame(Pixel = "", Cluster = "")

    # 2. Report descrittivo
    report_df <- data.frame(
      Pixel = c("=== CLUSTERING REPORT ===",
                "Cell size",
                "Clustering method:",
                "Distance method:",
                "Spatial penalization:",
                "Number of clusters:",
                if (input$algorithm == "K-means" && !(input$method %in% c("euclidean", "squared_euclidean"))) "MDS dimension:" else NULL,
                "Date:"),
      Cluster = c("",
                  input$cell_size,
                  input$algorithm,
                  input$method,
                  ifelse(input$weight_space, "Yes", "No"),
                  input$k,
                  if (input$algorithm == "K-means" && !(input$method %in% c("euclidean", "squared_euclidean"))) input$mds_var else NULL,
                  format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      stringsAsFactors = FALSE
    )

    # 3. Salva su file
    final_df <- rbind(cluster_df, empty_row, report_df)
    write.csv(final_df, file, row.names = FALSE)
  },
  contentType = "text/csv"
)

  })
  
  
}
    

# -------- AVVIO --------
shinyApp(ui, server)


