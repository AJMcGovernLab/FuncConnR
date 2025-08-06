# ConsensusConnectR - Multimethod Consensus-Based Preclinical Functional Connectivity Analysis
# Integrates five correlation methods with consensus approach for robust connectivity mapping

library(shiny)
library(shinydashboard)
library(DT)
library(shinyjs)
library(colourpicker)
# library(ggraph)  # Optional - will create alternative if not available
library(ggplot2)
# library(patchwork)  # Optional - will create alternative if not available
library(scales)

# Additional packages for enhanced correlation methods
suppressPackageStartupMessages({
  # Try to load optional packages for enhanced correlation methods
  tryCatch(library(psych), error = function(e) cat("Note: psych package not available - some correlation methods will use Pearson fallback\n"))
  tryCatch(library(corpcor), error = function(e) cat("Note: corpcor package not available - shrinkage correlation will use Pearson fallback\n"))
})

# Define %||% operator (null-coalescing operator)
`%||%` <- function(x, y) if(is.null(x)) y else x

# Source all original modules
source("modules/analysis_functions.R")
source("modules/visualization_functions.R") 
source("modules/ui_components.R")
source("modules/comprehensive_download_handler.R")

# ========================================================================
# ADVANCED ANALYSIS VISUALIZATION FUNCTIONS
# ========================================================================

# Render advanced MST metrics with central nodes identification
render_advanced_mst_metrics <- function(mst_results, group_colors = NULL) {
  if(is.null(mst_results) || length(mst_results) == 0) return()
  
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
  
  # 1. MST Basic Metrics
  groups <- names(mst_results)
  metrics_df <- data.frame(
    Group = groups,
    Nodes = sapply(mst_results, function(x) x$n_nodes),
    Edges = sapply(mst_results, function(x) x$n_edges),
    Diameter = sapply(mst_results, function(x) x$diameter),
    Avg_Path_Length = sapply(mst_results, function(x) x$avg_path_length),
    Max_Degree = sapply(mst_results, function(x) x$max_degree)
  )
  
  # Get group colors
  colors <- if(!is.null(group_colors)) {
    sapply(groups, function(g) group_colors[[g]] %||% "#3498DB")
  } else {
    rainbow(length(groups))
  }
  
  # Plot diameter comparison
  barplot(metrics_df$Diameter, names.arg = metrics_df$Group, 
          col = colors, main = "MST Diameter by Group",
          ylab = "Diameter", las = 2, cex.names = 0.8)
  
  # Plot max degree comparison
  barplot(metrics_df$Max_Degree, names.arg = metrics_df$Group,
          col = colors, main = "MST Max Degree by Group", 
          ylab = "Max Degree", las = 2, cex.names = 0.8)
  
  # Plot average path length
  valid_path_lengths <- metrics_df$Avg_Path_Length[is.finite(metrics_df$Avg_Path_Length)]
  if(length(valid_path_lengths) > 0) {
    barplot(ifelse(is.finite(metrics_df$Avg_Path_Length), metrics_df$Avg_Path_Length, 0),
            names.arg = metrics_df$Group, col = colors,
            main = "MST Average Path Length", ylab = "Avg Path Length", 
            las = 2, cex.names = 0.8)
  } else {
    plot(1, type = "n", main = "MST Average Path Length", axes = FALSE)
    text(1, 1, "No connected components", cex = 1.5)
  }
  
  # Summary statistics
  plot(1, type = "n", xlim = c(0, 2), ylim = c(0, 2), 
       main = "MST Summary", axes = FALSE, xlab = "", ylab = "")
  text(1, 1.7, paste("Groups analyzed:", length(groups)), cex = 1.2, font = 2)
  text(1, 1.4, paste("Total nodes:", sum(metrics_df$Nodes)), cex = 1.1)
  text(1, 1.1, paste("Total MST edges:", sum(metrics_df$Edges)), cex = 1.1)
  text(1, 0.8, paste("Avg diameter:", round(mean(metrics_df$Diameter), 2)), cex = 1.1)
}

# Render MST central nodes analysis
render_mst_central_nodes <- function(mst_results, brain_areas = NULL, area_colors = NULL, group_colors = NULL) {
  if(is.null(mst_results) || length(mst_results) == 0) return()
  
  par(mfrow = c(2, 2), mar = c(8, 4, 3, 2))
  
  groups <- names(mst_results)
  
  for (i in 1:min(4, length(groups))) {
    group <- groups[i]
    mst_data <- mst_results[[group]]
    
    if (!is.null(mst_data$central_nodes)) {
      # Get top central nodes by degree
      central_nodes <- mst_data$central_nodes$degree[1:min(5, length(mst_data$central_nodes$degree))]
      degree_values <- mst_data$degree_centrality[central_nodes]
      
      # Get colors for nodes based on brain areas
      node_colors <- rep("#3498DB", length(central_nodes))
      if (!is.null(brain_areas) && !is.null(area_colors)) {
        for (j in seq_along(central_nodes)) {
          node <- central_nodes[j]
          for (area_name in names(brain_areas)) {
            if (node %in% brain_areas[[area_name]]) {
              node_colors[j] <- area_colors[[area_name]]
              break
            }
          }
        }
      }
      
      barplot(degree_values, names.arg = central_nodes, col = node_colors,
              main = paste("Top MST Central Nodes:", group),
              ylab = "MST Degree", las = 2, cex.names = 0.7)
    } else {
      plot(1, type = "n", main = paste("MST Central Nodes:", group), axes = FALSE)
      text(1, 1, "No central nodes data", cex = 1.2)
    }
  }
}

# Render advanced MST networks visualization
render_advanced_mst_networks <- function(mst_results, correlations, brain_areas = NULL, area_colors = NULL, group_colors = NULL) {
  if(is.null(mst_results) || length(mst_results) == 0) return()
  
  groups <- names(mst_results)
  n_groups <- length(groups)
  
  if(n_groups == 1) {
    par(mfrow = c(1, 1), mar = c(2, 2, 3, 2))
  } else if(n_groups == 2) {
    par(mfrow = c(1, 2), mar = c(2, 2, 3, 2))
  } else {
    par(mfrow = c(2, 2), mar = c(2, 2, 3, 2))
  }
  
  for(group in groups[1:min(4, n_groups)]) {
    mst_data <- mst_results[[group]]
    
    if(!is.null(mst_data$mst_graph) && has_igraph) {
      mst_graph <- mst_data$mst_graph
      
      # Set node colors based on brain areas
      node_names <- V(mst_graph)$name
      if(is.null(node_names)) node_names <- paste0("Node", 1:vcount(mst_graph))
      
      node_colors <- rep("#95A5A6", vcount(mst_graph))
      if(!is.null(brain_areas) && !is.null(area_colors)) {
        for(area_name in names(brain_areas)) {
          matching_nodes <- which(node_names %in% brain_areas[[area_name]])
          if(length(matching_nodes) > 0) {
            node_colors[matching_nodes] <- area_colors[[area_name]]
          }
        }
      }
      
      # Set node sizes based on degree
      degrees <- degree(mst_graph)
      node_sizes <- scales::rescale(degrees, to = c(8, 20))
      
      # Highlight central nodes
      if(!is.null(mst_data$central_nodes$degree)) {
        central_indices <- which(node_names %in% mst_data$central_nodes$degree[1:3])
        if(length(central_indices) > 0) {
          node_sizes[central_indices] <- node_sizes[central_indices] * 1.3
        }
      }
      
      # Create layout
      layout <- layout_with_fr(mst_graph)
      
      # Plot MST
      plot(mst_graph, 
           layout = layout,
           vertex.size = node_sizes,
           vertex.color = node_colors,
           vertex.frame.color = "white",
           vertex.label.cex = 0.6,
           vertex.label.color = "black",
           edge.width = 2,
           edge.color = "gray60",
           main = paste("MST Network:", group))
      
      # Add legend for central nodes
      legend("topright", legend = c("Most Central", "Other Nodes"),
             pch = 21, pt.bg = c("red", "gray"), pt.cex = c(1.5, 1),
             cex = 0.7, bty = "n")
      
    } else {
      plot(1, type = "n", main = paste("MST Network:", group), axes = FALSE)
      text(1, 1, "MST not available", cex = 1.2)
    }
  }
}

# Render PCA analysis from precomputed results
render_pca_analysis_results <- function(pca_results, group_colors = NULL) {
  if(is.null(pca_results) || length(pca_results) == 0) {
    plot(1, type = "n", main = "PCA Analysis not available")
    return()
  }
  
  groups <- names(pca_results)
  n_groups <- length(groups)
  
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
  
  # Get group colors
  colors <- if(!is.null(group_colors)) {
    sapply(groups, function(g) group_colors[[g]] %||% "#3498DB")
  } else {
    rainbow(n_groups)
  }
  
  # Plot 1: PCA scores (PC1 vs PC2) for each group
  for(i in 1:min(4, n_groups)) {
    group <- groups[i]
    pca_data <- pca_results[[group]]
    
    if(!is.null(pca_data) && !is.null(pca_data$scores)) {
      scores <- pca_data$scores
      if(ncol(scores) >= 2) {
        plot(scores[,1], scores[,2], 
             col = colors[i], pch = 19, cex = 1.2,
             xlab = paste("PC1 (", round(pca_data$variance_explained[1]*100, 1), "%)"),
             ylab = paste("PC2 (", round(pca_data$variance_explained[2]*100, 1), "%)"),
             main = paste("PCA Scores:", group))
        abline(h = 0, v = 0, lty = 2, col = "gray")
        
        # Add region labels if available
        if(!is.null(rownames(scores))) {
          text(scores[,1], scores[,2], labels = rownames(scores), 
               pos = 3, cex = 0.6, col = colors[i])
        }
      } else {
        plot(1, type = "n", main = paste("PCA Scores:", group), axes = FALSE)
        text(1, 1, "Insufficient PCs", cex = 1.2)
      }
    } else {
      plot(1, type = "n", main = paste("PCA Scores:", group), axes = FALSE)
      text(1, 1, "PCA failed", cex = 1.2)
    }
  }
  
  # Fill remaining plots if less than 4 groups
  if(n_groups < 4) {
    for(i in (n_groups + 1):4) {
      plot.new()
    }
  }
}

# Render PCA loadings from precomputed results
render_pca_loadings_results <- function(pca_results, brain_areas = NULL, area_colors = NULL, group_colors = NULL) {
  if(is.null(pca_results) || length(pca_results) == 0) {
    plot(1, type = "n", main = "PCA Loadings not available")
    return()
  }
  
  groups <- names(pca_results)
  n_groups <- length(groups)
  
  par(mfrow = c(2, 2), mar = c(8, 4, 3, 2))
  
  # Plot loadings for PC1 for each group
  for(i in 1:min(4, n_groups)) {
    group <- groups[i]
    pca_data <- pca_results[[group]]
    
    if(!is.null(pca_data) && !is.null(pca_data$loadings)) {
      loadings_pc1 <- pca_data$loadings[,1]
      
      # Get colors for variables based on brain areas
      var_colors <- rep("#3498DB", length(loadings_pc1))
      if(!is.null(brain_areas) && !is.null(area_colors)) {
        var_names <- names(loadings_pc1)
        for(j in seq_along(var_names)) {
          var <- var_names[j]
          for(area_name in names(brain_areas)) {
            if(var %in% brain_areas[[area_name]]) {
              var_colors[j] <- area_colors[[area_name]]
              break
            }
          }
        }
      }
      
      # Sort by absolute loading values
      sorted_indices <- order(abs(loadings_pc1), decreasing = TRUE)
      top_indices <- sorted_indices[1:min(10, length(sorted_indices))]
      
      barplot(loadings_pc1[top_indices], 
              names.arg = names(loadings_pc1)[top_indices],
              col = var_colors[top_indices],
              main = paste("PC1 Loadings:", group),
              ylab = "Loading", las = 2, cex.names = 0.7)
      abline(h = 0, lty = 2)
      
    } else {
      plot(1, type = "n", main = paste("PC1 Loadings:", group), axes = FALSE)
      text(1, 1, "No loadings data", cex = 1.2)
    }
  }
}

# Render PCA variance explained from precomputed results
render_pca_variance_results <- function(pca_results, group_colors = NULL) {
  if(is.null(pca_results) || length(pca_results) == 0) {
    plot(1, type = "n", main = "PCA Variance not available")
    return()
  }
  
  groups <- names(pca_results)
  n_groups <- length(groups)
  
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 2))
  
  # Get group colors
  colors <- if(!is.null(group_colors)) {
    sapply(groups, function(g) group_colors[[g]] %||% "#3498DB")
  } else {
    rainbow(n_groups)
  }
  
  # Plot 1: Variance explained by each PC
  max_pcs <- max(sapply(pca_results, function(x) {
    if(!is.null(x) && !is.null(x$variance_explained)) length(x$variance_explained) else 0
  }))
  
  if(max_pcs > 0) {
    plot(1:max_pcs, rep(0, max_pcs), type = "n",
         ylim = c(0, max(sapply(pca_results, function(x) {
           if(!is.null(x) && !is.null(x$variance_explained)) max(x$variance_explained) else 0
         }))),
         xlab = "Principal Component", ylab = "Variance Explained",
         main = "Variance Explained by PC")
    
    for(i in seq_along(groups)) {
      group <- groups[i]
      pca_data <- pca_results[[group]]
      if(!is.null(pca_data) && !is.null(pca_data$variance_explained)) {
        lines(1:length(pca_data$variance_explained), pca_data$variance_explained,
              col = colors[i], lwd = 2, type = "b", pch = 19)
      }
    }
    
    legend("topright", legend = groups, col = colors, lwd = 2, cex = 0.8)
  }
  
  # Plot 2: Cumulative variance explained
  if(max_pcs > 0) {
    plot(1:max_pcs, rep(0, max_pcs), type = "n", ylim = c(0, 1),
         xlab = "Principal Component", ylab = "Cumulative Variance",
         main = "Cumulative Variance Explained")
    
    for(i in seq_along(groups)) {
      group <- groups[i]
      pca_data <- pca_results[[group]]
      if(!is.null(pca_data) && !is.null(pca_data$cumulative_variance)) {
        lines(1:length(pca_data$cumulative_variance), pca_data$cumulative_variance,
              col = colors[i], lwd = 2, type = "b", pch = 19)
      }
    }
    
    abline(h = c(0.8, 0.9), lty = 2, col = "gray")
    text(max_pcs/2, 0.82, "80%", cex = 0.8, col = "gray")
    text(max_pcs/2, 0.92, "90%", cex = 0.8, col = "gray")
    
    legend("bottomright", legend = groups, col = colors, lwd = 2, cex = 0.8)
  }
}

# Enhanced network plotting function using base igraph (compatible with existing setup)
create_enhanced_network_plots <- function(networks, brain_areas = NULL, area_colors = NULL, 
                                        group_colors = NULL, layout = "fr") {
  
  if (length(networks) == 0) {
    plot(1, type = "n", xlab = "", ylab = "", axes = FALSE)
    text(1, 1, "No networks available", cex = 1.5)
    return()
  }
  
  # Calculate grid layout for multiple networks
  n_networks <- length(networks)
  if (n_networks == 1) {
    par(mfrow = c(1, 1), mar = c(2, 2, 3, 2))
  } else if (n_networks == 2) {
    par(mfrow = c(1, 2), mar = c(2, 2, 3, 2))
  } else if (n_networks <= 4) {
    par(mfrow = c(2, 2), mar = c(2, 2, 3, 2))
  } else {
    par(mfrow = c(2, 3), mar = c(2, 2, 3, 2))
  }
  
  for (group_name in names(networks)) {
    network <- networks[[group_name]]
    
    if (!is.null(network) && vcount(network) > 0) {
      
      # Add brain area information to nodes
      if (!is.null(brain_areas)) {
        node_areas <- rep("Other", vcount(network))
        names(node_areas) <- V(network)$name
        
        for (area_name in names(brain_areas)) {
          regions <- brain_areas[[area_name]]
          matching_nodes <- intersect(regions, V(network)$name)
          if (length(matching_nodes) > 0) {
            node_areas[matching_nodes] <- area_name
          }
        }
        V(network)$brain_area <- node_areas
        
        # Set node colors based on brain areas
        if (!is.null(area_colors)) {
          node_colors <- rep("#808080", vcount(network))  # Default gray
          for (area_name in names(area_colors)) {
            area_nodes <- which(V(network)$brain_area == area_name)
            if (length(area_nodes) > 0) {
              node_colors[area_nodes] <- area_colors[area_name]
            }
          }
          V(network)$color <- node_colors
        }
      } else {
        V(network)$color <- "#1F78B4"  # Default blue
      }
      
      # Set node sizes based on degree
      node_degrees <- degree(network)
      if (max(node_degrees) > min(node_degrees)) {
        V(network)$size <- scales::rescale(node_degrees, to = c(5, 15))
      } else {
        V(network)$size <- 8
      }
      
      # Create layout
      if (layout == "circle") {
        graph_layout <- layout_in_circle(network)
      } else if (layout == "fr") {
        graph_layout <- layout_with_fr(network)  
      } else if (layout == "kk") {
        graph_layout <- layout_with_kk(network)
      } else if (layout == "grid") {
        graph_layout <- tryCatch({
          layout_on_grid(network)
        }, error = function(e) {
          layout_with_fr(network)
        })
      } else {
        graph_layout <- layout_with_fr(network)
      }
      
      # Get group color for border
      group_color <- if (!is.null(group_colors) && group_name %in% names(group_colors)) {
        group_colors[[group_name]]
      } else {
        "black"
      }
      
      # Plot the network
      plot(network,
           layout = graph_layout,
           vertex.color = V(network)$color,
           vertex.size = V(network)$size,
           vertex.label = V(network)$name,
           vertex.label.cex = 0.7,
           vertex.label.color = "black",
           vertex.label.dist = 1,
           vertex.frame.color = "white",
           edge.width = abs(E(network)$weight) * 3,
           edge.color = adjustcolor("gray60", alpha.f = 0.7),
           main = paste("Group:", group_name),
           sub = paste("Layout:", toupper(layout)))
      
      # Add border in group color
      box(col = group_color, lwd = 3)
      
    } else {
      # Empty plot for missing network
      plot(1, type = "n", xlab = "", ylab = "", axes = FALSE, main = paste("Group:", group_name))
      text(1, 1, "No network data", cex = 1.2)
    }
  }
  
  # Reset par
  par(mfrow = c(1, 1))
  
  # Add legend for brain areas (only if we have brain areas and colors)
  if (!is.null(brain_areas) && !is.null(area_colors) && length(area_colors) > 0) {
    # Create a simple legend
    legend("bottomright", 
           legend = names(area_colors),
           fill = area_colors,
           title = "Brain Areas",
           cex = 0.8,
           bg = "white")
  }
}

# Original results UI function from ui_components.R
create_results_ui <- function() {
  tagList(
    # Analysis not complete section
    conditionalPanel(
      condition = "!output.analysisComplete",
      fluidRow(
        column(
          width = 12,
          box(
            title = "⏳ Analysis Status", status = "warning", solidHeader = TRUE, width = NULL,
            h4("Analysis Not Started"),
            p("Please complete the data import and brain area assignment first."),
            actionButton("goToDataImportFromResults", "📊 Go to Data Import", class = "btn btn-primary")
          )
        )
      )
    ),
    
    # Analysis complete section
    conditionalPanel(
      condition = "output.analysisComplete",
      fluidRow(
        column(
          width = 12,
          box(
            title = "🎉 Analysis Complete!", status = "success", solidHeader = TRUE, width = NULL,
            verbatimTextOutput("analysisSummary")
          )
        )
      ),
      
      fluidRow(
        column(
          width = 12,
          tabBox(
            id = "results_tabs", width = NULL,
            
            # Step 1: Data Quality & Imputation
            tabPanel("1️⃣ Imputation",
              h4("🔍 Data Imputation Results"),
              p("Missing data handling and quality assessment before correlation analysis."),
              verbatimTextOutput("imputationSummary"),
              br(),
              plotOutput("imputationPlots", height = "600px")
            ),
            
            # Step 2: MultiMethod Correlation  
            tabPanel("2️⃣ MultiMethod Correlation",
              h4("🔬 5-Method Correlation Consensus"),
              p("Using 5-method correlation consensus (Pearson, Spearman, Biweight, Shrinkage, Partial) with median aggregation for optimal correlation values."),
              
              h5("📊 Consensus Quality Metrics"),
              DT::dataTableOutput("consensusMetricsTable"),
              br(),
              
              h4("📈 Consensus Correlation Matrices"),
              plotOutput("correlationPlots", height = "600px"),
              br(),
              h4("📊 Correlation Distributions"),
              plotOutput("correlationDistributionPlot", height = "600px")
            ),
            
            # Step 3: Topology Analysis
            tabPanel("3️⃣ Topology Analysis",
              tabsetPanel(
                # 3a: Percolation per group
                tabPanel("3a. Percolation Analysis",
                  h4("💧 Group-Specific Percolation Analysis"),
                  p("Data-driven threshold selection with percolation analysis performed per group (change from full data threshold)."),
                  plotOutput("groupPercolationPlot", height = "600px"),
                  br(),
                  h4("📊 Global vs Group-Specific Thresholds"),
                  plotOutput("thresholdComparisonPlot", height = "600px")
                ),
                
                # 3b: Networks, Network Analysis, Regional Connectivity, Individual Nodes, Edge Analysis
                tabPanel("3b. Networks",
                  fluidRow(
                    column(3,
                      h5("🎛️ Layout Options"),
                      selectInput("networkLayout", "Choose Layout:",
                                 choices = list(
                                   "Force-Directed (FR)" = "fr", 
                                   "Circular" = "circle", 
                                   "Kamada-Kawai" = "kk",
                                   "Grid" = "grid"
                                 ),
                                 selected = "fr")
                    ),
                    column(9,
                      h4("🕸️ Network Visualizations"),
                      plotOutput("networkPlots", height = "600px")
                    )
                  ),
                  br(),
                  h4("🎨 Network Gallery"),
                  plotOutput("networkGalleryPlot", height = "800px")
                ),
                
                tabPanel("3b. Network Analysis",
                  h4("📊 Network Topology Dashboard"),
                  plotOutput("networkDashboardPlot", height = "600px")
                ),
                
                tabPanel("3b. Regional Connectivity",
                  h4("🧠 Regional Functional Connectivity"),
                  p("Functional connectivity patterns grouped by anatomical regions"),
                  plotOutput("nodeCentralityPlot", height = "600px"),
                  br(),
                  h4("🌡️ Regional Connectivity Heatmap"),
                  plotOutput("nodeHeatmapPlot", height = "500px"),
                  br(),
                  h4("🔗 Inter-Regional Connectivity Patterns"),
                  plotOutput("brainAreaConnectivityPlot", height = "600px")
                ),
                
                tabPanel("3b. Individual Nodes",
                  h4("🎯 Top Individual Nodes"),
                  p("Node-level functional connectivity analysis with anatomical region mapping"),
                  plotOutput("individualNodeCentralityPlot", height = "600px"),
                  br(),
                  h4("🌡️ Individual Node Heatmap"),
                  plotOutput("individualNodeHeatmapPlot", height = "600px")
                ),
                
                tabPanel("3b. Edge Analysis",
                  h4("⚡ Edge Weight Distributions"),
                  plotOutput("edgeMetricsPlot", height = "600px")
                ),
                
                tabPanel("3b. Regional Connectivity",
                  h4("🧠 Regional Connectivity Analysis"),
                  p("Average eigenvector centrality and node strength values by anatomical region"),
                  plotOutput("regionalConnectivityPlot", height = "600px"),
                  br(),
                  h4("📊 Regional Summary Statistics"),
                  plotOutput("regionalSummaryPlot", height = "600px")
                )
              )
            ),
            
            # Step 4: Weighted Network Analysis
            tabPanel("4️⃣ Weighted Network Analysis",
              h4("🔍 Comprehensive Weighted Network Metrics"),
              p("Full correlation network analysis without threshold - using all correlation values as edge weights"),
              
              tabsetPanel(
                # 4a: Calculate weighted metrics for nodes in each group
                tabPanel("4a. Weighted Node Metrics",
                  h5("📊 Eigenvector Centrality & Node Strength Analysis"),
                  plotOutput("weightedNodeMetricsPlot", height = "600px"),
                  br(),
                  h5("📈 All Weighted Network Statistics"),
                  plotOutput("allWeightedStatsPlot", height = "600px")
                ),
                
                # 4b: Hub Node Comparison Across Groups
                tabPanel("4b. Hub Node Comparison",
                  h5("🌟 Hub Node Comparison Across Groups"),
                  plotOutput("weightedEigenvectorHubPlot", height = "600px"),
                  br(),
                  h5("📊 Top Weighted Eigenvector Nodes by Group"),
                  plotOutput("weightedEigenvectorComparison", height = "600px")
                ),
                
                # 4c: Weighted Eigenvector vs Nodestrength per group
                tabPanel("4c. Eigenvector vs Node Strength",
                  h5("🔗 Node Strength vs Weighted Eigenvector Centrality"),
                  plotOutput("strengthEigenvectorPlot", height = "600px"),
                  br(),
                  h5("⚖️ Comparison of Weighted vs Unweighted Eigenvector Centrality"),
                  plotOutput("weightedVsUnweightedPlot", height = "600px")
                ),
                
                # 4d: Stability analysis
                tabPanel("4d. Stability Analysis",
                  h5("📊 Cross-Group Eigenvector Stability"),
                  plotOutput("eigenvectorStabilityPlot", height = "600px"),
                  br(),
                  h5("📈 Rank Changes from Unweighted to Weighted Analysis"),
                  plotOutput("eigenvectorRankChangePlot", height = "600px")
                ),
                
                
                tabPanel("Data Table",
                  h5("📋 Weighted Network Analysis Results"),
                  br(),
                  br(), br(),
                  DT::dataTableOutput("weightedEigenvectorTable")
                )
              )
            ),
            
            # Step 5: Cross Method Comparison
            tabPanel("5️⃣ Cross Method Comparison",
              h4("⚖️ Weighted vs Percolated Network Metrics"),
              p("Comparison between weighted (threshold-free) and percolated (threshold-based) network metrics."),
              
              tabsetPanel(
                tabPanel("Weighted vs Percolation",
                  h5("🎯 Weighted vs Percolation Eigenvector by Group"),
                  plotOutput("weightedVsPercolationEigenvectorPlot", height = "600px"),
                  br(),
                  h5("💪 Weighted vs Percolation Node Strength by Group"),
                  plotOutput("weightedVsPercolationNodeStrengthPlot", height = "600px")
                ),
                
                tabPanel("Eigenvector + Node Strength",
                  h5("🎯 Eigenvector Centrality Comparison"),
                  plotOutput("crossMethodEigenvectorPlot", height = "600px"),
                  br(),
                  h5("💪 Node Strength Comparison"),
                  plotOutput("crossMethodNodeStrengthPlot", height = "600px")
                ),
                
                tabPanel("Method Correlation",
                  h5("📊 Cross-Method Metric Correlations"),
                  plotOutput("methodCorrelationPlot", height = "600px"),
                  br(),
                  h5("📈 Method Agreement Analysis"),
                  plotOutput("methodAgreementPlot", height = "600px")
                ),
                
                tabPanel("Hub Conservation",
                  h5("🔗 Hub Conservation Between Methods"),
                  plotOutput("hubConservationPlot", height = "600px"),
                  br(),
                  h5("📊 Method-Specific Hub Rankings"),
                  plotOutput("hubRankingComparisonPlot", height = "600px")
                )
              )
            ),
            
            # Step 5.5: Advanced Analysis
            tabPanel("5.5️⃣ Advanced Analysis",
              h4("🔬 Advanced Network Analysis"),
              p("MST analysis and PCA for deeper network understanding."),
              
              tabsetPanel(
                tabPanel("MST Analysis",
                  h5("🌳 Minimum Spanning Tree Metrics"),
                  plotOutput("advancedMstMetricsPlot", height = "600px"),
                  br(),
                  h5("🎯 MST Central Nodes"),
                  plotOutput("mstCentralNodesPlot", height = "600px"),
                  br(),
                  h5("🌐 MST Network Visualization"),
                  plotOutput("advancedMstNetworksPlot", height = "600px")
                ),
                
                tabPanel("PCA Analysis",
                  h5("📊 Principal Component Analysis"),
                  plotOutput("pcaAnalysisPlot", height = "600px"),
                  br(),
                  h5("🎯 PCA Loadings"),
                  plotOutput("pcaLoadingsPlot", height = "600px"),
                  br(),
                  h5("📈 PCA Variance Explained"),
                  plotOutput("pcaVariancePlot", height = "600px")
                )
              )
            ),
            
            # Step 6: Summary and Download
            tabPanel("6️⃣ Summary & Download",
              h4("🎯 Analysis Summary Dashboard"),
              plotOutput("summaryDashboardPlot", height = "600px"),
              
              hr(),
              
              h4("📊 Pipeline Overview"),
              plotOutput("pipelineOverviewPlot", height = "400px"),
              
              hr(),
              
              conditionalPanel(
                condition = "output.analysisComplete",
                div(
                  style = "text-align: center; padding: 20px;",
                  h4("📥 Download Complete Analysis Package"),
                  p("Download all results including plots with your custom colors, Excel/CSV data files, and JSON network files from the complete restructured pipeline"),
                  br(),
                  div(style = "text-align: center; margin: 20px;",
                    downloadUI("download_module")
                  )
                )
              )
            )
          )
        )
      )
    )
  )
}

# Enhanced UI with flexible import and brain area assignment
ui <- dashboardPage(
  dashboardHeader(
    title = "🧠 ConsensusConnectR - Multimethod Functional Connectivity Analysis",
    tags$li(
      class = "dropdown",
      actionButton("about_btn", "About", icon = icon("info-circle"))
    )
  ),
  
  dashboardSidebar(
    sidebarMenu(
      id = "sidebarMenu",
      menuItem("📁 Data Import", tabName = "import", icon = icon("upload")),
      menuItem("🧠 Anatomical Regions", tabName = "brain_areas", icon = icon("brain")),
      menuItem("📊 Analysis", tabName = "analysis", icon = icon("chart-bar"))
    ),
    
    hr(),
    
    # Progress indicator
    div(
      class = "progress-indicator",
      h5("Progress"),
      div(id = "progress_import", class = "progress-item progress-pending",
          icon("circle"), " Data Import"),
      div(id = "progress_brain_areas", class = "progress-item progress-pending",
          icon("circle"), " Anatomical Regions"),
      div(id = "progress_analysis", class = "progress-item progress-pending",
          icon("circle"), " Analysis")
    ),
    
    hr(),
    
    div(
      style = "padding: 15px;",
      p("Enhanced Interface", style = "font-weight: bold;"),
      p("Original Pipeline", style = "font-size: 90%;")
    )
  ),
  
  dashboardBody(
    useShinyjs(),
    dashboard_css(),
    
    # Citation Modal
    div(
      id = "citation-modal",
      class = "modal",
      style = "display: block;",
      div(
        class = "modal-content",
        div(
          class = "modal-header",
          h2("🧠 Welcome to ConsensusConnectR"),
          p("Multimethod Consensus-Based Preclinical Functional Connectivity Analysis Platform")
        ),
        div(
          class = "modal-body",
          style = "padding: 30px;",
          div(
            class = "alert alert-info",
            h4("📚 Citation Agreement"),
            p("This software implements a multimethod consensus approach for preclinical functional connectivity analysis."),
            p("By using ConsensusConnectR, you agree to cite our work in your publications."),
            tags$blockquote(
              style = "font-style: italic; background: #f8f9fa; padding: 15px; border-left: 4px solid #007bff;",
              p("Your citation reference will be provided here.")
            )
          ),
          checkboxInput("citationAgreement", "I agree to cite this software", value = FALSE),
          br(),
          div(
            style = "text-align: center;",
            actionButton("proceedToApp", "Proceed to Analysis", 
                        class = "btn-primary btn-lg", disabled = TRUE)
          )
        )
      )
    ),
    
    # Main app content
    div(
      id = "main-app",
      style = "display: none;",
      tabItems(
        # Enhanced Data Import Tab
        tabItem(
          tabName = "import",
          fluidRow(
            column(
              width = 4,
              box(
                title = "📁 Data Upload", status = "primary", solidHeader = TRUE, width = NULL,
                fileInput(
                  "datafile",
                  "Choose CSV File",
                  accept = c(".csv")
                ),
                
                hr(),
                
                downloadButton("download_template", "Download Template", class = "btn-info btn-sm"),
                
                hr(),
                
                div(
                  class = "alert alert-info",
                  h5("📋 Required Format"),
                  tags$ul(
                    tags$li("ID column (unique identifier for each subject/sample)"),
                    tags$li("Group columns (e.g., Sex, Treatment, Group)"),
                    tags$li("Regional activity/expression measurements"),
                    tags$li("Optional: Behavioral measures")
                  )
                )
              )
            ),
            
            column(
              width = 8,
              conditionalPanel(
                condition = "output.hasData",
                box(
                  title = "⚙️ Data Configuration", status = "success", solidHeader = TRUE, width = NULL,
                  
                  p("Select the appropriate columns from your dataset using the dropdown menus below:"),
                  
                  fluidRow(
                    column(
                      width = 6,
                      selectInput(
                        "id_column",
                        "ID Column:",
                        choices = c("Upload data first..." = "")
                      )
                    ),
                    column(
                      width = 6,
                      selectizeInput(
                        "group_columns",
                        "Group Columns:",
                        choices = NULL,
                        multiple = TRUE
                      )
                    )
                  ),
                  
                  fluidRow(
                    column(
                      width = 6,
                      selectizeInput(
                        "behavior_columns",
                        "Behavioral Columns (Optional):",
                        choices = NULL,
                        multiple = TRUE,
                        options = list(placeholder = "Optional")
                      )
                    ),
                    column(
                      width = 6,
                      br(),
                      checkboxInput("combine_groups", "Create Combined Groups", value = TRUE)
                    )
                  ),
                  
                  hr(),
                  
                  div(
                    style = "text-align: center;",
                    actionButton("configure_data", "Configure Data", 
                                icon = icon("check"), class = "btn-success")
                  )
                )
              )
            )
          ),
          
          fluidRow(
            column(
              width = 12,
              conditionalPanel(
                condition = "output.hasData",
                box(
                  title = "📊 Data Preview", status = "info", solidHeader = TRUE, width = NULL,
                  DTOutput("data_preview_table")
                )
              )
            )
          ),
          
          uiOutput("successPanel")
        ),
        
        # Anatomical Regions Assignment Tab
        tabItem(
          tabName = "brain_areas",
          conditionalPanel(
            condition = "!output.dataConfigured",
            div(
              class = "alert alert-warning",
              icon("exclamation-triangle"),
              " Please complete Data Import first."
            )
          ),
          
          conditionalPanel(
            condition = "output.dataConfigured",
            fluidRow(
              column(
                width = 6,
                box(
                  title = "🧠 Region & Group Selection", status = "primary", solidHeader = TRUE, width = NULL,
                  
                  selectizeInput(
                    "selected_regions",
                    "Measurement Variables to Include:",
                    choices = NULL,
                    multiple = TRUE,
                    options = list(plugins = list('remove_button'))
                  ),
                  
                  div(
                    class = "btn-group",
                    actionButton("select_all_regions", "Select All", class = "btn-sm btn-outline-primary"),
                    actionButton("clear_regions", "Clear", class = "btn-sm btn-outline-secondary")
                  ),
                  
                  hr(),
                  
                  selectizeInput(
                    "selected_groups",
                    "Groups to Include:",
                    choices = NULL,
                    multiple = TRUE,
                    options = list(plugins = list('remove_button'))
                  ),
                  
                  div(
                    class = "btn-group",
                    actionButton("select_all_groups", "Select All", class = "btn-sm btn-outline-primary"),
                    actionButton("clear_groups", "Clear", class = "btn-sm btn-outline-secondary")
                  ),
                  
                  hr(),
                  
                  checkboxInput("include_behavior_in_analysis", "Include behavioral data", value = FALSE)
                )
              ),
              
              column(
                width = 6,
                box(
                  title = "🎨 Anatomical Region Mapping", status = "success", solidHeader = TRUE, width = NULL,
                  
                  div(
                    style = "max-height: 300px; overflow-y: auto;",
                    uiOutput("brain_area_assignments")
                  ),
                  
                  hr(),
                  
                  fluidRow(
                    column(
                      width = 8,
                      textInput("new_brain_area", "New Anatomical Region:", placeholder = "e.g., Prefrontal Cortex, Hippocampus")
                    ),
                    column(
                      width = 4,
                      br(),
                      actionButton("add_brain_area", "Add", icon = icon("plus"), class = "btn-sm btn-primary")
                    )
                  )
                )
              )
            ),
            
            fluidRow(
              column(
                width = 12,
                box(
                  title = "🎨 Color Configuration", status = "warning", solidHeader = TRUE, width = NULL,
                  
                  tabsetPanel(
                    tabPanel(
                      "Anatomical Region Colors",
                      br(),
                      uiOutput("brain_area_colors")
                    ),
                    tabPanel(
                      "Group Colors", 
                      br(),
                      uiOutput("group_colors")
                    )
                  )
                )
              )
            ),
            
            
            fluidRow(
              column(
                width = 12,
                div(
                  style = "text-align: center; margin: 20px;",
                  actionButton("run_original_analysis", "🔬 Run Enhanced 5-Method Correlation Analysis", 
                              icon = icon("play"), class = "btn-success btn-lg")
                )
              )
            )
          )
        ),
        
        # Analysis Results Tab (Original)
        tabItem(
          tabName = "analysis",
          conditionalPanel(
            condition = "!output.analysisComplete",
            div(
              class = "alert alert-info",
              icon("info-circle"),
              " Complete data import and anatomical region mapping to run functional connectivity analysis."
            )
          ),
          
          conditionalPanel(
            condition = "output.analysisComplete",
            # Insert all original results UI here
            create_results_ui()
          )
        )
      )
    )
  )
)

# Enhanced server function 
server <- function(input, output, session) {
  
  # All original reactive values
  analysis_results <- reactiveValues(
    complete = FALSE,
    imputation = NULL,
    correlations = NULL,
    networks = NULL,
    adjacency_matrices = NULL,
    global_metrics = NULL,
    node_metrics = NULL,
    similarities = NULL,
    conservation = NULL,
    raw_data = NULL,
    imputed_data = NULL
  )
  
  # Enhanced reactive values
  ui_state <- reactiveValues(
    data_imported = FALSE,
    data_configured = FALSE,
    brain_areas_configured = FALSE,
    analysis_complete = FALSE,
    raw_data = NULL,
    processed_data = NULL,
    column_info = NULL,
    brain_areas = list(),
    area_colors = list(),
    group_colors = list()
  )
  
  # Default brain areas (from original)
  default_brain_areas <- list(
    "Dorsal HPC" = c("dDG", "dCA1", "dCA2", "dCA3"),
    "Ventral HPC" = c("vDG", "vCA1", "vCA3"),
    "Subiculum" = c("dSub", "vSub"),
    "Nucleus Accumbens" = c("NaC", "NaS"),
    "Frontal" = c("ACC", "IL", "PRL"),
    "Amygdala" = c("CeA", "BLA", "LA"),
    "Retrosplenial" = c("RSGab", "RSGc", "RSD")
  )
  
  default_area_colors <- c(
    "Dorsal HPC" = "#D3ADC4", "Ventral HPC" = "#C88AB1", "Subiculum" = "#9B59B6",
    "Nucleus Accumbens" = "#A3DFD7", "Frontal" = "#FAE9BD", 
    "Amygdala" = "#F0BC94", "Retrosplenial" = "#85C1E9"
  )
  
  # Citation modal handling
  observeEvent(input$citationAgreement, {
    if (input$citationAgreement) {
      runjs("$('#proceedToApp').prop('disabled', false);")
    } else {
      runjs("$('#proceedToApp').prop('disabled', true);")
    }
  })
  
  observeEvent(input$proceedToApp, {
    runjs("
      document.getElementById('citation-modal').style.display = 'none';
      document.getElementById('main-app').style.display = 'block';
    ")
    showNotification("✅ Welcome! Start by uploading your data.", duration = 3)
  })
  
  # === DATA IMPORT SECTION ===
  
  # File upload
  uploaded_data <- reactive({
    if (is.null(input$datafile)) return(NULL)
    
    tryCatch({
      # Robust CSV reading with encoding detection and error handling
      data <- tryCatch({
        # Try UTF-8 first (handles BOM)
        read.csv(input$datafile$datapath, stringsAsFactors = FALSE, encoding = "UTF-8")
      }, error = function(e) {
        # Fallback to default encoding
        tryCatch({
          read.csv(input$datafile$datapath, stringsAsFactors = FALSE)
        }, error = function(e2) {
          # Final fallback with different separators
          read.csv(input$datafile$datapath, stringsAsFactors = FALSE, sep = ";")
        })
      })
      
      # Clean column names comprehensively
      original_names <- names(data)
      # Remove BOM, leading/trailing spaces, and make valid R names
      names(data) <- trimws(names(data))  # Remove leading/trailing spaces
      names(data) <- gsub("^﻿", "", names(data))  # Remove BOM character
      names(data) <- make.names(names(data), unique = TRUE)  # Make valid R names
      
      ui_state$raw_data <- data  
      ui_state$data_imported <- TRUE
      
      # Update progress indicator
      runjs("
        $('#progress_import').removeClass('progress-pending').addClass('progress-complete');
        $('#progress_import i').removeClass('fa-circle').addClass('fa-check-circle');
      ")
      
      # Auto-detect columns and update choices
      col_names <- names(data)
      
      # Notify user if column names were cleaned
      cleaned_names <- which(original_names != names(data))
      if (length(cleaned_names) > 0) {
        showNotification(
          paste("Column names cleaned:", length(cleaned_names), "columns had spaces/special characters removed"), 
          type = "message", duration = 5
        )
      }
      
      # Update column selectors (user-driven selection only)
      updateSelectInput(session, "id_column", choices = c("Select ID column..." = "", col_names))
      
      # Update other selectors
      updateSelectizeInput(session, "group_columns", choices = col_names)
      updateSelectizeInput(session, "behavior_columns", choices = col_names)
      
      return(data)
    }, error = function(e) {
      showNotification(paste("Error reading file:", e$message), type = "error")
      return(NULL)
    })
  })
  
  output$hasData <- reactive({ !is.null(uploaded_data()) })
  outputOptions(output, "hasData", suspendWhenHidden = FALSE)
  
  # Data preview
  output$data_preview_table <- renderDT({
    req(uploaded_data())
    datatable(uploaded_data(), options = list(pageLength = 10, scrollX = TRUE))
  })
  
  # Configure data
  observeEvent(input$configure_data, {
    req(uploaded_data())
    
    # Validate inputs
    
    # More flexible validation
    if (is.null(input$id_column) || input$id_column == "") {
      showNotification("⚠️ Please select an ID column", type = "warning", duration = 10)
      return()
    }
    
    if (length(input$group_columns) == 0) {
      showNotification("⚠️ Please select at least one group column", type = "warning", duration = 10)
      return()
    }
    
    data <- uploaded_data()
    
    # Store column information
    all_cols <- names(data)
    metadata_cols <- c(input$id_column, input$group_columns, input$behavior_columns)
    region_cols <- setdiff(all_cols, metadata_cols)
    
    ui_state$column_info <- list(
      id_column = input$id_column,
      group_columns = input$group_columns,
      behavior_columns = input$behavior_columns,
      region_columns = region_cols
    )
    
    # Create combined group column if requested
    if (input$combine_groups && length(input$group_columns) > 1) {
      data$Group <- apply(data[input$group_columns], 1, paste, collapse = "_")
    } else if (length(input$group_columns) == 1) {
      data$Group <- data[[input$group_columns[1]]]
    }
    
    ui_state$processed_data <- data
    
    # Initialize brain areas immediately after data configuration
    region_columns <- ui_state$column_info$region_columns
    
    # Initialize brain areas with case-insensitive matching
    matched_areas <- list()
    unassigned_regions <- region_columns
    
    for (area_name in names(default_brain_areas)) {
      area_regions <- default_brain_areas[[area_name]]
      # Case-insensitive matching
      matched <- c()
      for (region in area_regions) {
        # Find case-insensitive matches
        matches <- region_columns[tolower(region_columns) == tolower(region)]
        matched <- c(matched, matches)
      }
      matched <- unique(matched)
      
      if (length(matched) > 0) {
        matched_areas[[area_name]] <- matched
        unassigned_regions <- setdiff(unassigned_regions, matched)
      }
    }
    
    # Always include unassigned regions in "Other" category
    if (length(unassigned_regions) > 0) {
      matched_areas[["Other"]] <- unassigned_regions
    }
    
    # If no brain areas were matched at all, put everything in "All Regions"
    if (length(matched_areas) == 0 || (length(matched_areas) == 1 && "Other" %in% names(matched_areas))) {
      matched_areas <- list("All Regions" = region_columns)
    }
    
    ui_state$brain_areas <- matched_areas
    # Set colors for matched areas
    ui_state$area_colors <- default_area_colors[names(matched_areas)]
    if ("Other" %in% names(matched_areas)) {
      ui_state$area_colors["Other"] <- "#95A5A6"
    }
    if ("All Regions" %in% names(matched_areas)) {
      ui_state$area_colors["All Regions"] <- "#3498DB"
    }
    
    # Set group colors
    groups <- unique(data$Group)
    group_palette <- c("#E74C3C", "#3498DB", "#2ECC71", "#F39C12", "#9B59B6")
    ui_state$group_colors <- setNames(group_palette[1:length(groups)], groups)
        
    ui_state$data_configured <- TRUE
    
    showNotification("✅ Data configured successfully!")
    
    # Automatically switch to brain areas tab after configuration
    updateTabItems(session, "sidebarMenu", "brain_areas")
    showNotification("📍 Navigating to Anatomical Regions tab...", type = "message", duration = 2)
  })
  
  # Render success panel when data is configured
  output$successPanel <- renderUI({
    if (isTRUE(ui_state$data_configured)) {
      # Render success panel with button
      fluidRow(
        column(
          width = 12,
          div(
            class = "alert alert-success",
            h4("✅ Data Successfully Configured!"),
            p("Proceed to Anatomical Regions tab to map your measurements to anatomical structures."),
            div(
              style = "text-align: center; margin-top: 15px;",
              actionButton("go_to_brain_areas", "Go to Anatomical Regions →", 
                          icon = icon("arrow-right"), class = "btn-primary")
            )
          )
        )
      )
    } else {
      # Show nothing when not configured
      NULL
    }
  })
  
  # Reactive outputs for conditionalPanels
  output$dataConfigured <- reactive({ ui_state$data_configured })
  outputOptions(output, "dataConfigured", suspendWhenHidden = FALSE)
  
  # Go to brain areas button
  observeEvent(input$go_to_brain_areas, {
    updateTabItems(session, "sidebarMenu", "brain_areas")
  })
  
  # === BRAIN AREAS SECTION ===
  
  # Initialize brain areas when data is configured
  observe({
    req(ui_state$data_configured)
    req(ui_state$column_info)
    
    if (ui_state$brain_areas_configured) {
      return()
    }
    
    region_columns <- ui_state$column_info$region_columns
    
    # Update region selection - ensure all regions are selected by default
    updateSelectizeInput(session, "selected_regions", 
                        choices = region_columns, selected = region_columns)
    
    showNotification("✅ All brain regions have been selected by default. You can modify the selection as needed.", 
                    type = "message", duration = 4)
    
    # Get unique groups
    if ("Group" %in% names(ui_state$processed_data)) {
      groups <- unique(ui_state$processed_data$Group)
    } else {
      group_col <- ui_state$column_info$group_columns[1]
      groups <- unique(ui_state$processed_data[[group_col]])
    }
    
    updateSelectizeInput(session, "selected_groups", 
                        choices = groups, selected = groups)
    
    showNotification(paste("✅ All", length(groups), "groups have been selected:", paste(groups, collapse = ", ")), 
                    type = "message", duration = 4)
    
    # Initialize brain areas with case-insensitive matching
    matched_areas <- list()
    unassigned_regions <- region_columns
    
    for (area_name in names(default_brain_areas)) {
      area_regions <- default_brain_areas[[area_name]]
      # Case-insensitive matching
      matched <- c()
      for (region in area_regions) {
        # Find case-insensitive matches
        matches <- region_columns[tolower(region_columns) == tolower(region)]
        matched <- c(matched, matches)
      }
      matched <- unique(matched)
      
      if (length(matched) > 0) {
        matched_areas[[area_name]] <- matched
        unassigned_regions <- setdiff(unassigned_regions, matched)
      }
    }
    
    # Always include unassigned regions in "Other" category
    if (length(unassigned_regions) > 0) {
      matched_areas[["Other"]] <- unassigned_regions
    }
    
    # If no brain areas were matched at all, put everything in "All Regions"
    if (length(matched_areas) == 0 || (length(matched_areas) == 1 && "Other" %in% names(matched_areas))) {
      matched_areas <- list("All Regions" = region_columns)
    }
    
    ui_state$brain_areas <- matched_areas
    # Set colors for matched areas
    ui_state$area_colors <- default_area_colors[names(matched_areas)]
    if ("Other" %in% names(matched_areas)) {
      ui_state$area_colors["Other"] <- "#95A5A6"
    }
    if ("All Regions" %in% names(matched_areas)) {
      ui_state$area_colors["All Regions"] <- "#3498DB"
    }
    
    # Set group colors
    group_palette <- c("#E74C3C", "#3498DB", "#2ECC71", "#F39C12", "#9B59B6")
    ui_state$group_colors <- setNames(group_palette[1:length(groups)], groups)
    
    # Update progress
    runjs("
      $('#progress_brain_areas').removeClass('progress-pending').addClass('progress-complete');
      $('#progress_brain_areas i').removeClass('fa-circle').addClass('fa-check-circle');
    ")
    
    # Mark brain areas as configured to prevent infinite loop
    ui_state$brain_areas_configured <- TRUE
  })
  
  # Region/group selection buttons
  observeEvent(input$select_all_regions, {
    req(ui_state$column_info)
    updateSelectizeInput(session, "selected_regions", selected = ui_state$column_info$region_columns)
  })
  
  observeEvent(input$clear_regions, {
    updateSelectizeInput(session, "selected_regions", selected = character(0))
  })
  
  observeEvent(input$select_all_groups, {
    req(ui_state$processed_data)
    if ("Group" %in% names(ui_state$processed_data)) {
      all_groups <- unique(ui_state$processed_data$Group)
    } else {
      group_col <- ui_state$column_info$group_columns[1]
      all_groups <- unique(ui_state$processed_data[[group_col]])
    }
    updateSelectizeInput(session, "selected_groups", selected = all_groups)
  })
  
  observeEvent(input$clear_groups, {
    updateSelectizeInput(session, "selected_groups", selected = character(0))
  })
  
  # Brain area assignments UI
  output$brain_area_assignments <- renderUI({
    req(ui_state$brain_areas)
    req(input$selected_regions)
    
    ui_list <- list()
    
    for (area_name in names(ui_state$brain_areas)) {
      area_id <- gsub("[^a-zA-Z0-9]", "_", area_name)
      current_regions <- ui_state$brain_areas[[area_name]]
      
      ui_list[[area_name]] <- div(
        style = "margin-bottom: 15px; padding: 10px; border: 1px solid #ddd; border-radius: 5px;",
        fluidRow(
          column(4, strong(area_name)),
          column(8, 
            selectizeInput(
              paste0("area_assign_", area_id),
              label = NULL,
              choices = input$selected_regions,
              selected = intersect(current_regions, input$selected_regions),
              multiple = TRUE,
              options = list(placeholder = "Select regions", plugins = list('remove_button'))
            )
          )
        )
      )
    }
    
    do.call(tagList, ui_list)
  })
  
  # Add new brain area
  observeEvent(input$add_brain_area, {
    req(input$new_brain_area)
    
    area_name <- trimws(input$new_brain_area)
    if (area_name != "" && !area_name %in% names(ui_state$brain_areas)) {
      ui_state$brain_areas[[area_name]] <- character(0)
      ui_state$area_colors[area_name] <- "#808080"
      updateTextInput(session, "new_brain_area", value = "")
      showNotification(paste("Added brain area:", area_name))
    }
  })
  
  # Brain area colors UI
  output$brain_area_colors <- renderUI({
    req(ui_state$brain_areas)
    
    ui_list <- list()
    for (area_name in names(ui_state$brain_areas)) {
      area_id <- gsub("[^a-zA-Z0-9]", "_", area_name)
      current_color <- ui_state$area_colors[area_name]
      if (is.na(current_color)) current_color <- "#808080"
      
      ui_list[[area_name]] <- fluidRow(
        style = "margin-bottom: 10px;",
        column(6, strong(area_name)),
        column(6, 
          colourInput(
            paste0("area_color_", area_id),
            label = NULL, value = current_color, showColour = "background"
          )
        )
      )
    }
    
    do.call(tagList, ui_list)
  })
  
  # Group colors UI
  output$group_colors <- renderUI({
    req(input$selected_groups)
    
    ui_list <- list()
    for (group_name in input$selected_groups) {
      group_id <- gsub("[^a-zA-Z0-9]", "_", group_name)
      current_color <- ui_state$group_colors[group_name]
      if (is.null(current_color) || is.na(current_color)) current_color <- "#808080"
      
      ui_list[[group_name]] <- fluidRow(
        style = "margin-bottom: 10px;",
        column(6, strong(group_name)),
        column(6, 
          colourInput(
            paste0("group_color_", group_id),
            label = NULL, value = current_color, showColour = "background"
          )
        )
      )
    }
    
    do.call(tagList, ui_list)
  })
  
  # === RUN ORIGINAL ANALYSIS ===
  
  observeEvent(input$run_original_analysis, {
    # Validate that regions and groups are selected
    if (is.null(input$selected_regions) || length(input$selected_regions) == 0) {
      showNotification("⚠️ Please select at least one brain region to analyze!", type = "error", duration = 5)
      return()
    }
    
    if (is.null(input$selected_groups) || length(input$selected_groups) == 0) {
      showNotification("⚠️ Please select at least one group to analyze!", type = "error", duration = 5)
      return()
    }
    
    # Update brain areas and colors from inputs
    for (area_name in names(ui_state$brain_areas)) {
      area_id <- gsub("[^a-zA-Z0-9]", "_", area_name)
      assign_input <- paste0("area_assign_", area_id)
      color_input <- paste0("area_color_", area_id)
      
      if (!is.null(input[[assign_input]])) {
        ui_state$brain_areas[[area_name]] <- input[[assign_input]]
      }
      if (!is.null(input[[color_input]])) {
        ui_state$area_colors[area_name] <- input[[color_input]]
      }
    }
    
    for (group_name in input$selected_groups) {
      group_id <- gsub("[^a-zA-Z0-9]", "_", group_name)
      color_input <- paste0("group_color_", group_id)
      if (!is.null(input[[color_input]])) {
        ui_state$group_colors[group_name] <- input[[color_input]]
      }
    }
    
    # Set the data for the original analysis
    analysis_results$raw_data <- ui_state$processed_data
    
    # Update progress and switch to analysis tab
    runjs("
      $('#progress_analysis').removeClass('progress-pending').addClass('progress-active');
      $('#progress_analysis i').removeClass('fa-circle').addClass('fa-spinner fa-spin');
    ")
    
    updateTabItems(session, "sidebarMenu", "analysis")
    
    showNotification("🔄 Running restructured analysis pipeline...", duration = NULL, id = "analysis_running")
    
    # === ORIGINAL ANALYSIS PIPELINE ===
    tryCatch({
      data <- analysis_results$raw_data
      
      # Filter data based on selected regions and groups
      region_cols <- input$selected_regions
      selected_groups <- input$selected_groups
      
      # Filter by groups
      if ("Group" %in% names(data)) {
        data <- data[data$Group %in% selected_groups, ]
      }
      
      # Include behavioral data if requested
      analysis_cols <- region_cols
      if (input$include_behavior_in_analysis && !is.null(ui_state$column_info$behavior_columns)) {
        analysis_cols <- c(analysis_cols, ui_state$column_info$behavior_columns)
      }
      
      # Step 1: Data imputation
      showNotification("Step 1/4: Imputing missing data...", duration = 2)
      imputation_result <- perform_mice_imputation(data[, analysis_cols, drop = FALSE])
      analysis_results$imputation <- imputation_result
      
      # Step 2: Correlation analysis
      showNotification("Step 2/4: Computing correlations...", duration = 2)
      
      # Create groups vector
      if ("Group" %in% names(data)) {
        groups <- data$Group
      } else {
        groups <- rep("All_Data", nrow(data))
      }
      
      correlations <- compute_correlations(imputation_result$imputed_data, groups)
      analysis_results$correlations <- correlations
      
      # Step 3: Topology Analysis - with group-specific percolation
      analysis_type <- "threshold_free"  # Always run comprehensive analysis
      
      showNotification("Step 3/5: Topology Analysis - Percolation per group...", duration = 2)
      
      # Calculate both global and group-specific percolation thresholds
      if(length(correlations) > 1) {
        pooled_consensus_matrix <- Reduce("+", correlations) / length(correlations)
      } else {
        pooled_consensus_matrix <- correlations[[1]]
      }
      global_threshold <- calculate_percolation_threshold(pooled_consensus_matrix)
      
      # NEW: Calculate group-specific thresholds (change from full data threshold)
      group_thresholds <- list()
      for(group_name in names(correlations)) {
        group_threshold <- calculate_percolation_threshold(correlations[[group_name]])
        group_thresholds[[group_name]] <- group_threshold
        cat(sprintf("  Group %s threshold: %.3f (vs global: %.3f)\n", 
                   group_name, group_threshold, global_threshold))
      }
      analysis_results$group_thresholds <- group_thresholds
      
      networks <- list()
      adjacency_matrices <- list()
      all_global_metrics <- data.frame()
      all_node_metrics <- data.frame()
      all_edge_metrics <- data.frame()
      all_brain_area_metrics <- data.frame()
      threshold_free_results <- NULL
      
      # Always run percolation-based network analysis (traditional approach)
      for(group_name in names(correlations)) {
        cor_matrix <- correlations[[group_name]]
        
        # NEW: Use group-specific percolation threshold instead of global
        group_threshold <- group_thresholds[[group_name]]
        network <- create_network(cor_matrix, threshold = group_threshold)
        
        if(!is.null(network)) {
          networks[[group_name]] <- network
          
          # Store adjacency matrix
          adj_matrix <- as_adjacency_matrix(network, attr = "weight", sparse = FALSE)
          adjacency_matrices[[group_name]] <- adj_matrix
          
          # Compute comprehensive metrics
          metrics <- compute_network_metrics(network, group_name, ui_state$brain_areas)
          if(!is.null(metrics)) {
            all_global_metrics <- rbind(all_global_metrics, metrics$global)
            all_node_metrics <- rbind(all_node_metrics, metrics$nodes)
            if(!is.null(metrics$edges)) {
              all_edge_metrics <- rbind(all_edge_metrics, metrics$edges)
            }
            if(!is.null(metrics$brain_area_metrics)) {
              all_brain_area_metrics <- rbind(all_brain_area_metrics, metrics$brain_area_metrics)
            }
          }
        }
      }
      
      # Always run threshold-free analysis with all metrics
      # Compute threshold-free metrics
      threshold_free_results <- perform_threshold_free_analysis(
        correlations, 
        names(correlations), 
        analysis_type = "all"  # Always compute all threshold-free metrics
      )
      
      analysis_results$networks <- networks
      analysis_results$adjacency_matrices <- adjacency_matrices
      analysis_results$global_metrics <- all_global_metrics
      analysis_results$node_metrics <- all_node_metrics
      analysis_results$edge_metrics <- all_edge_metrics
      analysis_results$brain_area_metrics <- all_brain_area_metrics
      analysis_results$global_threshold <- global_threshold
      analysis_results$analysis_type <- analysis_type
      analysis_results$threshold_free_results <- threshold_free_results
      
      # Step 4: Network conservation analysis (always run when multiple networks exist)
      showNotification("Step 4/4: Network conservation analysis...", duration = 2)
      
      if(length(networks) > 1) {
        # Network similarities using percolation thresholds
        similarities <- data.frame()
        group_names <- names(networks)
        
        for(i in 1:(length(group_names)-1)) {
          for(j in (i+1):length(group_names)) {
            g1 <- group_names[i]
            g2 <- group_names[j]
            
            cor1 <- correlations[[g1]]
            cor2 <- correlations[[g2]]
            
            sim <- compute_network_similarity(cor1, cor2, threshold = global_threshold)
          
          sim_row <- data.frame(
            Group1 = g1,
            Group2 = g2,
            Jaccard = sim$jaccard,
            Overlap = sim$overlap,
            Shared_Edges = sim$intersection,
            Edges_Group1 = sim$edges_net1,
            Edges_Group2 = sim$edges_net2,
            Edge_Preservation_1to2 = sim$edge_preservation_1to2,
            Edge_Preservation_2to1 = sim$edge_preservation_2to1
          )
          similarities <- rbind(similarities, sim_row)
          }
        }
        
        analysis_results$similarities <- similarities
        
        # Hub conservation
        hub_conservation <- compute_hub_conservation(networks)
        analysis_results$conservation <- list(
          hub_conservation = hub_conservation,
          conservation_possible = !is.null(hub_conservation)
        )
      } else {
        analysis_results$conservation <- list(conservation_possible = FALSE)
      }
      
      
      # Step 5: Weighted Network Analysis
      showNotification("Step 5/6: Weighted Network Analysis - All metrics...", duration = 2)
      
      # Compute weighted eigenvector centrality for full correlation networks
      weighted_eigen_results <- compute_weighted_eigenvector_centrality(correlations, groups)
      analysis_results$weighted_eigenvector <- weighted_eigen_results
      
      # Compare across groups if multiple groups exist
      if(length(networks) > 1 && !is.null(weighted_eigen_results)) {
        comparison_results <- compare_weighted_eigenvector_across_groups(weighted_eigen_results)
        analysis_results$weighted_eigenvector_comparison <- comparison_results
      }
      
      # Identify hub nodes
      if(!is.null(weighted_eigen_results)) {
        hub_results <- identify_weighted_eigenvector_hubs(weighted_eigen_results, top_n = 10)
        analysis_results$weighted_eigenvector_hubs <- hub_results
      }
      
      # NEW: Minimum Spanning Tree Analysis
      showNotification("Step 5/6: Computing MST metrics...", duration = 1)
      mst_results <- compute_mst_analysis(correlations)
      analysis_results$mst_results <- mst_results
      
      # NEW: PCA Analysis
      showNotification("Step 5/6: Computing PCA analysis...", duration = 1)
      pca_results <- compute_pca_analysis(correlations)
      analysis_results$pca_results <- pca_results
      
      # Step 6: Cross Method Comparison
      showNotification("Step 6/6: Cross Method Comparison - Weighted vs Percolated...", duration = 2)
      
      # Compare weighted vs percolated network metrics
      cross_method_comparison <- perform_cross_method_comparison(
        weighted_eigenvector = weighted_eigen_results,
        threshold_free_results = threshold_free_results,
        percolation_networks = networks,
        node_metrics = all_node_metrics
      )
      analysis_results$cross_method_comparison <- cross_method_comparison
      
      # Compute consensus eigenvector centrality (percolation + weighted)
      consensus_eigen_results <- compute_consensus_eigenvector_centrality(networks, weighted_eigen_results)
      analysis_results$consensus_eigenvector <- consensus_eigen_results
      
      analysis_results$complete <- TRUE
      ui_state$analysis_complete <- TRUE
      
      runjs("
        $('#progress_analysis').removeClass('progress-active').addClass('progress-complete');
        $('#progress_analysis i').removeClass('fa-spinner fa-spin').addClass('fa-check-circle');
      ")
      
      removeNotification(id = "analysis_running")
      showNotification("✅ Complete neural analysis pipeline finished successfully!", duration = 5)
      
    }, error = function(e) {
      runjs("
        $('#progress_analysis').removeClass('progress-active').addClass('progress-pending');
        $('#progress_analysis i').removeClass('fa-spinner fa-spin').addClass('fa-circle');
      ")
      
      removeNotification(id = "analysis_running")
      showNotification(paste("❌ Analysis error:", e$message), type = "error", duration = NULL)
      print(e)  # For debugging
    })
  })
  
  output$analysisComplete <- reactive({ ui_state$analysis_complete })
  outputOptions(output, "analysisComplete", suspendWhenHidden = FALSE)
  
  # === ALL ORIGINAL OUTPUT FUNCTIONS ===
  
  # Analysis Summary
  output$analysisSummary <- renderText({
    req(analysis_results$complete)
    
    total_nodes <- length(colnames(analysis_results$correlations[[1]]))
    total_groups <- length(analysis_results$correlations)
    threshold <- if(!is.null(analysis_results$global_threshold)) {
      round(analysis_results$global_threshold, 3)
    } else {
      "N/A"
    }
    
    analysis_type <- analysis_results$analysis_type %||% "percolation"
    
    summary_lines <- c(
      "🎯 Analysis Complete!\n",
      "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n",
      sprintf("📊 Total brain regions analyzed: %d\n", total_nodes),
      sprintf("👥 Groups compared: %d\n", total_groups),
      sprintf("🎚️ Global percolation threshold: %s\n", threshold),
      sprintf("🧠 Networks successfully created for all groups\n")
    )
    
    # Always show comprehensive analysis information since both analyses always run
    summary_lines <- c(summary_lines,
      sprintf("🔢 Threshold-free metrics computed\n"),
      sprintf("📈 Both percolation and threshold-free analyses available\n"),
      sprintf("🔗 Network conservation analysis completed\n")
    )
    
    summary_lines <- c(summary_lines,
      "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n",
      "✅ All visualizations and metrics are now available in the tabs below!"
    )
    
    paste(summary_lines, collapse = "")
  })
  
  # Imputation Summary
  output$imputationSummary <- renderText({
    req(analysis_results$imputation)
    
    imp_result <- analysis_results$imputation
    
    paste(
      "🔍 Data Imputation Summary\n",
      "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n",
      sprintf("📊 Original missing values: %d\n", imp_result$missing_count),
      sprintf("🎯 Imputation method: %s\n", imp_result$method),
      sprintf("✅ Missing data successfully handled\n"),
      sprintf("📈 Data quality improved for network analysis\n")
    )
  })
  
  # All the plot outputs
  output$imputationPlots <- renderPlot({
    req(analysis_results$imputation)
    render_imputation_plots(analysis_results, analysis_results$raw_data)
  })
  
  output$correlationPlots <- renderPlot({
    req(analysis_results$correlations)
    # Get experimental group colors
    exp_group_colors <- NULL
    if(!is.null(ui_state$group_colors) && length(ui_state$group_colors) > 0) {
      exp_group_colors <- ui_state$group_colors
    }
    render_correlation_plots(analysis_results$correlations, exp_group_colors)
  })
  
  output$correlationDistributionPlot <- renderPlot({
    req(analysis_results$correlations)
    render_correlation_distributions(analysis_results$correlations)
  })
  
  # Consensus metrics table
  output$consensusMetricsTable <- DT::renderDataTable({
    req(analysis_results$correlations)
    
    # Extract consensus metadata
    consensus_metadata <- attr(analysis_results$correlations, "consensus_metadata")
    if(is.null(consensus_metadata)) return(NULL)
    
    # Create summary table
    metrics_df <- data.frame()
    for(group_name in names(consensus_metadata)) {
      metadata <- consensus_metadata[[group_name]]
      metrics_df <- rbind(metrics_df, data.frame(
        Group = group_name,
        `Sample Size` = metadata$sample_size,
        `Brain Regions` = metadata$n_variables,
        `Methods Used` = metadata$n_methods,
        `Method Names` = paste(metadata$methods_used, collapse = ", "),
        `Mean |Correlation|` = round(metadata$mean_abs_correlation, 3),
        `Strong Correlations (≥0.4)` = metadata$strong_correlations,
        `Expected Network Density` = paste0(round(metadata$network_density_at_04 * 100, 1), "%"),
        check.names = FALSE
      ))
    }
    
    DT::datatable(metrics_df, 
                  options = list(scrollX = TRUE, pageLength = 10),
                  caption = "5-Method Correlation Consensus Quality Metrics")
  }, server = FALSE)
  
  output$globalPercolationPlot <- renderPlot({
    req(analysis_results$correlations)
    render_global_percolation_analysis(analysis_results$correlations, analysis_results$global_threshold)
  })
  
  output$thresholdAnalysisPlot <- renderPlot({
    req(analysis_results$correlations)
    render_threshold_analysis(analysis_results$correlations, analysis_results$global_threshold)
  })
  
  output$networkPlots <- renderPlot({
    req(analysis_results$networks)
    
    create_enhanced_network_plots(
      networks = analysis_results$networks,
      brain_areas = ui_state$brain_areas,
      area_colors = ui_state$area_colors,
      group_colors = ui_state$group_colors,
      layout = input$networkLayout %||% "fr"
    )
  })
  
  output$networkGalleryPlot <- renderPlot({
    req(analysis_results$networks)
    
    # Create a gallery view with different layouts in a 2x2 grid
    layout_types <- c("fr", "circle", "kk", "grid")
    par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
    
    for (layout_type in layout_types) {
      # For gallery, show just the first network with different layouts
      first_network_name <- names(analysis_results$networks)[1]
      first_network <- analysis_results$networks[[first_network_name]]
      
      if (!is.null(first_network) && vcount(first_network) > 0) {
        
        # Add brain area coloring
        if (!is.null(ui_state$brain_areas)) {
          node_areas <- rep("Other", vcount(first_network))
          names(node_areas) <- V(first_network)$name
          
          for (area_name in names(ui_state$brain_areas)) {
            regions <- ui_state$brain_areas[[area_name]]
            matching_nodes <- intersect(regions, V(first_network)$name)
            if (length(matching_nodes) > 0) {
              node_areas[matching_nodes] <- area_name
            }
          }
          V(first_network)$brain_area <- node_areas
          
          # Set node colors
          if (!is.null(ui_state$area_colors)) {
            node_colors <- rep("#808080", vcount(first_network))
            for (area_name in names(ui_state$area_colors)) {
              area_nodes <- which(V(first_network)$brain_area == area_name)
              if (length(area_nodes) > 0) {
                node_colors[area_nodes] <- ui_state$area_colors[area_name]
              }
            }
            V(first_network)$color <- node_colors
          }
        } else {
          V(first_network)$color <- "#1F78B4"
        }
        
        # Set node sizes
        node_degrees <- degree(first_network)
        if (max(node_degrees) > min(node_degrees)) {
          V(first_network)$size <- scales::rescale(node_degrees, to = c(5, 12))
        } else {
          V(first_network)$size <- 8
        }
        
        # Create layout
        if (layout_type == "circle") {
          graph_layout <- layout_in_circle(first_network)
        } else if (layout_type == "fr") {
          graph_layout <- layout_with_fr(first_network)  
        } else if (layout_type == "kk") {
          graph_layout <- layout_with_kk(first_network)
        } else if (layout_type == "grid") {
          graph_layout <- tryCatch({
            layout_on_grid(first_network)
          }, error = function(e) {
            layout_with_fr(first_network)
          })
        }
        
        # Plot the network
        plot(first_network,
             layout = graph_layout,
             vertex.color = V(first_network)$color,
             vertex.size = V(first_network)$size,
             vertex.label = V(first_network)$name,
             vertex.label.cex = 0.6,
             vertex.label.color = "black",
             vertex.label.dist = 0.8,
             vertex.frame.color = "white",
             edge.width = abs(E(first_network)$weight) * 2,
             edge.color = adjustcolor("gray60", alpha.f = 0.6),
             main = paste(toupper(layout_type), "Layout"))
      }
    }
    
    par(mfrow = c(1, 1))
    mtext("Network Layout Gallery", outer = TRUE, cex = 1.5, font = 2)
  })
  
  output$networkDashboardPlot <- renderPlot({
    req(analysis_results$global_metrics)
    # Get experimental group colors
    exp_group_colors <- NULL
    if(!is.null(ui_state$group_colors) && length(ui_state$group_colors) > 0) {
      exp_group_colors <- ui_state$group_colors
    }
    render_network_dashboard(analysis_results$global_metrics, exp_group_colors)
  })
  
  output$nodeCentralityPlot <- renderPlot({
    req(analysis_results$brain_area_metrics)
    render_node_centrality(analysis_results$brain_area_metrics, ui_state$brain_areas, ui_state$area_colors)
  })
  
  output$nodeHeatmapPlot <- renderPlot({
    req(analysis_results$brain_area_metrics)
    render_node_heatmap(analysis_results$brain_area_metrics, ui_state$brain_areas, ui_state$area_colors)
  })
  
  # Individual node analysis plots
  output$individualNodeCentralityPlot <- renderPlot({
    req(analysis_results$node_metrics)
    render_individual_node_centrality(analysis_results$node_metrics, ui_state$brain_areas, ui_state$area_colors)
  })
  
  output$individualNodeHeatmapPlot <- renderPlot({
    req(analysis_results$node_metrics)
    render_individual_node_heatmap(analysis_results$node_metrics, ui_state$brain_areas, ui_state$area_colors)
  })
  
  # Additional plots for enhanced metrics
  output$edgeMetricsPlot <- renderPlot({
    if(!is.null(analysis_results$edge_metrics)) {
      # Get experimental group colors
      exp_group_colors <- NULL
      if(!is.null(ui_state$group_colors) && length(ui_state$group_colors) > 0) {
        exp_group_colors <- ui_state$group_colors
      }
      render_edge_metrics(analysis_results$edge_metrics, exp_group_colors)
    }
  })
  
  output$brainAreaConnectivityPlot <- renderPlot({
    if(!is.null(analysis_results$brain_area_metrics)) {
      render_brain_area_connectivity(analysis_results$brain_area_metrics, ui_state$brain_areas, ui_state$area_colors)
    }
  })
  
  output$conservationStatsPlot <- renderPlot({
    req(analysis_results$similarities)
    render_conservation_stats(analysis_results$similarities)
  })
  
  output$networkSimilarityMatrixPlot <- renderPlot({
    req(analysis_results$similarities)
    render_network_similarity_matrix(analysis_results$similarities)
  })
  
  output$hubConservationPlot <- renderPlot({
    req(analysis_results$conservation)
    if(analysis_results$conservation$conservation_possible) {
      render_hub_conservation(analysis_results$conservation$hub_conservation)
    }
  })
  
  output$summaryDashboardPlot <- renderPlot({
    req(analysis_results$complete)
    render_summary_dashboard(analysis_results)
  })
  
  # NEW RESTRUCTURED PIPELINE OUTPUTS
  
  # Step 3: Topology Analysis - Group-specific percolation
  output$groupPercolationPlot <- renderPlot({
    req(analysis_results$correlations)
    if (!is.null(analysis_results$group_thresholds)) {
      render_group_specific_percolation(analysis_results$correlations, analysis_results$group_thresholds)
    } else {
      plot(1, type = "n", main = "Group-specific percolation analysis not available")
    }
  })
  
  output$thresholdComparisonPlot <- renderPlot({
    req(analysis_results$correlations)
    if (!is.null(analysis_results$group_thresholds)) {
      render_threshold_comparison(analysis_results$correlations, analysis_results$global_threshold, analysis_results$group_thresholds)
    } else {
      plot(1, type = "n", main = "Threshold comparison not available")
    }
  })
  
  # Step 4: Weighted Network Analysis
  output$weightedNodeMetricsPlot <- renderPlot({
    req(analysis_results$threshold_free_results)
    render_weighted_node_metrics(analysis_results$threshold_free_results, ui_state$group_colors)
  })
  
  output$allWeightedStatsPlot <- renderPlot({
    req(analysis_results$threshold_free_results)
    render_all_weighted_stats(analysis_results$threshold_free_results, ui_state$group_colors)
  })
  
  
  # Step 5: Cross Method Comparison
  output$crossMethodEigenvectorPlot <- renderPlot({
    req(analysis_results$weighted_eigenvector)
    req(analysis_results$node_metrics)
    render_cross_method_eigenvector_comparison(analysis_results$weighted_eigenvector, 
                                             analysis_results$node_metrics, ui_state$group_colors)
  })
  
  output$crossMethodNodeStrengthPlot <- renderPlot({
    req(analysis_results$threshold_free_results)
    req(analysis_results$node_metrics)
    render_cross_method_node_strength_comparison(analysis_results$threshold_free_results$node_strength,
                                                analysis_results$node_metrics, ui_state$group_colors)
  })
  
  output$methodCorrelationPlot <- renderPlot({
    if (!is.null(analysis_results$cross_method_comparison)) {
      render_method_correlation_analysis(analysis_results$cross_method_comparison)
    } else {
      plot(1, type = "n", main = "Cross-method correlation analysis not available")
    }
  })
  
  output$methodAgreementPlot <- renderPlot({
    if (!is.null(analysis_results$cross_method_comparison)) {
      render_method_agreement_analysis(analysis_results$cross_method_comparison)
    } else {
      plot(1, type = "n", main = "Method agreement analysis not available")
    }
  })
  
  output$hubRankingComparisonPlot <- renderPlot({
    if (!is.null(analysis_results$cross_method_comparison)) {
      render_hub_ranking_comparison(analysis_results$cross_method_comparison, ui_state$group_colors)
    } else {
      plot(1, type = "n", main = "Hub ranking comparison not available")
    }
  })
  
  # Step 6: Pipeline Overview
  output$pipelineOverviewPlot <- renderPlot({
    req(analysis_results$complete)
    render_pipeline_overview(analysis_results)
  })
  
  # Advanced Analysis Tab Outputs
  output$advancedMstMetricsPlot <- renderPlot({
    if (!is.null(analysis_results$mst_results)) {
      render_advanced_mst_metrics(analysis_results$mst_results, ui_state$group_colors)
    } else {
      plot(1, type = "n", main = "MST analysis not available")
    }
  })
  
  output$mstCentralNodesPlot <- renderPlot({
    if (!is.null(analysis_results$mst_results)) {
      render_mst_central_nodes(analysis_results$mst_results, ui_state$brain_areas, ui_state$area_colors, ui_state$group_colors)
    } else {
      plot(1, type = "n", main = "MST central nodes analysis not available")
    }
  })
  
  output$advancedMstNetworksPlot <- renderPlot({
    if (!is.null(analysis_results$mst_results)) {
      render_advanced_mst_networks(analysis_results$mst_results, analysis_results$correlations,
                                 ui_state$brain_areas, ui_state$area_colors, ui_state$group_colors)
    } else {
      plot(1, type = "n", main = "MST networks not available")
    }
  })
  
  output$pcaAnalysisPlot <- renderPlot({
    if (!is.null(analysis_results$pca_results)) {
      render_pca_analysis_results(analysis_results$pca_results, ui_state$group_colors)
    } else {
      plot(1, type = "n", main = "PCA analysis not available")
    }
  })
  
  output$pcaLoadingsPlot <- renderPlot({
    if (!is.null(analysis_results$pca_results)) {
      render_pca_loadings_results(analysis_results$pca_results, ui_state$brain_areas, ui_state$area_colors, ui_state$group_colors)
    } else {
      plot(1, type = "n", main = "PCA loadings not available")
    }
  })
  
  output$pcaVariancePlot <- renderPlot({
    if (!is.null(analysis_results$pca_results)) {
      render_pca_variance_results(analysis_results$pca_results, ui_state$group_colors)
    } else {
      plot(1, type = "n", main = "PCA variance analysis not available")
    }
  })
  
  # NEW ADDITIONAL OUTPUTS FOR REQUESTED FEATURES
  
  # Regional connectivity analysis
  output$regionalConnectivityPlot <- renderPlot({
    req(analysis_results$weighted_eigenvector)
    req(analysis_results$threshold_free_results)
    render_regional_connectivity_analysis(analysis_results$weighted_eigenvector, 
                                         analysis_results$threshold_free_results$node_strength,
                                         ui_state$brain_areas, ui_state$group_colors)
  })
  
  output$regionalSummaryPlot <- renderPlot({
    req(analysis_results$weighted_eigenvector)
    req(analysis_results$threshold_free_results)
    render_regional_summary_stats(analysis_results$weighted_eigenvector,
                                 analysis_results$threshold_free_results$node_strength,
                                 ui_state$brain_areas, ui_state$group_colors)
  })
  
  # MST network visualization
  
  # Weighted vs Percolation comparison by group
  output$weightedVsPercolationEigenvectorPlot <- renderPlot({
    req(analysis_results$weighted_eigenvector)
    req(analysis_results$node_metrics)
    render_weighted_vs_percolation_eigenvector_by_group(analysis_results$weighted_eigenvector,
                                                       analysis_results$node_metrics,
                                                       ui_state$group_colors)
  })
  
  output$weightedVsPercolationNodeStrengthPlot <- renderPlot({
    req(analysis_results$threshold_free_results)
    req(analysis_results$node_metrics)
    render_weighted_vs_percolation_node_strength_by_group(analysis_results$threshold_free_results$node_strength,
                                                         analysis_results$node_metrics,
                                                         ui_state$group_colors)
  })
  
  # Weighted Eigenvector Centrality Outputs
  output$weightedEigenvectorComparison <- renderPlot({
    req(analysis_results$weighted_eigenvector)
    render_weighted_eigenvector_comparison(analysis_results$weighted_eigenvector, 
                                          ui_state$group_colors)
  })
  
  output$weightedVsUnweightedPlot <- renderPlot({
    req(analysis_results$weighted_eigenvector)
    render_weighted_vs_unweighted_eigenvector(analysis_results$weighted_eigenvector,
                                              ui_state$group_colors)
  })
  
  output$eigenvectorRankChangePlot <- renderPlot({
    req(analysis_results$weighted_eigenvector)
    render_eigenvector_rank_change(analysis_results$weighted_eigenvector,
                                   ui_state$group_colors)
  })
  
  output$strengthEigenvectorPlot <- renderPlot({
    req(analysis_results$weighted_eigenvector)
    render_strength_eigenvector_relationship(analysis_results$weighted_eigenvector,
                                            ui_state$group_colors)
  })
  
  output$weightedEigenvectorHubPlot <- renderPlot({
    req(analysis_results$weighted_eigenvector_hubs)
    render_weighted_eigenvector_hub_comparison(analysis_results$weighted_eigenvector_hubs,
                                               ui_state$group_colors)
  })
  
  output$eigenvectorStabilityPlot <- renderPlot({
    req(analysis_results$weighted_eigenvector_comparison)
    render_eigenvector_stability(analysis_results$weighted_eigenvector_comparison)
  })
  
  output$weightedEigenvectorTable <- DT::renderDataTable({
    req(analysis_results$weighted_eigenvector)
    
    # Format the data for display
    display_data <- analysis_results$weighted_eigenvector
    
    # Helper function to safely round numeric columns
    safe_round <- function(x, digits = 4) {
      tryCatch({
        numeric_x <- as.numeric(x)
        if(all(is.na(numeric_x))) return(x)  # Return original if all NA
        round(numeric_x, digits)
      }, error = function(e) {
        return(x)  # Return original on error
      })
    }
    
    # Safely round numeric columns that exist
    if("Weighted_Eigenvector" %in% names(display_data)) {
      display_data$Weighted_Eigenvector <- safe_round(display_data$Weighted_Eigenvector, 4)
    }
    if("Node_Strength" %in% names(display_data)) {
      display_data$Node_Strength <- safe_round(display_data$Node_Strength, 4)
    }
    if("Avg_Edge_Weight" %in% names(display_data)) {
      display_data$Avg_Edge_Weight <- safe_round(display_data$Avg_Edge_Weight, 4)
    }
    if("Rank_Change" %in% names(display_data)) {
      display_data$Rank_Change <- safe_round(display_data$Rank_Change, 2)
    }
    
    DT::datatable(display_data,
                  options = list(pageLength = 25,
                                scrollX = TRUE,
                                order = list(list(2, 'desc'))),  # Sort by Weighted_Eigenvector
                  rownames = FALSE) %>%
      DT::formatStyle('Weighted_Eigenvector',
                      backgroundColor = styleInterval(c(0.5, 0.8), 
                                                     c('white', 'lightblue', 'darkblue')),
                      color = styleInterval(c(0.8), c('black', 'white'))) %>%
      if("Rank_Change" %in% names(display_data)) {
        DT::formatStyle('Rank_Change',
                        color = styleInterval(c(-0.5, 0.5), 
                                             c('red', 'gray', 'green')))
      }
  })
  
  
  # Navigation from results back to import
  observeEvent(input$goToDataImportFromResults, {
    updateTabItems(session, "sidebarMenu", "import")
  })
  
  # Download template (original)
  output$download_template <- downloadHandler(
    filename = function() { "neural_plasticity_template.csv" },
    content = function(file) {
      template_data <- data.frame(
        Subject = c(1, 2, 3, 4, 5),
        Group = c("Group_A", "Group_B", "Group_A", "Group_B", "Group_A"),
        Treatment = c("Control", "Treated", "Control", "Treated", "Control"),
        Behavior = c(5.054, 9.683, 0.302, 12.451, 3.782),
        ACC = c(0.0006015, 0.0002775, 0.000265, 0.0003892, 0.0004123),
        IL = c(0.0004275, 0.0001625, 0.0001615, 0.0002341, 0.0003456),
        PRL = c(0.000383, 0.000332, 0.0001585, 0.0002876, 0.0002934),
        NaC = c(0.00010715, 0.0000655, 0.00008915, 0.00009234, 0.00010123),
        NaS = c(0.000138, 0.0000899, 0.00009745, 0.00011234, 0.00012345),
        BLA = c(0.000889216, 0.00061946, 0.000785147, 0.00072345, 0.00081234),
        LA = c(0.00083405, 0.000779674, 0.000872365, 0.00089234, 0.00091234),
        CeA = c(0.001126815, 0.001035434, 0.001039471, 0.001112345, 0.001089234),
        dDG = c(0.0008234, 0.0007123, 0.0009234, 0.0008923, 0.0007892),
        dCA1 = c(0.0007234, 0.0006789, 0.0008123, 0.0007892, 0.0008234),
        vDG = c(0.0007892, 0.0006789, 0.0008456, 0.0008123, 0.0007234)
      )
      write.csv(template_data, file, row.names = FALSE)
    }
  )
  
  # About button
  observeEvent(input$about_btn, {
    showModal(modalDialog(
      title = "About ConsensusConnectR",
      size = "l",
      div(
        h4("Multimethod Consensus-Based Functional Connectivity Analysis"),
        p("ConsensusConnectR integrates five correlation methods to provide robust preclinical functional connectivity mapping through a consensus-based approach."),
        
        h5("🔬 Correlation Methods:"),
        tags$ul(
          tags$li("Pearson correlation - Linear relationships"),
          tags$li("Spearman correlation - Monotonic relationships"),
          tags$li("Partial correlation - Direct connections controlling for confounders"),
          tags$li("Shrinkage correlation - Regularized estimates for small samples"),
          tags$li("Polychoric correlation - Ordinal data relationships")
        ),
        
        h5("📊 Key Features:"),
        tags$ul(
          tags$li("Consensus-based integration of multiple correlation methods"),
          tags$li("Percolation-based network thresholding"),
          tags$li("Functional connectivity mapping by anatomical regions"),
          tags$li("Group comparison and conservation analysis"),
          tags$li("Comprehensive network topology metrics")
        ),
        
        h5("🎯 Designed for Preclinical Research:"),
        tags$ul(
          tags$li("Optimized for small sample sizes typical in animal studies"),
          tags$li("Handles missing data through multiple imputation"),
          tags$li("Flexible data import for various experimental designs"),
          tags$li("Interactive brain region assignment and visualization")
        )
      ),
      footer = modalButton("Close")
    ))
  })
  
  # Initialize download module
  downloadServer("download_module", analysis_results, ui_state)
}

# Run the application
shinyApp(ui = ui, server = server)