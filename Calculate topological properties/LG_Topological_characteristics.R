# === 1. Load required packages ===
library(igraph)
library(tidyverse)
library(ggraph)
library(tidygraph)
library(gridExtra)

# === 2. Define a function to read GraphML files and calculate network topology ===
analyze_network <- function(file_path) {
  # Extract the network name from the file path without the extension
  net_name <- tools::file_path_sans_ext(basename(file_path))
  
  g <- read_graph(file_path, format = "graphml")
  
  # Basic network information
  node_count <- vcount(g)
  edge_count <- ecount(g)
  is_connected_net <- is.connected(g)
  components_count <- components(g)$no
  graph_size <- node_count + edge_count
  
  # Calculate topological metrics
  degree_vals <- degree(g)
  avg_degree <- mean(degree_vals)
  avg_path <- if (is_connected_net) mean_distance(g) else NA
  clustering <- transitivity(g, type = "average")
  density_val <- edge_density(g)
  modularity_obj <- cluster_louvain(g)
  modularity_val <- modularity(modularity_obj)
  module_count <- length(unique(membership(modularity_obj)))
  betweenness_mean <- mean(betweenness(g))
  
  # Visualize the network structure
  plot <- ggraph(as_tbl_graph(g)) +
    geom_edge_link(alpha = 0.2) +
    geom_node_point(aes(size = degree_vals), color = "steelblue") +
    geom_node_text(aes(label = name), repel = TRUE, size = 3) +
    theme_void() +
    ggtitle(net_name)  # Use the file name as the plot title
  
  # Return the topological metrics and network plot
  return(list(
    metrics = data.frame(
      Network = net_name,  # Use the file name as the network name
      Nodes = node_count,
      Edges = edge_count,
      Size = graph_size,
      Connected = is_connected_net,
      Components = components_count,
      Avg_Degree = avg_degree,
      Avg_Path_Length = avg_path,
      Clustering = clustering,
      Density = density_val,
      Modularity = modularity_val,
      Module_Count = module_count,
      Betweenness_Mean = betweenness_mean
    ),
    plot = plot
  ))
}

# === 3. Batch-analyze GraphML files ===
files <- list.files(pattern = "\\.graphml$", full.names = TRUE)
results <- lapply(files, analyze_network)  # Pass file paths directly

# === 4. Summarize and display results ===
all_metrics <- do.call(rbind, lapply(results, function(x) x$metrics))
all_plots <- lapply(results, function(x) x$plot)

# Print the topological metrics table
print(all_metrics)

# Save the results as a CSV file
write.csv(all_metrics, "network_topology_summary.csv", row.names = FALSE)

# Display network structure plots using the first three examples
grid.arrange(grobs = all_plots[1:3], ncol = 3)