# === 1. Load required packages ===
library(tidyverse)
library(igraph)
library(tidygraph)
library(ggraph)
library(readr)

# === 2. Read data ===
df_raw <- read_csv("LG2005.csv", skip = 2)  # Skip the first two rows and read only the data
bin_info <- read_csv("LG2005.csv", n_max = 1, col_names = FALSE)  # The first row contains bin names

# === 3. Set column names ===
colnames(df_raw)[1:3] <- c("KO", "Gene", "Function")
bin_names <- as.character(bin_info[1, 4:ncol(bin_info)])
names(bin_names) <- colnames(df_raw)[4:ncol(df_raw)]

# Replace column names with bin names
colnames(df_raw)[4:ncol(df_raw)] <- bin_names
bin_cols <- colnames(df_raw)[4:ncol(df_raw)]

# === 3.5 Remove KO rows with total TPM equal to 0 ===
df_filtered <- df_raw %>%
  mutate(total_TPM = rowSums(across(all_of(bin_cols)))) %>%
  filter(total_TPM > 0) %>%
  select(-total_TPM)

# === 4. Build edge list (KO -> Bin) ===
edge_df <- df_filtered %>%
  pivot_longer(cols = all_of(bin_cols), names_to = "Bin", values_to = "TPM") %>%
  filter(TPM > 0) %>%
  select(KO, Bin, TPM, Function)

# === 5. Build nodes ===
nodes_ko <- df_filtered %>%
  select(KO, Gene, Function) %>%
  distinct() %>%
  rename(name = KO) %>%
  mutate(type = "KO", label = Gene)

nodes_bin <- tibble(name = bin_cols) %>%
  mutate(Function = NA, type = "Bin", label = name)

nodes_all <- bind_rows(nodes_ko, nodes_bin)

# === 6. Rebuild edge list for consistency ===
edges <- edge_df %>%
  select(from = KO, to = Bin, weight = TPM)

# === 7. Build network ===
g <- graph_from_data_frame(d = edges, vertices = nodes_all, directed = FALSE)

# === 8. Save the .graphml file ===
write_graph(g, file = "LG2005_bipartite_network.graphml", format = "graphml")

# === 9. Convert to tidygraph and plot ===
tg <- g %>%
  as_tbl_graph() %>%
  mutate(group = ifelse(type == "Bin", "Bin", Function))

# === 10. Export PDF figure ===
pdf("LG2005_bipartite_network.pdf", width = 12, height = 8)

ggraph(tg, layout = "fr") +
  geom_edge_link(aes(width = weight), alpha = 0.3) +
  geom_node_point(aes(color = group, shape = type), size = 4) +
  geom_node_text(aes(label = label), repel = TRUE, size = 4, max.overlaps = 100) +
  scale_shape_manual(values = c("KO" = 16, "Bin" = 15)) +
  theme_void() +
  theme(legend.position = "right") +
  labs(title = "LG2005 Gene-Bin Bipartite Network")

dev.off()