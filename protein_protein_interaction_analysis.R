rm(list = ls())
set.seed(101)
setwd('/Users/gwk/Desktop/Thesis Figures')

# Load required packages
library(STRINGdb)
library(igraph)
library(tidyverse)

# Load data to analyse 
mut_df <- read.csv('mutant_mutant/mutExposedwtUnexposed.csv', header = T)
head(mut_df)
attach(mut_df)

df <- mut_df|>
  filter(padj < 0.05) |>
  select(padj, log2FoldChange, symbols)

colnames(df) <- c('pvalue', 'logFC', 'gene')
  
string_db <- STRINGdb$new(
  version = "12.0",
  species = 7955,
  score_threshold = 400,
  network_type = "full" # Can be 'functional' or 'physical'
)

example1 <- string_db$map(
  df, 
  "gene",
  removeUnmappedRows = TRUE
)

head(example1)

top200 <- example1$STRING_id[1:200]

string_db$plot_network(top200)

# Remove non-DE genes and add `color` column
example1_filtered <- string_db$add_diff_exp_color(
  example1[example1$pvalue < 0.05, ],
  logFcColStr = "logFC"
)

head(example1_filtered)

# Post payload information to STRING server
pid <- string_db$post_payload(
  example1_filtered$STRING_id,
  colors = example1_filtered$color
)

pid

# Visualize network with borders ("halo")
string_db$plot_network(top200, payload_id = pid)

enrichment <- string_db$get_enrichment(top200)

head(enrichment)

# Using only genes in `example1` as background
bg <- unique(example1$STRING_id)

string_db$set_background(bg)

# Get functional annotation
annot <- string_db$get_annotations(top200)

head(annot)

# Find clusters using only the first 600 genes
clusters <- string_db$get_clusters(example1$STRING_id[1:600])

# Check first 4 clusters
clusters[1:4]

# Plotting first 4 clusters
par(mfrow = c(2,2))
for(x in seq_len(4)) {
  string_db$plot_network(clusters[[x]])
}

par(mfrow = c(1,1)) # back to original number of rows and cols
# Get a list of all proteins
all_proteins <- string_db$get_proteins()

head(all_proteins)

# Get STRING ID for proteins 'tp53' and 'atm'
tp53 <- string_db$mp("tp53")
atm <- string_db$mp("atm")

list(tp53, atm)

# Get neighbors of 'tp53'
string_db$get_neighbors(tp53) |> head()

