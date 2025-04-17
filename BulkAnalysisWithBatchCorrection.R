## Clear the environmental space here to start affresh 
rm(list = ls())

## Set up the analysis to run multiple cores, here we use 4 cores after removing lowely 
# expressed genes
library("BiocParallel")
register(MulticoreParam(4))

## Set the directory with the expression datasets here 
## Re-analysed samples with the new reference genome
setwd('/Users/gwk/Desktop/PhD/Data/PhD_data/Brain/GZ11_star_output/star')

## This is the output folder for the final analysis
output_folder = '/Users/gwk/Desktop/PhD/Data/PhD_data/March_03_25_Final_Analysis/DGE_Files/15-04-25/'
output_dir = '/Users/gwk/Desktop/PhD/Data/PhD_data/March_03_25_Final_Analysis/normalised_data_sets'

## Source user defined files containing some useful functions here
source("/Users/gwk/Desktop/Bioinformatics/bulk-rnaseq-zebrafish-mouse-human/scripts/functionsneededforanalysis.R")

## Set up variables here 
## Read the mapped and quantified files into r for further analysis here 
file1 <- "1_S1_ReadsPerGene.out.tab"
dr <- '/Users/gwk/Desktop/PhD/Data/PhD_data/Brain/GZ11_star_output/star'

## load the data into R for further analysis here 
# This is output from star - these are files with ReadsPerGene.out.tab
countMatrix <- staroutput_preprocessing(dr, file1)

## Read in the data from the sample informatio here 
# If no sample information exists, you can create to match the countMatrix
samples <- read.csv("sample_information.csv", row.names = 1)

## Explore the loaded data here
head(countMatrix)
head(samples, 15)

## Remove outliers from the dataset here
samples <- samples[-c(2,12),]
countMatrix <- countMatrix[,-c(11,12)]

## Fix the column names in the count matrix to match the sample information rownames
row_col_names <- as.vector(paste0("sample",rownames(samples)))
rownames(samples) <- row_col_names
colnames(countMatrix) <- row_col_names

## Makesure that sample rownames match the countmatrix colnames here
# The lines below must return true for all
# if false, please investigate the rownames match the countMatrix colnames
all(row.names(samples) %in% colnames(countMatrix))
all(row.names(samples) == colnames(countMatrix))

samples$Group <- as.factor(samples$Group)
samples$batch <- as.factor(samples$batch)

# Build a DESeq2 object here. The count data is given as a metrix, column names
# are the samples and the rownames are the Ensemble gene ids. This is important
# The design contains what makes up the model (negatve bionimial model in this case)
attach(countMatrix)
ddsE <- DESeqDataSetFromMatrix(countData = as.matrix(countMatrix),
                              colData = samples,
                              design = ~ Group
                              )

# Pre-filtration to remove genes that are lowely expressed across all samples
keep <- rowSums(counts(ddsE)) >= 10
ddsE <- dds[keep,]

## Differential gene expression
## probably the most important command :)
ddsE <- DESeq(ddsE)

## Need to correct for potential batch effects here
# Variance-stabilizing transformation (preserves design info)
vsd <- vst(ddsE, blind = FALSE)
# Remove batch effect using limma (for visualization only)
vsd_corrected <- removeBatchEffect(assay(vsd), batch = vsd$batch)

# Select top variable genes
var_genes <- apply(vsd_corrected, 1, var)
top_genes <- head(order(var_genes, decreasing = TRUE), 500)
tsne_input <- t(vsd_corrected[top_genes, ])

# Run t-SNE
set.seed(123)
tsne_result <- Rtsne(tsne_input, dims = 2, perplexity = 3)

# Create data frame for plotting
df_tsne <- data.frame(
  Sample = rownames(tsne_input),
  tSNE1 = tsne_result$Y[,1],
  tSNE2 = tsne_result$Y[,2],
  Condition = samples$Group
)

# Plot t-SNE
ggplot(df_tsne, aes(x = tSNE1, y = tSNE2, color = Condition)) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(title = "t-SNE on Batch-Corrected Expression", color = "Group")

# Heatmap of batch-corrected expression for same genes
top_matrix <- vsd_corrected[top_genes, ]

# Annotation for columns
annotation_col <- data.frame(
  Condition = samples$Group,
  Batch = samples$batch
)
rownames(annotation_col) <- rownames(samples)

# Plot heatmap
pheatmap(top_matrix,
         annotation_col = annotation_col,
         show_rownames = FALSE,
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         scale = "row",
         main = "Heatmap of Top 1000 Variable Genes (Batch Corrected)")

