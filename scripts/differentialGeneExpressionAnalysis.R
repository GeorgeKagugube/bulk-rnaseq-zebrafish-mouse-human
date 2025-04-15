## Clear the environmental space here to start affresh 
rm(list = ls())

## Set the directory with the expression datasets here 
## Re-analysed samples with the new reference genome
setwd('/Users/gwk/Desktop/PhD/Data/PhD_data/Brain/GZ11_star_output/star')

## This is the output folder for the final analysis
output_folder = '/Users/gwk/Desktop/PhD/Data/PhD_data/March_03_25_Final_Analysis'
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
all(row.names(samples) == colnames(countMatrix[,-c(11,12)]))

# Split the data so that mutant and wt 
# Mutant sample information
mut <- samples %>%
  filter(Genotype == "mutant") %>%
  dplyr::select(Group)

mutt <- c(3,4,8,9,10,11)
# Wild type sample information
wt <- samples %>%
  filter(Genotype == "wt") %>%
  dplyr::select(Group) 

## Split the matrix here to analyse each genotype seperately here
# Mutant count matrix
countMut <- countMatrix[,mutt]
names(countMut)

# WT countMatrix
countWT <- countMatrix[,-mutt]

# Check that all the rows and columns still match (sample rowname/count colnames)
all(row.names(mut) %in% colnames(countMut))
all(row.names(wt) %in% colnames(countWT))

## Change conditions to factor. This is the column or columns to be used by the
## model when buidling the object. Not Important, but if not done will return a
## warning message in the next step
mut$Group <- as.factor(mut$Group)
wt$Group <- as.factor(wt$Group)
samples$Group <- as.factor(samples$Group)
samples$batch <- as.factor(samples$batch)

# Build a DESeq2 object here. The count data is given as a metrix, column names
# are the samples and the rownames are the Ensemble gene ids. This is important
# The design contains what makes up the model (negatve bionimial model in this case)
dds <- DESeqDataSetFromMatrix(countData = as.matrix(countMatrix),
                              colData = samples,
                              design = ~ Group)

dds_mut <- DESeqDataSetFromMatrix(countData = countMut,
                                  colData = mut,
                                  design = ~ Group)

dds_wt <- DESeqDataSetFromMatrix(countData = as.matrix(countWT),
                                  colData = wt,
                                  design = ~ Group)

##Although filtering with DESeq is not reccomended, here we remove all genes whose
## row sum is less than 10
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

## Perform some quality checks here. Use the developed function here. Input is the
## dds object and in some cases you need the sample_id as the function inputs expect
## check out the individual functions below in the environmnetal space
## calculate the pca values here
vsdmut <- vst(dds_mut)
vsdwt <- vst(dds_wt)
vsd <- vst(dds)

## Run tsne analysis here
expr_matrix <- assay(vsd)

# Keep genes with high variance across samples
var_genes <- apply(expr_matrix, 1, var)
top_genes <- head(order(var_genes, decreasing = TRUE), 1000)
expr_matrix_filtered <- expr_matrix[top_genes, ]

tsne_input <- t(expr_matrix_filtered)

set.seed(123)  # For reproducibility
tsne_result <- Rtsne(tsne_input, dims = 2, perplexity = 1, verbose = TRUE)

df <- data.frame(
  Sample = rownames(tsne_input),
  tSNE1 = tsne_result$Y[,1],
  tSNE2 = tsne_result$Y[,2],
  Group = samples$Group
)

ggplot(df, aes(x = tSNE1, y = tSNE2, color = Group)) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(title = "t-SNE of Bulk RNA-seq Data", color = "Group")

## Plot visualise the PCA
plotPCA(vsd, intgroup="Group")
plotPCA(vsdmut, intgroup = 'Group')
plotPCA(vsdwt, intgroup = 'Group')

## Normalise the count data
all_norm_data <- normalisation_func(dds)
mut_exposed <- normalisation_func(dds_mut)
wt_exposed <- normalisation_func(dds_wt)

## calculate the principle component analysis
pca_mut <- principle_component(dds_mut)
pca_wt <- principle_component(dds_wt)
pca_all <- principle_component(dds)
  
# Create a data frame that can used going forward from here on 
count_mut_df <- cbind(mut,pca_mut$x)
count_wt_df <- cbind(wt, pca_wt$x)
count_all <- cbind(samples, pca_all$x)

## Visualise the data here 
ggplot(data = count_all) +
  geom_point(aes(x=PC2, y=PC8, colour = Group), size=5) +
  theme_minimal() +
  labs(x = 'PC1: 68% variance',
       y = 'PC2: 22% varience') +
  theme(
    axis.text = element_text(size = 20),
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    axis.title.x = element_text(size = 15, vjust = 0.5),
    axis.text.y = element_text(angle = 90, vjust = 0.5),
    axis.title.y = element_text(size = 15, vjust = 0.5)
  )
  
## Perform hierarchical clustering 
### Compute pairwise correlation values
rld_cor <- cor(rld_mat)    ## cor() is a base R function
rld_wt_corr <- cor(rld_wt_mat)

## Generat a heatmap plot here 
pheatmap(rld_cor)
pheatmap(rld_wt_corr)

## Differential gene expression
## probably the most important command :)
dds <- DESeq(dds)

## Need to correct for potential batch effects here
# Variance-stabilizing transformation (preserves design info)
vsd <- vst(dds, blind = FALSE)
# Remove batch effect using limma (for visualization only)
vsd_corrected <- removeBatchEffect(assay(vsd), batch = vsd$batch)

# Select top variable genes
var_genes <- apply(vsd_corrected, 1, var)
top_genes <- head(order(var_genes, decreasing = TRUE), 100)
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

dds_mut <- DESeq(dds_mut)
dds_wt <- DESeq(dds_wt)

res <- results(dds)
ress_mut <- results(dds_mut)
res_wt <- results(dds_wt)

resultsNames(dds)
resultsNames(dds_mut)
resultsNames(dds_wt)

res <- results(dds, name = "Group_wt_unexposed_vs_mut_exposed")
ress_mut <- results(dds_mut, name = "Group_mut_unexposed_vs_mut_exposed")
resLFC_mutComb <- results(dds, name = "Group_mut_unexposed_vs_mut_exposed")
ress_wt <- results(dds_wt, name = "Group_wt_unexposed_vs_wt_exposed")
ress_wtVSmutexposed <- results(dds, name = "Group_wt_exposed_vs_mut_exposed")

resLFC <- lfcShrink(dds, coef="Group_wt_unexposed_vs_mut_exposed", 
                    type="apeglm")
resLFC_mut <- lfcShrink(dds_mut, coef="Group_mut_unexposed_vs_mut_exposed", 
                    type="apeglm")
resLFC_mutComb <- lfcShrink(dds, coef = "Group_mut_unexposed_vs_mut_exposed",
                            type = 'apeglm')
resLFC_wt <- lfcShrink(dds_wt, coef="Group_wt_unexposed_vs_wt_exposed", 
                    type="apeglm")
ress_wtVSmutexposed <- lfcShrink(dds, coef = "Group_wt_exposed_vs_mut_exposed")
summary(resLFC)
resLFC_mut <- na.omit(resLFC_mut)
resLFC_wt <- na.omit(resLFC_wt)

dim(resLFC_wt %>%
  as.data.frame() %>%
  filter(padj < 0.05))
lfc <- na.omit(resLFC)

createObjectSet <- function(df, pvalue=0.05){
  if (is.data.frame(df)){
    set1 <- df %>%
      filter(padj < pvalue) %>%
      rownames()
  }else{
    set1 <- df %>%
      as.data.frame() %>%
      filter(padj < pvalue) |>
      rownames()
  }
  
  return(set1)
}

mutOnly <- resLFC_mut %>%
  as.data.frame() %>%
  filter(padj < 0.05) %>%
  rownames()

mutComb <- resLFC_mutComb %>%
  as.data.frame() %>%
  filter(padj < 0.05) %>%
  rownames()

mut_wt <- resLFC %>%
  as.data.frame() %>%
  filter(padj < 0.05) %>%
  rownames()

mute_wte <- ress_wtVSmutexposed %>%
  as.data.frame() %>%
  filter(padj < 0.05) %>%
  rownames()

wt <- resLFC_wt %>%
  as.data.frame() %>%
  filter(padj < 0.05) %>%
  rownames()

dge1 <- list('wt exposed' = wt,
             'mut Exposed' = mutOnly)
ggvenn(dge1)

dge2 <- list('wt exposed' = wt,
             'mut exposed' = mutComb)
ggvenn(dge2)
dge <- list('wt_mut' = rownames(lfc[lfc$padj < 0.05,]),
            'mut_exp' = rownames(resLFC_mut[resLFC_mut$padj < 0.05,]),
            'wt_expo' = rownames(resLFC_wt[resLFC_wt$padj < 0.05,]))

wt[wt %in% mutComb]
ggvenn(dge)

res05 <- results(dds_wt, alpha=0.05)
summary(res05)
res05[res05$padj < 0.05]
res05 %>%
  filter(padj < 0.05)
diff_expression_func = function(data_obj, condition_ = 'mut_unexposed',
                                test_comparision = 'Group_mut_exposed_vs_mut_unexposed',
                                coefficient = 2){ 
  data_obj$Group <- relevel(data_obj$Group, ref = condition_)
  data_obj <- DESeq(data_obj)
  
  ## 
  res_data_condition <- results(data_obj,
                                name = test_comparision,
                                alpha = 0.05)
  ## Shrink the data object to extract log fold changes 
  res_data_shrink <- lfcShrink(dds = data_obj, 
                               coef = coefficient, 
                               res = res_data_condition)
  
  ## Annotate the dat with entrezi and Ensembl gene IDs here
  annotated_data <- annot_data(res_data_shrink)
  
  ## Return the annotated data set for downstreem analysis
  return(annotated_data)
}


## Annotat the data before exporting the data 
annot_data <- function(data, org = org.Dr.eg.db) {
  data$entrezid <- mapIds(org,
                          keys = row.names(data),
                          column = c("ENTREZID"),
                          keytype = "SYMBOL",
                          multiVals = "first")
  data$esembl_id <- mapIds(org,
                           keys = row.names(data),
                           column = c("ENSEMBL"),
                           keytype = "SYMBOL",
                           multiVals = "first")
  
  return(data)
}

res_all_wt_exposed_shrink <- annot_data(res_all_wt_exposed_shrink)

## Export data to a csv file for further analysis and inspection later
## Export the data 
write.csv(res_all_wt_exposed_shrink,
          file = paste0(output_folder,'/wt_exposed.csv'))

################################################################################
############. Heatmap Plot of the DEG list here ################################
################################################################################
sig_gene <- function(datafframe, basemean = 100, lfc = 1.0) {
  sigs <- datafframe %>%
    as.data.frame() %>%
    filter(padj < 0.05)
  
  ## sig
  sig_df <- sigs[(sigs$baseMean > basemean) & abs((sigs$log2FoldChange >= lfc)),]
  sig_df <- na.omit(sig_df)
  
  ## Return value from the data frame here
  return(sig_df)
}

sig_df <- sig_gene(res_wt_shrink)

data_heatmap <- function(normalised_data, sample){
  sig_df <- sig_gene(res_wt_shrink)
  ## Return the z-score of the matrix here
  normalised_data <- normalised_data[rownames(sig_df),]
  mut_mat.z <- t(apply(normalised_data, 1, scale))
  colnames(mut_mat.z) <- sample$Group
  
  return(mut_mat.z)
  
}
data_heatmap(normalised_counts_wt, wt)
