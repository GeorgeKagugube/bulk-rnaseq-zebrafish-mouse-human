## Clear the environmental space here to start affresh 
## Clear environmental space
rm(list = ls())

set.seed(101)

## Set the directory with the expression datasets here 
setwd('/Users/gwk/Desktop/Bioinformatics/bulk-rnaseq-zebrafish-mouse-human/data')

# This is the output folder for the final analysis
output_dir = '/Users/gwk/Desktop/Thesis Figures'

## Source user defined files containing some useful functions here
source("/Users/gwk/Desktop/Bioinformatics/bulk-rnaseq-zebrafish-mouse-human/scripts/functionsneededforanalysis.R")

## Load the datasets to be analysed here 
countMatrix <- read.csv('raw_count_matrix.csv', row.names = 1)
samples <- read.csv("sample_information.csv", row.names = 1)

## Split the datasets into WT and Home
wt_sample_info <- samples |>
  filter(Genotype == 'wt')
wtcountmtx <- countMatrix |>
  select(row.names(wt_sample_info))

mut_sample_info <- samples |>
  filter(Genotype == 'mutant')
mutcountmtx <- countMatrix |>
  select(row.names(mut_sample_info))

# Build a DESeq2 object here. The count data is given as a metrix, column names
# are the samples and the rownames are the Ensemble gene ids. This is important
# The design contains what makes up the model (negatve bionimial model in this case)
dds <- DESeqDataSetFromMatrix(countData = as.matrix(countMatrix),
                              colData = samples,
                              design = ~ Group)

## Wild type only datasets, see PCA sample clustering for more insight into this choice
ddswt <- DESeqDataSetFromMatrix(countData = as.matrix(wtcountmtx),
                                colData = wt_sample_info,
                                design = ~ Group)

## Mutant only datasets, see PCA sample clustering for more insight into this choice
ddsmut <- DESeqDataSetFromMatrix(countData = as.matrix(mutcountmtx),
                                colData = mut_sample_info,
                                design = ~ Group)

##Although filtering with DESeq is not reccomended, here we remove all genes whose
## row sum is less than 10
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

keep <- rowSums(counts(ddswt)) >= 10
ddswt <- ddswt[keep,]

keep <- rowSums(counts(ddsmut)) >= 10
ddsmut <- ddsmut[keep,]

## Relevel the conditions here
dds$Group <- relevel(dds$Group, ref = 'wt_unexposed')
#dds$Group <- relevel(dds$Group, ref = 'mut_unexposed')

# Relevel the Wild type data here so the unexposed is the reference group
ddswt$Group <- relevel(ddswt$Group, ref = 'wt_unexposed')
ddsmut$Group <- relevel(ddsmut$Group, ref = 'mut_unexposed')

## Perform some quality checks here. Use the developed function here. Input is the
## dds object and in some cases you need the sample_id as the function inputs expect
## check out the individual functions below in the environmnetal space
## calculate the pca values here
vsd <- vst(ddsmut)

## Plot visualise the PCA
plotPCA(vsd, intgroup="Group")


## Normalise the count data
all_norm_data <- normalisation_func(dds)
wt_norm_data <- normalisation_func(ddswt)
mut_norm_data <- normalisation_func(ddsmut)

## calculate the principle component analysis
pca_all <- principle_component(dds)
  
# Create a data frame that can used going forward from here on 
count_all <- cbind(samples, pca_all$x)

## Visualise the data here 
ggplot(data = count_all) +
  geom_point(aes(x=PC1, y=PC2, colour = Treatment), size=5) +
  theme_minimal() +
  labs(x = 'PC1: 29% variance',
       y = 'PC2: 23% varience') +
  theme(
    axis.text = element_text(size = 20),
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    axis.title.x = element_text(size = 15, vjust = 0.5),
    axis.text.y = element_text(angle = 90, vjust = 0.5),
    axis.title.y = element_text(size = 15, vjust = 0.5)
  )
  
## Perform differenPC10## Perform differential gene expression analysis here
dds <- DESeq(dds)
ddsmut <- DESeq(ddsmut)
ddswt <- DESeq(ddswt)

res <- results(dds)
res
summary(res)
# ress_mut <- results(dds_mut)
# res_wt <- results(dds_wt)

resultsNames(dds)
# resultsNames(dds_mut)
# resultsNames(dds_wt)

## Genotype specific changes
## These are changes seen just be the mutation without the 
## Only the unexposed are explored thereby picking up changes araising from the mutation without the exposure

## Maganese responsiveness changes
## These are changes seen when WT are exposed to Mn
treatment <- results(ddswt, name = "Group_wt_exposed_vs_wt_unexposed")
treatment <- lfcShrink(dds = dds, coef = 2, type = 'apeglm')
treatment <- addDirectionlabel(treatment)
treatment <- annot_data(treatment)
treatment <- treatment[!grepl('LOC', rownames(treatment)), ]
#write.csv(treatment,
#          file = paste0(output_dir,'/mnExposureEffect/mnExposuredge.csv'))
volcanoPlot(treatment)#, xlimlimit = c(-2.5, 8), ylimlimit = c(0, 15))

## Effect of mutation (wt_unexposed vs mutant unexposed)
mutation_specific_changes <- results(ddsmut, name = "Group_mut_exposed_vs_mut_unexposed")
mutation_specific_changes <- lfcShrink(dds = ddsmut, coef = 2, type = 'apeglm')
mutation_specific_changes <- addDirectionlabel(mutation_specific_changes)
mutation_specific_changes <- annot_data(mutation_specific_changes)
mutation_specific_changes <- mutation_specific_changes[!grepl('LOC', rownames(mutation_specific_changes)), ]
write.csv(mutation_specific_changes,
          file = paste0(output_dir,'/mutant_mutant/mnExposuredge.csv'))
volcanoPlot(mutation_specific_changes, xlimlimit = c(-2.5, 8), ylimlimit = c(0, 15))

## Interaction between treatment and genotype
## These are comparing unexposed WT with exposed mutations
interactions <- results(dds, name = "Group_mut_exposed_vs_wt_unexposed")
interactions <- lfcShrink(dds = dds, coef = 2, type = 'apeglm')
interactions <- addDirectionlabel(interactions)
interactions <- annot_data(interactions)
interactions <- interactions[!grepl('LOC', rownames(interactions)), ]
#write.csv(interactions,
#          file = paste0(output_dir,'/mnmutInteraction/mutExposedwtUnexposed.csv'))
volcanoPlot(interactions)#, xlimlimit = c(-2.5, 8), ylimlimit = c(0, 15))

## These are comparing unexposed WT with exposed mutations
mutationOnly <- results(dds, name = "Group_mut_unexposed_vs_wt_unexposed")
mutationOnly <- lfcShrink(dds = dds, coef = 3, type = 'apeglm')
mutationOnly <- addDirectionlabel(mutationOnly)
mutationOnly <- annot_data(mutationOnly)
mutationOnly <- mutationOnly[!grepl('LOC', rownames(mutationOnly)), ]
#write.csv(mutationOnly,
#          file = paste0(output_dir,'/mutantEffects/mutunexposedwtUnexposed.csv'))
volcanoPlot(mutationOnly, xlimlimit = c(-5, 8), ylimlimit = c(0, 14))


df = as.data.frame(mutation_specific_changes)
# Filter significant genes
sig_genes <- df %>%
  filter(padj < 0.05 & diffExpression != "NO")

# Define top genes
top10_up <- sig_genes %>%
  arrange(desc(log2FoldChange)) %>%
  head(10)

top20_up <- sig_genes %>%
  arrange(desc(log2FoldChange)) %>%
  head(20)

ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = diffExpression)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("NS" = "grey", "Upregulated" = "firebrick", "Downregulated" = "steelblue")) +
  geom_text_repel(data = top10_up, aes(label = symbols), size = 3, box.padding = 0.5, max.overlaps = 10) +
  labs(title = "Figure 3.1: Volcano Plot",
       x = "log2 Fold Change",
       y = "-log10 Adjusted P-value") +
  theme_minimal()

ggplot(top10_up, aes(x = reorder(symbols, log2FoldChange), y = log2FoldChange)) +
  geom_col(fill = "firebrick") +
  coord_flip() +
  labs(title = "Figure 3.2: Top 10 Upregulated Genes",
       x = "Gene",
       y = "log2 Fold Change") +
  theme_minimal()

mat <- matrix(top20_up$log2FoldChange, nrow = 1)
colnames(mat) <- top20_up$symbols
rownames(mat) <- "log2FC"

pheatmap(mat,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         main = "Figure 3.3: Heatmap of Top 20 Upregulated Genes",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))


## Venn Diagram 
set1 = sig_gene_names(treatment)
set2 <- sig_gene_names(mutation_specific_changes)
set3 <- sig_gene_names(mutationOnly)
set4 <- sig_gene_names(interactions)
dgset <- list(#'WT+Mn' = set1,
              'Mutation' = set3,
              'Interaction' = set4)

venPlot(dgeset = dgset)

dge_mn <- mutation_specific_changes |>
  as.data.frame() |>
  filter(padj < 0.05) |>
  arrange(desc(log2FoldChange))

write.csv(dge_mn,
          file = paste0(output_dir,'/signifcant_genes/mutExpo_mutUnexpo_dge.csv'))

## Generate a heatmap here 
signs <- sig_gene(treatment)
hmap <- data_heatmap(as.data.frame(wt_norm_data),wt_sample_info)
figure_heatmap(hmap, signs)

## Generate a heatmap here 
signs <- sig_gene(mutationOnly)
hmap <- data_heatmap(as.data.frame(mut_norm_data),wt_sample_info)
figure_heatmap(hmap, signs)

## Generate a heatmap here 
signs <- sig_gene(treatment)
hmap <- data_heatmap(as.data.frame(wt_norm_data),wt_sample_info)
figure_heatmap(hmap, signs)

## Generate a heatmap here 
signs <- sig_gene(treatment)
hmap <- data_heatmap(as.data.frame(wt_norm_data),wt_sample_info)
figure_heatmap(hmap, signs)

mutation_specific_changes |>
  as.data.frame() |>
  filter(padj <= 0.05)
