## Clear the environmental space here to start affresh 
## Clear environmental space
rm(list = ls())

set.seed(101)

## Set the directory with the expression datasets here 
setwd('/Users/gwk/Desktop/Bioinformatics/bulk-rnaseq-zebrafish-mouse-human/data')

# This is the output folder for the final analysis
output_dir = '/Users/gwk/Desktop/Thesis Figures/DifferentialGeneExpression/'
output_dir = '/Users/gwk/Desktop/Thesis Figures/'

## Source user defined files containing some useful functions here
source("/Users/gwk/Desktop/Bioinformatics/bulk-rnaseq-zebrafish-mouse-human/scripts/functionsneededforanalysis.R")

## Load the datasets to be analysed here 
countMatrix <- read.csv('raw_count_matrix.csv', row.names = 1)
samples <- read.csv("./samplieinformation.csv", row.names = 1)

## Split the datasets into WT and Home
wt_sample_info <- samples |>
  filter(Genotype == 'wt')
wtcountmtx <- countMatrix |>
  select(row.names(wt_sample_info))

mut_sample_info <- samples |>
  filter(Genotype == 'mutant')
mutcountmtx <- countMatrix |>
  select(row.names(mut_sample_info))

## Filter before creating the DESeq object
##Although filtering with DESeq is not reccomended, here we remove all genes whose
## row sum is less than 10
keep <- rowSums(countMatrix > 10) >= 6
countMatrix <- countMatrix[keep,]

keep_wt <- rowSums(wtcountmtx > 10) >= 3
wtcountmtx <- wtcountmtx[keep_wt,]

keep_mut <- rowSums(mutcountmtx > 10) >= 3
mutcountmtx <- mutcountmtx[keep_mut,]

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

## Relevel the conditions here
dds$Group <- relevel(dds$Group, ref = 'wt_unexposed')
#dds$Group <- relevel(dds$Group, ref = 'mut_unexposed')

# Relevel the Wild type data here so the unexposed is the reference group
ddswt$Group <- relevel(ddswt$Group, ref = 'wt_unexposed')
ddsmut$Group <- relevel(ddsmut$Group, ref = 'mut_unexposed')

## Perform some quality checks here. Use the developed function here. Input is the
## dds object and in some cases you need the sample_id as the function inputs expect
# ## check out the individual functions below in the environmnetal space
# ## calculate the pca values here
# vsd <- vst(ddsmut)
# 
# ## Plot visualise the PCA
# plotPCA(vsd, intgroup="Group")
# 
# ## Normalise the count data
# all_norm_data <- normalisation_func(dds)
# wt_norm_data <- normalisation_func(ddswt)
# mut_norm_data <- normalisation_func(ddsmut)
# 
# ## calculate the principle component analysis
# pca_all <- principle_component(dds)
#   
# # Create a data frame that can used going forward from here on 
# count_all <- cbind(samples, pca_all$x)
# 
# ## Visualise the data here 
# ggplot(data = count_all) +
#   geom_point(aes(x=PC1, y=PC2, colour = Group), size=5) +
#   theme_minimal() +
#   labs(x = 'PC1: 29% variance',
#        y = 'PC2: 23% varience') +
#   theme(
#     axis.text = element_text(size = 20),
#     axis.text.x = element_text(angle = 90, vjust = 0.5),
#     axis.title.x = element_text(size = 15, vjust = 0.5),
#     axis.text.y = element_text(angle = 90, vjust = 0.5),
#     axis.title.y = element_text(size = 15, vjust = 0.5)
#   )
#   
# ## Perform differenPC10## Perform differential gene expression analysis here
dds <- DESeq(dds)
ddsmut <- DESeq(ddsmut)
ddswt <- DESeq(ddswt)


mutantExposed |>
  as.data.frame() |>
  filter(padj < 0.05) |>
  nrow()
## Extract Group comparisons here
## ========================= Mutant,Exposed vs Mutant unexposed ===============
resMut <- results(ddsmut)
summary(resMut)
resultsNames(ddsmut)
mutantExposed <- results(ddsmut, name = "Group_mut_exposed_vs_mut_unexposed")
mutantExposed <- lfcShrink(ddsmut, coef = 2, type = 'apeglm')
mutantExposed <- addDirectionlabel(mutantExposed)
mutantExposed <- annot_data(mutantExposed)
summary(mutantExposed)
# Remove all genes that are not annotated here
mutantExposed <- mutantExposed[!grepl('LOC', rownames(mutantExposed)), ]
# visualise the data using a volcano plot here
volcanoPlot(mutantExposed)

## Export the the DGE here
write.csv(mutantExposed,
          file = paste0(output_dir,'/MutantExposed/mutExposed_vs_mutUnexposed.csv'))
## =======================Wild type Exposed vs Unexposed =======================
resWT <- results(ddswt)
resWT
summary(resWT)
resultsNames(ddswt)
WTExposed <- results(ddswt, name = "Group_wt_exposed_vs_wt_unexposed")
WTExposed <- lfcShrink(ddswt, coef = 2, type = 'apeglm')
WTExposed <- addDirectionlabel(WTExposed)
WTExposed <- annot_data(WTExposed)
summary(WTExposed)
# Remove all genes that are not annotated here
WTExposed <- WTExposed[!grepl('LOC', rownames(WTExposed)), ]
# visualise the data using a volcano plot here
volcanoPlot(WTExposed)

## Export the the DGE here
write.csv(WTExposed,
          file = paste0(output_dir,'/WildtypeExposed/WTExposed_vs_WTUnexposed.csv'))

## ============= Combined data analysis for interactions and mutant effect =====
res <- results(dds)
res
summary(res)
resultsNames(dds)

## ================================= Mutant only effects =======================
mutantUnexposed <- results(dds, name = "Group_mut_unexposed_vs_wt_unexposed")
mutantUnexposed <- lfcShrink(dds = dds, coef = 3, type = 'apeglm')
mutantUnexposed <- addDirectionlabel(mutantUnexposed)
mutantUnexposed <- annot_data(mutantUnexposed)
summary(mutantUnexposed)

## Remove all genes that are not annotated in the reference genome
volcanoPlot(mutantUnexposed)

## Export the datset here
write.csv(mutantUnexposed,
          file = paste0(output_dir,'/MutantUnexposed/mutUnexposed_vs_wtUnexposed.csv'))


## ===================== Interaction between genotype and exposure =============
interaction <- results(dds, name = "Group_mut_exposed_vs_wt_unexposed")
interaction <- lfcShrink(dds = dds, coef = 4, type = 'apeglm')
interaction <- addDirectionlabel(interaction)
interaction <- annot_data(interaction)
summary(interaction)

## Remove all genes that are not annotated in the reference genome
volcanoPlot(interaction)

## Export the datset here
write.csv(interaction,
          file = paste0(output_dir,'/Interactions/mutExposed_vs_mutUnexposed.csv'))


## ========================= Extracting key genes from each dataset ============
significant_genes<- mutantExposed %>%
  as.data.frame() %>%
  #dplyr::filter(padj <=0.01, abs(log2FoldChange) >= 2) %>% 
  dplyr::filter(padj <=0.01) %>% 
  rownames()


significant_genes
## ====================== Over representation analysis =======================
significant_genes_map<- clusterProfiler::bitr(geneID = significant_genes,
                                              fromType="SYMBOL", toType="ENTREZID",
                                              OrgDb="org.Dr.eg.db")

head(significant_genes_map)

## =============== Prepare background gene list ================================
## background genes are genes that are detected in the RNAseq experiment 
background_genes<- mutantExposed %>% 
  as.data.frame() %>% 
  filter(baseMean != 0) %>%
  tibble::rownames_to_column(var = "gene") %>%
  pull(gene)

background_genes_map<- bitr(geneID = background_genes, 
                            fromType="SYMBOL", 
                            toType="ENTREZID",
                            OrgDb="org.Dr.eg.db")

## ========================== GO TERM ANALYSIS =================================
ego <- enrichGO(gene          = significant_genes_map$ENTREZID,
                universe      = background_genes_map$ENTREZID,
                OrgDb         = org.Dr.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)

## ================ Visualise the genes here ============
barplot(ego, showCategory=10)
dotplot(ego)


### =================== MsigDB Hallmark gene sets ===============
# install.packages("msigdbr")
library(msigdbr)

m_df <- msigdbr(species = "zebrafish")
head(m_df)

## ========= 
# Note, the new version of msigdbr changed the arguments: category vs collection
# welcome to Bioinformatics!..

m_t2g <- msigdbr(species = "zebrafish", collection = "C4") %>% 
  dplyr::select(gs_name, ncbi_gene)


head(m_t2g)

## ===========
em <- enricher(significant_genes_map$ENTREZID, TERM2GENE=m_t2g, 
               universe = background_genes_map$ENTREZID )
head(em)

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
set1 = sig_gene_names(as.data.frame(mutantExposed))
set2 <- sig_gene_names(WTExposed)
set3 <- sig_gene_names(mutantUnexposed)
set4 <- sig_gene_names(interaction)
dgset <- list('mutantsExposed' = set1,
              'WTexposed' = set2,
              'MutantUnexposed' = set3,
              'Interaction' = set4
              )

venPlot(dgeset = dgset)

dge_mn <- interaction |>
  as.data.frame() |>
  filter(padj < 0.05) |>
  arrange(desc(log2FoldChange))

write.csv(dge_mn,
          file = paste0(output_dir,'/signifcant_genes/interaction_dge.csv'))

## Generate a heatmap here 
signs <- sig_gene(WTExposed)
wthmap <- data_heatmap(as.data.frame(wt_norm_data),wt_sample_info)
figure_heatmap(wthmap, signs)

## Generate a heatmap here 
signs <- sig_gene(mutantExposed)
hmap <- data_heatmap(as.data.frame(mut_norm_data),mut_sample_info)
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
