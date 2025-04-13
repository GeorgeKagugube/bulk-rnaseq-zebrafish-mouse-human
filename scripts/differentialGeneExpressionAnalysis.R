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

## Fix the column names in the count matrix to match the sample information rownames
row_col_names <- as.vector(paste0("sample",rownames(samples)))
rownames(samples) <- row_col_names
colnames(countMatrix) <- row_col_names

## Makesure that sample rownames match the countmatrix colnames here
# The lines below must return true for all
# if false, please investigate the rownames match the countMatrix colnames
all(row.names(samples) %in% colnames(countMatrix))
all(row.names(samples) == colnames(countMatrix))

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

## Perform some quality checks here. Use the developed function here. Input is the
## dds object and in some cases you need the sample_id as the function inputs expect
## check out the individual functions below in the environmnetal space
pca_plot(dds_mut, mut)
pca_plot(dds_wt, wt)
pca_plot(dds, samples)

vsdmut <- vst(dds_mut)

vsdwt <- vst(dds_wt)

## calculate the pca values here
vsd <- vst(dds)
plotPCA(vsd, intgroup="Group")
plotPCA(vsdmut, intgroup = 'Group')
plotPCA(vsdwt, intgroup = 'Group')

all_norm_data <- normalisation_func(dds)
mut_exposed <- normalisation_func(dds_mut)
wt_exposed <- normalisation_func(dds_wt)

## Export the csv file for further analysis and interogation
write.csv(wt_exposed,
          file = paste0(output_dir,'/wt_exposed.csv'))

pca_mut <- principle_component(dds_mut)
pca_wt <- principle_component(dds_wt)
pca_all <- principle_component(dds)
  
# Create a data frame that can used going forward from here on 
count_mut_df <- cbind(mut,pca_mut$x)
count_wt_df <- cbind(wt, pca_wt$x)
count_all <- cbind(samples, pca_all$x)

## Visualise the data here 
ggplot(data = count_all) +
  geom_point(aes(x=PC1, y=PC2, colour = Group), size=5) +
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
