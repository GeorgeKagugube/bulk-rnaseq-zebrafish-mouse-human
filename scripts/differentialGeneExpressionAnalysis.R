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


## Optional drop the outliers samples here
#samples <- samples[-c(2,12),]
#countMatrix <- countMatrix[,-c(11,12)]

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

## Relevel the conditions here
dds$Group <- relevel(dds$Group, ref = 'wt_unexposed')
## Perform some quality checks here. Use the developed function here. Input is the
## dds object and in some cases you need the sample_id as the function inputs expect
## check out the individual functions below in the environmnetal space
## calculate the pca values here
vsdmut <- vst(dds_mut)
vsdwt <- vst(dds_wt)
vsd <- vst(dds)

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
  geom_point(aes(x=PC1, y=PC2, colour = Group), size=5) +
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
  
## Perform differential gene expression analysis here
dds <- DESeq(dds)
dds_mut <- DESeq(dds_mut)
dds_wt <- DESeq(dds_wt)

res <- results(dds)
ress_mut <- results(dds_mut)
res_wt <- results(dds_wt)

resultsNames(dds)
resultsNames(dds_mut)
resultsNames(dds_wt)

## Set up group comparisions here 
res_wt_mn <- results(dds, contrast = c('Group','wt_unexposed','wt_exposed'))
res_wt_mn <- lfcShrink(dds, contrast = c('Group','wt_unexposed','wt_exposed'),  type = 'ashr')

## Baseline mutant effect on transcriptomics
res_controls <- results(dds, contrast = c('Group', 'wt_unexposed','mut_unexposed'))
res_controls <- lfcShrink(dds, contrast = c('Group', 'wt_unexposed','mut_unexposed'),  type = 'ashr')

## Manganese effect in mutant
mn_effect_mut <- results(dds, contrast = c('Group','mut_unexposed', 'mut_exposed'))
mn_effect_mut <- lfcShrink(dds, contrast = c('Group','mut_unexposed', 'mut_exposed'),  type = 'ashr')

## Interaction between genotype and manganese
mn_genotype_interaction <- results(dds, contrast = c('Group', 'wt_exposed', 'mut_exposed'))
mn_genotype_interaction <- lfcShrink(dds, contrast = c('Group','wt_exposed', 'mut_exposed'),  type = 'ashr')


dge1 <- list('wt exposed' = wt,
             'mut Exposed' = mut_wt)
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

## Export data to a csv file for further analysis and inspection later
## Export the data 
write.csv(resLFC,
          file = paste0(output_folder,'Group_mut_exposed_vs_mut_unexposed.csv'))

