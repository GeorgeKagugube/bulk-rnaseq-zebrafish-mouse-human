## Clear the environmental space here to start affresh 
## Clear environmental space
rm(list = ls())

## Set the directory with the expression datasets here 
setwd('/Users/gwk/Desktop/PhD/Data/PhD_data/Brain/GZ11_star_output/star')

# This is the output folder for the final analysis
output_folder = '/Users/gwk/Desktop/PhD/Data/PhD_data/March_03_25_Final_Analysis/DGE_Files/'
output_dir = '/Users/gwk/Desktop/PhD/Data/PhD_data/March_03_25_Final_Analysis/normalised_data_sets'

## Source user defined files containing some useful functions here
source("/Users/gwk/Desktop/Bioinformatics/bulk-rnaseq-zebrafish-mouse-human/scripts/functionsneededforanalysis.R")

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

# Wild type sample information
wt <- samples %>%
  filter(Genotype == "wt") %>%
  dplyr::select(Group) 

## Split the matrix here to analyse each genotype seperately here
# Mutant count matrix
countMut <- countMatrix[,rownames(mut)]
names(countMut)

# WT countMatrix
countWT <- countMatrix[,rownames(wt)]

# Check that all the rows and columns still match (sample rowname/count colnames)
all(row.names(mut) %in% colnames(countMut))
all(row.names(wt) %in% colnames(countWT))

## Change conditions to factor. This is the column or columns to be used by the
## model when buidling the object. Not Important, but if not done will return a
## warning message in the next step
mut$Group <- as.factor(mut$Group)
wt$Group <- as.factor(wt$Group)
samples$Group <- as.factor(samples$Group)

write.csv(samples,
          '/Users/gwk/Desktop/PhD/Data/PhD_data/March_03_25_Final_Analysis/normalised_data_sets/sample_information')

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
#dds$Group <- relevel(dds$Group, ref = 'wt_unexposed')
dds$Group <- relevel(dds$Group, ref = 'mut_unexposed')

## Perform some quality checks here. Use the developed function here. Input is the
## dds object and in some cases you need the sample_id as the function inputs expect
## check out the individual functions below in the environmnetal space
## calculate the pca values here
vsd <- vst(dds)

## Plot visualise the PCA
plotPCA(vsd, intgroup="Group")


## Normalise the count data
all_norm_data <- normalisation_func(dds)

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
#dds_mut <- DESeq(dds_mut)
#dds_wt <- DESeq(dds_wt)

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
mutation_specific <- results(dds, name =  "Group_mut_unexposed_vs_wt_unexposed")
mutation_specific <- lfcShrink(dds = dds, coef = 3, type = 'apeglm')

## Maganese responsiveness changes
## These are changes seen when WT are exposed to Mn
treatment <- results(dds, name = "Group_wt_exposed_vs_wt_unexposed")
treatment <- lfcShrink(dds = dds, coef = 4, type = 'apeglm')

## Interaction between treatment and genotype
## These are comparing unexposed WT with exposed mutations
interactions <- results(dds, name = "Group_mut_exposed_vs_wt_unexposed")
interactions <- lfcShrink(dds = dds, coef = 2, type = 'apeglm')

## Effect of mutation (wt_unexposed vs mutant unexposed)
mutation_specific_changes <- results(dds, name = "Group_mut_exposed_vs_mut_unexposed")
mutation_specific_changes <- lfcShrink(dds = dds, coef = 2, type = 'apeglm')

## Annotate dataframes with artificial columns here
mutation_specific <- addDirectionlabel(mutation_specific)
treatment <- addDirectionlabel(treatment)
interactions <- addDirectionlabel(interactions)
mutation_specific_changes <- addDirectionlabel(mutation_specific_changes)


## Annotate the datasets with entrezid IDs
mutation_specific <- annot_data(mutation_specific)
treatment <- annot_data(treatment)
interactions <- annot_data(interactions)
mutation_specific_changes <- annot_data(mutation_specific_changes)

## Remove all genes/transcripts that are not annotated in the reference genome 
## Remove genes without annotation in zdfin
mutation_specific <- mutation_specific[!grepl('LOC', rownames(mutation_specific)), ]
treatment <- treatment[!grepl('LOC', rownames(treatment)), ]
interactions <- interactions[!grepl('LOC', rownames(interactions)), ]
mutation_specific_changes <- mutation_specific_changes[!grepl('LOC', rownames(mutation_specific_changes)), ]


## Export the data for sharing and further analysis/visualisation
write.csv(mutation_specific,
          file = paste0(output_folder,'mut_unexposed_vs_wt_unexposed.csv'))

write.csv(treatment,
          file = paste0(output_folder,'wt_exposed_vs_wt_unexposed.csv'))

write.csv(interactions,
          file = paste0(output_folder,'wt_unexposed_vs_mut_exposed.csv'))

write.csv(mutation_specific_changes,
          file = paste0(output_folder,'mut_exposed_vs_mut_unexposed.csv'))
