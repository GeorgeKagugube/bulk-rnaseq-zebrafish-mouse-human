## Clear the environmental space here to start affresh 
## Clear environmental space
rm(list = ls())

## Set the directory with the expression datasets here 
setwd('/Users/gwk/Desktop/PhD/Data/PhD_data/Brain/GZ11_star_output/star')

# This is the output folder for the final analysis
output_folder = '/Users/gwk/Desktop/PhD/Data/PhD_data/March_03_25_Final_Analysis/DGE_Files/21_04_25/'
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

# Build a DESeq2 object here. The count data is given as a metrix, column names
# are the samples and the rownames are the Ensemble gene ids. This is important
# The design contains what makes up the model (negatve bionimial model in this case)
dds <- DESeqDataSetFromMatrix(countData = as.matrix(countMatrix),
                              colData = samples,
                              design = ~ Genotype + Treatment + Genotype:Treatment)

dds_mut <- DESeqDataSetFromMatrix(countData = countMut,
                                 colData = mut,
                                 design = ~ Group)

dds_wt <- DESeqDataSetFromMatrix(countData = as.matrix(countWT),
                                 colData = wt,
                                 design = ~ Group)

##Although filtering with DESeq is not reccomended, here we remove all genes whose
## row sum is less than 10
#keep <- rowSums(counts(dds)) >= 10
#dds <- dds[keep,]

## Filter mutants here
#keep_mut <- rowSums(counts(dds_mut)) >= 10
#dds_mut <- dds_mut[keep_mut,]

## Filter wild type data here
#keep_wt <- rowSums(counts(dds_wt)) >= 10
#dds_wt <- dds_wt[keep_wt,]

## Relevel the conditions here
#dds$Group <- relevel(dds$Group, ref = 'wt_unexposed')
#dds_mut$Group <- relevel(dds_mut$Group, ref = 'mut_unexposed')
#dds_wt$Group <- relevel(dds_wt$Group, ref = 'wt_unexposed')

## Perform some quality checks here. Use the developed function here. Input is the
## dds object and in some cases you need the sample_id as the function inputs expect
## check out the individual functions below in the environmnetal space
## calculate the pca values here
#vsdmut <- vst(dds_mut)
#vsdwt <- vst(dds_wt)
vsd <- vst(dds)
vsdmut <- vst(dds_mut)
vsdwt <- vst(dds_wt)

## Plot visualise the PCA
plotPCA(vsd, intgroup="Treatment")
#plotPCA(vsdmut, intgroup = 'Group')
#plotPCA(vsdwt, intgroup = 'Group')

## Normalise the count data
all_norm_data <- normalisation_func(dds)
mut_exposed <- normalisation_func(dds_mut)
wt_exposed <- normalisation_func(dds_wt)

## calculate the principle component analysis
#pca_mut <- principle_component(dds_mut)
#pca_wt <- principle_component(dds_wt)
pca_all <- principle_component(dds)
  
# Create a data frame that can used going forward from here on 
#count_mut_df <- cbind(mut,pca_mut$x)
#count_wt_df <- cbind(wt, pca_wt$x)
count_all <- cbind(samples, pca_all$x)

## Visualise the data here 
ggplot(data = count_all) +
  geom_point(aes(x=PC7, y=PC8, colour = Group), size=5) +
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
#dds_mut <- DESeq(dds_mut)
#dds_wt <- DESeq(dds_wt)

res <- results(dds)
# ress_mut <- results(dds_mut)
# res_wt <- results(dds_wt)

resultsNames(dds)
# resultsNames(dds_mut)
# resultsNames(dds_wt)

## Genotype specific changes
genotype <- results(dds, name =  "Genotype_wt_vs_mutant")
genotype <- lfcShrink(dds = dds, coef = 2, type = 'apeglm')

## Maganese responsiveness changes
treatment <- results(dds, name = "Treatment_unexposed_vs_exposed")
treatment <- lfcShrink(dds = dds, coef = 3, type = 'apeglm')

## Interaction between treatment and genotype
interactions <- results(dds, name = "Genotypewt.Treatmentunexposed")
interactions <- lfcShrink(dds = dds, coef = 4, type = "apeglm")

summary(genotype)


## Annotate dataframes with artificial columns here
geno <- addDirectionlabel(genotype)
treat <- addDirectionlabel(treatment)
inter <- addDirectionlabel(interactions)

## Export the data 
write.csv(geno,
          file = paste0(output_folder,'genotype_specific_changes.csv'))

write.csv(treat,
          file = paste0(output_folder,'exposure_specific_changes.csv'))

write.csv(inter,
          file = paste0(output_folder,'genotype_treament_interactions.csv'))


## Annotate dataframes with entrezi_ids
anotated_interactions <- annot_data(geno)
anotated_treatment <- annot_data(treatment)
anotated_genotype <- annot_data(genotype)

## Generate Volcano plots here 
vPlot(geno)
vPlot(treat)
vPlot(inter)

## Generate Ven daigrams

mn_treatment <- sig_gene_names(treat)
genotype_effect <- sig_gene_names(geno)
inte <- sig_gene_names(inter) 

dge1 <- list('Exposed' = mn_treatment,
             'Genotype' = genotype_effect)
ggvenn(dge1)

dge2 <- list('Exposed' = mn_treatment,
             'Genotype' = genotype_effect,
             'inter' = inte)
ggvenn(dge1)
dge <- list('wt_mut' = rownames(lfc[lfc$padj < 0.05,]),
            'mut_exp' = rownames(resLFC_mut[resLFC_mut$padj < 0.05,]),
            'wt_expo' = rownames(resLFC_wt[resLFC_wt$padj < 0.05,]))
levels(samples$Group)
samples$Group <- ordered(samples$Group,
                   levels = c('wt_unexposed','wt_exposed',
                              'mut_unexposed','mut_exposed'))
## Overlapping genes between
overlap_genes <- mn_treatment[mn_treatment %in% genotype_effect]
over_genes <- treat[overlap_genes,]
over_genes <- over_genes |>
  arrange(padj)
## Generate heatmaps here 
max_genes = 40
## Signifacnt genes 
signifcant_genes <- sig_gene(inter)[1:max_genes,]
heatmap_mtx <- data_heatmap(mut_exposed, sample = wt, lfc = inter)[1:max_genes,]

## Figure drwan here
figure_heatmap(heatmap_mtx, signifcant_genes)

## Annotate dataframes with entrezi_ids
anotated_interactions <- annot_data(geno)
anotated_treatment <- annot_data(treatment)
anotated_genotype <- annot_data(genotype)

## GO AND KEGG pathway analysis
inter_gene_list <- creategenelist(anotated_interactions)
geno_gene_list <- creategenelist(anotated_genotype, analysis = 'other')
treat_gene_list <- creategenelist(anotated_treatment, analysis = 'other')

##
go_inter <- rungseondata(geno_gene_list)
go_inter |>
  as.data.frame() |>
  filter(ONTOLOGY == 'BP') |>
  select('ONTOLOGY', 'Description', 'NES','qvalue') |>
  arrange(NES) |>
  head(40)
  
endoplasmic %in% go_inter$Description 
## Visualise these data
dotplot(go_inter, showCategory=20)

wt[wt %in% mutComb]
ggvenn(dge)

res05 <- results(dds_wt, alpha=0.05)
summary(res05)
res05[res05$padj < 0.05]
res05 %>%
  filter(padj < 0.05)

## Export data to a csv file for further analysis and inspection later
## Export the data 
#write.csv(resLFC,
#          file = paste0(output_folder,'Group_mut_exposed_vs_mut_unexposed.csv'))


# ## sET UP GROUP COMPARSIONS HERE 
# res <- results(dds_mut, name = 'Group_mut_exposed_vs_mut_unexposed')
# res <- lfcShrink(dds = dds_mut, coef = 2, type = 'apeglm')
# 
# ## Set up group comparisions here 
# res_wt_mn <- results(dds, contrast = c('Group','wt_unexposed','wt_exposed'))
# res_wt_mn <- lfcShrink(dds, contrast = c('Group','wt_unexposed','wt_exposed'),  type = 'ashr')
# 
# ## Baseline mutant effect on transcriptomics
# res_controls <- results(dds, contrast = c('Group', 'wt_unexposed','mut_unexposed'))
# res_controls <- lfcShrink(dds, contrast = c('Group', 'wt_unexposed','mut_unexposed'),  type = 'ashr')
# 
# ## Manganese effect in mutant
# mn_effect_mut <- results(dds, contrast = c('Group','mut_unexposed', 'mut_exposed'))
# mn_effect_mut <- lfcShrink(dds, contrast = c('Group','mut_unexposed', 'mut_exposed'),  type = 'ashr')
# 
# ## Interaction between genotype and manganese
# mn_genotype_interaction <- results(dds, contrast = c('Group', 'wt_exposed', 'mut_exposed'))
# mn_genotype_interaction <- lfcShrink(dds, contrast = c('Group','wt_exposed', 'mut_exposed'),  type = 'ashr')
# 
# ## annotate the file before exporting them
# mn_genotype_interaction <- annot_data(mn_genotype_interaction)
# mn_effect_mut <- annot_data(mn_effect_mut)
# res_controls <- annot_data(res_controls)
# res_wt_mn <- annot_data(res_wt_mn)

