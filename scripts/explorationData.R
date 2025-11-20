## Clear the workspce here 
rm(list = ls())

## Set a global seed here for reproducibility
set.seed(101)

## Set the working directory here 
#setwd('/Users/gwk/OneDrive - University College London/PhD/Data/PhD_data/March_03_25_Final_Analysis/Figure/GO Terms analysis/gsea')
setwd('/Users/gwk/Desktop/Thesis Figures/DifferentialGeneExpression/')

## Load the required libraries here
library(tidyverse)
library(ggplot2)
library(forcats)
library(gtsummary)
library(ggvenn)
library(RColorBrewer)


# Testing the palette with three colors
display.brewer.pal(3, "YlOrRd")

# Define a palette
mypalette <- brewer.pal(3, "YlOrRd")

# how are the colors represented in the mypalette vector?
mypalette

## ================= Differential gene expression ==============================
wtexpo <- read.csv('./WildtypeExposed/WTExposed_vs_WTUnexposed.csv', row.names = 1)
mutunexpo <- read.csv('./MutantUnexposed/mutUnexposed_vs_wtUnexposed.csv', row.names = 1)
mutexpo <- read.csv('./MutantExposed/mutExposed_vs_mutUnexposed.csv', row.names = 1)
mutMn <- read.csv('./Interactions/mutExposed_vs_mutUnexposed.csv', row.names = 1)

## ================= Load GSEA datasets ========================================
wtGSEA <- read.csv('./WildtypeExposed/gseaanalysis.csv', row.names = 1)
mutunexpoGSEA <- read.csv('./MutantUnexposed/gseaanalysis.csv', row.names = 1)
mutexpoGSEA <- read.csv('./MutantExposed/gseaanalysis.csv', row.names = 1)
mutMn_interactionGSEA <- read.csv('./Interactions/gseaanalysis.csv', row.names = 1)

## ============= Load KEGG pathway =============================================
wtKEGG <- read.csv('./WildtypeExposed/kegg_pathway_analysis.csv', row.names = 1)
mutunexpoKEGG <- read.csv('./MutantUnexposed/kegg_mutation_only.csv', row.names = 1)
mutexpoKEGG <- read.csv('./MutantExposed/kegg_mn_mutant.csv', row.names = 1)
mutMn_interactionKEGG <- read.csv('./Interactions/kegg_pathway_analysis.csv', row.names = 1)

## ============= Generate data for venn diagrams here =========================
## DGE Venn diagrams
set1 <- mutexpo|> filter(padj <0.05)|>select(symbols)
set2 <- mutunexpo|> filter(padj <0.05)|>select(symbols)
set3 <- mutMn|> filter(padj <0.05)|>select(symbols)
set4 <- wtexpo|> filter(padj <0.05)|>select(symbols)

## GSEA Venn diagrams
set4 <- mutexpoGSEA$ID
set5 <- mutunexpoGSEA$ID
set6 <- mutMn_interactionGSEA$ID
set7 <- wtGSEA$ID
  
## KEGG Venn diagrams 
set8 <- mutexpoKEGG$ID
set9 <- mutunexpoKEGG$ID
set10 <- mutMn_interactionKEGG$ID
set11 <- wtKEGG$ID

## The set object
dgeset <- list('Homs' = set1$symbols,
               'WT' = set4$symbols,
               'Homs_only' = set2$symbols)

dgsetGSEA <- list('MnEffects' = set4,
                  'Mut_only' = set5,
                  'WTExpo' = set7)

dgsetKEGG <- list('MnEffects' = set8,
                  'Mut_only' = set9,
                  'WTExpo' = set11)

## Draw a venn diagram here
ggvenn(dgeset,
       set_name_size = 10, 
       text_size = 6.5,
       fill_alpha = 0.25,
       stroke_size = 1.5,
       show_outside = "auto")

## Shared terms and pathways
commonGO <- intersect(intersect(as.vector(mn_effect$ID), as.vector(mutExpo$ID)),
                    as.vector(interaction$ID))
commonKEGG <- intersect(intersect(as.vector(mnEffectKEGG$ID), as.vector(mutExpoKEGG$ID)),
                        as.vector(interactionKEGG$ID))

commonkegg_mutation <- interaction(intersect(as.vector(mnEffectKEGG$ID), as.vector(mutationKEGG$ID)),
                                   as.vector(interactionKEGG$ID))
## Shared GO Terms by of the biological test conditions
d1 <- mnEffectKEGG[commonKEGG,] |>
  select(1,2,3,5,7,11) |>
  filter(NES < 0)

d2 <- mutExpo[commonGO,] |>
  select(2,3,6,1) |>
  Filter(NES < 0)
d3 <- interaction[commonGO,] |>
  select(2,3,6,1) |>
  filter(NES < 0)
list(d1,d2,d3) |>
  reduce(full_join, by='ID')
d4 <- merge(d1, d2, by = 'ID') |>
  select()

d1$group = rep('MN_Exposure',27)
d2$group = rep('mut_Exposure',27)
d3$group = rep('Interaction',27)

rbind(d1,d2, d3) |>
  select(2,5,3) |>
  as.matrix() |>
  ggplot(aes(x = group, y = Description, fill = NES)) +
  geom_tile(color = 'black') +
  #scale_fill_gradientn(colors = hcl.colors(20, "RdYlGn")) +
  geom_text(aes(label = NES), color = "white", size = 4) +
  coord_fixed()

d6 <- merge(D5, d3, by = 'Description')
D5 <- d6 |>
  select(1,2,3,5) |>
  rename(NES = 'interaction') |>
  arrange(desc(WTexposed)) |>
  gtable::gtable()
rownames(D5) <- D5$Description
heatmap(as.matrix(D5[,c(2,3,4)]), xlab = '')
  tidyheatmap(df = D5,
              rows = )
## Shared KEGG Pathways
mnEffectKEGG[commonKEGG,]

wtKEGG |>
  filter(NES < 0)|>
  ggplot() +
  geom_point(aes(x = NES, y = Description, color = -log10(p.adjust), 
                 size = setSize)) +
  theme_bw() +
  theme(axis.text.x = element_text(size=rel(2.15), face = 'bold'),
        axis.title = element_text(size=rel(2.15), face = 'bold'),
        axis.text.y = element_text(size=rel(1.5), face = 'bold')) +
  xlab("Normalised Enrichment Score") +
  #ylab('Enriched pathways') +
  scale_color_gradientn(name = "-log10(padj)", colors = mypalette) +
  theme(legend.title = element_text(size=rel(1.15), hjust=0.5, face="bold")) #+
  facet_wrap(~ONTOLOGY)


## Seperate the mitochondrial and cytosolic ribosomes
## Seperate mitochondrial DNA encoded genes versus nuclear subunits
## Nuclear encoded
## Make fresh buffer
## Protease inhibitor (Ask Gareth abot this)
## Protein quantification

## Check how many rows are present in the dataset
nrow(mn_effect)
nrow(mutation)
nrow(mutExpo)
nrow(interaction)

## Set a copy of the dataset that is being analysed here
df <- interaction

## Chec the actual terms in the dataset here

dd <- df |>
  filter(ONTOLOGY == 'CC' & NES > 0) |>
  select(1,3,4,6,8,12) 

dd

df1 <- mnEffectKEGG[commonKEGG,]
lst <- as.vector(unlist(strsplit(d1['dre04020','core_enrichment'],'/')))
mnEffectKEGG[mnEffectKEGG$entrezid == lst,]
mnEffect[mnEffect$entrezid == lst,]

## Explore more
extractGenes <- function(df.data, ontology = 'dre04020 '){
  coreGenes <- as.vector(unlist(strsplit(df.data[ontology,'core_enrichment'],
                                         '/')))
  return(coreGenes)
}

## Create a table that can be used in the report here
## Managese effect in WT exposed
mnEffect[extractGenes(mn_effect),c(7,2,5,6)] |>
  #filter(diffExpression == 'UP' | diffExpression == 'DOWN') |>
  arrange(desc(symbols)) |>
  #rename(symbols = 'Gene', log2FoldChange = 'Log2FC',padj = 'p.adjusted',
  #       diffExpression = 'DGE') |>
  gt::gt() |>
  gt::fmt_number(suffixing = T, n_sigfig = 3) #|>

# Effect of Mn on mutants (MutExposed vs WT unexposed)
mut_mnEffect[extractGenes(d1),c(7,2,5,6)] |>
  filter(diffExpression == 'UP' | diffExpression == 'DOWN') |>
  arrange(desc(diffExpression)) |>
  rename(symbols = 'Gene', log2FoldChange = 'Log2FC',padj = 'p.adjusted',
         diffExpression = 'DGE') |>
  gt() |>
  fmt_number(suffixing = T, n_sigfig = 3) #|>
  #gt::gtsave(filename = ".") # use extensions .png, .html, .docx, .rtf

  ## Mn effects on exposed mutatnts vs mutant Unexposed  
mutEffect[extractGenes(d1),c(7,2,5,6)] |>
    #filter(diffExpression == 'UP' | diffExpression == 'DOWN') |>
  arrange(desc(diffExpression)) |>
    gt() |>
    fmt_number(suffixing = T, n_sigfig = 3) 

mutMut[extractGenes(d1),c(7,2,5,6)] |>
  filter(diffExpression == 'UP' | diffExpression == 'DOWN') |>
  arrange(desc(diffExpression)) |>
  rename(symbols = 'Gene', log2FoldChange = 'Log2FC',padj = 'p.adjusted',
         diffExpression = 'DGE') |>
  gt() |>
  fmt_number(suffixing = T, n_sigfig = 3) 

## Visualisation
wtGSEA |>
  #filter(ONTOLOGY == 'CC' & NES > 0) |>
  select(1,3,4,6,8) |>
ggplot() +
  geom_point(aes(x = NES, y = Description, color = -log10(p.adjust), 
                 size = setSize)) +
  theme_bw() +
  theme(axis.text.x = element_text(size=rel(2.15), face = 'bold'),
        axis.title = element_text(size=rel(2.15)),
        axis.text.y = element_text(size=rel(1.5), face = 'bold')) +
  xlab("Normalised Enrichment Score") +
  ylab("Top 30 significant GO terms") +
  ggtitle("Dotplot of top 30 significant GO terms") +
  theme(plot.title = element_text(hjust=0.5,face = "bold")) +
  scale_color_gradientn(name = "-log10(p.adjust)", colors = mypalette) +
  theme(legend.title = element_text(size=rel(1.15), hjust=0.5, face="bold")) 



  


