## Clear the workspce here 
rm(list = ls())

## Set a global seed here for reproducibility
set.seed(101)

## Set the working directory here 
#setwd('/Users/gwk/OneDrive - University College London/PhD/Data/PhD_data/March_03_25_Final_Analysis/Figure/GO Terms analysis/gsea')
setwd('/Users/gwk/Desktop/Thesis Figures')

## Load the required libraries here
library(tidyverse)
library(ggplot2)
library(forcats)
library(gtsummary)

## Differential gene expression
mnEffect <- read.csv('./mnExposureEffect/mnExposuredge.csv', row.names = 1)
mutEffect <- read.csv('./mutantEffects/mutunexposedwtUnexposed.csv', row.names = 1)
mut_mnEffect <- read.csv('./mnmutInteraction/mutExposedwtUnexposed.csv', row.names = 1)
mutMut <- read.csv('./mutant_mutant/mutExposedwtUnexposed.csv', row.names = 1)

## Load datasets for analysis 
mn_effect <- read.csv('./mnExposureEffect/gseaanalysis.csv', row.names = 1)
mutation <- read.csv('./mutantEffects/gseaanalysis.csv', row.names = 1)
mutExpo <- read.csv('./mutant_mutant/gseaanalysis.csv', row.names = 1)
interaction <- read.csv('./mnmutInteraction/gseaanalysis.csv', row.names = 1)

## Load the KEGG pathway data here 
mnEffectKEGG <- read.csv('./mnExposureEffect/kegg_pathway_analysis.csv', 
                         row.names = 1)
mutationKEGG <- read.csv('./mutantEffects/kegg_pathway_analysis.csv',
                         row.names = 1)
mutExpoKEGG <- read.csv('./mutant_mutant/kegg_pathway_analysis.csv', 
                        row.names = 1)
interactionKEGG <- read.csv('./mnmutInteraction/kegg_pathway_analysis.csv',
                            row.names = 1)
## Explore the dataste here 
head(mn_effect)
names(mn_effect)
head(mutantion)

## 
## Set up the terms for a venn diagrams
set1 <- mn_effect$ID
set2 <- mutExpo$ID
set3 <- interaction$ID

## Pathways
set4 <- mnEffectKEGG$ID
set5 <- mutExpoKEGG$ID
set6 <- interactionKEGG$ID

## The set object
dgeset <- list('MnEffect' = set1,
               'Mn_MutMut' = set2,
               'Interaction' = set3)

dgsetPath <- list('MnEffects' = set4,
                  'Mn_MutMut' = set5,
                  'interaction' = set6)

## Draw a venn diagram here
ggvenn(dgsetPath,
       set_name_size = 10, 
       text_size = 6.5,
       fill_alpha = 0.25,
       stroke_size = 1.5,
       show_outside = "auto")



## Shared terms and pathways
commonGO <- intersect(intersect(as.vector(mn_effect$ID), as.vector(mutExpo$ID)),
                    as.vector(interaction$ID))
commonKEGG <- intersect(intersect(as.vector(mnEffectKEGG$ID), as.vector(mutationKEGG$ID)),
                        as.vector(interactionKEGG$ID))
## Shared GO Terms
d1 <- mn_effect[commonGO,] |>
  select(2,3,6,8,1,12)
#mutation[commonGO,] |>
#  select(3,6,1)
d2 <- mutExpo[commonGO,] |>
  select(2,3,6,1)
d3 <- interaction[commonGO,] |>
  select(2,3,6,1)
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
  gt()
rownames(D5) <- D5$Description
heatmap(as.matrix(D5[,c(2,3,4)]), xlab = '')
  tidyheatmap(df = D5,
              rows = )
## Shared KEGG Pathways
mnEffectKEGG[commonKEGG,]

mutExpo[commonGO,] |>
  filter(NES < 0)|>
  ggplot() +
  geom_point(aes(x = NES, y = Description, color = -log10(p.adjust), 
                 size = setSize)) +
  theme_bw() +
  theme(axis.text.x = element_text(size=rel(2.15), face = 'bold'),
        axis.title = element_text(size=rel(2.15), face = 'bold'),
        axis.text.y = element_text(size=rel(1.5), face = 'bold')) +
  xlab("Normalised Enrichment Score") +
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
lst <- as.vector(unlist(strsplit(df1['dre03008','core_enrichment'],'/')))
mnEffectKEGG[mnEffectKEGG$entrezid == lst,]
mnEffect[mnEffect$entrezid == lst,]

## Explore more
extractGenes <- function(df.data, ontology = 'GO:0099537'){
  coreGenes <- as.vector(unlist(strsplit(df.data[ontology,'core_enrichment'],
                                         '/')))
  return(coreGenes)
}

## Create a table that can be used in the report here
## Managese effect in WT exposed
mnEffect[extractGenes(d1),c(7,2,5,6)] |>
  filter(diffExpression == 'UP' | diffExpression == 'DOWN') |>
  arrange(desc(diffExpression)) |>
  rename(symbols = 'Gene', log2FoldChange = 'Log2FC',padj = 'p.adjusted',
         diffExpression = 'DGE') |>
  gt() |>
  fmt_number(suffixing = T, n_sigfig = 3) #|>

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
df |>
  filter(ONTOLOGY == 'CC' & NES < 0) |>
  select(1,3,4,6,8) |>
ggplot() +
  geom_point(aes(x = NES, y = Description, color = -log10(p.adjust), 
                 size = setSize)) +
  theme_bw() +
  theme(axis.text.x = element_text(size=rel(2.15), face = 'bold'),
        axis.title = element_text(size=rel(2.15)),
        axis.text.y = element_text(size=rel(1.5), face = 'bold')) +
  xlab("Gene ratios") +
  ylab("Top 30 significant GO terms") +
  ggtitle("Dotplot of top 30 significant GO terms") +
  theme(plot.title = element_text(hjust=0.5,face = "bold")) +
  scale_color_gradientn(name = "-log10(p.adjust)", colors = mypalette) +
  theme(legend.title = element_text(size=rel(1.15), hjust=0.5, face="bold"))



  


