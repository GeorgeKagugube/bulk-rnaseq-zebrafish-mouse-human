## Clear the global environment space 
rm(list = ls())

set.seed(101)

## Load the file with all the required libraries here
source("/Users/gwk/Desktop/Bioinformatics/BulkRNA_Analysis/functionsneededforanalysis.R")
source("/Users/gwk/Desktop/Bioinformatics/BulkRNA_Analysis/visualisationFunctions.R")

### Set the working directory here. This is the path to the GO term/pathway files
setwd('/Users/gwk/Desktop/PhD/Data/PhD_data/March_03_25_Final_Analysis/DGE_Files')
#setwd('/Users/gwk/Desktop/PhD/Data/PhD_data/March_03_25_Final_Analysis/gene_ontology_analysis')

## Set the file output directory here 
output_dir <- '/Users/gwk/Desktop/PhD/Data/PhD_data/March_03_25_Final_Analysis/gene_ontology_analysis'

## Load the datasets to be analysed here (these are dge files)


## Zfin data
zfin <- read.delim("https://zfin.org/downloads/wildtype-expression_fish.txt",
                   sep = '\t', header = FALSE)
colnames(zfin) <- c('zfinid','X','genotype','v4','tissue','v6','v7',
                    'stage','stage_1','transcript','v11','v12','v13','v14','v15')
zfin <- zfin %>%
  select(2,3,5,8)

## Merge the two dataframes here 
df <- merge(x=df, y=zfin, by = 'X')

df %>%
  filter(padj < 0.05)

## Exploratroy analysis starts from here
head(mut_expo)
str(mut_expo)

mut_expo[mut_expo$X == 'bcl2a',]
treat <- treatment |>
  as.data.frame() |>
  select(8,7,2) |>
  drop_na()

# Define a palette
mypalette <- brewer.pal(3, "YlOrRd")

## Mn exposure effectes
## GSEA Analysis
df = treatment
gene_list <- creategenelist(df)#  creatgenelist(treatment, analysis = 'other')
gse <- rungseondata(geneList = gene_list)
write.csv(gse,
          file = paste0(output_dir,'/mutantEffects/gseaanalysis.csv'))

## Visualisation
gdd <- gse |>
  as.data.frame() |>
  filter(ONTOLOGY == 'MF' & NES > 0) |>
  arrange(desc(NES))

## Visualisation
ggplot(gdd) +
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
  scale_color_gradientn(name = "-log10(padj)", colors = mypalette) +
  theme(legend.title = element_text(size=rel(1.15), hjust=0.5, face="bold"))

## Run KEGG Pathway analysis here 
df2 <- df |>
  select(7,8,2) |>
  drop_na()
kegg_list <- creategenelist(df2, analysis = 'other')
kegg_organism = "dre"
kk2 <- gseKEGG(geneList     = kegg_list,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")
write.csv(kk2,
          file = paste0(output_dir,'/mutantEffects/kegg_pathway_analysis.csv'))
bp_plot <- kk2 |>
  as.data.frame() |>
  filter(NES < 0) |>
  arrange(desc(NES))  

ggplot(bp_plot) +
  geom_point(aes(x = NES, y = Description, color = -log10(p.adjust), 
             size = setSize)) +
  theme_bw() +
  theme(axis.text.x = element_text(size=rel(1.15)),
        axis.title = element_text(size=rel(1.15))) +
  xlab("Gene ratios") +
  ylab("Top 30 significant GO terms") +
  ggtitle("Dotplot of top 30 significant GO terms") +
  theme(axis.text.x = element_text(size=rel(2.15), face = 'bold'),
        axis.title = element_text(size=rel(2.15)),
        axis.text.y = element_text(size=rel(1.5), face = 'bold')) +
  scale_color_gradientn(name = "Significance \n (-log10(padj))", colors = mypalette) +
  theme(legend.title = element_text(size=rel(1.15),
                                    hjust=0.5, 
                                    face="bold"))


dotplot(kk2, showCategory = 10)#, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
emapplot(kk2)
gse_dataframe <- as.data.frame(gse)

gseGO_gene_name_decoding(as.data.frame(pont), 'Calcium')
#write.csv(gse_dataframe,
#          file = paste0(output_dir,'/mut_exposed_vs_wt_unexposed_all_unsimplified.csv'))

## To eliminate redundant terms, one can run simplify here
gse_simplifeid <- simplify(gse, cutoff=0.7,by="p.adjust",select_fun=min)
gse_simplifeid_dataframe <- as.data.frame(gse_simplifeid)

## Export the csv file for further analysis and interogation
#write.csv(gse_dataframe,
#          file = paste0(output_dir,'/mut_exposed_vs_wt_unexposed_all_simplifeid.csv'))

## Generate a dot plot here 
dotplot(gse_simplifeid, showCategory = 10) +#, split = '.sign') +
  facet_grid(.~.sign) 
  theme(axis.text.x =element_text(angle =45, hjust = 1),
        axis.text.y = element_text(angle = 45, size = 10))

## Create the cnet plot here
cnetplot(gse, categorySize = 'pvalue',
         foldChange = gene_list, showCategory = 2)

## Create a heatmap 
p1 <- heatplot(gse_simplifeid, showCategory = 2)
p2 <- heatplot(gse_simplifeid, foldChange = gene_list[1:10],
               showCategory = 2)

## Combine the plots using cowplot grid function here
cowplot::plot_grid(p1, p2, ncol = 1, labels = LETTERS[1:2])

## Create a ridge plot interpretation of up and down regulated pathways 
ridgeplot(gse_simplifeid) +
  labs(x = 'Enrichment distribution') +
  theme(axis.text.y = element_text(size = 8))

## create a GSEA plot
gseaplot(gse, by = 'all',title = gse$Description[1], geneSetID = 1)

## We can create a gseplot2 here
gseaplot2(gse,title = gse$Description[1], geneSetID = 1)

## We can create a gseplot2 here
gseaplot2(gse,geneSetID = 5:10)#,
          #color = c('blue','red','black'))

###############################################################################
### Gene Universe
anndf <- na.omit(transcriptome_annotation(df4))
background <-anndf$log2FoldChange
names(background) <- anndf$ENTREZID

################################################################################
### Brain data and figure generation 
gse_all <-gseFunc(gse_gene_list,item = "ALL")
gse_all <- simplify(gse_all, cutoff=0.7,by="p.adjust",select_fun=min)

## 
terms <- gse$Description[1:3]
pmcplot(terms, 2011:2024, proportion = FALSE)
pmcplot(terms)
#### Visualisation 
dotplot(gse, showCategory = 15, 
        title = "Gene Ontology Terms: Eye WT Exposed" , 
        split=".sign") + facet_grid(.~.sign)

treeplot(pairwise_termsim(gse), nCluster=12)


## Save the data here
df <- as.data.frame(setReadable(gse_all, 
                                OrgDb = "org.Dr.eg.db",
                                keyType = 'auto'))

write.csv(df, file = "Eyes_WT_goterms.csv")

### kegg pathway analysis
pathway_analysis <- pathway_ont(genelist = gene_list)

## Visualise the data 
dotplot(pathway_analysis, showCategory = 20)

## Export the data here
df22 <- as.data.frame(setReadable(pathway_analysis, 
                                  OrgDb = "org.Dr.eg.db",
                                  keyType = "ENTREZID"))
write.csv(df22, file = "eyes_WT_keggPathways.csv")
################################################################################
## Pathway visualisation
dotplot(pathway_analysis, showCategory=10)
treeplot(pairwise_termsim(pathway_analysis))

fgsea_react <- gsePathway(geneList = gse_kegg_gene_list, 
                          organism = "zebrafish",
                          minGSSize = 15, 
                          maxGSSize = 500, 
                          eps = 0, 
                          nPermSimple = 10000, 
                          seed = TRUE)
dotplot(fgsea_react)
as.data.frame(setReadable(pathway_analysis, 
                          OrgDb = "org.Dr.eg.db",
                          keyType = "ENTREZID"))
## Pathway analysis
############ KEGG PATHWAYS ####################################
# dre04020 ==> Calcium Signaling pathway
# dre04744 ==> Phototransduction
# dre04210 ==> Aptoptosis
# dre00190 ==> Oxidative phosphorylation
# dre04130 ==> SNARE interactions in vesicular transport
# dre04216 ==> Ferroptosis
############################################################

keggPathway <- function(data_frame, organism = 'dre', pathWay = 'dre04137'){
  # Load the required library
  library(pathview)
  
  ## Generate the gene list for the analysis here 
  gene_list <- creatgenelist(data_frame, analysis = 'other')
  
  # Produce the native KEGG plot (PNG)
  dme <- pathview(gene.data=gene_list, pathway.id=pathWay, species = organism)
  
  # Produce a different plot (PDF) (not displayed here)
  dme <- pathview(gene.data=gene_list, pathway.id=pathWay, species = 'dre', kegg.native = F)
  
  return (knitr::include_graphics(paste0(pathWay,".pathview.png")))
}  
 
keggPathway(mut_expo, pathWay = 'dre03410')

#####
makeSPIAdata(kgml.path=system.file("extdata/keggxml/dre",package="SPIA"),organism="dre",out.path="./")
makeSPIAdata(kgml.path="./dre",organism="dre",out.path=".")

spia_result <- spia(de=gse_kegg_gene_list, all=background, organism="dre")


#### Visualisation of the data
c.chaperones <- c('hspa5', 'hsp90b1','calr3a','canx')
head(all_norm_data)
er.stress <- all_norm_data[c.chaperones,]
colnames(er.stress) <- samples$Group

er.stress <- er.stress |>
  t() |>
  as.data.frame() 
er.stress['Samples'] <- rownames(er.stress)  
er.stress$Samples <- gsub("\\.\\d+$", "", er.stress$Samples)
attach(er.stress)
er.stress$Samples <- as.factor(er.stress$Samples)

my_comparisons <- list(c('wt_unexposed', 'wt_exposed'), c('mut_unexposed', 'mut_exposed'), 
                       c('wt_unexposed', 'mut_exposed'))

ggplot(er.stress, aes(x=Samples, y=calr3a, colour = Samples)) + 
    geom_dotplot(binaxis='y', stackdir='center') +
    theme_classic() +
    theme(axis.text.y = element_text(size = 15, face = 'bold')) +
    stat_compare_means(comparisons = my_comparisons, label = "p.signif")

gse_simplifeid  
