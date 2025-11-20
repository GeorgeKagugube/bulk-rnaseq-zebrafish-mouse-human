## Clear the global environment space 
rm(list = ls())

set.seed(101)

## Load the file with all the required libraries here
source("/Users/gwk/Desktop/Bioinformatics/bulk-rnaseq-zebrafish-mouse-human/scripts/functionsneededforanalysis.R")
source("/Users/gwk/Desktop/Bioinformatics/bulk-rnaseq-zebrafish-mouse-human/scripts/visualisationFunctions.R")

### Set the working directory here. This is the path to the GO term/pathway files
output_dir = '/Users/gwk/Desktop/Thesis Figures/DifferentialGeneExpression'#/MutantUnexposed/Figures'
setwd('/Users/gwk/Desktop/Thesis Figures')

## Load the datasets to be analysed here (these are dge files)
mutExposed <-read.csv('mutant_mutant/mnExposuredge.csv', row.names = 1)
## Filter out the un-annotated genes
mutExposed <- mutExposed[!grepl('si:', rownames(mutExposed)), ]
mutExposed <- mutExposed[!grepl('zgc:', rownames(mutExposed)), ]

mutUnexpo <- read.csv('mutantEffects/mutunexposedwtUnexposed.csv', row.names = 1)
## Filter out the un-annotated genes
mutUnexpo <- mutUnexpo[!grepl('si:', rownames(mutUnexpo)), ]
mutUnexpo <- mutUnexpo[!grepl('zgc:', rownames(mutUnexpo)), ]

wtExposed <- read.csv('mnExposureEffect/mnExposuredge.csv', row.names = 1)
## Filter out the un-annotated genes
wtExposed <- wtExposed[!grepl('si:', rownames(wtExposed)), ]
wtExposed <- wtExposed[!grepl('zgc:', rownames(wtExposed)), ]

## Create the table needed for this work
wtExposed |>
  as.data.frame() |>
  filter(padj < 0.05 & log2FoldChange > 0) |>
  select(c(7,2,5,6)) |>
  mutate('log2FoldChange'=round(log2FoldChange,2),'padj' = round(padj,2)) |>
  arrange(desc(log2FoldChange)) |>
  head(10) |>
  gt(auto_align = T) |>
  gt::gtsave(filename = "./Tables/Top10DGE_DOWNgenesWtExpxp.rtf")
## Perform GSEA analysis here 
mutExposed_gsea <- simplify(rungseondata(geneList = creategenelist(mutExposed)),
         cutoff = 0.5, by = 'p.adjust', measure = 'Wang')

write_tsv(mutExposed, paste0(output_dir,'/mutExposedGOTERMS.tsv'))

mutUnexo_gsea <- simplify(rungseondata(geneList = creategenelist(mutUnexpo)),
         cutoff = 0.5, by = 'p.adjust', measure = 'Wang')

wtExposed_gsea <- simplify(rungseondata(geneList = creategenelist(wtExposed)),
         cutoff = 0.5, by = 'p.adjust', measure = 'Wang')

mutExposed_gsea |>
  as.data.frame() |>
  select(c(1,3,4,6,8)) |>
  mutate('NES'=round(NES,2),'p.adjust' = round(p.adjust,3)) |>
  arrange(desc(NES)) |>
  gt(auto_align = T) |>
  gt::gtsave(filename = "./Tables/goTerm_wtExposed.rtf")

## Dotplots
dotplot(mutExposed_gsea, showCategory = mutupcat,
        font.size = 15)

dotplot(mutUnexo_gsea, x = 'NES',
        font.size = 15)

dotplot(wtExposed_gsea, showCategory = 5,
        split = 'ONTOLOGY')

## ===================== GSEA KEGG Pathway analysis here =======================
mutExposed_kegg <- gseKEGG(geneList     = creategenelist(na.omit(mutExposed[,c(7,8,2)]), 'other'),
                           organism     = 'dre',
                           nPerm        = 10000,
                           minGSSize    = 3,
                           maxGSSize    = 800,
                           pvalueCutoff = 0.05,
                           pAdjustMethod = "none",
                           keyType       = "ncbi-geneid")

mutExposed_kegg |>
  as.data.frame()

mutUnexpo_kegg <- gseKEGG(geneList     = creategenelist(na.omit(mutUnexpo[,c(7,8,2)]), 'other'),
                           organism     = 'dre',
                           nPerm        = 10000,
                           minGSSize    = 3,
                           maxGSSize    = 800,
                           pvalueCutoff = 0.05,
                           pAdjustMethod = "none",
                           keyType       = "ncbi-geneid")

mutUnexpo_kegg |>
  as.data.frame()
wtExposed_kegg <- gseKEGG(geneList     = creategenelist(na.omit(wtExposed[,c(7,8,2)]), 'other'),
                           organism     = 'dre',
                           nPerm        = 10000,
                           minGSSize    = 3,
                           maxGSSize    = 800,
                           pvalueCutoff = 0.05,
                           pAdjustMethod = "none",
                           keyType       = "ncbi-geneid")
wtExposed_gsea |>
  as.data.frame() |>
  nrow()
write.csv(wtExposed_kegg,
         file = paste0(output_dir,'/WildtypeExposed/kegg_pathway_analysis.csv'))
## Perform some basic visualisaton here 
## doplot
mutCategories <- c('synaptic membrane','chemical synaptic transmission',
                   'synapse organization','xon guidance', 
                   'central nervous system neuron differentiation', 'receptor localization to synapse',
                   'calcium-ion regulated exocytosis')

mutupcat <- c('oxidoreductase complex','proteasome complex','endoplasmic reticulum lumen',
              'ribosome biogenesis','protein maturation','protein folding',
              'respiratory electron transport chain','structural constituent of ribosome',
              'endopeptidase regulator activity')

mutExposed_gsea[mutCategories %in% mutExposed_gsea$Description,]
dotplot(wtExposed_gsea, showCategory = mutCategories,
        x = 'NES', font.size = 15)

dotplot(mutUnexo_gsea, x = 'NES',
        font.size = 15)

dotplot(wtExposed_gsea, showCategory = 5,
        split = 'ONTOLOGY')

## overlapping terms
mn <- mutExposed_kegg |>
  as.data.frame() |>
  select(2) |>
  unlist()

mut <- mutUnexpo_kegg |>
  as.data.frame() |>
  select(2) |>
  unlist()

wt <- wtExposed_kegg |>
  as.data.frame() |>
  select(2) |>
  unlist()

dgsea <- list(mutExposed = mn,
              mutUnexposed = mut,
              wtExposed = wt)
ggvenn(dgsea)
membership <- intersect(mn, wt)

mutExposed_gsea[membership,]

inter = rep(0, length(unlist(membership)))
for (item in 1:length(unlist(membership))){
  inter[item] = unlist(membership)[item]
}

appendixFigure <- mutExposed_gsea |>
  as.data.frame() |>
  select(c(1,2,3,4,6,8,12))

write.csv(appendixFigure,
          file = paste0(output_dir,'/mutExposedGOTERMS.csv'))

dotplot(mutExposed_gsea, showCategory = test)
test <- as.vector(gsub(' ', ',', inter))

## Check the known transporters in each comparision
mnTransporters <- c('slc39a14', 'slc39a8', 'slc30a10', 'slc11a2')

## Check for the expression of these transporter 
mutExposed[mnTransporters,] |>
  select(c(7,2,4,5)) |>
  mutate(`Fold Change` = logratio2foldchange(log2FoldChange, base = 2)) |>
  select(c(1,5,4)) |>
  rename('symbols'='Gene') |>
  gt() #|>
  gt::gtsave(filename = "./Tables/mnSpecificTransporters_in_mutExposed.png")

mutUnexpo[mnTransporters,] |>
  select(c(7,2,4,5)) |>
  mutate(`Fold Change` = logratio2foldchange(log2FoldChange, base = 2)) |>
  select(c(1,5,4)) |>
  rename('symbols'='Gene') |>
  gt() |>
  gt::gtsave(filename = "./Tables/mnSpecificTransporters_in_mutUnexposed.png")

wtExposed[mnTransporters,] |>
  select(c(7,2,4,5)) |>
  mutate(`Fold Change` = logratio2foldchange(log2FoldChange, base = 2)) |>
  select(c(1,5,4)) |>
  rename('symbols'='Gene') |>
  gt() |>
  gt::gtsave(filename = "./Tables/mnSpecificTransporters_in_WTExposed.png")


# Define a palette
mypalette <- brewer.pal(3, "YlOrRd")

## Load the data to be analyised here 
#df <- read.csv('/Users/gwk/Desktop/Thesis Figures/mnExposureEffect/gseaanalysis.csv')
df <- read.csv('/Users/gwk/Desktop/Thesis Figures/mutant_mutant/mutExposedwtUnexposed.csv',  row.names = 1)
head(df)

volcanoPlot(wtexpo, xlimlimit = c(-1.5, 5), ylimlimit = c(0, 20))
## Mn exposure effectes
## ========================= GSEA Analysis: GO Term analysis ===================
gene_list <- creategenelist(mutExposed)#  creatgenelist(treatment, analysis = 'other')
gse <- rungseondata(geneList = gene_list)

gse_simplified <- simplify(gse, cutoff=0.5, by = 'p.adjust', measure = 'Wang')
simplify(rungseondata(geneList = creategenelist()),
         cutoff = 0.5, by = 'p.adjust', measure = 'Wang')

write.csv(gse,
          file = paste0(output_dir,'./Interactions/gseaanalysis.csv'))

## Visualisation
dotplot(gse_simplified, showCategory=10)#, split=".sign") + facet_grid(.~.sign)

## Network plot
xx <- compareCluster(gse, fun="enrichGO", OrgDb="org.Dr.eg.db")
x2 <- pairwise_termsim(gse_simplified)
emapplot(gse_simplified, showCategory = 10)

for (item in 1:nrow(mutUnexo_gsea)){
  # plot1 <- gseaplot(mutExposed_gsea, by = "all", title = mutExposed_gsea$Description[item], geneSetID = item)
  # 
  # print (plot1)
  fileName <- file.path(output_dir,paste0(mutUnexo_gsea$Description[item],'.jpeg'))
  
  ## Open jpeg file
  #jpeg(fileName,width=550, height=400)
  plot1 <- gseaplot(mutUnexo_gsea, by = "all", title = mutUnexo_gsea$Description[item], geneSetID = item)
  ggsave(filename = fileName,
         plot = plot1,width = 5.5, height = 4, dpi = 300, limitsize = F)

}
# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
gseaplot(mutExposed_gsea, by = "all", title = mutExposed_gsea$Description[2], geneSetID = 2)

## Poolout the genes that contributing to a GO Term 
go_term <- mutExposed_gsea |>
  as.data.frame() |>
  filter(Description == 'protein folding') |>
  select(core_enrichment) |>
  unlist() |>
  strsplit('/') |>
  unlist() |>
  as.vector() 

# ## Selected gene lists to explore
erProteostasis <- c('canx','calr3a','p4hb','pdia4','txndc5','fkbp9','hsp90b1',
                     'hspa5','hsp90aa1.1','cct2')
# caBuffer <- c('calub','calua','calr3a','hyou1','p4hb','pdia4','txndc5','casq1b',
#               'casq2','hspa5','pdia3')
# amp <- c('gria2b','gria3a','gria3b','gria4a','grin1a','grin1b','grin2aa',
#          'grin2bb','grik1a','grid1b')
# aux <- c('cacng2a', 'cacng3b', 'cacng4b', 'cacng5b')
# 
# activity <- c('bdnf','npas4a','ntrk2b','camk2d1','camk2d2','prkch')
# 
# upreg <- c('chia.5','chia.6','calr3a','calua', 'calub', 'canx', 'fkbp9', 
#            'si:dkey-22i16.3','vtnb','fjx1')
# 
# downreg <- c('slc4a5','irbp','oatx','rtn4rl2a', 'rtn4rl2b', 'vgf','adgrb1b',
#              'scg2a','si:dkey-175g6.2')
# psdGenes <- c('cadps2', 'cadpsa', 'cadpsb', 'baiap3')

t1 <- mutExposed[go_term,] |>
  filter(padj < 0.05)
  select(c(7,2,5)) |>
  mutate('log2FoldChange' = round(log2FoldChange, 2), 
         padj = round(padj, 3)) |>
  rename('log2FoldChange'='mutExposed') 

t2 <- wtExposed[erProteostasis,] |>
  select(c(7,2,5)) |>
  mutate('log2FoldChange' = round(log2FoldChange, 2), 
         padj = round(padj, 3)) |>
  rename('log2FoldChange' = 'wtExposed')

t3 <- mutUnexpo[erProteostasis,] |>
  select(c(7,2,5)) |>
  mutate('log2FoldChange' = round(log2FoldChange, 2), 
         padj = round(padj, 3)) |>
  rename('log2FoldChange' = 'MutantionOnly')

cbind(t1, t2,t3) |>
  select(c(1,2,5,8)) |>
  rename(symbols = 'Genes') |>
  arrange(Genes) |>
  gt(auto_align = F) |>
  gt::gtsave(filename = "./subcellularTables/endoplasmic reticulum lumen.rtf")

# mutExposed |>
#   filter(padj < 0.05 & log2FoldChange > 0) |>
#   select(c(7,2,5)) |>
#   mutate(`Fold Change` = logratio2foldchange(log2FoldChange, base = 2)) |>
#   select(c(1,4),3) |>
#   arrange(`Fold Change`) |>
#   gt(auto_align = F) |>
#   gt::gtsave(filename = "./Tables/potassium ion transmembrane transport_mutExposed.rtf")

gdd <- wtExposed_gsea |>
  as.data.frame() |>
  filter(NES < -2) |>
  arrange(desc(NES))

## Visualisation
ggplot(gdd) +
  geom_point(aes(x = NES, y = Description, color = -log10(p.adjust), 
                 size = setSize)) +
  facet_wrap(~ONTOLOGY) +
  theme_bw() +
  theme(axis.text.x = element_text(size=rel(2.00), face = 'bold', angle = 90),
        axis.title = element_text(size=rel(2.15)),
        axis.text.y = element_text(size=rel(1.5), face = 'bold')) +
  xlab("Normalised Enrichment Score") +
  #ylab("Top 30 significant GO terms") +
  ggtitle("Dotplot of top 30 significant GO terms") +
  theme(plot.title = element_text(hjust=0.5,face = "bold")) +
  scale_color_gradientn(name = "-log10(padj)", colors = mypalette) +
  theme(legend.title = element_text(size=rel(1.15), hjust=0.5, face="bold")) +
  scale_x_continuous(n.breaks = 4)

## ===================== GSEA KEGG Pathway analysis here =======================
mutExposed_kegg <- gseKEGG(geneList     = creategenelist(na.omit(mutExposed[,c(7,8,2)]), 'other'),
               organism     = 'dre',
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")
write.csv(kk2,
          file = paste0(output_dir,'/Interactions/kegg_pathway_analysis.csv'))

## Visualise KEGG pathways here in a figure ====================================
bp_plot <- kk2 |>
  as.data.frame() |>
  mutate(direction = ifelse(NES>0, 'UP', 'DOWN')) |>
  arrange(NES) |>
  filter(NES < -1)

## Generate a visualisation here 
ggplot(bp_plot) +
  geom_point(aes(x = NES, y = Description, color = -log10(p.adjust), 
             size = setSize)) +
  facet_wrap(~direction) +
  theme_bw() +
  theme(axis.text.x = element_text(size=rel(1.15)),
        axis.title = element_text(size=rel(1.15))) +
  xlab("Normalised Enrichment Score") +
  ylab("Top 30 significant GO terms") +
  ggtitle("Dotplot of top 30 significant GO terms") +
  theme(axis.text.x = element_text(size=rel(2.15), face = 'bold'),
        axis.title = element_text(size=rel(2.15)),
        axis.text.y = element_text(size=rel(1.5), face = 'bold')) +
  scale_color_gradientn(name = "Significance \n (-log10(padj))", colors = mypalette) +
  theme(legend.title = element_text(size=rel(1.15),
                                    hjust=0.5, 
                                    face="bold"))


dotplot(gse, showCategory = 10)#, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
emapplot(gse_simplifeid)
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
cnetplot(gse_simplified, categorySize = 'pvalue',
         foldChange = gene_list, showCategory = 5)

## Create a heatmap 
p1 <- heatplot(gse_simplified, showCategory = 2)
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
