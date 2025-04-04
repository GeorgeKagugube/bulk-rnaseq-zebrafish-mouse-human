## Function for preparing gene lists for gsea or over representation analysis
## Developed by George William Kagugube 
## Date: 29 Jan 2024 
## Usage, read teh file into R preferably using source (file_name)
## prepare an object from DESeq of the differentially expressed genes
## This is given as an  data frame to the function here

# Load the libraries that will be needed for the analysis herein 
requiredPackages <- c('ggvenn', 'dplyr', 'tidyr', 'data.table', 'ggplot2',
                      'pheatmap', 'RColorBrewer', 'ComplexHeatmap', 'circlize',
                      'EnhancedVolcano', 'gridExtra', 'grid', 'enrichplot',
                      'clusterProfiler', 'org.Dr.eg.db', 'GOSemSim','DOSE',
                      'ggupset','AnnotationDbi', 'AnnotationHub', 'DESeq2', 
                      "pathview","ReactomePA", "gtExtras", "docstring")

##################### Loading the required libraries here ######################
lapply(requiredPackages, library, character.only=TRUE)

################################################################################
################### Functions start from here. #################################
################################################################################
## Function to perform data normalisation 
normalisation_func = function(data_obj){
  #' @description Computes normalised count matrix of the input data
  #' 
  #' Parameters
  #' _________________________
  #' @param data_obj: A dataframe that is output from function x from STAR aligner.
  
  norm_data = counts(estimateSizeFactors(data_obj), normalized = T)
  
  return(norm_data)
}

## Function that creates a gene list for gse analysis
creategenelist <- function(dge_data_set, analysis = 'gsea'){
  #' @title Go or KEGG Gene list  
  #' 
  #' @description Creates a gene list that can be passed on for GO/Pathway analysis
  #' 
  #' Parameters
  #' _________________________
  #' @param dge_data_set: A DESeq log shrink output, annotated with entrezid and esemblid
  #' 
  #' @param analysis: Choice between gsea or others: str. gsea and genelist name will be

  ## Prepare the list to use in the analysis here                      
  gse_gene_list <- dge_data_set$log2FoldChange 

  ## Generate names depending on the type of anlysis
  if (analysis == 'gsea'){
    names(gse_gene_list) <- dge_data_set$X 
    
  } else if (analysis == 'other'){
    names(gse_gene_list) <- dge_data_set$entrezid
  }                                               
  
  ## Remove any missing values in the gene list
  gene_list <- na.omit(gse_gene_list)###########
  
  ## Gene set enrichment analysis needs a sorted list of genes
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  return(gene_list)
}

## Create a list of genes for ORA 
oragenelist <- function(dge_data_set, padj = 0.05, basemean = 50){
  
  #' @title ORA Gene list
  #' _____________________
  #' 
  #' @description
    #' AThis function creates a gene list baseb on a user defined basemean and 
    #' adjusted p-value...
    #' 

  #' @param dge_data_set:  DESeq or EdgR differential gene expression output
  #' 
  #' @param padj: adjusted p-value
  #' @param basemean: descbasemean from the expression values, defualt value is
  #' 50. 
  #' 
  #' ______________________
  #' @example oregenelist.R.
  
  ora_gene_list <- na.omit(dge_data_set[(dge_data_set$baseMean > basemean & dge_data_set$padj < padj),])
  gene_list <- rownames(ora_gene_list)

  return(gene_list)
}

## Create a gene universe to use as the background for the GO term analysis 
gene_universe <- function (dge_data_set, key = 'symbol'){
  if (key == 'symbol'){
    gene_uni <- na.omit(dge_data_set$x)
  } else {
    gene_uni <- na.omit(dge_data_set$entrezid)
  }
  return(gene_uni)
}

## This function runs fgsea using clusterprofiler. By defualt, it looks at all the ontology
## Terms, but this can be modified in the function call via the ONT parameter
## Where tha above function has run successfully, its output should be right input for
## this function
rungseondata <- function(geneList, ONT = 'all'){
  gse <- gseGO(
    geneList = geneList,
    ont = ONT,
    keyType = 'SYMBOL',
    minGSSize = 3,
    maxGSSize = 250,
    eps = 0,
    pvalueCutoff = 0.05,
    verbose = TRUE,
    OrgDb = 'org.Dr.eg.db',
    pAdjustMethod = 'fdr'
  )
  
}

## The difference partly between Over representation analysis and gene set enrichment
## analysis is the 
oraFunc <- function (gene_list, ONT = "ALL") {
  go_enrich <- enrichGO(gene = gene_list,
                        OrgDb = "org.Dr.eg.db", 
                        keyType = 'ENTREZID',
                        ont = ONT,
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "fdr",
                        qvalueCutoff = 0.10)
  return(go_enrich)
}


## KEGG and Reactome Pathway analyses
## KEGG
pathway_ont <- function(genelist, mode = "gse"){
    kk <- gseKEGG(geneList = gene_list,
                  keyType = 'kegg',
                  organism = 'dre',
                  minGSSize = 3,
                  maxGSSize = 1000,
                  pAdjustMethod = 'fdr')
  return(kk)
}

## Reactome analysis
reactome_analysis <- function(geneList, mode = 'gse', lfc = 1.5){

  if (mode == 'gse'){
    print ('Performing ORA analysis')
    pathways <- gsePathway(geneList,
      pvalueCutoff = 0.2,
      pAdjustMethod = 'fdr',
      verbose = F)
  } else {
    print ('Performing ORA analysis')
    pathways <- enrichPathway(names(geneList)[abs(geneList) > lfc])
  }
return(pathways)
}

## Gene Symbol extraction function
gseGO_gene_name_decoding <- function(df, word = "synap"){
  
  # Return the gene names
  gene_symbols <- gsub('/', ',', df$Description[grep(word, df$Description)])
  
  # Return the genes
  return(gene_symbols)
}

## Visualise the keeg pathway
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

zfin <- read.delim("https://zfin.org/downloads/wildtype-expression_fish.txt",
                   sep = '\t', header = FALSE)
