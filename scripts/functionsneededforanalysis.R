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
#' Importing the data for analysis
#' @description
#' This function reads data from STAR mapping software
#' 
#' Parameters
#' @param path_to_star files directory containing all the files 
#' @param one_star_file Import one of the files to be imported in for analysis
#' 
staroutput_preprocessing <- function(path_to_star_files, one_star_file){
  # Function takes in two inputs
  # 1.path to directory with the star out put files
  # 2. On of these files to read and continue to create/build the data matrix
  # I would define these outside of the function as variables that I supply to
  # the function
  ## Set the working directory/project directory here
  setwd(path_to_star_files)
  
  ## Loading data into R
  raw_data <- read.table(file = one_star_file,sep = "\t",header = F,as.is = T)
  
  ## Gather all the data files in the directory of interest here
  data_files <- dir(pattern = "*ReadsPerGene.out.tab")
  
  # Create an empty vetor here
  countmatrix <- c()
  
  ## Loop through the files and extract the counts from star
  for (sample in seq_along(data_files)) {
    input_data <- read.table(file = data_files[sample],sep = "\t",header = F, 
                             as.is = T)
    
    ## Combine teh samples into a count table here 
    countmatrix <- cbind(countmatrix, input_data[,2])
  }
  
  ## Convert the countmatrix to a data frame here
  CountMatrix <- as.data.frame(countmatrix)
  
  ## Add rw names (ENSEMBL Gene IdS)
  rownames(CountMatrix) <- raw_data[,1]
  
  ## Rename the columns with the sample names 
  colnames(CountMatrix) <- sub("ReadsPerGene.out.tab","",data_files)
  
  CountMatrix <- CountMatrix[-c(1:4),]
  ## Return the count matrix here
  return(CountMatrix)
}

## Function to perform data normalisation 
#' Normalise transcript level data
#' 
#' @description Computes normalised count matrix of the input data
#' Parameters
#' _________________________
#' @param data_obj: A dataframe that is output from function x from STAR aligner.
#' @export
normalisation_func = function(data_obj){
  norm_data = counts(estimateSizeFactors(data_obj), normalized = T)
  
  return(norm_data)
}

## Transform normalised data for visualisation here 
principle_component = function (data){
  rld_data <- rlog(data, blind = T)
  pca = prcomp(t(assay(rld_data)))
  return(pca)
}
sig_gene <- function(datafframe) {
  sigs <- datafframe %>%
    as.data.frame() %>%
    filter(padj < 0.05)
  ## Return value from the data frame here
  return(sigs)
}

sig_gene_names <- function(datafframe) {
  sigs <- datafframe %>%
    as.data.frame() %>%
    filter(padj < 0.05) %>%
    rownames()
  ## Return value from the data frame here
  return(sigs)
}

#' Generate a matrix for constructing a heatmap of selected data
#' @description
#' Function that computes a matrix that can be used to heatmap
#' 
#' @param normalised_data A dataframe of the normalised count data ()
#' @param sample data frame containing sample information used to build DESeq2 object
#' @param lfc log2foldchange (shrunk) dataset, needed to extract significant genes
#' 
data_heatmap <- function(normalised_data, sample){
  mut_mat.z <- t(apply(normalised_data, 1, scale))
  colnames(mut_mat.z) <- sample$Group
  
  return(mut_mat.z)
}

## Plot a heatMap using this function
#' Plot Heatmap
#' 
#' @description
#' This function plots a heatmap
#' 
#' @params: data_frame: a matrix of normalised counts with the rownames set to the gene names
#'          signs: a dataframe of selected gens from the DGE analysis
#' 
figure_heatmap <- function(data_frame, sigs){
  valid_rows <- intersect(rownames(sigs), rownames(data_frame))
  plot1 <- Heatmap(data_frame[valid_rows,], name = "Z - Score", row_km = 2, column_km = 2,row_labels = valid_rows,
                   column_labels = colnames(data_frame),
                   #top_annotation = HeatmapAnnotation(data_frame = 1:dim(data_frame)[2]) #, bar1 = anno_points(runif(dim(data_frame)[2]))),
                   #right_annotation = rowAnnotation(data_frame = dim(data_frame[rownames(signs),])[1]:1, 
                   #                                  bar2 = anno_barplot(runif(dim(data_frame[rownames(signs),])[1])))
  )
  
  return(plot1)
}

#' Adding gene symbols as a coulmn and another column for to include the direction of change of expression
#' Create a new dataframe for the existing lfcshrunk dataframe 
#' @description
#' Function that takes in a dataset from DESeq2 and adds two columns one
#' 
#' @param dataset An annotated dataframe from DESeq with gene symbols as rows 

addDirectionlabel <- function(dataset) {
  mutated_dataset <- dataset %>%
    as.data.frame() %>%
    mutate(diffExpression = case_when(
      padj < 0.05 & log2FoldChange > 0 ~ "UP",
      padj < 0.05 & log2FoldChange < -0 ~ "DOWN",
      TRUE ~ "NO"
    ),
    symbols = row.names(dataset))
  
  return(mutated_dataset)
}


####### valcano plots
volcanoPlot <- function(df, noGenes, xlimlimit = c(-8, 8), ylimlimit = c(0, 25)) {
  df$delabel <- ifelse(df$symbols %in% head(df[order(df$padj), "symbols"], 50), df$symbols, NA)
  
  myvolcanoplot <- ggplot(data = df, aes(x = log2FoldChange, y = -log10(padj), col = diffExpression, label = delabel)) +
    geom_vline(xintercept = c(-0.5, 0.5), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
    geom_point(size = 3) + 
    scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"))+#, # to set the colours of our variable  
                      # labels = c("Downregulated", "Not significant", "Upregulated")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
    coord_cartesian(ylim = ylimlimit, xlim = xlimlimit) + # since some genes can have minuslog10padj of inf, we set these limits
    labs(color = 'Severe', #legend_title, 
         x = expression("Log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
    scale_x_continuous(breaks = seq(-10, 10, 2)) + # to customise the breaks in the x axis
    #geom_text_repel(max.overlaps = Inf) + # To show all labels 
    geom_label_repel(label.size = 0.5) +
    theme_classic(base_size = 30) +
    theme(
      axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.5), color = 'black'),
      axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.5), color = 'black'),
      axis.text.x = element_text(face = 'bold', size = rel(1.0)),
      axis.text.y = element_text(face = 'bold', size = rel(1.0)),
      plot.title = element_text(hjust = 0.5)
    )
  
  return(myvolcanoplot)
}

## Venn Diagram of any number of sets to be drawn here
venPlot <- function(dgeset){
  p3 <- ggvenn(dgeset, 
               set_name_size = 10, 
               text_size = 7.5,
               fill_alpha = 0.25,
               stroke_size = 1.5,
               show_outside = "auto")
  
  return(p3)
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
    names(gse_gene_list) <- dge_data_set$symbols 
    
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
oragenelist <- function(dge_data_set, padj = 0.05, logFC = 1){
  
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
  
  ora_gene_list <- na.omit(dge_data_set[(abs(dge_data_set$log2FoldChange) > logFC & dge_data_set$padj < padj),])
  gene_list <- rownames(ora_gene_list)

  return(gene_list)
}

## Create a gene universe to use as the background for the GO term analysis 
gene_universe <- function (dge_data_set, key = 'symbol'){
  if (key == 'symbol'){
    gene_uni <- na.omit(rownames(dge_data_set))
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
oraFunc <- function (gene_list, background, ONT = "ALL") {
  go_enrich <- enrichGO(gene = gene_list,
                        OrgDb = "org.Dr.eg.db", 
                        keyType = 'SYMBOL',
                        ont = ONT,
                        universe = background,
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "fdr",
                        qvalueCutoff = 0.05)
  return(go_enrich)
}


## KEGG and Reactome Pathway analyses
## KEGG
pathway_ont <- function(gene_list){
    kk <- gseKEGG(geneList = gene_list,
                  keyType = 'ncbi-geneid',
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


# Load data from zfin for ontology level analysis
#zfin <- read.delim("https://zfin.org/downloads/wildtype-expression_fish.txt",
#                   sep = '\t', header = FALSE)
