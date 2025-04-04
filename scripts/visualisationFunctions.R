sample_dist <- function(df, sample_id){
  
  # To plot the distribution of the data of any chosen sample
  # Inputs: 
          # 1. df, a normalised data frame containing the read count data
          # 2. sample_id, any chosen column of the count data frame here
  # Outputs: Visualisation of the sample distribution 
  plot1 <- ggplot(data = df) +
    geom_histogram(aes(x=df$sample_id), stat = "bin", bins = 200) +
    # zoom in to see the distribution more clearly
    xlim(-2, 50) +
    xlab("Raw expression counts") +
    ylab("Number of genes") +
    theme_classic() +
    theme(
      text = element_text(family = "Arial", face = "bold", size = 30),
      axis.text = element_text(family = "Arial" ,face = "bold",size = 30),
      axis.text.y = element_text(family = "Arial", face = "bold", size = 30),
      axis.title.y = element_text(family = "Arial", face = "bold", size = 30),
      axis.title.x = element_text(family = "Arial", face = "bold", size = 30)
    )
  
  return(plot1)
}

mean_var_plot <- function(df){
  mean_counts <- apply(countData[, 1:length(colnames(df))], 1, mean)
  variance_counts <- apply(df[, 1:length(colnames(df))], 1, var)
  df_count <- data.frame(mean_counts, variance_counts)
  
  ## Plot here
  plot2 <- ggplot(df_count) +
    geom_point(aes(x=mean_counts, y=variance_counts)) + 
    geom_line(aes(x=mean_counts, y=mean_counts, color="red")) +
    scale_y_log10() +
    scale_x_log10() +
    theme_classic() +
    theme(
      text = element_text(family = "Arial", face = "bold", size = 30),
      axis.text = element_text(family = "Arial" ,face = "bold",size = 30),
      axis.text.y = element_text(family = "Arial", face = "bold", size = 30),
      axis.title.y = element_text(family = "Arial", face = "bold", size = 30),
      axis.title.x = element_text(family = "Arial", face = "bold", size = 30)
    )
  return(plot2)
}

pca_plot <- function(dds, samples){
  # Transform counts for data visualization
  rld <- rlog(dds, blind=TRUE)
  # Extract the rlog matrix from the object
  rld_mat <- assay(rld)
  colnames(rld_mat) <- samples$Group
  pca <- prcomp(t(rld_mat))
  
  z <- plotPCA(rld, intgroup="Group")
  z + geom_label(aes(label = samples)) +
    theme_classic() +
    theme(
      #text = element_text(family = "Arial", face = "bold", size = 30),
      #axis.text = element_text(family = "Arial" ,face = "bold",size = 30),
      axis.text.y = element_text(family = "Arial", face = "bold", size = 30),
      axis.text.x = element_text(size = 12, face = 'bold'),
      axis.title.y = element_text(family = "Arial", face = "bold", size = 30),
      axis.title.x = element_text(family = "Arial", face = "bold", size = 30)
    )
  
  return(z)
}

################################################
############# Gene Counts plots   #############
################################################
gene_count <- function(df, gene) {
  plot_gene <- df %>%
    filter(Name == gene) %>%
    dplyr::select(starts_with("Sa")) %>%
    gather(key = "Samples", value = "Counts")  %>%
    mutate(Treatment = c("wt_unexposed","wt_unexposed","wt_unexposed",
                         "mut_unexposed","mut_unexposed","mut_unexposed",
                         "wt_exposed","wt_exposed","wt_exposed",
                         "mut_exposed","mut_exposed","mut_exposed")) %>%
    ggplot(aes(x = Treatment, y = Counts)) +
    geom_jitter(position = position_dodge(0.2)) +
    geom_boxplot(alpha=0.2) +
    scale_x_discrete(limits = c("wt_unexposed","wt_unexposed","wt_unexposed",
                                "wt_exposed","wt_exposed","wt_exposed",
                                "mut_unexposed","mut_unexposed","mut_unexposed",
                                "mut_exposed","mut_exposed","mut_exposed")) +
    labs(
      title = gene,
      y = "Normalised Counts",
      x = "Treatment (µM)"
    ) +
    theme_classic() +
    theme(
      axis.text.y = element_text(family = "Times-Roman", size = 15),
      axis.text.x = element_text(family = "Times-Roman",size = 15),
      axis.title.y = element_text(family = "Times-Roman",size = 20),
      axis.title.x = element_text(family = "Times-Roman",size = 20),
      legend.text = element_text(family = "Times-Roman",size = 15),
      legend.title = element_text(family = "Times-Roman",size = 20)
    )
  
  return(plot_gene)
}

####### valcano plots
volcanoPlot <- function(df){
  plot_1 <- EnhancedVolcano(df,
                            lab = rownames(df),
                            x = 'log2FoldChange',
                            y = 'pvalue',
                            axisLabSize = 30,
                            xlim = c(-3, 3),
                            ylim = c(0, 30),
                            pointSize = 2.0,
                            labSize = 8.5,
                            labFace = "plain",
                            title = 'DESeq2 results',
                            subtitle = 'Differential expression',
                            caption = 'FC cutoff, 1.333; p-value cutoff, 0.05',
                            legendPosition = "right",
                            legendLabSize = 30,
                            pCutoff = 0.05,
                            FCcutoff = 1.0)
  
  return(plot_1)
}

venPlot <- function(dgeset){
  p3 <- ggvenn(dgeset, 
               set_name_size = 14, 
               text_size = 7,
               fill_alpha = 0.25,
               stroke_size = 1.5,
               show_outside = "auto")
  
  return(p3)
}

### Visualisation functions
plotting_GO <- function(obj, mode = "gsea") {
  
  if (mode == "gsea"){
    plot <- dotplot(obj, 
                    color = "p.adjust",
                    showCategory=10, 
                    font.size = 20,
                    split=".sign") + facet_grid(.~.sign)
    
    return(plot)
  }else{
    plot <- dotplot(obj, 
                    color = "p.adjust",
                    showCategory=10, 
                    font.size = 20)
    
    return(plot)
  }
  plot <- dotplot(obj, 
                  color = "p.adjust",
                  showCategory=10, 
                  font.size = 17,
                  split=".sign") + facet_grid(.~.sign)
  
  return(plot)
}

emap_plotting <- function(obj){
  ego2 <- pairwise_termsim(obj,
                           method = "JC")
  plot <- emapplot(ego2)
  
  return(plot)
}

