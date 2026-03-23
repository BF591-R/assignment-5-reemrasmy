library('tidyverse')
library('SummarizedExperiment')
library('DESeq2')
library('biomaRt')
library('testthat')
library('fgsea')

#' Function to generate a SummarizedExperiment object with counts and coldata
#' to use in DESeq2
#'
#' @param csv_path (str): path to the file verse_counts.tsv
#' @param metafile (str): path to the metadata sample_metadata.csv
#' @param selected_times (list): list of sample timepoints to use
#' 
#'
#' @return SummarizedExperiment object with subsetted counts matrix
#'   and sample data. Ensure that the timepoints column used as input 
#'   to the model design has 'vP0' set as the reference factor level. Your 
#'   colData dataframe should have columns named samplename and timepoint.
#' @export
#'
#' @examples se <- make_se('verse_counts.tsv', 'sample_metadata.csv', c('vP0', 'vAd'))
make_se <- function(counts_csv, metafile_csv, selected_times) {
  
  counts <- read_tsv(counts_csv)
  metadata <- read_csv(metafile_csv)
  
  # Subet the metadata for the portion we are using 
  subset_meta <- metadata %>%
    filter(timepoint %in% selected_times) %>%
    mutate(timepoint = factor(timepoint, levels=c("vP0", "vAd"))) %>%
    dplyr::select(samplename, timepoint)
  
  # Subset the counts table using the subsetted metadata 
  # all_of() -- ; 1 is the gene column, instead of hardcoding "gene"
  subset_counts <- counts %>%
    dplyr::select(1, all_of(subset_meta$samplename))
  
  # Create the counts matrix
  counts_matrix <- subset_counts %>%
    column_to_rownames("gene") %>%
    as.matrix()
  
  # Create SummarizedExperiment object 
  se <- SummarizedExperiment(
    assays = list(counts = counts_matrix),
    colData = subset_meta
  )
  
  return(se)
}

#' Function that runs DESeq2 and returns a named list containing the DESeq2
#' results as a dataframe and the dds object returned by DESeq2
#'
#' @param se (obj): SummarizedExperiment object containing counts matrix and
#' coldata
#' @param design: the design formula to be used in DESeq2
#'
#' @return list with DESeqDataSet object after running DESeq2 and results from
#'   DESeq2 as a dataframe
#' @export
#'
#' @examples results <- return_deseq_res(se, ~ timepoint)
return_deseq_res <- function(se, design) {
  # Create the dds object by converting se to DESeqDataSet object
  dds <- DESeqDataSet(se, design = design) %>%
    DESeq()
  # Extract results table as a df
  res_df <- results(dds) %>%
    as.data.frame()
  
  # Return results and DESeqDataSet as list
  deseq_out <- list(dds = dds, results = res_df)
  
  return(deseq_out)
}

#' Function that takes the DESeq2 results dataframe, converts it to a tibble and
#' adds a column to denote plotting status in volcano plot. Column should denote
#' whether gene is either 1. Significant at padj < .10 and has a positive log
#' fold change, 2. Significant at padj < .10 and has a negative log fold change,
#' 3. Not significant at padj < .10. Have the values for these labels be UP,
#' DOWN, NS, respectively. The column should be named `volc_plot_status`. Ensure
#' that the column name for your rownames is called "genes". 
#'
#' @param deseq2_res (df): results from DESeq2 
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return Tibble with all columns from DESeq2 results and one additional column
#'   labeling genes by significant and up-regulated, significant and
#'   downregulated, and not significant at padj < .10.
#'   
#' @export
#'
#' @examples labeled_results <- label_res(res, .10)
#' 
label_res <- function(deseq2_res, padj_threshold) {
  deseq_res_tib <- deseq2_res %>%
    tibble::rownames_to_column("genes") %>%
    as_tibble() %>%
    
    mutate(
      volc_plot_status = case_when(
        is.na(padj) | is.na(log2FoldChange) ~ "NS", 
        padj < padj_threshold & log2FoldChange > 0 ~ "UP",
        padj < padj_threshold & log2FoldChange < 0 ~ "DOWN",
        TRUE ~ "NS"
      )
    )
  return(deseq_res_tib)
}

#' Function to plot the unadjusted p-values as a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#'
#' @return ggplot: a histogram of the raw p-values from the DESeq2 results
#' @export
#'
#' @examples pval_plot <- plot_pvals(labeled_results)
plot_pvals <- function(labeled_results) {
  
  pval_plot <- labeled_results %>%
    filter(!is.na(pvalue)) %>%
    ggplot(aes(x = pvalue)) +
    geom_histogram(bins = 50, fill="lightblue") +
    labs(
      title = "Histogram of raw p-values",
      x = "Raw p-value",
      y = "Gene count"
    ) +
    
    theme_minimal()
  
  return(pval_plot)
}

#' Function to plot the log2foldchange from DESeq2 results in a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return ggplot: a histogram of log2FC values from genes significant at padj 
#' threshold of 0.1
#' @export
#'
#' @examples log2fc_plot <- plot_log2fc(labeled_results, .10)
plot_log2fc <- function(labeled_results, padj_threshold) {
    
  log2fc_plot <- labeled_results %>%
    filter(!is.na(log2FoldChange), !is.na(padj), padj < padj_threshold) %>%
    ggplot(aes(x = log2FoldChange)) +
    geom_histogram(bins = 50, fill="lightcoral") +
    labs(
      title = "Distribution of log2 fold changes for significant genes",
      x = "log2FoldChange",
      y = "Gene count"
    ) +
    
    theme_minimal()
  
  
  return(log2fc_plot)
}

#' Function to make scatter plot of normalized counts for top ten genes ranked
#' by ascending padj
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param dds_obj (obj): The object returned by running DESeq (dds) containing
#' the updated DESeqDataSet object with test results
#' @param num_genes (int): Number of genes to plot
#'
#' @return ggplot: a scatter plot with the normalized counts for each sample for
#' each of the top ten genes ranked by ascending padj
#' @export
#'
#' @examples norm_counts_plot <- scatter_norm_counts(labeled_results, dds, 10)
scatter_norm_counts <- function(labeled_results, dds_obj, num_genes){
  
  # extract sample metadata
  sample_metadata <- as.data.frame(colData(dds_obj))
  
  # find top genes by smallest padj
  top_genes <- labeled_results %>%
    filter(!is.na(padj)) %>%
    arrange(padj) %>%
    slice_head(n = num_genes) %>%
    pull(genes)
  
  # extract normalized counts from dds
  norm_counts <- counts(dds_obj, normalized = TRUE)
  
  # subset counts to top genes
  top_gene_counts <- norm_counts[top_genes, ] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("genes") %>%
    as_tibble() %>%
    pivot_longer(
      cols= -genes,
      names_to="samplename",
      values_to="counts"
    ) %>%
    left_join(sample_metadata, by = "samplename")
  
  
  # make scatter plot
  scatt_plot <- ggplot(top_gene_counts, aes(x = samplename, y = counts, color = timepoint)) +
    geom_point() +
    facet_wrap(vars(genes)) +
    labs(
      title = "Normalized counts for top differentially expressed genes",
      x = "Sample name",
      y = "Normalized counts"
    )
  
  return(scatt_plot)
}

#' Function to generate volcano plot from DESeq2 results
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#'
#' @return ggplot: a scatterplot (volcano plot) that displays log2foldchange vs
#'   -log10(padj) and labeled by status
#' @export
#'
#' @examples volcano_plot <- plot_volcano(labeled_results)
#' 
plot_volcano <- function(labeled_results) {
  
  volcano_res <- labeled_results %>%
    filter(!is.na(padj)) %>%
    mutate(neg_log10_padj = -log10(padj))
  
  volcano_plot <- ggplot(volcano_res, aes(x = log2FoldChange, y = neg_log10_padj, color = volc_plot_status)) +
    geom_point() +
    labs(
      title = "Volcano plot of differential expression results",
      x = "Log2 fold change",
      y = "-log10(adjusted p-value)",
      color = "Status"
    ) +
    theme_minimal()
  
  return(volcano_plot)
}

#' Function to generate a named vector ranked by log2FC descending
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param id2gene_path (str): Path to the file containing the mapping of
#' ensembl IDs to MGI symbols
#'
#' @return Named vector with gene symbols as names, and log2FoldChange as values
#' ranked in descending order
#' @export
#'
#' @examples rnk_list <- make_ranked_log2fc(labeled_results, 'data/id2gene.txt')

make_ranked_log2fc <- function(labeled_results, id2gene_path) {
  gene_names <- read_tsv(id2gene_path,col_names = c("genes", "symbols"))
  
  join_symbols <- left_join(labeled_results, gene_names, by = "genes") %>%
    filter(!is.na(log2FoldChange), !is.na(symbols))
  
  gene_fc_vec <- setNames(join_symbols$log2FoldChange, join_symbols$symbols)
  ranked_vec <- sort(gene_fc_vec, decreasing=TRUE)
  
  return(ranked_vec)
}

#' Function to run fgsea with arguments for min and max gene set size
#'
#' @param gmt_file_path (str): Path to the gene sets of interest in GMT format
#' @param rnk_list (named vector): Named vector generated previously with gene 
#' symbols and log2Fold Change values in descending order
#' @param min_size (int): Minimum number of genes in gene sets to be allowed
#' @param max_size (int): Maximum number of genes in gene sets to be allowed
#'
#' @return Tibble of results from running fgsea
#' @export
#'
#' @examples fgsea_results <- run_fgsea('data/m2.cp.v2023.1.Mm.symbols.gmt', rnk_list, 15, 500)
run_fgsea <- function(gmt_file_path, rnk_list, min_size, max_size) {
  
  # returns a list of pathways from the gmt file (input) 
  pathways <- gmtPathways(gmt_file_path)
  
  # removing the duplicate gene names (failed in previous test)
  rnk_list <- rnk_list[!duplicated(names(rnk_list))]
  
  # run fgsea and convert results to a tibble
  fgsea_tib <- fgsea(pathways, rnk_list, minSize = min_size, maxSize = max_size) %>%
    as_tibble()
  
  return(fgsea_tib)
}

#' Function to plot top ten positive NES and top ten negative NES pathways
#' in a barchart
#'
#' @param fgsea_results (tibble): the fgsea results in tibble format returned by
#'   the previous function
#' @param num_paths (int): the number of pathways for each direction (top or
#'   down) to include in the plot. Set this at 10.
#'
#' @return ggplot with a barchart showing the top twenty pathways ranked by positive
#' and negative NES
#' @export
#'
#' @examples fgsea_plot <- top_pathways(fgsea_results, 10)
top_pathways <- function(fgsea_results, num_paths){
  
  top_pos_paths <- fgsea_results %>% 
    arrange(padj) %>%  
    filter(NES > 0) %>%
    slice_max(NES, n=num_paths)
  
  top_neg_paths <- fgsea_results %>%
    arrange(padj) %>%
    filter(NES < 0) %>%
    slice_max(NES, n=num_paths)
  
  top_paths <- bind_rows(top_pos_paths, top_neg_paths)
  
  fgsea_plot <- top_paths %>%
    mutate(pathway = forcats::fct_reorder(pathway, NES)) %>%
    ggplot() +
    geom_bar(aes(x=pathway, y=NES, fill = padj < .25), stat='identity') +
    scale_fill_manual(values = c('TRUE' = 'red', 'FALSE' = 'blue')) + 
    theme_minimal() +
    ggtitle('Top Enriched Pathways from FGSEA (vP0 vs vAd Heart Samples)') +
    ylab('Normalized Enrichment Score (NES)') +
    xlab('') +
    coord_flip()
  
  return(fgsea_plot)
}

