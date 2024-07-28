library(DESeq2)
library(tidyverse) # tibbleのために使用

# DESeq2 with 3 groups
dds <- DESeqDataSetFromTximport(txi, meta, design = ~ Group)
  dds
    # class: DESeqDataSet 
    # dim: 25434 10 
    # metadata(1): version
    # assays(2): counts avgTxLength
    #   rownames(25434): EGFP ENSDARG00000000001 ... ENSDARG00000117204
    # mScarletI-C1-2A-epNTR
    # rowData names(0):
    #   colnames(10): Cont-1 Cont-2 ... MTZ7d-2 MTZ7d-3
    # colData names(2): id Group
  
  dds$Group <- relevel(dds$Group, ref = "Ctrl")
  
  
  
dds_lrt <- DESeq(dds, test="LRT", reduced=~1)
res_lrt <- results(dds_lrt)
  head(res_lrt)
  summary(res_lrt) #summary of results
  
  #pvalueの低い順に
  res_lrt <- res_lrt[order(res_lrt$padj),]
  




##### QC ##########
# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)
### Plot PCA 
plotPCA(rld, intgroup="Group", add.labels = TRUE)


rld_mat <- assay(rld)
### Compute pairwise correlation values
rld_cor <- cor(rld_mat)    ## cor() is a base R function
pheatmap(rld_cor)

head(rld_cor)   ## check the output of cor(), make note of the rownames and colnames


# ref) https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html


# Subset the LRT results to return genes with padj < 0.05
sig_res_lrt <- res_lrt %>%
  data.frame() %>%
    rownames_to_column(var="gene") %>% 
      as_tibble() %>% 
        filter(padj < 0.05)

# Get sig gene lists
siglrt_genes <- sig_res_lrt %>% 
  pull(gene)
  
  length(siglrt_genes)

# Compare to numbers we had from Wald test
  nrow(sigOE)
  nrow(sigKD)
  

# Subset results for faster cluster finding (for classroom demo purposes)
  clustering_sig_genes <- sig_res_lrt %>% arrange(padj)
  
  
  # Obtain rlog values for those significant genes
  cluster_rlog <- rld_mat[clustering_sig_genes$gene, ]
  
  # Use the `degPatterns` function from the 'DEGreport' package to show gene clusters across sample groups
  clusters <- degPatterns(cluster_rlog, metadata = meta, time = "Group", col=NULL)
