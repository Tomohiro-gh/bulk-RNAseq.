library(dittoSeq)

## Cellcycle geneについて，Heatmapを書いてみる

# CellCycleのGOにannotateされている遺伝子を読み込む
GO0007049 <- read.xlsx("/Users/tomohiro/Dropbox/FukuharaLab_Res/Database/GO_Genes_zf/GO0007049_zf_cellcycle.xlsx", rowNames = FALSE, startRow = 2)　
  head(GO0007049)
GO0007049 <- GO0007049 %>% 
    distinct(bioentity_label, .keep_all=TRUE) %>%
    pull(bioentity_label)
  length(GO0007049) #695
  
# geneに変換

  
# GO0007049を並び替える
res_wld <- results(dds_wld)
res_wld <- res_wld[order(res_wld$padj), ]
 
vsd <- vst(dds_wld, blind=FALSE)
vst <- assay(vst)
vst <- as.data.frame(vst)
vst_sig <- vst[rownames(vst) %in% significant_gene_names,]





  
dds_wld <- readRDS("FlowvsStop_DESeq2_wald.rds")




#############################333
d <- dittoHeatmap(dds, 
                  genes=,
                  annot.by = "Group",
                  complex = TRUE,
                  column_split = dds_wld$Group)

png("heatmap_GO0007049_cellcycle.png",res=300, height=2000, width=2000)
d
dev.off()  


# Heatmpa
select <- order(rowMeans(counts(dds_wld, normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds_wld)[,c("Group","sample.names")])

png("heatmap_GO0007049_cellcycle.png", height = 600, width = 600)
pheatmap(assay(rld)[select,],
         cluster_rows=FALSE,
         show_rownames=FALSE,
         cluster_cols=FALSE,
         annotation_col=df,
         cellwidth = 40,
         cellheight = 20,
         fontsize = 20)
dev.off()




# https://support.bioconductor.org/p/133313/
normalized_data <- subset(counts(dds,normalized=T), rownames(counts(dds,normalized=T)) %in% significant_gene_names)

scaled_data <- t(scale(t(normalized_data)))

         