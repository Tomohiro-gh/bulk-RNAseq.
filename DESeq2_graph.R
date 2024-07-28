library(DESeq2)


##########  plotCount arrange version  ############
# Script : 220508_DESeq_lrt.R

# DESeq plot countの自前
# required
genes.df # gene symbol <-> ensembl id 対応表
#marker geneはSYMBOLに名前が入っていないので手動で入れる
  genes.df.modif <- genes.df
  genes.df.modif$SYMBOL[genes.df.modif$ENSEMBL=="EGFP"] <- "EGFP"
  genes.df.modif$SYMBOL[genes.df.modif$ENSEMBL=="mScarletI-C1-2A-epNTR"] <- "mScarletI-C1-2A-epNTR"

### Function
myCountPlot <- function(deseqobject, GOI){
  # GOIが１つの時用に genesymbolのvectorへ入れておく
  genesymbol <- GOI
  for (i in 1:length(genesymbol)){
    #genesymbolで与えたものをensemblへ変換する
    ensemblid <- genes.df.modif[genes.df.modif$SYMBOL %in% genesymbol[i], ] %>%
      pull(ENSEMBL)
    #変換した geneが含まれていれば，グラフを角
    if (length(ensemblid) == 1){
      # DESeq内の plotCounts関数： returnData Tでggplotへわたす
      p <- plotCounts(deseqobject, gene=ensemblid, intgroup="Group", returnData=TRUE)
      
      # ggplotで整形
      ggplot(p, aes(x=Group, y=count, color=Group)) + 
        geom_point(size=4, position=position_jitter(w=0.12,h=0)) +
        # color はjco palletを参照
        scale_color_manual(values=c(Ctrl="#646970",MTZ2d="#EFC000",MTZ7d="#0073C2")) +
        ylim(0, max(p$count)*1.1) +
        # scale_y_log10(breaks=c(25,100,400))
        theme_minimal() +
        labs(title = genesymbol[i],
             x=NULL,
             y="Nornalized Count") +
        theme(plot.title = element_text(hjust=0.5, size=32, face="bold"),
              axis.text.x = element_text(size=24, face="bold"),
              axis.text.y = element_text(size=24, face="bold"),
              axis.title.x = element_text(hjust=0.5, size= 16, face="bold"), 
              axis.title.y = element_text(hjust=0.5, size= 24, face="bold"),
              legend.position = "none")
      
      ggsave(paste0("Count(DESeq2)_",genesymbol[i], ".png"), dpi=300, width=6, height=6)
      
    } else print(paste0("Warning: the data does not include ", genesymbol[i]))
    
  }
  
}


# ex) 
  ECs <- c("kdrl","cdh5","fli1a","pecam1", "tek")
  myCountPlot(dds_lrt, ECs) # OK
