# DEseq2 tutorial

-------------------
# 1) 2 group with biolocical replicates

## Step1 dataの読み込み
#### count data by rsem
#### data import by tximport

```r
library(DESeq2)
library(apeglm)
library(ashr)
library(tidyverse) # tibbleのために使用
library(ggrepel)
library(pheatmap)
library(AnnotationDbi)
library(org.Dr.eg.db)
library(rtracklayer)　# gene id をgene symbolに変えるために使用する
library(SummarizedExperiment)
library(openxlsx)
library(tximport)
library(RColorBrewer)
library(vsn)  # meanSdPlot で使用
library(patchwork)

wd="path/to/mydirectory"
setwd(wd)

##rsemからのcount dataの読み込み
countdir = "/Path_to_CountData" # Directoryを指定
# sampleの名前
sample.names <- c("Ctrl1","Ctrl2","Ctrl3","Drug1","Drug2","Drug3")
group <- c("Ctrl","Ctrl","Ctrl","Drug","Drug","Drug")
# create a named vector pointing to the rsem files
files <- file.path(countdir, paste0(sample.names, "_rsem.genes.results"))
  names(files) <- sample.names
  all(file.exists(files)) # TRUE　ならOK

## (optionnal ) >>>>>>>>>>>
## 特定のサンプルのみにしたい場合（dataを除く場合）
FOI <- c(1,3,4,6) #例2,5番目のサンプルを除外
files <- files[FOI]
sample.names <- sample.names[FOI]
(group <- group[FOI])
## <<<<<<<<<<<<<<<<<<<

# groupの順序を決めておく
group <- factor(x=group,
                levels = c("Ctrl", "Drug"))

# Reading rsem files with tximport
txi <- tximport(files, type="rsem", txIn=F, txOut=F)
  head(txi$counts)

# length=0を1に置換：これがないと下記のようにエラーが出る ref2参照
    ## https://support.bioconductor.org/p/92763/
txi$length[txi$length==0] <- 1
```

## Step2 DESeq objectの準備，　解析
```r 
# DESEQ用のmetadataの作成: 幾つの情報を持たせることができる
id <- sample.names
Group <- group
Trial <- rep(c("1st", "2nd", "3rd"),2)

meta <- data.frame(sample.names, Group, Trial)

## creating DEseq object
dds <- DESeqDataSetFromTximport(txi, 
                                colData = meta,
                                design = ~ Group)
      dds #確認１
      resultsNames(dds)　#確認2

# low count geneを削除: count10以上が少なくとも２サンプルで見られているものに限定
smallestGroupSize <- 2
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]


## -----------------------------------------------------------------------------
## Differential expression analysis
dds <- DESeq(dds)
  dds
  resultsNames(dds)

res <- results(dds)

```

## Ensembl Gene IDを symbolへ変換
・ Zebrafishの場合， ensembl IDしか付与されていないことがある

・ ensemblにしかないものはensembl　IDを付与


```r
genes.df <- 
  AnnotationDbi::select(x = org.Dr.eg.db,
                        keys = rownames(dds),
                        keytype = "ENSEMBL",
                        columns = "SYMBOL") 
  colnames(genes.df) #[1] "ENSEMBL" "SYMBOL" 
  dim(genes.df) #[1] 14166     2
  genes.df[which(duplicated(genes.df[[1]])),]
  ## ENSEMBL1に対してSYMBOLが複数の場合が多数ある -> 1:1になるうにロールアップする
  ## -> geneをロールアップする
  genes.df <- genes.df %>% 
    group_by(ENSEMBL) %>% 
    arrange(ENSEMBL, .by_group = TRUE) %>%
    summarise(
      across(everything(),
             ~paste0(na.omit(.x), collapse = ': ')))
  
    dim(gene.df)  #genes.df
  
res1 <- as.data.frame(res) %>% mutate(ENSEMBL=rownames(.))
    head(res1)
    dim(res1)

## mergeする: merege(x, y) xの行数に合うようにmergeする
res.merge <- merge(genes.df, res1, by="ENSEMBL")
  head(res.merge)
  dim(res.merge) #[1] 14006     8
  colnames(res.merge)
## 列順の変更  
res.merge <- res.merge %>% dplyr::select(ENSEMBL, SYMBOL, everything())
  colnames(res.merge) # [1] "ENSEMBL"        "SYMBOL"         "baseMean"       "log2FoldChange" "lfcSE"          "stat"           "pvalue"        "padj -> OK

# 
# # rownames(genes.df) <- genes.df$ENSEMBL # <- これはduplicateで怒られる．．．
# # ensemblの重複行を削除する
# genes.df <- genes.df %>% dplyr::distinct(., ENSEMBL, .keep_all=TRUE)
#   rownames(genes.df) <- genes.df$ENSEMBL
#   dim(genes.df)# [1] 14523     2

which(duplicated(res.merge[[2]])) %>% length() #SYMBOLの重複 -> 826
## duplicateしたsymbolをensemblへ置き換え
res.merge[[2]][which(duplicated(res.merge[[2]]))] <- 
  res.merge[[1]][which(duplicated(res.merge[[2]]))]
      which(duplicated(res.merge[[2]])) %>% length() #確認 -> ０ならOK

# NAはないか？
which(is.na(res.merge[[2]])) %>% length() #1 -> EGFPがNAなので NAにENSEMBLを入れておく
res.merge[[2]][which(is.na(res.merge[[2]]))] <- 
  res.merge[[1]][which(is.na(res.merge[[2]]))]
      which(is.na(res.merge[[2]])) %>% length() #0 #確認 -> ０ならOK
res.merge[1,2] <- "EGFP"
      
## ddsのrownamesを変更
  head(res.merge)
rownames(dds) <- res.merge$SYMBOL
```

### Dounstream analysis 1: MA plot
```r

```


### Dounstream analysis 2: Sample Distance & PCA
```r


```

### Dounstream analysis 3: Volcano plot
