

## Sample dataの読み込み
```r
# edgeR
library(edgeR)
# gene id conversion
library(AnnotationDbi)
library(org.Dr.eg.db)
library(rtracklayer)
## Graphics 
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(patchwork)
library(RColorBrewer)
## data frame
library(openxlsx)
library(dplyr)

wd="path/to/mydirectory"
setwd(wd)

## rsemからのcount dataの読み込み
countdir = "/Path_to_CountData" # Directoryを指定
# sampleの名前
sample.names <- c("Ctrl1", "Ctrl2", "Ctrl3", "Drug1", "Drug2", "Drug3")
group <- c("Ctrl", "Ctrl", "Ctrl", "Drug", "Drug", "Drug")
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


# Reading rsem files with tximport
txi <- tximport(files_selected, type="rsem", txIn=F, txOut=F)
  head(txi$counts)
cts <- txi$counts
```

## 発現解析 case 1 : 2群比較 - biological replicateあり
```r
## edgeR test
y <- DGEList(counts = cts, group = group)
  y$samples
    levels(y$samples$group) 

# filter out lowly expressed genes using the following commands:
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
  # a CPM of 1 corresponds to a count of 6-7 in the smallest sample.
  # A requirement for expression in two or more libraries is used as the minimum number of samples in each group is two.

## Noramlization
y <- calcNormFactors(y)
design <- model.matrix(~group_selected)

# どちらかを実行
y <- estimateDisp(y, design)
y <- estimateCommonDisp(y, design)

# to estimate tagwise dispersions:
y <- estimateTagwiseDisp(y)

et <- exactTest(y)

```

## 発現解析　 case 2 : Glm
"glmQLFit() fits the negative binomial GLM for each tag and produces an object of class DGEGLM with some new components. This DGEGLM object can then be passed to glmQLFTest() to carry out the QL F-test."
"Alternatively, one can perform likelihood ratio test to test for differential expression. The testing can be done by using the functions glmFit() and glmLRT()."
```r
y <- estimateDisp(y, design)

y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

group <- factor(c(1,1,2,2,3,3))
design <- model.matrix(~group)
fit <- glmQLFit(y, design)

# To compare 2 vs 1:
qlf.2vs1 <- glmQLFTest(fit, coef=32)
# To compare 3 vs 1:
qlf.3vs1 <- glmQLFTest(fit, coef=3)
# To compare 3 vs 2:
qlf.3vs2 <- glmQLFTest(fit, contrast=c(0,-1,1))

qlf <- glmQLFTest(fit, coef=2:3)
topTags(qlf)

```
## 発現解析 case 3：　２群比較 - no biological replicate　(edgeR のmanualから)
シンプルにこの４行だけを通す．　bcvはmanualより抜粋

"Typical values for the common BCV (square-root- dispersion) for datasets arising from well-controlled experiments are 0.4 for human data, 0.1 for data on genetically identical model organisms or 0.01 for technical replicates. Here is a toy example with simulated data:"


```r
bcv <- 0.2
counts <- matrix(rnbinom(40,size=1/bcv^2,mu=10), 20,2)
y <- DGEList(counts=counts, group=1:2)
et <- exactTest(y, dispersion=bcv^2)
```

```r
## ENSEMBL geneをconversionする ------------------------------------------------
genes.df <- 
  AnnotationDbi::select(x = org.Dr.eg.db,
                        keys = rownames(DEtable),
                        keytype = "ENSEMBL",
                        columns = "SYMBOL") 
  colnames(genes.df) #[1] "ENSEMBL" "SYMBOL" 
  dim(genes.df) #[1] 14673     2

DEtable <- DEtable %>% mutate(ENSEMBL=rownames(.))
  head(DEtable)
  dim(DEtable)

## mergeする: merege(x, y) xの行数に合うようにmergeする
  DEtable.merge <- merge(DEtable, genes.df, by="ENSEMBL")
    dim(DEtable.merge) #[1] 13836     7
    colnames(DEtable.merge)

## 列順の変更  
DEtable.merge <- DEtable.merge %>% 
  dplyr::select(ENSEMBL, SYMBOL, everything())
    colnames(DEtable.merge) # [1] "ENSEMBL"        "SYMBOL"         "baseMean"       "log2FoldChange" "lfcSE"          "stat"           "pvalue"        "padj -> OK

# # rownames(genes.df) <- genes.df$ENSEMBL # <- これはduplicateで怒られる．．．
# # ensemblの重複行を削除する
# genes.df <- genes.df %>% dplyr::distinct(., ENSEMBL, .keep_all=TRUE)
#   rownames(genes.df) <- genes.df$ENSEMBL
#   dim(genes.df)# [1] 14523     2

    which(duplicated(DEtable.merge[[2]])) %>% length() #885の重複
## duplicateしたsymbolをensemblへ置き換え
DEtable.merge[[2]][which(duplicated(DEtable.merge[[2]]))] <- 
  DEtable.merge[[1]][which(duplicated(DEtable.merge[[2]]))]
    which(duplicated(DEtable.merge[[2]])) %>% length() #確認 -> ０ならOK

# NAはないか？
    which(is.na(DEtable.merge[[2]]))  %>% length() #1 -> EGFPがNAなので NAにENSEMBLを入れておく
DEtable.merge[[2]][which(is.na(DEtable.merge[[2]]))] <- 
  DEtable.merge[[1]][which(is.na(DEtable.merge[[2]]))]
    which(is.na(DEtable.merge[[2]])) %>% length() #0 #確認 -> ０ならOK
  #確認
  head(DEtable.merge)

##  ------------------------------------------------

```
