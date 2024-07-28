#tixmport
library(tximport)
# ref1) https://bioconductor.statistik.tu-dortmund.de/packages/3.3/bioc/vignettes/tximport/inst/doc/tximport.html

# ref2) https://nakazakimasahito.wordpress.com/2020/06/02/%E7%99%BA%E7%8F%BE%E8%A7%A3%E6%9E%906%E3%80%80deseq2/

# RSEM -> DESeq2へ渡すための

# ex) Exp132のcountdata
dir = "/Users/tomohiro/Dropbox/FukuharaLab_Res/Experiment/Exp132_SkinEC_PCablation_bulkRNAseq/2_GeneCount_v2_rsem"

sample.names <- c("Cont-1","Cont-2","Cont-3","MTZ2d-1","MTZ2d-2","MTZ2d-3","MTZ2d-4","MTZ7d-1","MTZ7d-2","MTZ7d-3")


# create a named vector pointing to the rsem files
files <- file.path(dir, paste0(sample.names, ".rsem.genes.results"))
#file.path関数
#,で区切ると，ディレクトリの階層になる．paste0は，separotorなしで使用する
names(files) <- sample.names

all(file.exists(files)) # fileが存在するか確認．TRUEが返って来ればOK

# filesを元に tximportでファイルを読み込む読み込む

txi <- tximport(files, type="rsem", txIn=F, txOut=F)
# gene countの場合は，txIn txOutどちらもFalse
# list型の物が生成される
# 確認
head(txi$counts)

# length=0を1に置換：これがないと下記のようにエラーが出る ref2参照
# error :    all(lengths > 0) は TRUE ではありません 
txi$length[txi$length==0] <- 1


# DESEQ用のmetadataの作成
id <- as.factor(colnames(txi$counts)) 
Group <- factor(c("Ctrl","Ctrl","Ctrl","MTZ2d","MTZ2d","MTZ2d","MTZ2d","MTZ7d","MTZ7d","MTZ7d"))
meta <- data.frame(id, Group)



#### ここからはDESeq2で
#DESeqへわたす
dds <- DESeqDataSetFromTximport(txi, meta, ~Group)

#Wald検定（wt）
dds_wt <- DESeq(dds)
res_wt <- results(dds_wt)
res_wt_naomit <- na.omit(res_wt) # NA を除外


# 尤度比検定（lrt）
dds_lrt <- DESeq(dds, test="LRT", reduced=~1)
res_lrt <- results(dds_lrt)
res_lrt_naomit <- na.omit(res_lrt)


