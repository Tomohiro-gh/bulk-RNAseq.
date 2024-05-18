# bulk-RNAseq. work flow


## Step1. Check fastq file - fastp & multiQC

```sh
SEQLIBS=(Ctrl1 Ctrl2 Ctrl3 Drug1 Drug2 Drug3)

mkdir ./rawfastq_QCreport

for seqlib in ${SEQLIBS[@]}; do
  fastp -i ./raw_fastq/${seqlib}_L001_R1_001.fastq.gz \
        -I ./raw_fastq/${seqlib}_L001_R2_001.fastq.gz \
        -h ./rawfastq_QCreport/${seqlib}_fatsp_report.html
done
```

-----------------
## Step2. Quality control - fastp & multiQC
```sh
mkdir -p ./afterQC

DIR_fatsq=$HOME/...../raw_fastq

for seqlib in ${SEQLIBS[@]}; do
  fastp -i $DIR_fatsq/${seqlib}_L001_R1_001.fastq.gz \
        -I $DIR_fatsq/${seqlib}_L001_R2_001.fastq.gz \
        -o ./afterQC/${seqlib}_R1_fastp.fastq.gz \
        -O ./afterQC/${seqlib}_R2_fastp.fastq.gz \
        -h ./afterQC/${seqlib}_fastp_report.html \
        -j ./afterQC/${seqlib}_fastp_report.json \
        -w 4 \
        -q 30 \
        -n 30 \
        -x 30 \
        -t 3 \
        -T 3 \
        -l 70\
        --correction \
          --overlap_len_require 30 \
          --overlap_diff_limit 5 \
          --overlap_diff_percent_limit 20
done

multiqc ./afterQC
```
##### fastp options
- -j, -h: report file, html and json format
- -w: worker thread number, default is 3
- -q: the quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified.
- -t: trimming how many q bases in tail for read1, default is 0 (int [=0])
- -x: enable polyX trimming in 3' ends.
- -n: nの数だけ同じ塩基が出てきたらそれを除去する
- -l: これ以下の長さのリードは捨てる
- --correction でオーバーラップしたペアエンドリード領域のエラーコレクションのフラグを立てる


-----------------
## Step3. Mapping & Counts - STAR/RSEM
##### ペアエンドの場合
```sh
for seqlib in ${SEQLIBS[@]}; do

  STAR  --runThreadN 12 \
        --genomeDir $HOME/../STARindex \
        --readFilesIn $DIR_fastq/${seqlib}_R1_fastp.fastq.gz $DIR_fastq/${seqlib}_R2_fastp.fastq.gz \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode TranscriptomeSAM \
        --outFileNamePrefix $DIR_fastq_aligned/${seqlib}_STAR_      
          # STAR option
            # outSAMtype BAM SortedByCoordinate : sorted BAM files
            # readFilesCommand zcat:  to uncompress .gz files
            # quantMode TranscriptomeSAM : output SAM/BAM alignments to transcriptome into a separate file
            # 発現定量するには--quantModeの指定は必須

  rsem-calculate-expression -p 8 \
      --no-bam-output \
      --alignments \
      --paired-end \
      $DIR_fastq_aligned/${seqlib}_STAR_Aligned.toTranscriptome.out.bam \
      $HOME/../rsemindex \
      ${seqlib}_rsem
          # RSEM option
            # --estimate-rspd: リードがトランスクリプト全体にどのように分布しているか推定する。rspd (read start position distribution)。
done

# Counts matrixを作成する (遺伝子のカウントだけだったら.genes.resultsを使用する)
rsem-generate-data-matrix \
  Ctrl1_rsem.genes.results \
  Ctrl2_rsem.genes.results \
  Ctrl3_rsem.genes.results \
  Drug1_rsem.genes.results \
  Drug2_rsem.genes.results \
  Drug3_rsem.genes.results \
> Counts_alignment.tsv
```
-----------------
## Step4. Creating Expression matrix (TPM) - (Trinity)

```sh
$HOME/full_path_to/abundance_estimates_to_matrix.pl \
    --est_method RSEM \
    --gene_trans_map none\
    --out_prefix TPMmatrix \
    --cross_sample_norm TMM \
    $wd/Ctrl1_rsem.genes.results \
    $wd/Ctrl2_rsem.genes.results \
    $wd/Ctrl3_rsem.genes.results \
    $wd/Drug1_rsem.genes.results \
    $wd/Drug2_rsem.genes.results \
    $wd/Drug3_rsem.genes.results   
```



bulk-RNAseq. analysis
