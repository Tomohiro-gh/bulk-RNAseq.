## 公共データベースからSRRをダウンロードするまで


###1 SRA toolkitのダウンロード




--------------
####2 SRRのダウンロード
####2-1 : fasterq-dumpを使用 [参考：京橋インフォマティシャン](http://xn--sskume-h43e.net/entry/2022/01/09/205032)

```sh
fasterq-dump  SRRxxxxxx
fasterq-dump --split-files SRRxxxxxx
```
これでIndexリードがうまくわかれない場合がある-> [解決策](https://github.com/ncbi/sra-tools/issues/399)

例えばscRNAseqの場合の３つのファイルが存在する場合は下記のようにするとうまくいく 25/01/08
```sh
fasterq-dump SRRxxxx --include-technical -S

  ## -S|--split-files:  write reads into different files
  ## --include-technical:   include technical reads 
```

・　注意：fasterq-dumpは .gzの圧縮では別途やる必要がある
実例のコード
```sh
SEQLIBS=(SRR23344352 SRR23344353)

cd $HOME

for seqlib in ${SEQLIBS[@]}; do

    fasterq-dump ${seqlib} --split-files  --include-technical | gzip

done

```


