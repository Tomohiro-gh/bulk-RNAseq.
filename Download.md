## 公共データベースからSRRをダウンロードするまで


###1 SRA toolkitのダウンロード




--------------
####2 SRRのダウンロード
####2-1 : fasterq-dumpを使用

```sh
fasterq-dump  SRRxxxxxx
fasterq-dump --split-files SRRxxxxxx
```
これでうまくいかない場合がある-> [解決策](https://github.com/ncbi/sra-tools/issues/399)

例えばscRNAseqの場合の３つのファイルが存在する場合は下記のようにするとうまくいく 25/01/08
```sh
fasterq-dump SRRxxxx --include-technical -S
```
