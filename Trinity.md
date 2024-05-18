# Trinityの使い方


元のサイト：
https://github.com/trinityrnaseq/trinityrnaseq/tree/master/util


trinityrnaseq ON NIsp


Install と環境構築 on NIG super computer

	Condaで環境を作って動かす
	https://kazumaxneo.hatenablog.com/entry/2020/12/04/182413

### Installation & Setup
```sh
conda create -n trinity python=3.9
conda activate trinity
# Trinityのインストール： versionを指定したほうがいいらしい
mamba install -c conda-forge -c bioconda -y trinity=2.15

align_and_estimate_abundance.pl　-h 
helpが出ればOK
```
 [オリジナルファイル](https://github.com/trinityrnaseq/trinityrnaseq/releases) からのダウンロード，インストールする場合：

-----------------------------------
### 参考

## TPMの算出
```sh
$HOME/path/to/abundance_estimates_to_matrix.pl \
    --est_method RSEM \
    --gene_trans_map none\
    --out_prefix TPMmatrix \
    --cross_sample_norm TMM \
    $wd/Ctrl_2_S5_rsem.genes.results \
    $wd/Ctrl_3_S7_rsem.genes.results \
    $wd/Ctrl_4_S9_rsem.genes.results \
    $wd/MTZ_1_S4_rsem.genes.results \
    $wd/MTZ_2_S6_rsem.genes.results \
    $wd/MTZ_3_S8_rsem.genes.results   
```
abundance_estimates_to_matrix.plはfull pathを指定

Trinity 付属スクリプト：align_and_estimate_abundance.pl　を使う




#### abundance_estimates_to_matrix.pl

スクリプトで発現行列を出力、filter_low_expr_transcripts.plスクリプトで低発現転写産物をフィルタリングする [坂一馬先生の記事](https://kazumaxneo.hatenablog.com/entry/2021/12/25/212853
)
https://github.com/trinityrnaseq/trinityrnaseq/blob/master/util/abundance_estimates_to_matrix.pl

うまく走らない:
abundance_estimates_to_matrix.pl line 153.
->  abundance_estimates_to_matrix.plをフルパスで指定する



#### align_and_estimate_abundance.pl
[坂一馬先生の記事](https://kazumaxneo.hatenablog.com/?page=1607097062)


7. リードカウントのマトリクスを作成する by Kanako Bessho-Uehara, PhD

https://kbessho512.wixsite.com/kanabu/7-%E3%83%AA%E3%83%BC%E3%83%89%E3%82%AB%E3%82%A6%E3%83%B3%E3%83%88%E3%81%AE%E3%83%9E%E3%83%88%E3%83%AA%E3%82%AF%E3%82%B9%E3%82%92%E4%BD%9C%E6%88%90%E3%81%99%E3%82%8B




Trinity 2.6.6 run in miniconda, line 196.error message #634

https://github.com/trinityrnaseq/trinityrnaseq/issues/634

abundance_estimates_to_matrix.pl Error #345
https://github.com/trinityrnaseq/trinityrnaseq/issues/345


使い方は，下記HPがわかりやすい
Kanako Bessho-Uehara, PhD
https://kbessho512.wixsite.com/kanabu/7-%E3%83%AA%E3%83%BC%E3%83%89%E3%82%AB%E3%82%A6%E3%83%B3%E3%83%88%E3%81%AE%E3%83%9E%E3%83%88%E3%83%AA%E3%82%AF%E3%82%B9%E3%82%92%E4%BD%9C%E6%88%90%E3%81%99%E3%82%8B
