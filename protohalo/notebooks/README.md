# Protohalo notebooks

IllustrisTNG300-1 の既存解析を読むための索引です。ファイル名の先頭番号は、Git履歴で確認できる初回追加順です（実際の作成日時とは限りません）。TNG100 の学習・試行用ノートブックは `examples/` に分け、別の連番にしています。


## 作成順の索引

| 番号 | Notebook | Git初回追加日（日本時間） |
|---|---|---|
| 01 | [01_tng300_z8_descendant_mass.ipynb](01_tng300_z8_descendant_mass.ipynb) | 2026-09-11 |
| 02 | [02_tng300_descendant_with_galaxy_properties.ipynb](02_tng300_descendant_with_galaxy_properties.ipynb) | 2026-09-11 |
| 03 | [03_tng300_descendant_with_environment.ipynb](03_tng300_descendant_with_environment.ipynb) | 2026-09-25 |
| 04 | [04_tng300_environment_robustness.ipynb](04_tng300_environment_robustness.ipynb) | 2026-09-25 |
| 05 | [05_tng300_protohalo_lagrangian_regions.ipynb](05_tng300_protohalo_lagrangian_regions.ipynb) | 2026-09-25 |
| 06 | [06_tng300_z10_descendants_colibre_comparison.ipynb](06_tng300_z10_descendants_colibre_comparison.ipynb) | 2026-10-07 |
| 07 | [07_tng300_main_results_summary.ipynb](07_tng300_main_results_summary.ipynb) | 2026-10-08 |
| 08 | [08_tng300_descendant_predictive_performance.ipynb](08_tng300_descendant_predictive_performance.ipynb) | 2026-10-08 |
| 09 | [09_tng300_latent_host_posterior.ipynb](09_tng300_latent_host_posterior.ipynb) | 2026-10-14 |

`09` は外部の作業フォルダで作成され、2026-10-13にこのフォルダへ移動、2026-10-14にGitへ追加されました。
`examples/` の2冊はともに2026-09-02の同一コミットで初回追加され、作成の前後関係は不明です。便宜上、tutorialを01、mergertree_testを02としています。

## 最初に読む

[07_tng300_main_results_summary.ipynb](07_tng300_main_results_summary.ipynb) に、主要4冊の結果表と、z≈8 のハロー質量選択に情報を追加した場合の z=0 descendant host mass PDF をまとめています。図は実行番号 `[6]` のセルにあります。保存済みCSVだけで再実行でき、シミュレーション本体は不要です。

詳しい手法・図の読み方は [解析まとめ](../analysis_summary/tng300_analysis_summary_ja.md)、COLIBRE 論文との違いは [比較まとめ](../analysis_summary/tng300_vs_colibre_ja.md) を参照してください。

## z≈8 の主解析

| 順序 | Notebook | 役割 | 保存データ／前提 |
|---|---|---|---|
| 1 | [01_tng300_z8_descendant_mass.ipynb](01_tng300_z8_descendant_mass.ipynb) | 初期ハローの質量選択、中心subhaloから z=0 FoF host への対応付け、基準の質量分布 | `data/tng300_z8_M11_descendants.csv` を生成 |
| 2 | [02_tng300_descendant_with_galaxy_properties.ipynb](02_tng300_descendant_with_galaxy_properties.ipynb) | 中心銀河の恒星質量・SFR、5 cMpc の銀河数と descendant mass の関係 | 1のCSVとTNGカタログを使用。銀河特性のみの独立CSVは保存していない |
| 3 | [03_tng300_descendant_with_environment.ipynb](03_tng300_descendant_with_environment.ipynb) | 1–20 cMpc の球内の銀河数・ハロー数・質量和、条件付きPDFと回帰 | 1のCSVを読み、`data/tng300_z8_M11_descendants_environment.csv` を生成 |
| 4 | [04_tng300_environment_robustness.ipynb](04_tng300_environment_robustness.ipynb) | tracerの質量閾値、円柱投影、解析的Lagrangian scaleによる頑健性確認 | 1のCSVを読み、`data/tng300_environment_*_summary.csv` と `data/tng300_environment_cylinder_correlations.csv` を生成 |

1、3、4が保存データの生成経路です。2と3は1から分岐し、2は3の実行前提ではありません。主要4冊を再計算するには TNG300-1 本体と `illustris_python` が必要です。

## 補足解析・比較

| Notebook | 役割 | 実行時の入力 |
|---|---|---|
| [05_tng300_protohalo_lagrangian_regions.ipynb](05_tng300_protohalo_lagrangian_regions.ipynb) | z=0 のDM粒子IDを z≈8 へ照合し、protohaloの空間的広がりを調べる | 主解析のCSV、TNG300-1の粒子スナップショット。キャッシュは `cache/protohalo_particle_cache/` |
| [06_tng300_z10_descendants_colibre_comparison.ipynb](06_tng300_z10_descendants_colibre_comparison.ipynb) | z≈10 のサンプルを追跡し、COLIBRE Fig.7/8に対応する図と分類を比較 | TNG300-1、必要に応じて z≈8 の基準CSV。出力は `data/tng300_z10_colibre/` |

z≈10 の比較は z≈8 の環境解析とは別の枝です。図の保存・PDF描画の例は、比較ノートブックの z≈10 / z≈8 halo-mass PDF 比較セルを参照してください。セルの実行番号は再実行で変わります。

## 学習・試行用

- [examples/01_tutorial.ipynb](examples/01_tutorial.ipynb)：IllustrisTNG のカタログ、粒子、merger tree を扱う TNG100 チュートリアル。
- [examples/02_mergertree_test.ipynb](examples/02_mergertree_test.ipynb)：TNG100 の merger tree 操作の試行。

これらは主解析の依存先ではありません。元のコード・出力を保存しており、実行する場合は各セルの `basePath` を使用環境の TNG100 のパスに設定してください。

## 実行場所とデータ保存先

主要なTNG300ノートブックは、このリポジトリ内の作業ディレクトリから `protohalo` を検出します。主解析のCSVは `protohalo/data/`、粒子・treeのキャッシュは `protohalo/cache/` に保存します。TNG300-1本体が別の場所にある場合は、カーネル起動前に環境変数 `TNG_OUTPUT_DIR` をその `output` ディレクトリへ設定してください。

保存済みのセル出力は過去の計算結果です。整理時には再計算していないため、一部の出力に旧保存パスが表示されていても、現在のコードの保存先は上記のディレクトリです。CSVを生成するセルを再実行すると、同名の既存CSVを上書きします。

## 後続の予測性能検証で再利用するもの

新しい予測性能検証ノートブックはまだ作成していません。入力の候補は以下の既存データです。

| 保存ファイル（`../data/`） | 内容 | 用途 |
|---|---|---|
| `tng300_z8_M11_descendants.csv` | 3,099行、初期ハローと最終ホストのID・質量・追跡状態 | 基準サンプル、初期質量、目的変数、照合 |
| `tng300_z8_M11_descendants_environment.csv` | 3,099行・52列、同じサンプルの中心銀河特性と8半径の環境量 | 初期質量だけのモデルと環境を追加したモデルの比較 |
| `tng300_environment_robustness_summary.csv` | 16行、tracer閾値ごとの最適半径・相関 | 既存の探索結果の参照 |
| `tng300_environment_cylinder_correlations.csv` | 24行、円柱サイズごとの相関 | 投影効果の参照。各ハローの円柱内環境量は含まない |
| `tng300_environment_lagrangian_summary.csv` | 5行、最終質量binごとの解析的スケール | 環境半径の物理解釈 |

個体データの照合キーは `GroupID_z8`、目的変数は `logM200c_z0`、初期質量は `logM200c_z8` です。同じ `GroupID_z0` を持つ複数の初期ハローがあるため、学習・評価を分ける際には最終ホストの共有を考慮する必要があります。環境半径やモデルの選択も学習側で行い、評価側からの情報流入を避けます。

現在の要約PDFは環境などで選んだ部分集団の分布であり、未使用データでの予測精度を直接測る図ではありません。保存カタログから始めれば、球内環境量の予測性能検証のために粒子照合やSubLink追跡をやり直す必要はありません。

## 潜在ホスト事後分布の保存先

`09_tng300_latent_host_posterior.ipynb` の標本CSV・追跡SQLite・推論結果・図は、すべて `protohalo/cache/tng300_bayesian_latent_results/` 以下へ保存します。リポジトリ直下・`protohalo/`・`notebooks/` のどこから実行しても同じ保存先を検出します。`protohalo/cache/` はGitの追跡対象外です。
