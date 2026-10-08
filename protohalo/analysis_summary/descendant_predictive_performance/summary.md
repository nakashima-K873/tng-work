### 自動計算summary（halo-level、95% paired host-bootstrap CI）

A baseline: RMSE **0.629 dex**, Brier **0.2172**。

- B Galaxy: RMSE 0.629 dex; 改善率 0.0% [-0.5, 0.5]。Brier 0.2160; 改善率 0.6% [-0.4, 1.4]。回帰/分類半径 対象外/対象外 cMpc。区間幅 1.379 dex [1.371, 1.388]、coverage 0.692 [0.635, 0.755]。

- C Galaxy environment: RMSE 0.531 dex; 改善率 15.7% [10.1, 21.1]。Brier 0.1537; 改善率 29.2% [18.8, 39.4]。回帰/分類半径 15/15 cMpc。区間幅 0.996 dex [0.978, 1.014]、coverage 0.671 [0.617, 0.723]。

- D Halo environment: RMSE 0.416 dex; 改善率 34.0% [28.8, 38.9]。Brier 0.1174; 改善率 46.0% [36.2, 55.0]。回帰/分類半径 10/10 cMpc。区間幅 0.710 dex [0.696, 0.723]、coverage 0.694 [0.648, 0.740]。

- E Galaxy + environment: RMSE 0.530 dex; 改善率 15.8% [10.5, 20.9]。Brier 0.1530; 改善率 29.6% [19.0, 39.7]。回帰/分類半径 15/15 cMpc。区間幅 0.995 dex [0.977, 1.012]、coverage 0.662 [0.607, 0.714]。

- F Full simulation: RMSE 0.417 dex; 改善率 33.8% [28.7, 38.6]。Brier 0.1167; 改善率 46.3% [36.6, 55.6]。回帰/分類半径 10/10 cMpc。区間幅 0.709 dex [0.695, 0.722]、coverage 0.690 [0.647, 0.738]。

Aの区間幅 1.395 dex、coverage 0.715。名目coverageは0.68で、幅の縮小だけを成功とみなさない。

C: test RMSEについて改善のCIは正、Brierについて改善のCIは正。

D: test RMSEについて改善のCIは正、Brierについて改善のCIは正。

training/testの最終ホスト共有はassertで0を確認した。unique-descendantとの比較はmodel_comparisonとFig.6を参照。異なるホスト間の空間相関、単一simulation volume、モデル選択の変動は残る。

A calibration: binned ECE=0.063。N≥30の2 bin中、平均予測確率が観測率のhost-bootstrap 95% CI外にあるbinは0。これは多重比較を補正しない診断で、完全なcalibrationを証明する検定ではない。

C calibration: binned ECE=0.052。N≥30の6 bin中、平均予測確率が観測率のhost-bootstrap 95% CI外にあるbinは0。これは多重比較を補正しない診断で、完全なcalibrationを証明する検定ではない。

D calibration: binned ECE=0.026。N≥30の5 bin中、平均予測確率が観測率のhost-bootstrap 95% CI外にあるbinは0。これは多重比較を補正しない診断で、完全なcalibrationを証明する検定ではない。

F calibration: binned ECE=0.019。N≥30の5 bin中、平均予測確率が観測率のhost-bootstrap 95% CI外にあるbinは0。これは多重比較を補正しない診断で、完全なcalibrationを証明する検定ではない。

unique-descendantのcoverageは主sampleと別に確認する。観測selectionの変化によって区間の較正も変わり得る。

Cluster確率の信頼性はFig.3とcalibration_bins.csvの観測率・CI・標本数から評価する。Brierの改善は完全なcalibrationを保証しない。Mock JWSTへの拡張には、halo mass proxyの誤差、観測selection、不完全性、恒星質量/SFRの推定誤差、投影・photo-z効果を取り入れ、Nhalo benchmarkを観測可能量と区別した再検証が必要。