"""Regenerate the standalone TNG notebook from the tested implementation.

Only the generated .ipynb is needed on TNG JupyterLab. Run from any directory:
    python protohalo/tools/build_observational_notebook.py
"""
from pathlib import Path
import hashlib
import json
import textwrap

ROOT = Path(__file__).resolve().parents[1]
SOURCE = (ROOT/'observational_environment.py').read_text()
SOURCE_SHA = hashlib.sha256(SOURCE.encode()).hexdigest()
cells = []


def md(value):
    cells.append(dict(cell_type='markdown', metadata={}, source=textwrap.dedent(value).strip()+'\n'))


def code(value):
    cells.append(dict(cell_type='code', metadata={}, source=textwrap.dedent(value).strip()+'\n',
                      execution_count=None, outputs=[]))


def implementation(start=None, end=None):
    left = 0 if start is None else SOURCE.index(start)
    right = len(SOURCE) if end is None else SOURCE.index(end)
    code(SOURCE[left:right])


md(r"""
# TNG300：観測条件を通した環境情報と descendant mass の直接予測

**実行場所：IllustrisTNG公式JupyterLab。TNG300-1のgroup catalogとSubLinkを使用する。**
この1冊に読み込み・追跡・環境計算・推論の全関数を含む。追加のローカルPythonファイル、
09/10の実行済みkernelは不要。上から順に実行する。

研究の問いは、銀河の観測に対応する情報から $z=0$ のhost質量をどこまで予測できるか、
特に広い環境の情報がhalo massだけを使う予測へ何を加えるか、である。
ノートブック10の $M_\star$ と3D球内銀河数の比較を、以下の実験へ拡張する。

| 解析 | 目的 | 比較する情報 |
|---|---|---|
| 1. 模擬観測と完全性 | 3Dの環境を観測条件でどこまで回収できるか | 真の位置・質量、測光赤方偏移、分光、深い選択 |
| 2. 多スケール環境 | コアと外縁が異なる情報を持つか | 1,2,5,10,15 cMpcの銀河数・質量和・形状 |
| 3. Intermediate | 真の $M_z$ を与えても環境が予測に役立つか | $P(M_0\mid M_z,z)$ と $P(M_0\mid M_z,D,z)$ |
| 4. Finalと二段階法 | 同じ観測入力から直接予測する価値を調べる | $P(M_0\mid D,z)$ と $\int P(M_0\mid M_z,z)P(M_z\mid D,z)dM_z$ |
| 5. 領域の行先監査 | 銀河中心のラベルが領域全体を代表するか | 同一領域の銀河がどの $z=0$ hostへ配分されるか |

**対象とラベルを固定する。** 最初の実装は「観測で選んだanchor銀河」を単位とし、
ラベルはその銀河のSubLink descendantが属する $z=0$ FoFの $M_{200c}$。
各環境・追観測実験でanchorとラベルを変えない。対象は中心・衛星の両方で、
初期halo massや最終cluster massによる事前選択をしない。
銀河の密度ピークを探索・重複除去して作るprotocluster候補カタログとは異なる。
解析5は将来の「領域全体を入力する」モデルへ進むためのラベル監査であり、主解析のラベルは変更しない。

**観測モデルの範囲。** この実装は、誤差を加えたSED恒星質量、赤方偏移、質量依存の検出率、
分光抽出率による制御実験である。実際のJWST flux/EW選択、UVの放射輸送、
ALMAの線・塵flux、Lyαや中性ガスの尤度は含まない。SFR追加実験も瞬間SFRに誤差を加えたproxy。
これらを実観測の較正済み予測、あるいは観測時間あたりの最適化結果として解釈しない。

背景：[文献調査](../analysis_summary/highz_protocluster_observables_literature_ja.md)、
[Witten et al. 2026, §5.1](https://arxiv.org/html/2605.28930v1#S5.SS1)、
[Helton et al. 2024](https://arxiv.org/abs/2311.04270)、
[Baxter et al., updated v3](https://arxiv.org/abs/2504.03836v3)。
文献調査ファイルは、単独アップロード時にはなくても実行に影響しない。
""")

md(r"""
## 実行設定と保存先

既定は `MODE='tng'`、snapshot 8（実際の $z\simeq8.012$）。まずこの時刻で全処理を実行する。
`Config(snapshots=(6,7,...,16))` とすれば複数時刻にも対応するが、各時刻に独立した予測器を学習し、
CVの平均損失は時刻ごとに等重みとする。時刻間の補間・lightconeの生成は行わない。

`BASE_PATH=None` は `TNG_OUTPUT_DIR` と公式JupyterLabの標準パスを探索する。
見つからなければ `'/home/tnguser/sims.TNG/TNG300-1/output'` など実際の場所を指定する。
必要なpackageは `numpy,pandas,scipy,scikit-learn,h5py,matplotlib,illustris_python`。
インストールを自動実行しない。

`MODE='synthetic_demo'` はソフトウェアの動作確認専用。銀河位置・質量・行先を人工的に作り、
すべての表と図にdemo表示を付ける。TNG解析の代替・科学的な結果ではない。

出力は `protohalo/cache/tng300_observational_environment/` 以下。リポジトリから実行すれば既存の
`protohalo` を検出し、単独アップロードなら作業ディレクトリの下に作成する。
設定・入力ファイルのstat・実装のhashごとに別のrunを作り、09/10の結果を上書きしない。
""")
implementation(None, 'def load_tng_snapshot(')
code(f"""
from dataclasses import replace
from IPython.display import display
import matplotlib.pyplot as plt
import platform
import scipy
import sklearn

MODE = 'tng'  # 'synthetic_demo' is software verification only
BASE_PATH = None
# Optional: notebook-09 galaxies_all.csv with its manifest.json alongside.
# None means directly trace selected anchors; notebook 09 is not a prerequisite.
LEGACY_CSV = None
CFG = Config(snapshots=(8,))
IMPLEMENTATION_SHA256 = '{SOURCE_SHA}'

cwd = Path.cwd().resolve()
PROJECT_ROOT = next((p for p in [cwd, *cwd.parents] if (p/'protohalo'/'notebooks').is_dir()), cwd)
CACHE_ROOT = PROJECT_ROOT/'protohalo'/'cache'/'tng300_observational_environment'
if MODE not in ('tng', 'synthetic_demo'):
    raise ValueError('Unknown MODE')
DEMO = MODE == 'synthetic_demo'
LABEL_PREFIX = '[SYNTHETIC DEMO] ' if DEMO else 'TNG300-1 mock observing experiment: '
if DEMO:
    CFG = replace(CFG, n_trees=8, leaf_grid=(4,), ywidth_grid=(.2,), mwidth_grid=(.15,),
                  max_anchors_per_snapshot=400, cv_cap_per_snapshot=15, test_cap_per_snapshot=20,
                  min_train=15, min_validation=3, min_hosts=3, bootstrap_draws=30,
                  region_audit_n=5, bridge_grid_size=61)
CFG.validate()
print(LABEL_PREFIX+str(asdict(CFG)))
print('Cache:', CACHE_ROOT)
""")

md(r"""
## 入力：物理的な銀河集団と、独立した子孫ラベル

**目的。** 周辺銀河の観測可能性と、子孫を追えるかどうかを混同しない。
tracerは同時刻の `SubhaloFlag=1`、真の $M_\star\geq10^8M_\odot$ の全保存銀河。
環境の計算にtracerの将来質量・追跡成功・train/test所属を使わない。
この保存下限は観測完全性でも分解能収束の保証でもない。

**手法。** 恒星質量は `SubhaloMassType[:,4]`、SFRは `SubhaloSFR`。
位置はckpc/hからcMpc、`SubhaloVel`は物理的なpeculiar velocity [km/s]。
anchorの真の $M_z$ は現在のFoFの `Group_M_Crit200`。最終ラベルも $M_{200c}$ で統一する。
中心・衛星を区別する列を保存し、最初は同じ銀河重み付き母集団で比較する。

子孫はSubLink `RootDescendantID` の終点を読み、snapshot・SubfindID・tree IDを照合する。
`not_in_tree` と $z=0$ 以前の正常な終端を記録し、I/O errorを未追跡に置き換えない。
追跡はSQLiteへ100天体ごとにcommitし、通常のInterrupt時にもfinallyで保存する。
元のtreeファイルを更新した場合は別cacheになる。
任意指定の09 CSVはmanifestと元ファイルのstatが一致する場合だけ再利用する。

出典：[TNG field specifications](https://www.tng-project.org/data/docs/specifications/)、
[SubLink documentation](https://www.tng-project.org/data/docs/scripts/#sublink)。
""")
implementation('def load_tng_snapshot(', 'def mock_observe(')

md(r"""
## 解析1：同じ銀河集団を、浅い／深い観測と測光／分光へ写す

**目的。** 環境が持つ理想的な情報のうち、赤方偏移誤差・検出の不完全性・分光抽出で失われる情報を調べる。

**手法。** 各銀河に恒星質量のGaussian誤差（既定0.2 dex）を加える。
推定恒星質量に依存するlogistic検出率を用い、浅い／深い観測で同じ一様乱数を共有する。
分光tracerも同じ乱数で抽出し、deepをshallowの上位集合にする。
既定のtracer閾値は浅い観測が $\log M_{\star,obs}\geq9$、深い観測が8.5。
深い実験では分光抽出率も0.5から0.8へ増やすので、両方が変わる合同実験である。

anchorは浅い撮像で検出され、$\log M_{\star,obs}\geq9$ を満たす銀河から、ラベルを見る前に無作為抽出する。
最大2000天体/時刻のcapは計算量の制御で、質量や環境による上位選択ではない。
**すべてのanchorには別途systemic分光がある**と仮定し、誤差は既定50 km/s。
周辺tracerの測光誤差は $\sigma_z/(1+z)=0.02$、分光誤差は50 km/sとする。
catastrophic photo-z、flux依存の線検出、誤差間の相関は未実装。

SFR追加実験の量は $\log_{10}(\mathrm{SFR}+0.01)$ に0.3 dexの誤差を加えたproxy。
0.01はゼロSFRを有限値へ写す固定のfloorであり、観測上限値の尤度ではない。
""")
implementation('def mock_observe(', 'def environment_features(')

md(r"""
## 解析2：多スケールの銀河数、質量分布、空間形状

**目的。** 中心の情報を、外縁・周辺の銀河の情報から分けて評価する。
半径は1,2,5,10,15 cMpcを事前固定する。testの相関を見て最適半径を選ばない。

**3D参照。** 真の位置の $r<R$ 球。`SubhaloID`で自分だけを除き、同じFoFの他銀河を含む。
**観測に対応する環境。** 視線を既定z軸に固定し、$R_\perp<R$ の全LOSを検索する。
赤方偏移空間の座標は $u=H(z)x_\parallel/(1+z)+v_{pec}$ に測定誤差を加えたもの。
分光環境では $|\Delta u|<1200$ km/sのtracerを数える。
測光環境ではGaussian誤差の下で、この窓に入る確率を重みにする：

$$w_i=\Phi[(v_{max}-\Delta u_i)/\sigma_{pair}]-\Phi[(-v_{max}-\Delta u_i)/\sigma_{pair}],
\qquad N_{eff}(R)=\sum_{R_{\perp,i}<R}w_i.$$

箱は周期境界なので、Gaussianの周期像の寄与も足す。測光のLOSを10/15 cMpcで事前切断しない。
これは平坦な赤方偏移priorでの模擬membershipであり、局所密度priorや実測 $p(z)$ で較正した値ではない。
宇宙膨張とpeculiar velocityを含むため、3D球と円柱は同じ体積・母集団ではない。
snapshotを周期的に使う実験であり、lightcone上の進化やbox外の独立構造は含まない。

保存する量は $\ln(1+N_{eff})$、$\log_{10}(1+\sum_i w_iM_{\star,obs,i})$、
最近傍距離、anchorと最も重い近傍のmass gap、重み付き2D位置共分散の楕円率。
形状は有効近傍数が3未満なら0とし、`has_shape`を併記する。近傍なしも専用flagを持つ。
proxy SFRの平均には`has_sfr_neighbors`を付ける。真の／推定されたhalo massや子孫情報を
これらの観測特徴量へ混ぜない。観測数は逆完全性補正せず、その観測条件のまま学習する。

$z\simeq8$で10 cMpcは半径約3.8分角。既定の15 cMpcまでの円形視野は
広域mosaicを前提とする理想化した一様視野で、実際のJWST footprintやmaskはまだ含まない。
""")
implementation('def environment_features(', 'def spatial_split(')

md(r"""
## 評価設計：空間と最終hostを分離する

**目的。** 同じ最終hostや、重なった周辺領域がtrainとtestに入り、予測性能を高く見せることを避ける。

**手法。** LOSと直交する軸に、周期箱の20%幅のtest slabを一つ固定する。
境界の両側に $2R_{max}=30$ cMpcのembargoを置き、testと同じ $z=0$ hostへ行くtraining銀河も除く。
直交軸で分けるため、測光環境のLOSが長くてもtrain/testの円柱は重ならない。
残ったtraining slabを3分割し、CVでも同じ空間embargoと最終host purgeを行う。
複数時刻の場合は全時刻を通してhostをpurgeする。

既定の最小学習数は40、validationは8、学習側の最終host数は10。
これらは計算を成立させるための下限で、統計的精度を保証する数ではない。
候補leaf sizeは10/20/50とし、最も少ないCV学習標本でも2つのleafを作れる候補だけを
**全fold・全モデルに共通の規則で**使用する。例えば最小foldが53天体ならleaf=10/20を使い、
leaf=50は除外する。fold自体を抜いたり、空間embargoやhost purgeを緩めたりしない。

候補bandwidthや木のleaf sizeはtraining CVだけで選ぶ。CVが全foldで成立しない時刻は、
全モデルで共通に除外し、理由と天体数を保存する。疎な時刻へ無検証でパラメータを移さない。
testは一度固定して、全モデルを同じ天体上で採点する。
このembargoでも大スケールの相関・宇宙分散がなくなるわけではない。
""")
implementation('def spatial_split(', 'class GaussianMixture:')

md(r"""
## 条件付き確率分布の推定

**目的。** 中央値だけでなく、幅・非対称性・clusterになる確率を比較する。
高次元の環境特徴量を扱うため、10のproduct KDEからleaf分布を使うdensity estimatorへ変更する。

`ExtraTreesRegressor`をbootstrapなしで学習し、queryと同じleafに入ったtraining銀河へ
leaf内等重みを与え、全treeで平均する。応答は $y=\log_{10}(M_0/M_\odot)$。

$$\hat p(y\mid x)=\sum_j w_j(x)\,\mathcal N(y;y_j,h_y^2),\quad
w_j(x)=\frac1T\sum_t\frac{\mathbf1[l_t(x)=l_t(x_j)]}{n_{t,l_t(x)}}.$$

Gaussian混合は全実数上で正規化され、密度の単位はdex$^{-1}$。
標本数による候補判定を先に行い、その後leaf sizeと出力幅 $h_y$ をCVのNLLで選ぶ。
標本数より大きいleafを指定して、すべてのqueryが同じ分布になる候補を比較に混ぜない。
平滑化幅は観測誤差やモデル不確実性そのものではない。
これはシミュレーションの母集団に基づく経験的な条件付き分布であり、木のパラメータのBayesian posteriorではない。

**二段階法。** 同じ観測特徴量から $p(m_z\mid D)$ を学び、別の $p(y\mid m_z)$ と積分する。
$m_z=\log_{10}M_z$ のGaussian成分を101個のbinに解析的に積分し、各bin中心で第2段階を評価する。
点推定 $\hat M_z$ の代入にはしない。第1・第2段階の幅も、最終 $y$ のCV NLLで同時に選ぶ。
grid外の確率の切り落としを検査し、trainingの $M_z$ 範囲外へ積分された確率も保存する。
有限grid近似の感度は `bridge_grid_size` を変えたtraining側の再実験で確認できる。
""")
implementation('class GaussianMixture:', 'def tune_models(')

md(r"""
## 解析3・4：共通標本でのモデル比較

| モデル | 入力 | 意味 |
|---|---|---|
| `halo_oracle` | 真の $M_z$ | 観測では得られない基準 |
| `halo_plus_environment` | 真の $M_z$＋分光環境の銀河数 | Intermediate：halo massに対する環境の追加価値 |
| `core` | 推定恒星質量 | 現在の10に対応する出発点 |
| `photo_environment` | core＋測光環境の多半径銀河数 | 周辺の測光情報の価値 |
| `spec_narrow / spec_wide` | core＋分光環境、最大5 / 15 cMpc | 広い視野の情報の価値 |
| `deep_spec` | core＋深い選択での分光環境 | 深さと分光抽出率を増やす価値 |
| `geometry` | spec_wide＋質量和、近傍、形状 | 個数以外の空間・質量分布 |
| `population_proxy` | geometry＋SFR proxy | 恒星集団情報の探索的な追加 |
| `truth3d_reference` | 真の恒星質量＋3D環境 | 理想的な参照。性能の数学的上限とは限らない |
| `two_stage_population` | population_proxyと全く同じ観測入力 | 直接法と同条件の二段階control |

`halo_oracle`と直接法の違いを、同じ観測情報量での優劣とは解釈しない。
二段階controlは同じ標本で両factorを学び、WittenのSMHM relationを数値的に再現するものではない。
学習母集団・ラベルが同じcontrolを置くことで、halo massへ圧縮する近似を検証する。

主指標は未使用testのNLL、$M_0\geq10^{14}M_\odot$のBrier score、68/95%包含率。
区間幅だけが小さくなったことを改善としない。featureのtraining範囲外、または混合のhost ESS<5を
support不足として表示するが、全体の採点から除かない。これは簡易診断で、joint supportを保証しない。

平均は各時刻内で銀河を等重みとする。同じhostへ行く複数銀河は残し、同じ時刻を繰り返すと
独立標本数が増えるとは主張しない。差の区間はtestの最終hostを単位にしたpaired bootstrap。
host間の空間相関・宇宙分散・シミュレーション系統誤差はこの区間に含まれない。
""")
implementation('def tune_models(', 'def region_audit(')

md(r"""
## 解析5：一つの領域は、一つのdescendantになるか

**目的。** anchorの将来hostを、領域全体の将来と同一視できるかを調べる。
これは観測可能な特徴量ではなく、シミュレーションの正解側の監査である。

**手法。** labelを見る前に無作為抽出した30 anchorについて、真の位置の半径10 cMpc球を固定する。
親集団は保存下限以上の全銀河で、anchor自身も含む。個数重みと真の恒星質量重みを別々に集計し、
それぞれの $z=0$ hostに入る割合を求める。最も大きい割合を受け取るhostをdominantとする。
「最も大きいdescendant mass」を選ぶラベルは作らない。

未追跡銀河を分母から落とさず、dominant割合のlower/upper boundを示す。
未追跡の全重みが他のhostへ入っても順位が変わらない場合だけ、dominantの同定がstableと記録する。
この監査の真の3D領域は観測で定義した候補そのものではない。主解析のanchorラベルを
dominantラベルへ切り替えて比較するには、別に共通標本とsplitを再定義する必要がある。
""")
implementation('def region_audit(', 'class MemoryStore:')

md("""
## 実装の確認用データと、共通標本の組立

次の人工データ関数は `MODE='synthetic_demo'` でのみ呼ばれる。
TNG modeでは人工データや人工の最終hostを使用しない。
""")
implementation('class MemoryStore:')

md("""
## 実行1：入力の読み込みと自己検査

最初に周期境界・membership確率・混合分布の正規化を確認する。
TNG modeでは、読んだheaderとDM粒子数からTNG300-1であることも検査する。
入力の読み込みに失敗した場合は停止し、demoへ自動的に切り替えない。
""")
code("""
print(analytic_self_checks())
CACHE_ROOT.mkdir(parents=True, exist_ok=True)
if DEMO:
    raw_catalogs, metadata, store = synthetic_inputs(CFG, n_groups=180, per_group=20)
    make_store = lambda: MemoryStore(store.labels)
else:
    BASE = discover_base(BASE_PATH)
    print('Reading TNG300-1:', BASE)
    raw_catalogs, metadata = {}, {}
    for sn in CFG.snapshots:
        raw_catalogs[sn], metadata[sn] = load_tng_snapshot(BASE, sn, CFG, CACHE_ROOT/'inputs')
    make_store = lambda: DescendantStore(BASE, CACHE_ROOT/'inputs', legacy_csv=LEGACY_CSV)
    store = make_store()

manifest = experiment_manifest(CFG, metadata, store, MODE, IMPLEMENTATION_SHA256)
manifest['packages'] = dict(python=platform.python_version(), numpy=np.__version__,
                           pandas=pd.__version__, scipy=scipy.__version__, sklearn=sklearn.__version__)
RUN_DIR = CACHE_ROOT/'runs'/MODE/signature(manifest)[:16]
RUN_DIR.mkdir(parents=True, exist_ok=True)
(RUN_DIR/'manifest.json').write_text(json.dumps(manifest, indent=2, ensure_ascii=False))
display(pd.DataFrame([dict(snapshot=sn, redshift=metadata[sn]['redshift'], N_tracers=len(raw_catalogs[sn]),
                           box_cMpc=metadata[sn]['boxsize_cMpc']) for sn in CFG.snapshots]))
print(LABEL_PREFIX+'run directory:', RUN_DIR)

def save_figure(fig, name):
    fig.suptitle(LABEL_PREFIX+name.replace('_', ' '), fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, .94))
    fig.savefig(RUN_DIR/(name+'.png'), dpi=160, bbox_inches='tight')
    plt.show()
""")

md(r"""
## 実行2：共通anchor、追跡完全性、多半径環境の構築

各観測条件で同じanchorを用い、feature計算後に有効な共通ラベルを持つ行で推論標本を作る。
以下の表で、選択数・有効ラベル数・星粒子数を確認する。
推論は追跡成功と有限な現在host質量に条件づけられ、未追跡の母集団への補正はしない。
`Nstar`は保存されたPartType4要素数で、windを除いたformed-star数とは限らず、収束の証明には使わない。
同じrunではfeature CSVを再利用する。別の誤差seedや深さは別runになる。
""")
code("""
try:
    if (RUN_DIR/'features.csv').exists() and (RUN_DIR/'tracking_quality.csv').exists():
        data = pd.read_csv(RUN_DIR/'features.csv', float_precision='round_trip')
        tracking_quality = pd.read_csv(RUN_DIR/'tracking_quality.csv')
        observed_by_snap = {sn: mock_observe(raw_catalogs[sn], metadata[sn], CFG) for sn in CFG.snapshots}
        parent_anchors = {sn: select_anchors(observed_by_snap[sn], CFG) for sn in CFG.snapshots}
        print('Reused features from this fingerprinted run')
    else:
        data, tracking_quality, observed_by_snap, parent_anchors = build_experiment(
            raw_catalogs, metadata, CFG, store)
        atomic_csv(data, RUN_DIR/'features.csv')
        atomic_csv(tracking_quality, RUN_DIR/'tracking_quality.csv')
finally:
    store.close()
if data.duplicated(['snapshot', 'SubhaloID']).any():
    raise ValueError('Invalid feature-cache keys')
display(tracking_quality)
""")

md("""
## 実行3：分割とCVの成立条件の確認

spatial/host purgeのため、多くのanchorが学習・testの両方から外れることがある。
その数を隠さず保存する。`cv_fold_audit.csv`には各条件の必要数と不成立の理由も記録する。
CV不成立なら停止し、testを使って分割やパラメータを調整しない。
対象時刻・計算上限を変える再設計は、新しいrunとして扱う。
""")
code("""
boxes = [metadata[sn]['boxsize_cMpc'] for sn in CFG.snapshots]
if not np.allclose(boxes, boxes[0]):
    raise ValueError('All epochs must use the same periodic box')
partitioned, split_spec = spatial_split(data, metadata[CFG.snapshots[0]], CFG)
train = partitioned.loc[partitioned.partition.eq('train')].copy()
test = partitioned.loc[partitioned.partition.eq('test')].copy()
folds, cv_audit, EPOCHS = cv_splits(train, split_spec, CFG)
atomic_csv(partitioned[['snapshot', 'SubhaloID', 'GroupID_z0', 'split_coordinate', 'partition', 'cv_fold']],
           RUN_DIR/'split.csv')
atomic_csv(cv_audit, RUN_DIR/'cv_fold_audit.csv')
(RUN_DIR/'split_spec.json').write_text(json.dumps(split_spec, indent=2))
display(partitioned.groupby(['snapshot', 'partition']).size().rename('N').reset_index())
display(cv_audit)
if not EPOCHS:
    details = cv_audit.loc[~cv_audit.eligible, ['snapshot', 'fold', 'reason']].to_string(index=False)
    raise ValueError('No epoch supports every CV fold. Inspect cv_fold_audit.csv; '
                     'the spatial/host separation has been retained.\\n'+details)
if set(EPOCHS) != set(CFG.snapshots):
    print('Excluded from ALL model comparisons:', sorted(set(CFG.snapshots)-set(EPOCHS)))
print('Epochs used by ALL models:', EPOCHS)
""")

md(r"""
## 実行4：trainingで環境の回収を可視化する

各半径について、真の3D球内銀河数と、測光・分光・深い分光環境での有効銀河数の中央値を比較する。
3Dのtracerは真の質量閾値、観測tracerは推定質量と検出で選ぶため、比は単なる測定biasではない。
異なる奥行き・閾値・分光率による違いを含む。ここではtrainingだけを描き、半径は事前固定のまま。
""")
code("""
recovery_rows = []
for sn in EPOCHS:
    frame = train.loc[train.snapshot.eq(sn)]
    for scenario in ('truth3d', 'photo', 'spec', 'deep_spec'):
        for radius in CFG.radii_cMpc:
            values = np.expm1(frame[scenario+'__lnN_R'+rtag(radius)])
            q16, median, q84 = np.quantile(values, [.16, .5, .84])
            recovery_rows.append(dict(snapshot=sn, scenario=scenario, radius_cMpc=radius,
                                      q16=q16, median=median, q84=q84))
recovery = pd.DataFrame(recovery_rows)
atomic_csv(recovery, RUN_DIR/'environment_recovery_training.csv')
fig, axes = plt.subplots(1, len(EPOCHS), figsize=(6*len(EPOCHS), 4), squeeze=False)
for ax, sn in zip(axes[0], EPOCHS):
    for scenario, g in recovery.loc[recovery.snapshot.eq(sn)].groupby('scenario', sort=False):
        ax.plot(g.radius_cMpc, g['median'], marker='o', label=scenario)
    ax.set(xlabel='Radius [cMpc]', ylabel='Median neighbor count / effective count', title='snapshot '+str(sn))
    ax.legend()
save_figure(fig, 'environment_recovery_training')
""")

md("""
## 実行5：training CVでモデルを固定する

全モデルで同じvalidation銀河を使用する。計算量を抑えるcapもモデル間で共有する。
`cv_candidate_audit.csv`で、最小foldの学習数と候補ごとの必要数・採用可否を確認する。
候補の採用可否にはvalidation/testの質量や予測損失を用いない。
CVの1モデル・1候補・1fold・1時刻ごとにCSVへ保存し、通常の中断後には完了分を再利用する。
`selected_models.json`はtestの評価前に確定する。
""")
code("""
candidate_audit = density_candidate_audit(folds, EPOCHS, CFG)
atomic_csv(candidate_audit, RUN_DIR/'cv_candidate_audit.csv')
display(candidate_audit.groupby(['leaf', 'minimum_cv_train', 'required_train', 'eligible']).size().rename('N_candidates').reset_index())
selected, cv_scores = tune_models(train, folds, EPOCHS, CFG, checkpoint=RUN_DIR/'cv_scores.csv')
(RUN_DIR/'selected_models.json').write_text(json.dumps(selected, indent=2, ensure_ascii=False))
display(pd.DataFrame([dict(model=k, cv_nll=v['cv_nll'], n_features=len(v['features']), **v['params'])
                      for k, v in selected.items()]))
""")

md(r"""
## 実行6：未使用testで分布を評価する

NLLは $-\log p(y_{true}\mid D)$、Brierは $(P[M_0\geq10^{14}]-\mathbf1[M_0\geq10^{14}])^2$。
PITから中心68/95%包含率を求める。分布の狭さ、包含率、support率を一緒に読む。
サポート不足や外挿queryも採点に含め、良い点だけを選ばない。
""")
code("""
predictions, validation, fitted_models = evaluate_models(train, test, selected, EPOCHS, CFG)
atomic_csv(predictions, RUN_DIR/'test_predictions.csv')
atomic_csv(validation, RUN_DIR/'test_validation.csv')
differences = paired_differences(predictions, CFG)
atomic_csv(differences, RUN_DIR/'paired_differences.csv')
calibration = calibration_table(predictions)
atomic_csv(calibration, RUN_DIR/'calibration_bins.csv')
display(validation)
display(differences.loc[differences.metric.eq('nll')])
""")

md("""
## 実行7：改善量、確率calibration、予測区間を読む

paired差は「新しい情報を加えたモデル − 比較元」。NLL/Brierは負なら改善する。
区間幅の差が負でも、包含率が悪化していれば改善とはしない。
同じtestを何度も見てモデルを選び直した場合、このtestは未使用評価ではなくなる。
""")
code("""
delta = differences.loc[differences.metric.eq('nll')].reset_index(drop=True)
fig, ax = plt.subplots(figsize=(9, 5))
ypos = np.arange(len(delta))
for i, row in delta.iterrows():
    ax.plot([row.lo95, row.hi95], [i, i], color='tab:blue')
    ax.plot(row.delta, i, 'o', color='tab:blue')
ax.axvline(0, color='black', ls='--', lw=1)
ax.set_yticks(ypos)
ax.set_yticklabels(delta.comparison)
ax.set_xlabel('Delta NLL: new minus reference (host bootstrap interval)')
save_figure(fig, 'paired_information_gain')

chosen_models = ['core', 'spec_wide', 'population_proxy', 'two_stage_population']
fig, axes = plt.subplots(1, 3, figsize=(15, 4))
for name in chosen_models:
    p = predictions.loc[predictions.model.eq(name)]
    axes[0].hist(p.pit, bins=np.linspace(0, 1, 11), density=True, histtype='step', label=name)
    # Bins are shown per epoch; do not pool incompatible epoch calibrations.
    for sn in EPOCHS:
        g = calibration.loc[calibration.model.eq(name) & calibration.snapshot.eq(sn)]
        axes[1].plot(g.mean_probability, g.cluster_fraction, 'o-', label=name+' sn='+str(sn))
    g = validation.loc[validation.model.eq(name)]
    axes[2].scatter(g.W68_dex, g.coverage68, label=name)
axes[0].axhline(1, color='gray', ls='--')
axes[0].set(xlabel='PIT', ylabel='Density')
axes[1].plot([0, 1], [0, 1], 'k--')
axes[1].set(xlabel='Predicted P(cluster)', ylabel='Observed cluster fraction', xlim=(0, 1), ylim=(0, 1))
axes[2].axhline(.68, color='gray', ls='--')
axes[2].set(xlabel='Mean 68% interval width [dex]', ylabel='68% coverage', ylim=(0, 1))
axes[0].legend(fontsize=7)
axes[1].legend(fontsize=6)
axes[2].legend(fontsize=7)
save_figure(fig, 'probability_calibration')
""")

md(r"""
## 実行8：真の $M_z$ に条件づけた比較の診断

halo mass baselineへの環境の追加価値は、`halo_plus_environment − halo_oracle`の未使用NLLで判定する。
以下は真の $M_z$ の0.5 dex binに分けた診断。各binの天体数も示す。
このbin分割を使ってtestに最適な閾値を選ぶことはしない。
改善は、この母集団・観測条件・推定器の範囲での予測的な追加情報を示す。
有限標本の結果だけから、数学的な条件付き独立性の一般的な証明・反証とはしない。
""")
code("""
halo_test = predictions.loc[predictions.model.isin(['halo_oracle', 'halo_plus_environment'])].merge(
    test[['snapshot', 'SubhaloID', 'logMz']], on=['snapshot', 'SubhaloID'], validate='many_to_one')
halo_test['logMz_bin'] = np.floor(halo_test.logMz*2)/2
halo_bins = halo_test.groupby(['snapshot', 'model', 'logMz_bin']).agg(
    N=('nll', 'size'), N_hosts=('GroupID_z0', 'nunique'), NLL=('nll', 'mean'), Brier=('brier', 'mean')).reset_index()
atomic_csv(halo_bins, RUN_DIR/'halo_conditional_bins.csv')
display(halo_bins)
""")

md(r"""
## 実行9：同じtest銀河の直接法・二段階法のPDF

最初の採点対象銀河を機械的に選び、真の結果に合わせて見栄えの良い例を探さない。
二段階法は中間halo massの分布を積分している。線は有限シミュレーションの条件付き分布で、
観測装置・銀河形成モデル・宇宙論の系統誤差まで含む事後分布ではない。
""")
code("""
example = predictions.loc[predictions.model.eq('population_proxy')].iloc[0]
query = test.loc[test.snapshot.eq(example.snapshot) & test.SubhaloID.eq(example.SubhaloID)].iloc[[0]]
grid = np.linspace(min(train.logM0.min(), query.logM0.iloc[0])-1,
                   max(train.logM0.max(), query.logM0.iloc[0])+1, 400)
fig, ax = plt.subplots(figsize=(8, 4))
example_rows = []
for name in ('core', 'spec_wide', 'population_proxy', 'two_stage_population'):
    model = fitted_models[(int(example.snapshot), name)]
    post, outside = model.predict(query)
    p = post[0]
    ax.plot(grid, np.exp(p.logpdf(grid)), label=name)
    example_rows.append(dict(model=name, outside_range=bool(outside[0]),
        median=p.quantile(.5), q16=p.quantile(.16), q84=p.quantile(.84),
        P_cluster=1-p.cdf(14), **p.diagnostics()))
ax.axvline(query.logM0.iloc[0], color='black', ls='--', label='TNG label' if not DEMO else 'Synthetic label')
ax.set(xlabel='log10(M0 / Msun)', ylabel='Conditional density [dex^-1]')
ax.legend(fontsize=8)
save_figure(fig, 'example_conditional_densities')
example_table = pd.DataFrame(example_rows)
atomic_csv(example_table, RUN_DIR/'example_posterior_summary.csv')
display(example_table)
""")

md("""
## 実行10：領域の行先とanchorの代表性

子孫追跡を領域内のtracerにも拡張するため、追加I/Oが必要になる。
最初は30 anchorで監査する。`region_audit_n=0`にすれば明示的に省略できる。
追跡済みの天体は同じSQLiteから再利用し、主モデルの入力や採点にはこの表を用いない。
""")
code("""
region_tables = []
if CFG.region_audit_n > 0:
    audit_store = make_store()
    try:
        for sn in EPOCHS:
            region_tables.append(region_audit(observed_by_snap[sn], parent_anchors[sn],
                                              metadata[sn], CFG, audit_store))
    finally:
        audit_store.close()
if region_tables:
    region_fates = pd.concat(region_tables, ignore_index=True)
    atomic_csv(region_fates, RUN_DIR/'region_descendant_audit.csv')
    display(region_fates)
    fig, ax = plt.subplots(figsize=(7, 4))
    for weighting, g in region_fates.groupby('weighting'):
        ax.scatter(g.f_main_lower, g.f_anchor_lower, label=weighting, alpha=.6)
    ax.plot([0, 1], [0, 1], 'k--')
    ax.set(xlabel='Dominant host fraction: lower bound', ylabel='Anchor host fraction: lower bound',
           xlim=(0, 1), ylim=(0, 1))
    ax.legend()
    save_figure(fig, 'region_descendant_allocation')
else:
    print('Region audit explicitly disabled')
""")

md(r"""
## 結果の保存と、次の観測モデルへ進める際の判断

`test_validation.csv`と`paired_differences.csv`で、広い視野、分光、深い選択、形状、SFR proxyの追加を比べる。
優先順位は、NLL/Brierが改善し、包含率・supportを維持するかで判断する。
同じtestで最善の追観測設定を選んだ後の性能を確定するには、別の未使用領域が必要になる。

この段階で実装していない、実測のfilter flux・line EW選択、dust、金属量、ALMA、Lyαは、
実データまたは妥当性を確認したsynthetic observable catalogとselection functionを用意してから追加する。
真の物理量をそのまま観測特徴量と呼んで置き換えない。
観測時間あたりの情報量の比較には、装置の感度・pointing・target配置のモデルも別途必要である。

保存される主なファイル：`manifest.json`、`tracking_quality.csv`、`features.csv`、`split.csv`、
`cv_fold_audit.csv`、`cv_candidate_audit.csv`、`cv_scores.csv`、`selected_models.json`、`test_predictions.csv`、
`test_validation.csv`、`paired_differences.csv`、`calibration_bins.csv`、各PNG。
region監査を有効にすると`region_descendant_audit.csv`も保存する。
モデルそのもののpickleは保存しない。固定した設定とcatalogから再学習する。

誤差・深さ・seed・半径・モデルgridを変えた場合は、新しいrunになる。
古いrunと比較する際は、標本・選択・ラベル・splitが一致するかmanifestを先に確認する。
""")
code("""
completion = dict(mode=MODE, scientific_TNG_run=not DEMO, completed_epochs=[int(s) for s in EPOCHS],
                  n_test_anchors=int(predictions.loc[predictions.model.eq('core')].shape[0]),
                  result_directory=str(RUN_DIR), implementation_sha256=IMPLEMENTATION_SHA256)
(RUN_DIR/'completion.json').write_text(json.dumps(completion, indent=2))
print(LABEL_PREFIX+'completed')
print('Results:', RUN_DIR)
if DEMO:
    print('Artificial software-verification results only. Run MODE=tng for scientific analysis.')
""")

notebook = dict(cells=cells, metadata=dict(
    kernelspec=dict(display_name='Python 3', language='python', name='python3'),
    language_info=dict(name='python', file_extension='.py', mimetype='text/x-python',
                       pygments_lexer='ipython3'),
    implementation_sha256=SOURCE_SHA,
    intended_runtime='IllustrisTNG JupyterLab; standalone notebook'),
    nbformat=4, nbformat_minor=4)
# Match Jupyter's line-array serialization so regenerating code does not
# replace every unchanged source cell after a user has saved an execution.
for cell in cells:
    cell['source'] = cell['source'].splitlines(keepends=True)
destination = ROOT/'notebooks'/'11_tng300_observational_environment_descendant_prediction.ipynb'
destination.write_text(json.dumps(notebook, ensure_ascii=False, indent=1)+'\n')
print('Wrote', destination, 'cells=', len(cells))
