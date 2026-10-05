# Rainbow：角度分布の可視化と Mie 比較

対応する C++ 基準コミット：`b1778ca041214466edae9b0c1059932afbbd89d5`
（`Add: patch evaluation on GPU`）。対応 CSV：`rainbow_patch_optics_v1`。

**C++、CUDA、OptiX、CMake は今回変更しません。** `PatchOptics::write_csv()` が既に保存する
線形の光学量・角度格子・診断情報を直接使用します。ZIP の `vis/` をリポジトリ直下に置いてください。
`src/optix_context.cpp` に対するローカルの変更には触れません。

## 到達点と限界

現行 C++ の出力は、評価できる正則パッチの部分和です。焦線位相、回折、欠損セルの補完、
非正則セル・共有境界の光学的処理は未完了です。Exit code 0 は現在の検査を通って完了した意味で、
完成版の位相関数が検証済みという意味ではありません。

Python は欠けている物理を追加・代用・調整しません。Mie に合うように係数をフィットしたり、
各曲線を別々に最大値や積分で正規化したりもしません。Mie 比較は、この段階では期待されるモデル差も
含む診断です。入力の `optical_complete=false` / `source_coverage_certified=false` に対しては、
`--allow-partial-model` なしの比較を拒否します。この指定後も数値エラー・既知の未処理交差は評価から除外します。

## ファイル

- `optics_io.py`：CSV・メタデータ検査、格子への対応付け、単位換算、二入射状態の強度平均。
- `analysis_tools.py`：カラーマップ、状態図、断面図、立体角重み付き指標。
- `mie_reference.py`：実際の `miepython` を使う球形の参照計算。代替ソルバはありません。
- `plot_optics.py`：可視化 CLI。
- `compare_mie.py`：Mie 比較 CLI。
- `run_sphere_pair.cmd`：既存 C++ executable を使う球形 x/y 二入射状態の実行例。
- `tests/test_vis.py`：形式・単位・偏光・マスク・比較処理・描画の検査と、実 Mie の検査。

## 導入（Windows / cmd.exe）

リポジトリのルートで実行します。Python 3.10 以上を使用してください。

```bat
py -3 -m venv .venv
.venv\Scripts\python.exe -m pip install -r vis\requirements.txt
.venv\Scripts\python.exe -m unittest discover -s vis\tests -v
```

Mie の API は `miepython==3.3.0` に固定しています。使用した Python / NumPy / Matplotlib / miepython の
版は実行結果の JSON に記録されます。実験環境を厳密に再現する際は、さらに
`.venv\Scripts\python.exe -m pip freeze > outputs\python_environment.txt` を保存してください。

生成物を Git に入れない場合は、既存 `.gitignore` の末尾に次を追加してください。
既存の `.gitignore` 自体を置換する必要はありません。

```gitignore
/.venv/
/outputs/
__pycache__/
*.py[cod]
```

## 1. いま得られている非球形雨粒を表示

前回の C++ 実行で生成した **`drop_1mm_x_optics.csv`** を使用します。
`drop_1mm_x.csv`（頂点）や `*_queries.csv`（交差件数）ではありません。

```bat
.venv\Scripts\python.exe vis\plot_optics.py outputs\drop_1mm_x_optics.csv --out outputs\vis_drop --stage path --component total
```

生成物：

- `intensity.png`：対数色スケールの 2D 強度分布。
- `status.png`：既知交差の評価完了・真のゼロ・未処理・数値エラーを区別する図。
- `summary.json`：条件、単位、マスク、立体角範囲、部分積分、入力ファイルの SHA-256、実行環境。
- `angular_data.npz`：**線形値のまま**の `density_mm2_per_sr[Ntheta,Nphi]`、角度、重み、マスク。

`np.load(..., allow_pickle=False)` で読み込めます。`source_metadata_json` も文字列として含みます。
これは C++ の生 CSV を勝手に差し替えるキャッシュではなく、今回選択した stage/component の保存です。

虹の領域を拡大：

```bat
.venv\Scripts\python.exe vis\plot_optics.py outputs\drop_1mm_x_optics.csv --out outputs\vis_drop_bow --theta-range 120 150
```

非干渉の強度和も別に表示：

```bat
.venv\Scripts\python.exe vis\plot_optics.py outputs\drop_1mm_x_optics.csv --out outputs\vis_drop_incoherent --stage incoherent
```

`--component s` / `p` は **出力の散乱面基底の成分**です。`s` は垂直、`p` は出射方向と s の外積方向です。
C++ の `--polarization x` / `y` は **固定された入射基底の成分**なので、同じ名称の基底ではありません。

`--linear` は色だけを線形表示にします。`--vmin` / `--vmax` は対数表示でも**線形の物理単位値**を指定します。
複数条件を比較するときは同じ値を指定してください。既定の min/max は選択領域の有効値から決定します。

`--show-partial` は未処理交差を含む方向の有限な部分和も**表示する**ための明示オプションです。
数値エラーの NaN は表示しません。定量比較ではこのオプションを用意していません。

## 2. 球形データを新たに生成

半径 1 mm の既定形状は球ではないので、Mie 比較用のデータには必ず `--sphere` を指定します。
スクリプトは係数 8 個が**すべて厳密にゼロ**であることを検査し、非球形を Mie の参照対象にしません。

```bat
vis\run_sphere_pair.cmd
```

これは、既存の `build\Debug\rainbow_trace.exe` と module を使い、
半径 0.4 mm、真空波長 700 nm、内側屈折率 1.3314、傾斜 0 度、
入射格子 513×513、出射格子 1800×64 に対して x/y の二回を実行します。
二つのログと出力は `outputs\sphere_validation\` へ保存します。
**同じ名前の過去の結果を上書きします。比較を残す場合は OUT の値を変えてください。**

この解像度は初期比較用で、収束を保証する値ではありません。入射格子と出射角度分解能は別々に変化させて
収束を調べてください。90×180 の出射格子は theta 間隔 2 度であり、補間で高精細表示しても計算情報は増えません。
全周の大規模格子へ先に上げず、球形では theta を細かくして phi の数を抑える方法が適しています。
x/y 一方向の入力は偏光による phi 依存性があるので、単一ファイルで軸対称になると期待してはいけません。

この CMD は作成済み C++ 実行ファイルを起動するだけです。ビルド・validation 設定・キャッシュ・精度オプションは変更しません。
GPU 実機でこの CMD を実行する試験は提供側では未実施です。

## 3. Mie と比較

単一の入射 x 状態を対応する偏光の Mie と比較する例：

```bat
.venv\Scripts\python.exe vis\compare_mie.py outputs\sphere_validation\sphere_x_optics.csv --out outputs\compare_sphere_x --theta-range 120 150 --allow-partial-model
```

x/y の**強度平均**を無偏光として比較する例：

```bat
.venv\Scripts\python.exe vis\compare_mie.py outputs\sphere_validation\sphere_x_optics.csv --orthogonal outputs\sphere_validation\sphere_y_optics.csv --out outputs\compare_sphere_unpolarized --theta-range 120 150 --allow-partial-model
```

`--orthogonal` は入力 Jones ベクトルの内積を検査します。同じ x ファイルを二回渡すと失敗します。
粒径・波長・屈折率・形状・入射基底・角度格子・近似条件も照合します。異なる格子への再サンプリングは行いません。
二つの入射光量が異なる場合も、それぞれを入射 Jones ノルムで割ってから強度を等重みで平均します。
電場を状態間で合成しないので、意図せず斜め直線偏光を作ることはありません。

可視化側も同じ `--orthogonal` に対応します。

### 比較の生成物

- `model.png` / `mie.png`：同じ単位・同じ有効画素マスク・**共通の色範囲**。
- `log10_ratio.png`：log10(model/Mie)。正で有限なペアのみ。0 は一致、1 は10倍。
- `theta_slice.png`：同じ phi の theta 断面。既定 `--phi-deg 0` に最も近い**実際のサンプル**を使用。
  選んだ角度を図と JSON に記録し、phi=0 を偽って表示しません。
- `theta_azimuth_mean.png`：無偏光二入力の場合のみ。両者に共通の有効 phi 上で平均。
  phi が欠けた行もその有効数を JSON に記録します。
- `status.png`：元のモデルの状態。
- `comparison.csv` / `comparison.npz`：線形量、参照、比較可否、比。
- `metrics.json`：立体角重み付き相対 L1/L2、部分積分比、正値ペアの log10 RMSE、ゼロ件数、
  マスクの被覆立体角、理論 Mie 全散乱断面積と格子上の数値積分との差、使用ライブラリ版。

**`exit=0` は評価スクリプトが完走した意味です。物理的妥当性の PASS/FAIL は自動付与しません。**
現在の近似不足・未処理被覆・粗い格子の影響があるので、出力指標を実装誤差だけに帰属させないでください。

## 座標と保存形式

使用するのはメルカトルではなく、theta/phi が等間隔の正距円筒的な角度配置です。
画像は横が phi（-180〜180度）、縦が theta（0〜180度、0が上）。theta は入射**進行方向**との角度で、
0 が前方散乱、180 度が後方散乱です。theta は仰角や緯度ではありません。
C++ の direction_id = theta_index*Nphi + phi_index を検査して [Ntheta,Nphi] の行列へ置きます。
座標はすでに CSV で保存されているので、保存用の新しい C++ バッファは不要です。

実メルカトルなら緯度 beta=pi/2-theta に対して y=log(tan(pi/4+beta/2)) で、極が無限遠になります。
今回の全方向の角度評価には使いません。密度の色表示時に sin(theta) を掛けることもしません。
積分・指標にのみ、CSV の DeltaOmega=DeltaPhi*(cos(theta_lo)-cos(theta_hi)) を使用します。

CSV は監査しやすい一方、最終解像度では巨大になります。今回のローダーは必要列だけ読みますが、
全角度バッファを RAM に置くので、大規模格子では相応のメモリを要します。NPZ は再利用用の圧縮保存です。
C++ の直接バイナリ出力とチャンク読み込みは今後の拡張であり、今回は実装していません。

## 物理単位と偏光の対応（重要）

C++ 保存値を D、入射 Jones ベクトルを (Ex,Ey)、等価半径を a_mm とすると、比較量は

```
model_mm2_per_sr = D * a_mm**2 / (abs(Ex)**2 + abs(Ey)**2)
```

です。これはモデルを微分散乱断面積の単位にそろえる操作であって、物理的な完全性の保証ではありません。

Mie 側は高水準 API `intensities(n_sphere, diameter_nm, wavelength_nm, mu,
 n_env=n_environment, norm='qsca', n_pole=0)` を使用します。
半径 mm → 直径 nm の換算は **2*a_mm*10^6**。屈折率は絶対値を渡し、媒質補正の二重適用を避けます。
コードから読み取れる波長規約は真空波長（光路は n*l/lambda）です。

`norm='qsca'` の無偏光強度は球面積分が Qsca なので、`pi*a_mm**2` を掛けると微分散乱断面積になります。
戻り値の順序は **(parallel, perpendicular)** です。既定の albedo 正規化は使用しません。

固定入射基底 e0,e1 に対する phi ごとの散乱面基底は

```
s    = -sin(phi)*e0 + cos(phi)*e1
p_in = -cos(phi)*e0 - sin(phi)*e1
```

です。入射 Jones を単位ノルムにして、この s と p_in に投影した係数の絶対値二乗で Mie の二成分を重み付けします。
したがって単一の x/y 入射も、正しく対応する phi 依存分布と比較できます。
無偏光二入力の場合は各成分の重みが1/2です。
Mie は CSV の公称 theta セル中心で評価します（C++ のレイ方向は FP32 で格納されるため、
その丸め差まで含めたビット単位の同一方向計算ではありません）。

複素振幅の絶対位相は Mie 実装との規約差があるため比較しません。比較対象は s/p/total 強度です。
`n_pole=4` などとすると四つの光路族を選べる、という解釈は誤りです。ここでは全 multipole の Mie を使用します。
Mie は高次の内部反射や前方回折を含み、R/TT/TRT/TRRT と同じ打ち切りモデルではありません。

## 対数表示・ゼロ・未処理値

`LogNorm` は**色変換だけ**です。保存値に log をかけたり epsilon を足したりしません。
真のゼロは対数表示不能なので強度図では非表示になり、状態図で独立に確認できます。
数値エラーの NaN と未処理パッチの部分和も状態図・マスク・JSON で分けます。

既定では `known_hits_complete=0` の画素を強度図と指標から除きます。
それでも GAS に存在しないセルの寄与は既知交差マスクでは検出できません。
したがって全画像に `source_coverage_certified=false` という制限が残ります。
部分的な coherent intensity は全体の下界でもありません。未処理の波が破壊的干渉を起こす可能性があるためです。

## 指標の定義と解釈

同じ有効マスク M、セル立体角 w で

```
rel_L1 = sum_M w*abs(model-mie) / sum_M w*mie
rel_L2 = sqrt(sum_M w*(model-mie)**2 / sum_M w*mie**2)
```

を計算します。比の log10 RMSE は両値が正の集合だけで計算し、その件数とゼロ対の件数を別記録します。
`--theta-range` はセル中心による選択で、そのセル全体の立体角を使います。
全方向の数値積分も公称セル中心の点サンプルによる求積です。狭い前方ピークを解像できないと、
Mie 側だけでも理論全散乱断面積と一致しません。その誤差を報告し、モデル側を無理に規格化しません。

現行 `path` 出力は光路位相のみ。論文の焦線位相は補間後の位相に必要で、回折近似は別の後処理です。
論文の Fig.9/Fig.11 と同じ完成版比較になるのは、それらと境界処理を実装し、解像度収束を確認した後です。
現在でも分布の向き、偏光の対応、スケール、ピークのおおまかな位置、修正による変化を確認できます。

## 参照

- 対象論文：Sadeghi et al., *Physically-Based Simulation of Rainbows* (2012), §4.1.1–4.1.2, Eq.(3)–(4), Fig.9/11.
- CSV：上記基準コミットの `src/patch_optics.cpp`, `include/rainbow/query_direction_grid.hpp`, `patch_optics_data.hpp`。
- miepython normalization：<https://miepython.readthedocs.io/en/latest/03a_normalization.html>
- miepython intensities API：<https://miepython.readthedocs.io/en/stable/api/miepython.core.intensities.html>
- Matplotlib LogNorm：<https://matplotlib.org/stable/api/_as_gen/matplotlib.colors.LogNorm.html>
- PROJ equidistant cylindrical / Mercator：<https://proj.org/en/stable/operations/projections/eqc.html>、<https://proj.org/en/stable/operations/projections/merc.html>

提供側での実行範囲と未検証項目は `TEST_REPORT.md` を参照してください。
