# 六角氷柱の配向合成と CDF 出力 — M2

## 適用対象と到達点

`a9538941dcb9df2de6a4130fa66f625f0133c949` に **hex-trace-m1 差分を適用済み**の
`ice-crystal-phase` が対象です。M1 の追跡・偏光・交差処理と既存の雨滴ソルバは変更しません。
追加は `hex_phase/` とルート CMake の追加ブロックです。GitHub への push は行っていません。

この段階では、固定の入射進行方向・波長・形状に対して、任意の重み付き配向列を合成し、
既存 `PhaseRecord` で読み込める4ファイルを生成します。

```
metadata.json
phi_cdf.npy                 [Nphi+1]
theta_given_phi_cdf.npy     [Nphi, Ntheta+1]
u_edges.npy                [Ntheta+1]
```

FP64 / little endian / C order、phi 周辺 → theta 条件付き、
`u=(1-cos(theta))/2`、セル内立体角密度一定という契約を維持します。
**最終レコードに点群 CSV は追加しません。配向 CSV は再現用に別途保管してください。**
メタデータには入力 CSV の SHA256 を保存します。巨大な配向列を metadata に埋め込まず、
既存 reader の 1 MiB 上限も維持します。

### 物理・数値モデル

- M1 と同じ非吸収・等方媒質の幾何光学。二入力の光束 Jones 応答を伝播し、無偏光入射の強度を合成します。
- 内部全反射の s/p 相対位相は保持しますが、異なる経路・粒子姿勢の複素振幅は足しません。
- `--ior` は **指定波長で呼び出し側が評価した実屈折率**です。自動の氷の分散式ではありません。
  波長ラベルだけを変えて同じ `--ior` を使っても分散は生じません。
- 吸収・複屈折・有限粒径回折・経路間干渉・表面粗さを追加していません。
- 雨滴の Table II フィルタも、追加の球面 Gaussian も適用しません。
- 対称性による入射方向や姿勢の削減はまだ行いません。全ての指定姿勢をそのまま処理する比較基準です。
- 配向 CSV は有限の離散測度です。観測分布の選定や連続姿勢積分の収束認定は、この出力成功だけでは行いません。

## ビルド・テスト

### CPU 参照経路（CUDA/OptiX 不要）

リポジトリのルートから実行します。M1 と既存の NPY writer / reader が存在する必要があります。

```powershell
cmake -S hex_phase -B build-phase-host -DICE_PHASE_ENABLE_CUDA=OFF -DICE_PHASE_BUILD_TESTS=ON
cmake --build build-phase-host --config Release
ctest --test-dir build-phase-host -C Release --output-on-failure
```

CLI 相互運用テストは NumPy がある場合だけ登録します。`ctest -N` で
`ice_phase_cli_numpy` が含まれることを確認してください。NumPy がない場合、C++ 側テストを
通しただけで Python reader の相互運用を検証したことにはなりません。

### CUDA 追跡・GPU 集計・既存 GPU CDF

**CUDA M2 はリポジトリルートで構成します。** M1 の standalone CUDA ビルドと異なり、
既存 `rainbow_phase_dataset` をリンクするためです。既存の CUDA / OptiX / GPU アーキテクチャ設定を維持してください。

```powershell
cmake -S . -B build -DBUILD_TESTING=ON -DICE_BUILD_HEX_TRACE=ON -DICE_HEX_ENABLE_CUDA=ON -DICE_BUILD_HEX_PHASE=ON -DICE_PHASE_ENABLE_CUDA=ON -DICE_PHASE_BUILD_TESTS=ON
cmake --build build --config Release --target ice_generate_phase ice_generate_phase_cpu ice_phase_host_tests ice_phase_cuda_tests
ctest --test-dir build -C Release -R "^ice_phase_" --output-on-failure
```

`ice_generate_phase.exe` は CUDA 専用、`ice_generate_phase_cpu.exe` は明示的な CPU 参照経路です。
CUDA の失敗による CPU への自動切り替えはありません。
CUDA binary の隣に `modules/hex_trace.fatbin`, `phase_accumulate.fatbin`, `phase_cdf.fatbin` を配置する
CMake target を追加しています。`--modules DIR` で配置先を明示することもできます。

**提供環境では CPU の GCC / Clang / ASan+UBSan を検証しました。CUDA のコンパイル・実行と
Windows/MSVC は未検証です。ユーザー側で通った M1 のテストは、新規 M2 CUDA コードの検証ではありません。**

## 実行例

まず、付属の2配向・異なるレイ数の例を使用します。Visual Studio generator の配置を想定しています。

```powershell
.\build\hex_phase\Release\ice_generate_phase.exe --out out\phase_m2 --orientations hex_phase\configs\two_orientations.csv --radius-mm 0.1 --length-mm 0.2 --ki 1 -0.4 0.3 --ior 1.31 --wavelength-nm 550 --theta 72 --phi 144 --cdf-theta 36 --cdf-phi 72
```

CPU standalone の場合は実行ファイルを次に変更し、それ以外の引数は同じにします。

```text
.\build-phase-host\Release\ice_generate_phase_cpu.exe
```

この例の `1.31` は動作検証用の明示的な値であり、氷の 550 nm 光学定数を自動取得した結果ではありません。
`--radius-mm` は六角形の **外接円半径**、`--length-mm` は柱の **全長**です。
粒子は原点中心で、body 座標系の柱軸は +y です。`--ki` は ensemble 基準座標系での光の進行方向です。
正規化した実際の方向を metadata に保存します。

`--theta`, `--phi` は集計する細格子です。`--cdf-theta`, `--cdf-phi` を省略すると同じ解像度で保存します。
縮約を行う場合、保存解像度は細格子解像度の整数約数に限ります。
質量保存は形状保存の保証ではないため、`storage_processing.coarsening_tv` を確認してください。
既定の `--max-coarsening-tv 1` は構造的な上限であって、研究上の許容誤差ではありません。

### 揺らぎを持つテスト用配向列

例えば柱軸が水平面の上下3度以内に分布し、軸まわりの回転を一様にする列は次のように作成できます。

```powershell
python hex_phase\tools\make_orientations.py --out configs\horizontal_band_test.csv --mode axis-band --half-angle-deg 3 --count 128 --samples-per-orientation 128 --seed 7
```

生成器の組込み分布は、検証・モデル例に限定した次の4種類です。

| mode | 分布 |
|---|---|
| `fixed` | identity 一姿勢。`--count 1` を要求 |
| `isotropic` | SO(3) の一様分布から独立サンプル |
| `axis-cone` | body +y 軸が基準 +y から指定半角の円錐内に立体角一様、spin 一様 |
| `axis-band` | body +y 軸が水平面の上下指定半角の帯内に立体角一様、spin 一様 |

cone / band は **気象観測にフィットした Gaussian の代替ではありません**。
文献に基づく分布を使用する場合は、その定義に基づく姿勢・確率質量を同じ CSV 契約で与えます。
連続分布の姿勢数による収束と、各姿勢内の入射位置サンプル数による収束は別々に検証してください。

## 配向 CSV の厳密な契約

UTF-8（BOM なし）、ヘッダは次の7列です。CRLF と行先頭 `#` のコメントを扱います。

```csv
qw,qx,qy,qz,weight,samples,seed
1,0,0,0,0.5,1024,12345
0.7071067811865475244,0,0,0.7071067811865475244,0.5,2048,54321
```

`qw,qx,qy,qz` は scalar-first Hamilton quaternion による **body → ensemble** の能動回転です。
単位長を要求し、二乗ノルムと1との差が2e-12以内の丸め誤差だけを許します。任意長の4ベクトルは受け付けません。
行の順序・`samples`・`seed` が有限サンプル計算を定義します。

`weight` は結晶 **個数分布の確率質量**です。正値だけを受け付け、総和1を2e-12以内で要求します。
既に正規化された値の丸め差だけを総和で割って除き、その係数を保存します。
ゼロ質量の姿勢は行を省略します。

- 姿勢 PDF の点評価値そのものを `weight` に入れません。
- 求積の場合は、SO(3) の測度・必要な Jacobian・求積重みまで含めた確率質量を与えます。
- 配向分布自体から独立サンプルする場合は `weight=1/M` です。さらに同じ配向 PDF を掛けません。
- `samples` が各姿勢で違ってもよく、多い姿勢を余分に重くすることはありません。
- 結晶への衝突確率に面積を反映済みの姿勢分布を、この「個数分布」入力にそのまま渡しません。

各姿勢 R_j の入射方向は `R_j^T ki`、出射方向は `R_j ko_body` です。
その出射サンプルの球面セルへの寄与は

```
weight_j * projected_area(R_j^T ki) / samples_j * power_fraction
```

です。初回入射点は M1 と同じく可視面の投影面積に比例して直接生成するため、
ヒット率をさらに掛けることはありません。

## 数値処理の順序

```
配向CSV
  -> 各姿勢の body 入射方向と投影面積
  -> M1 の追跡（CPU または CUDA）
  -> レイ別監査を検査
  -> 出射方向を ensemble 基準へ戻して、未正規化光束を加算
  -> 極の成分を表現規約に従って割り当て
  -> 密度 D = セル光束 W / セル立体角
  -> CDF の構築・監査・保存
```

GPU では出射点群を CPU にダウンロードしません。出射点の生成は M1 の GPU CSR バッファ、
質量集計は GPU FP64 `atomicAdd`、密度バッファは既存の
`rainbow::PhaseCdf::build(PhaseDensityView, ...)` の `Scalar` 入口へ渡します。
`Scalar` は後段へ渡す密度の形式であり、追跡中の偏光を無視する指定ではありません。
追跡から CDF 構築までの途中で CPU に回収するのは、レイ監査情報と小さい集計統計です。
最終 CDF 配列は既存 writer が分割して読み戻し、ファイルへ保存します。

**既存 CDF は入力密度にセル立体角を掛けます。セル光束をそのまま渡すと二重計上になります。**
この誤接続が一様球面分布を変えてしまうこともテストしています。

CPU 参照経路は、セル質量を補償加算し、保守的に縮約してから順序付き FP64 prefix sum を作ります。
こちらは既存 GPU CDF の CPU 実行ではなく、同じ出力測度・セル契約に対する独立した実装です。
GPU テストで、既存 CDF に通した結果とのセル確率・g・縮約 TV を比較します。

### 丸め誤差と極・セル境界

GPU の atomic 加算順序は固定していないため、ビット一致は保証しません。
独立なレイ光束監査から得た総光束と、ヒストグラム積分を比較します。
`atomic_absorbed_addend_mass_mm2` は加算で消えた微小項の診断値で、厳密な誤差上界ではありません。

球面セル境界のごく近傍では、閉開区間の割当てを CPU/CUDA 間で安定化する丸め規則を使います。
u 境界は64 epsの座標スケール、phi は64 epsの格子スケールを上限とし、さらに **セル幅の1e-7以下**に制限します。
変更したサンプル数を `boundary_snapped_samples` に記録します。物理的な平滑化や確率 floor ではありません。

幾何光学の前方・後方デルタは削除しません。進行軸まわりの横成分が `256*DBL_EPSILON` 以下のものを
極の丸め域として別に集計し、細格子の最初／最後の極冠内に方位角一様で配分します。
undefined な phi の任意値による人工的な一本の子午線ピークを作らないための、**有限セル表現の規約**です。
これは回折による光学的な広がりではありません。細格子を変えるとこの表現幅も変わります。

`hg.g` は保存した CDF からの厳密なセル内平均余弦
`1-U[i]-U[i+1]` で計算します。元の点群の一次モーメント `histogram.g_rays` とは区別します。
また `sampling_frame_columns` の列は `[e0,e1,ki]`、JSON の外側配列は行です。
その行列はサンプリング局所座標から **ensemble 基準座標**へ写します。
一つの結晶の body frame と取り違えず、読む側で別の ONB に置き換えないでください。

## 光束収支と失敗時の扱い

各姿勢の平均投影面積を個数重みで平均し、非吸収・幾何学的衝突モデルとして
`Csca=Cext=mean_Aproj`, `Cabs=0` を保存します。Maxwell の波動散乱断面積ではありません。
実際に回収した光束積分と未回収光束積分は `transport_audit` に別々に保存します。

CDF は全姿勢の **回収済み光束を合成した後**に一度だけ正規化します。その前に必ず次を検査します。

- M1 の全レイが許容終了状態で、光束収支と各界面の診断が有限かつ許容範囲か。
- 全ての計画姿勢・レイ・出射点が処理されたか。
- 面積重み付きの未回収比率が `--max-unresolved-fraction` 以下か。
- ヒストグラム積分が独立なレイ監査の回収光束と一致するか。
- CDF 丸め誤差・失われた確率質量・縮約 TV が指定値以内か。

既定値は、最大内部相互作用1024、レイ別残存閾値1e-12、全体未回収比率1e-10、
光束収支1e-10、ヒストグラム積分相対誤差1e-10、CDF L1誤差1e-10、CDF消失質量1e-12です。
1024回で収束する保証はありません。上限到達は失敗とし、自動でレイ数・格子を減らしません。

生成中は `<out>.part` へ書きます。成功時だけ最終ディレクトリへ改名します。
生成時の失敗は通常終了コード2と `failure.json`、引数・事前検査の失敗は終了コード1です。
ディスク自体の書き込み失敗では診断ファイルも書けないことがあります。
既存出力・既存 `.part` は上書き・削除しません。

**`complete=true` は実行と構造監査の完了です。物理モデル誤差、配向積分誤差、レイ位置積分誤差、
角度離散化、保存時の縮約、NF近似の誤差は別々に評価する必要があります。**

## 主なファイル

| ファイル | 役割 |
|---|---|
| `orientation_plan.cpp` | CSV・quaternion・個数重みの検証、SHA256 |
| `phase_math.hpp` | CPU/CUDA 共通の回転・セル対応付け |
| `phase_table.cpp` | CPU 質量集計と独立した CDF 参照実装 |
| `phase_accumulate.cu` / `phase_cuda.cpp` | GPU 集計と密度バッファ化 |
| `phase_generate.cpp` | 合成の監査・metadata |
| `generate_hex_phase.cpp` | 全体の実行と既存 CDF 接続 |
| `phase_tests.cpp` | 個数重み・面積・回転・測度・CDF・実追跡のホスト検証 |
| `phase_cuda_tests.cpp` | CPU/GPU 集計と既存 CDF の比較、極・HDR・不正入力 |
| `phase_cli_test.py` | 無変更の `PhaseRecord` による相互運用と失敗時の非公開化 |

## 続く検証・拡張

この差分で「追跡 → 配向合成 → 既存形式の参照 CDF」は接続されます。
残る主要項目は、採用する氷の光学定数と配向分布の文献準拠入力、条件ごとの一括生成、
独立バッチ・配向数・入射レイ数・角度解像度に対する収束評価です。
現在の均一円錐／帯のテスト用配向を、観測に基づく配向モデルとして扱わないでください。
