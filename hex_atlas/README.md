> **M4 spectral/NF extension:** 波長から氷の屈折率を自動評価し，最終レコードを
> 実際の `RainbowReference` と同じ契約で保存する．一括生成・本番前確認は
> [`SPECTRAL_NF.md`](SPECTRAL_NF.md) を参照．以下の M3 の固定屈折率例は
> `--material constant` を明示した回帰試験であり，分散を含む通常生成ではない．

# M3: 固定入射応答の事前計算と決定論的な配向合成

## この差分の範囲

M1/M2/hotfix1 を変更せず、`hex_atlas/` とルート CMake の追加だけで実装する。
基準リポジトリは `Somatic-Cell/ice-crystal-phase`、確認した main は
`36477f4bed099b2c067f4cf57739dc4d1704b5f3`。

- **事前計算:** CPU または CUDA。固定六角柱の入射方向別の非正規化セル光束を保存する。
- **配向合成:** CPU/NumPy。既存応答の補間・対称変換・循環シフト・重み付き和のみ。
  レイトレーサを呼ばず、配向モデルを変更して再利用できる。
- **最終出力:** 各入射方向・波長に対する `rainbow.phase_cdf.numpy.v2` の4ファイル。
- **ランダム性:** 新しい経路は、姿勢・入射位置・光路分岐のいずれにも乱数を使わない。
  旧 M1/M2 の乱数サンプラーは削除せず、M3 からは呼ばない。
- **GPU の加算:** floating atomicAdd を使わない。レイごとに分離したセル列を固定順で集約する。
  異なる CPU/GPU/コンパイラ間でのビット一致までは主張しない。

これは有限格子・数値求積モデルの実装であり、連続配向積分の閉形式解ではない。
特に **入射方向間の補間** と **有限の姿勢求積点による積分誤差** は残る。
CUDA カーネルと CPU/CUDA 比較テストを同梱するが、配布側環境には CUDA Toolkit と GPU がなく、
CUDA のコンパイル・実行と Windows/MSVC は未検証。CPU の試験結果は別紙参照。

## 1. 高価な計算と配向モデルを分離する

```text
六角柱 + 波長に対して評価済みの屈折率
  -> 対称基本領域の固定入射格子
  -> 固定面積求積 + M1 trace_from_entry
  -> mass_mm2.npy (再利用する単粒子応答)

単粒子応答 + 配向分布 + 複数の世界入射方向
  -> 固定姿勢求積と既存応答の参照
  -> phi 反転 / phi 循環シフトの重みを事前集約
  -> 非負の直接畳み込み
  -> 入射条件ごとに一度だけ正規化
  -> metadata.json と3個の CDF/座標 NPY
```

単粒子応答は CDF ではない。出射セルへ入る断面積次元の光束を保持する。
配向を変更しても単粒子の光路を再計算しない。形状・屈折率を変更する場合は事前計算をやり直す。
現在は自動的な氷の分散モデルを含まない。`--ior` は指定波長で呼び出し側が与える実屈折率。
波長だけ変更し同じ `--ior` を与えた場合、非分散の同一光学条件になる。

## 2. 決定論的な入射位置と光路

側面6枚を長方形、底面2枚を各6個の三角形に分ける。
各面/三角形の面積一様な2次元パラメータ領域を S x S セルに分け、
各セルに Gauss-Legendre 2 x 3 点を配置する。三角形は面積保存の sqrt マッピングを使う。
各点の重みは対応する面の面積、投影余弦、セル面積、Gauss 重みの積。
**面ごとの本数が同じでも光束を均等配分しない。** 投影面積を重みに含める。

各正の投影面積を持つ領域を必ず被覆する。ヒット率推定や棄却法は使わない。
入力面積重みの総和と解析的投影面積を照合する。

光路は既存の `trace_from_entry()` による R, TT, TRT, TRRT, ... の決定論的分岐。
二つの入力偏光に対する JonesFlux、全反射の相対位相、光束収支をそのまま使用する。
辺・頂点の曖昧さ、上限到達、数値異常は失敗とする。乱数 jitter や黙った光路除去で回避しない。
`--max-internal-hits` は今回の既定値2048、`--tail-tolerance` は1e-12。
未回収光束は独立に保存し、正規化で隠さない。

## 3. 対称領域と入射補間（重要）

柱軸は body +y。入射極角 t はこの軸から測る。入射方位角 p は x-z 平面。
事前計算の領域は t in [0,pi/2], p in [0,pi/6]。
`--incident-theta-intervals 6 --incident-phi-intervals 3` は、端点を含む **7 x 4 = 28節点**。
領域の立体角は全体の1/24。ただし端点重複があるため節点数の削減率が常に24倍になるわけではない。

```text
k(t,p)  = (sin(t)cos(p), cos(t), sin(t)sin(p))
e0(t,p) = (cos(t)cos(p), -sin(t), cos(t)sin(p))
e1       = cross(k,e0)
```

実際の相対入射 s=R^T ki を D6h の24個の直交対称操作で基本領域へ折り畳む。
**O(3) の鏡映を含む。鏡映 Q を物理的姿勢 R in SO(3) と誤認しない。**
鏡映では出射方位角の符号を反転する。姿勢 R と Q の合成からシフト量を求める。

入射格子間では、上記局所フレームを同一視し、
**投影面積で割った応答 M_a/A_a を t,p に関して双線形補間**する。
これを、実際の照会方向 s で解析計算した A(s) で重み付けする。
割るのは投影面積 A_a であり、回収済み光束 E_a ではない。
したがって未回収比率は残り、同じ補間係数で未回収光束も伝播する。

この補間は有限入力格子の明示的な数値モデルであり、光学応答が本当に双線形という仮定の証明ではない。
ピーク位置や光路の切り替わりを格子間で厳密に復元できるわけではない。
補間誤差の収束認定は行わないが、既知の入力範囲外への外挿や最寄り点への黙った置換もしない。

出射 t_scatter は入射軸からの散乱角。回転しても変わらない。
出射 phi は、鏡映/循環シフトだけで共通フレームへ移せる。
非整数セルのシフトは区分一定セル同士の重なり幅によって質量を分配する。
これはその有限セルモデルに対する保守的な再配置であり、画像補間/ガウシアンフィルタではない。
FFT の負の丸め値や floor を避けるため、現在は非負の直接シフト和を使う。

## 4. 配向モデルと数値積分

Forster, Weber & Mayer (2025), **CrystalTrace**, Atmos. Meas. Tech. 18, 7853–7863,
DOI: 10.5194/amt-18-7853-2025, Section 2, Step 2 の singly oriented columns 型を参照する。
同論文の MC アルゴリズムを移植するのではなく、記載された Gaussian 傾斜モデルを決定論的に積分する。

- `horizontal_gaussian`: body +y 軸の鉛直からの極角 beta が N(pi/2,sigma^2)。
  密度は **d beta に対するもの**。物理的範囲 [0,pi] に明示的に切り詰めて正規化する。
  追加の sin(beta) を掛けない。sigma=0 は厳密水平というデルタモデル。
- `isotropic`: SO(3) 一様。cos(beta) に対する一様な固定求積を使う。
- `fixed`: 明示的な unit quaternion [w,x,y,z] と個数比。

R = Ry(alpha) Rz(beta) Ry(gamma)。alpha は [0,2pi) 一様。
上の二つの連続モデルでは gamma は一様であり、六回対称性により [0,pi/3) に削減する。
Parry 型の軸まわり拘束はこの実装に含めない。分布名を変えずに代用しない。

Gaussian では CDF 上の等確率区間を固定し、その中央確率の逆 CDF を64回の決定論的二分法で求める。
各 beta 区間の確率は定義上1/N_beta。alpha,gamma も等幅区間の中点で評価する。
配向平均はこの固定テンソル積求積であり、乱数・seed は存在しない。
少ない求積点で全ての細いピークが解像するという意味ではない。

混合係数 `fraction` は **粒子個数比**。各照会方向/姿勢での投影面積を含めて光束を混合するため、
散乱光束の混合比が fraction と同じになるとは限らない。総和1の正の比率を要求する。

付属設定の sigma=10度、水平成分0.9などは動作確認用のパラメータであり、観測から推定した標準値ではない。

## 5. ビルド

CPU のみ（ルートではなく standalone）：

```bat
cmake -S hex_atlas -B build-atlas-host -DICE_ATLAS_ENABLE_CUDA=OFF -DICE_ATLAS_BUILD_TESTS=ON
cmake --build build-atlas-host --config Release
ctest --test-dir build-atlas-host -C Release --output-on-failure
```

ルートから CUDA 版を追加（既存 CUDA/OptiX/architecture 設定を維持）：

```bat
cmake -S . -B build -DICE_BUILD_HEX_ATLAS=ON -DICE_ATLAS_ENABLE_CUDA=ON -DICE_ATLAS_BUILD_TESTS=ON
cmake --build build --config Release --target ice_precompute_atlas ice_precompute_atlas_cpu ice_atlas_host_tests ice_atlas_cuda_tests
ctest --test-dir build -C Release -R "^ice_atlas_" --output-on-failure
```

NumPy が必要：`python -m pip install -r hex_atlas\python\requirements.txt`。
NumPy を検出できない場合 Python/CLI テストは登録されないので `ctest -N` で確認する。
新規の `ice_atlas_` は CPUでは3件、CUDAを有効にしたルートでは5件。
Standalone の全 CTest には依存する M1 の2件も含まれる。

## 6. 事前計算は一度だけ

Visual Studio の出力配置の例：

```bat
.\build\hex_atlas\Release\ice_precompute_atlas.exe --out out\atlas_550 --radius-mm 0.1 --length-mm 0.2 --material constant --ior 1.31 --wavelength-nm 550 --incident-theta-intervals 6 --incident-phi-intervals 3 --theta 36 --phi 72 --area-subdivisions 2
```

CPU standalone の場合、実行ファイルを `build-atlas-host\Release\ice_precompute_atlas_cpu.exe` に置き換える。
CMake の単一構成 generator では `Release\` を外す。

- 外接円半径0.1 mm、全長0.2 mm。雨滴版の等体積球半径ではない。
- n=1.31は指定された実屈折率。自動分散、吸収、複屈折、回折、経路間干渉は含まない。
- 出射格子36 x 72、格納時の縮約/追加の平滑化なし。
- これは低解像度の統合確認例で、収束済みGT設定ではない。
- 追加メモリの中心は GPU per-ray histogram、16*batch_rays*Ntheta*Nphi バイト。
  `--batch-rays 64`、`--max-work-mib 512` が既定。上限超過は拒否し、勝手に解像度を変更しない。

中間出力：

```text
out/atlas_550/
    metadata.json       # atlasの別schema、ノードの面積・残差・フレーム・監査
    mass_mm2.npy        # [incident_theta,incident_phi,theta,phi], FP64
    u_edges.npy
```

## 7. 複数入射条件へ合成

```bat
python hex_atlas\python\compose_atlas.py --atlas out\atlas_550 --config hex_atlas\configs\horizontal_gaussian.json --out out\deterministic_phase
```

この設定は入射傾斜 [0,15,30,45,60,75,90] 度の7条件。
`ki=(cos(alpha),-sin(alpha),0)` という rainbow と同じ入射進行方向の規約を使う。
config の `directions` に三次元ベクトルの列を与えることも可能（inclinations_deg と同時指定は拒否）。
入射条件同士は平均しない。各条件内で姿勢を平均し、別々の CDF を作る。

```text
out/deterministic_phase/
    batch_plan.json                # 完了したレコードのインデックスとハッシュ
    generation_plan.json           # 開始時の計画（complete=false の計画スナップショット）
    orientation_model.json
    orientation_rotations.npy
    orientation_weights.npy
    records/i0000/w0000/
        metadata.json
        phi_cdf.npy
        theta_given_phi_cdf.npy
        u_edges.npy
    records/i0001/w0000/ ...
```

1レコードの形式は既存互換。最上位 `batch_plan.json` の schema は M3 独自の明示的インデックス。
NF 側が既存の一括生成 manifest の特定キーをハードコードしている場合、そのインデックス読み取りは別途対応が必要。
`PhaseRecord` に渡すのは各 `records/iXXXX/wXXXX`。単粒子atlasディレクトリではない。

配向幅を3度へ変更して再利用：

```bat
python hex_atlas\python\compose_atlas.py --atlas out\atlas_550 --config hex_atlas\configs\horizontal_gaussian_narrow.json --out out\deterministic_phase_narrow
```

ここには事前計算実行ファイルへの呼び出しがない。光路を再追跡しない。
複数波長は、各波長の評価済み屈折率でatlasを作り、
`--atlas out\atlas_550 out\atlas_650` のように渡す。同じ形状/格子を要求し、
1個の配向結合重みを波長間でも再利用する。波長を黙って補間しない。

## 8. 保存・監査・限界

出射セルは立体角密度一定。phi周辺CDF、theta|phi条件付きCDFは FP64 の順序付き走査で構築する。
CDF復元質量の L1 と lost-mass を監査し、陰に補修しない。HG g は保存CDFのセル平均余弦から求める。
合成器は NumPy CPU 実装を使い、既存GPU `PhaseCdf` を実行する方式ではない。
その代わり、既存の **無変更 reader** との end-to-end 回帰テストを含める。

前方デルタ状成分は除去しない。有限極冠セルでの表現であり、物理的回折とは区別する。
未回収光束も既存応答と同じ入射補間で混合し、総光束の監査を通してから正規化する。
全ての出力は `.part` に生成し、完了時だけ確定する。既存出力の上書き/自動再開/失敗のskipは実装しない。
途中失敗の `.part` は勝手に削除しない。

入力格子補間、入射面積求積、姿勢積分、有限出射セル、幾何光学モデルの誤差は別である。
`complete=true` は全レコードの生成/数値監査完了のみを意味し、物理GT認定ではない。
