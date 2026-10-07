# 六角柱トレーサ M1

対象基準: `Somatic-Cell/ice-crystal-phase` の
`a9538941dcb9df2de6a4130fa66f625f0133c949`。

## 範囲と未検証事項

この差分は、**固定配向の正六角柱に対する、偏光を保持した幾何光学トレーサ**です。
交差判定、投影面積、投影領域一様の入射点、反射・屈折・全反射、出射点群、
各入射レイの光束監査を実装しています。CPU と CUDA の追跡コアは同じ FP64 コードです。
**CUDA ホスト側とカーネルも実装を含みますが、提供環境に Toolkit/GPU がないため、
CUDA コンパイル・実行および Windows/MSVC は未検証です。CPU テストの合格を GPU の合格とはしません。**
CPU 専用ビルドは明示的な参照・検証経路です。CUDA の失敗で自動的に CPU へ切り替えません。

配向平均、球面セル集計、CDF 接続、氷の分散データの読み込みはこの M1 の範囲外です。
`phase_cdf.*`、既存 NPY 形式、既存 `PhaseRecord`、雨滴の数値実装は変更していません。
M1 の CSV は診断・レビュー用です。将来の本生成処理が CSV を経由する設計ではありません。
本段階だけで NF 用の完成した参照位相関数を出力するものではありません。

採用モデルは非吸収・非磁性・等方誘電体・理想平面です。複屈折、粗面、回折、
異なる光路間の干渉は含みません。全反射の s/p 相対位相は偏光状態に影響するため保持します。
単一光路の全成分に共通する伝播位相は、経路間で振幅合成しないので計算しません。
`--ior` は指定波長で呼び出し側が評価した実屈折率であり、既定値はありません。
`--wavelength-nm` は必須ラベルですが、この段階で氷の光学定数を自動評価する機能はありません。
例の `--ior 1.31` は数値検証条件です。氷の分散を実装済みと主張する値ではありません。

## ファイルと読む順番

1. `include/ice_crystal/hex_prism.hpp`: 形状、投影面積、入射点、半空間交差。
2. `include/ice_crystal/flux_interface.hpp`: FP64 の光束用 Fresnel/Jones 応答。
3. `include/ice_crystal/hex_trace_core.hpp`: 分岐追跡、終了条件、出射レコード、監査。
4. `src/hex_trace.cpp`: ホスト検証、CPU バッチ、集計。
5. `device/hex_trace.cu` / `src/hex_trace_cuda.cpp`: CUDA の count → prefix → replay/write。
6. `apps/trace_hex.cpp`: バッチ実行、診断出力、失敗時の `.part` 保持。

幾何・光学のヘッダ内実装は CPU/CUDA で同じ数式をコンパイルするためです。
ホスト処理、所有権、I/O は `.cpp` に分離しています。
FP64 の型と関数は `iceCrystal` 名前空間に追加しました。
既存の FP32 `DielectricInterface` / `JonesResponse32` の意味や精度を変更していません。
既存の `CudaContext`, `CudaModule`, `DeviceBuffer` を GPU の所有権管理に利用します。

## 形状・座標

中心原点、右手系、柱軸 `+y`。`R=circumradius_mm` は外接円半径、
`L=length_mm` は全長です。雨滴の等体積球半径とは異なります。

- 側面 0..5: 法線 `(cos(j*pi/3),0,sin(j*pi/3))`、平面距離 `sqrt(3)*R/2`。
- 面 6: `+y` の底面。面 7: `-y` の底面。
- 底面頂点角は `pi/6 + j*pi/3`。面0の正面は `+x`、頂点方向とは30度違います。
- `ki` と `ko` は光の実際の進行方向。方向は `make_trace_settings()` で明示的に正規化。
- `Frame` は `e1 = cross(k,e0)`。各出射レコードは自身の偏光基底を保持します。
- 幾何コアの位置・距離は `R` で割った無次元値。`OutgoingSample` の位置と内部長は mm。

8半空間の凸多面体を直接扱い、BVH/OptiX を交差に使いません。
これは三角形近似を避けた解析的平面交差の選択であり、処理速度の比較結果ではありません。
CUDA は Driver API のみで起動します。

投影面積は `A = sum_f area(f)*max(0,-dot(n_f,ki))`。
この項に比例して入射面を選択し、その面上を一様にサンプルします。
結果は投影シルエット内一様です。初回ヒット率を推定する必要はありません。
方位角の対称性による削減はまだ適用しません。任意の三次元 `ki` を受け取ります。

## 光束・偏光・経路

`JonesFlux` の二列は二つの直交入力偏光への応答です。
無偏光を単一の Jones ベクトルで表すものではありません。
単位無偏光入射に対する残存光束比は

```
power_fraction = 0.5 * sum(abs(J_ab)^2)
```

です。透過振幅は

```
t_flux = sqrt(n_t*cos(theta_t)/(n_i*cos(theta_i))) * t_E
```

に対応する係数です。コードでは同値な対称形で計算し、斜入射の大きな前置係数と
小さな電場透過係数を別々に計算することを避けています。
s/p 基底を毎回回転し、全反射では複素位相を掛けます。
各面で偏光を無偏光にリセットしません。異なる出射点の複素振幅は合成しません。

凸な一粒子なので、外へ出た枝はその粒子へ再入射しません。
最初の外部反射を出力し、内部では各界面からの透過枝を出力、反射枝だけを継続します。
このため保持する次数までの全枝を、指数的な木ではなく線形の内部反射列として列挙できます。
枝をランダムに一つ選ぶ処理も Russian roulette も使いません。
外部反射、直接前方透過、高次内部反射経路を、方向を理由に除外しません。

入射レイ k の出射 b の面積重みは

```
w_kb = A / total_incident_samples * power_fraction_kb
```

です。バッチサイズでは割りません。次段の配向合成で初めて姿勢の重みを掛けます。
`total_incident_samples` は全バッチ共通です。seed とグローバル sample ID により
バッチの区切りを変えても同じ入射点を生成します。

## 終了・失敗の契約

デフォルトは `residual_power_tolerance=1e-12`, `max_internal_hits=256`。
いずれも全レイの収束認定値ではありません。

- 光束が正確に0なら `exhausted`。
- 正の残存光束が設定閾値以下なら `tail_tolerance_reached`。残存量自体も保存。
- 上限までに閾値に届かなければ `interaction_limit`。不合格。
- 辺・頂点で法線を一意に選べない場合は `ambiguous_boundary`。不合格。
- 無効な入力、算術異常、出力容量不足、count/write不一致も不合格。

交差の許容幅は無次元座標で `128*DBL_EPSILON*max(1,half_length,|p|,|t|)`。
これは厳密な区間誤差保証ではありません。狭い形状やほぼ辺上の経路では
不合格になることがあります。未知の法線で続行したり、失敗レイを捨てて正常扱いしたりしません。
レイ原点に固定 epsilon を足して自交差を避けることもしていません。
初回面と次回面を識別して処理します。

光束収支は `escaped + unresolved - 1` を監査します。CLI の許容値は `1e-10`。
失敗・打ち切りで逃した光束を、出射光束の再正規化で隠しません。
`complete=true` は全予定レイがこの数値的受理条件を満たした診断出力の完了を示します。
角度分布や入射サンプル数の収束、Maxwell 解との一致を認定するものではありません。

テストで見つけた長い全反射列では、256回で止めると約0.9813の光束が未回収です。
このケースは意図どおり拒否され、上限を1024に設定すると366回で閾値に到達します。
この回帰試験を `long_TIR_residual` として含めています。上限を増やすことはユーザーの明示操作です。

## ビルド

### 既存リポジトリのルートから（CUDA）

既存の CUDA/OptiX 設定を維持して、追加対象だけをビルドします。
`CMakeLists.txt` への変更は末尾の `add_subdirectory(hex_trace)` 用5行だけです。

```powershell
cmake -S . -B build -DBUILD_TESTING=ON -DICE_BUILD_HEX_TRACE=ON -DICE_HEX_ENABLE_CUDA=ON -DICE_HEX_BUILD_TESTS=ON
cmake --build build --config Release --target ice_trace_hex ice_trace_hex_cpu ice_hex_host_tests ice_hex_cuda_tests
ctest --test-dir build -C Release -R "^ice_hex_" --output-on-failure
```

Visual Studio generator では `build/hex_trace/Release/` に実行ファイルが置かれます。
fatbin は対応する実行ファイルの `modules/hex_trace.fatbin` に配置します。
Ninja の場合は `Release/` を挟みません。
CUDA の構造体サイズと ABI識別子をロード時に照合します。
実行時には CUDA が利用できなければ失敗します。成功やスキップとして扱いません。

### CPU参照・ホスト単体テスト（CUDA/OptiX不要）

```powershell
cmake -S hex_trace -B build-hex-host -DICE_HEX_ENABLE_CUDA=OFF -DICE_HEX_BUILD_TESTS=ON
cmake --build build-hex-host --config Release
ctest --test-dir build-hex-host -C Release --output-on-failure
```

### 六角柱サブディレクトリだけを CUDA ビルド

```powershell
cmake -S hex_trace -B build-hex-cuda -DICE_HEX_ENABLE_CUDA=ON -DICE_HEX_BUILD_TESTS=ON
cmake --build build-hex-cuda --config Release
ctest --test-dir build-hex-cuda -C Release --output-on-failure
```

この経路は既存リポジトリの Driver API ラッパー4ソースをリンクするので、
`hex_trace` だけを他所へ移動して CUDA ビルドする使い方ではありません。
OptiX はこの standalone CUDA ビルドでは要求しません。
アーキテクチャは既存 root の設定を継承し、standalone の未指定時は86です。
使用環境に合わせる場合は `-DCMAKE_CUDA_ARCHITECTURES=...` を明示します。

## 実行例（数値検証用の入力屈折率）

```powershell
.\build\hex_trace\Release\ice_trace_hex.exe --out out\hex_m1 --radius-mm 0.1 --length-mm 0.2 --ki 1 -0.4 0.3 --ior 1.31 --wavelength-nm 550 --samples 4096 --seed 12345 --max-internal-hits 1024 --tail-tolerance 1e-12
```

CPU 参照ビルドでは実行ファイルを `build-hex-host/Release/ice_trace_hex_cpu.exe` に変更します。
上書き禁止です。再実行には別の出力先を指定してください。

```
out/hex_m1/
  trace_metadata.json   # ice.hex_trace.points.v1。phase_cdfではない
  outgoing.csv         # 出射方向・位置・Jones応答・光束重み・経路情報
  rays.csv             # 各入射レイの状態・残存光束・収支
```

終了コード0: 全予定レイが受理され、`.part` を最終ディレクトリへ変更。
終了コード2: 数値的品質不合格。`.part` に診断出力を残し、それ以降のバッチは実行しない。
終了コード1: 無効入力、I/O、GPU、メモリ等のエラー。成功レコードは作らない。
プロセス中断時やI/O失敗時は `.part` が部分状態で残ることがあります。

GPU API は `HexPrismTracer::outgoing()`, `offsets()`, `audits()` を公開し、
CUDA上のデータを次段の集計器へ渡せます。`has_result()` はバッファが得られたことだけを示し、
物理・数値品質の受理ではありません。利用前に audits を検証する責任があります。
出射データは sample-major CSR。各 ray i の範囲は `offsets[i]..offsets[i+1]`。

## テスト

ホスト11群: 投影面積、独立した三角形meshとの交差比較、入射面分布と底面の二次モーメント、
Fresnel/Snell/Brewster/TIR、偏光履歴、平行平板の解析光束、斜入射とバッチ不変性、
60度回転と一様スケール、300条件×100レイ、長い全反射列、異常入力/容量/辺上。
CLI テスト: JSON、CSV、重み、上書き禁止、不合格時の `.part`、CUDAへの偽のフォールバック防止。
GPUテスト: ホストとのCSR・方向・位置・光束・Jones全要素・異常系の比較（提供時は未実行）。
同じコアの CPU/GPU 一致は移植検証であり、独立の物理ソルバによる検証ではありません。

## 参照

- Mitsuba 3, Polarization / Fresnel equations:
  https://mitsuba.readthedocs.io/en/stable/src/key_topics/polarization.html
- NVIDIA, CUDA Driver API / modules:
  https://docs.nvidia.com/cuda/cuda-programming-guide/03-advanced/driver-api.html
- CMake, CUDA_FATBIN_COMPILATION:
  https://cmake.org/cmake/help/latest/prop_tgt/CUDA_FATBIN_COMPILATION.html

本実装はこれらのソースコードをコピーしたものではなく、上記の光束・座標契約で記述したものです。
Sadeghi2012 の雨滴用パッチ、焦線、回折フィルタを六角柱に流用していません。
