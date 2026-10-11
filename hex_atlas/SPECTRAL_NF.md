# M4: 単一粒子・波長/入射角データセット

基準: `1cf6758d28473678ed141654f8c0c9b4f21197de` (2026-10-11 に確認した main)。

## 今回の範囲

- 氷 Ih の実屈折率を真空波長から自動評価する (Warren–Brandt 2008)。
- 同じ形状・指定した配向分布で、複数波長と複数入射角の CDF を作る。
- 最終レコードを `spherical-hg-flow` の**無変更の RainbowReference**で検査する。
- 完了した単粒子 atlas を同一ビルド・同一条件で再利用する簡単なバッチ実行を追加する。

形状・粒径・乱流から配向幅を予測するモデルは**追加していない**。
柱/板/Parry を網羅する拡張、回折・吸収・複屈折・経路間干渉も追加しない。
M3 の固定面積求積、R/TT/TRT/... の全保持次数の追跡、決定論的配向求積、
投影面積の重み、非負の直接循環和、鉛直入射の解析的方位積分を変更していない。
これは現在の幾何光学/有限格子モデルを学習へ渡すための実装であり、Maxwell 基準解ではない。

## ビルド

リポジトリのルートから、既に動作した CUDA/OptiX/architecture 設定を維持して実行する。
`ICE_NF_SOURCE_DIR` は実際の学習プロジェクトのディレクトリに合わせる。

```bat
cmake -S . -B build -DBUILD_TESTING=ON -DICE_BUILD_HEX_ATLAS=ON -DICE_ATLAS_ENABLE_CUDA=ON -DICE_ATLAS_BUILD_TESTS=ON -DICE_NF_SOURCE_DIR=C:/Users/sy415/workspace/spherical-hg-flow
cmake --build build --config Release --target ice_precompute_atlas ice_precompute_atlas_cpu ice_optics_host_tests ice_atlas_host_tests ice_atlas_cuda_tests
ctest --test-dir build -C Release -R "^ice_atlas_" --output-on-failure
```

CPU のみなら以下。CUDA/OptiX 不要。

```bat
cmake -S hex_atlas -B build-atlas-host -DICE_ATLAS_ENABLE_CUDA=OFF -DICE_ATLAS_BUILD_TESTS=ON -DICE_NF_SOURCE_DIR=C:/Users/sy415/workspace/spherical-hg-flow
cmake --build build-atlas-host --config Release
ctest --test-dir build-atlas-host -C Release --output-on-failure
```

外部 NF プロジェクトを指定しなければ、その integration test だけ登録しない。
その場合に他のテストが通っても NF integration の合格とは言わない。
手動の `check_nf_compatibility.py` は生成後に使用できる。Python/NumPy が必要。
**この差分の配布側検証は CPU のみ。CUDA・Windows/MSVC の新しい実機試験は未実施。**
CUDA カーネル、Driver API ラッパー、実行時ライブラリの方式には変更を加えていない。

## 1. 先に12条件だけ生成して log PDF を確認

`configs/spectral_preview.json` は次を明示する (収束済みの認定条件ではない)。

- 外接円半径 R=0.1mm、全長 L=0.2mm。これまでの例と同じ形状。L/(2R)=1。
  細長い柱や薄板の観測代表形状と主張するものではない。
- c 軸の水平付近 Gaussian、sigma=3度、個数比0.9、SO(3)一様の個数比0.1。
  測定フィットではなく固定した基準分布。密度は d beta、追加 sin(beta) なし。
- 固定配向求積: azimuth=72、tilt=12、spin=24。乱数なし。
- 相対入射格子: theta/phi とも2度 (45区間と15区間、46x16節点)。
- 出射格子: 180x360 (1度)、縮約・Gaussian 後処理なし。
- 面積求積: subdivisions=2。GPU batch-rays=32、max-work-mib=512。
- 波長: 450/550/650nm。入射傾斜: 0/30/60/90度。

```bat
python -m pip install -r hex_atlas\python\requirements.txt
python hex_atlas\tools\run_spectrum.py --exe build\hex_atlas\Release\ice_precompute_atlas.exe --plan hex_atlas\configs\spectral_preview.json --atlas-root out\m4_atlases --out out\m4_preview --dry-run
python hex_atlas\tools\run_spectrum.py --exe build\hex_atlas\Release\ice_precompute_atlas.exe --plan hex_atlas\configs\spectral_preview.json --atlas-root out\m4_atlases --out out\m4_preview
```

`--dry-run` は容量目安と実行予定を表示するだけで、ファイルを作らない。
通常実行は前半の atlas 作成だけを指定バイナリで行い、後半は CPU/NumPy で合成する。
CUDA に失敗しても CPU へ切り替えない。子プロセスの出力はバイト列のログとして保存する。

**--material ice が通常モード。--ior は指定しない。**
従来の一定値試験だけ `--material constant --ior 1.31` を使う。
波長から CPU 上で n(lambda) を一度計算し、その値を CPU/CUDA トレーサ双方に渡す。
450/550/650nm の n はそれぞれ 1.3157/1.3110/1.3080。
200..1100nm の範囲をサポートし、範囲外は拒否する。温度補正式は追加していない。
虚部は参照値を metadata に残すが、solver の虚部は0。外部媒質は通常1。
詳細な出典・補間規則は `data/README.md`。

```bat
python -m pip install -r tools\phase_analysis_requirements.txt
python tools\phase_analysis.py inspect out\m4_preview\records\i0001\w0001 --out analysis\m4_30deg_550nm --label "M4: incidence 30 deg, 550 nm"
start "" "analysis\m4_30deg_550nm\pdf_log.png"
```

`i0001/w0001` はこの**preview設定**では30度/550nm。別設定ではindexの意味が変わるため、
必ず batch_plan/metadata の実際の条件値を読む。画像の色範囲は各図で自動設定される。
ゼロはゼロのまま保存し、log可視化では空白。NaN/Inf は拒否する。
既存可視化ツールの「two orientations」「CSV」注意文は旧 M2 前提で、M3/M4 に合わない場合がある。
配向の実設定と固定節点はデータセット直下に保存されている。

目視で判断する対象: 主要ピーク/アークの形と入射条件依存、波長間の変化、
穴/周期模様/格子段差の有無。鉛直入射では一様な方位角成分が横一様になること。
**滑らかに見えることだけで物理的正しさは認定しない。**
固定求積由来の細かい模様が残る場合、配向求積を増やして同じ atlas を再合成できる。
出射格子だけを増やすと、相対入射補間/面積求積の誤差を詳細表示するだけの場合がある。
floor や画像ぼかしで穴を埋めて学習データを作る処理は実装していない。

## 2. 無変更の学習側 reader で最終確認

```bat
python hex_atlas\tools\check_nf_compatibility.py --dataset out\m4_preview --nf-root C:\Users\sy415\workspace\spherical-hg-flow --out out\m4_preview_nf_check.json
```

PyTorchや学習は起動しない。学習側 checkout の `src/phaseflow/rainbow.py` を直接ロードし、
コードの SHA256 と Git blob SHA1 を記録する。全レコードについて、条件、質量、g、
各セル中心の固体角 PDF と NF フレームへの変換、CDF 逆変換と log_prob の一致を検査する。
データの書換え/修復/再正規化はしない。out JSON も上書きしない。

- 条件: `[wavelength_nm, -ki_y]` (= `[wavelength_nm, sin(inclination)]`)。
- 配列: `<f8`, C-order、phi_cdf `(Nphi+1,)`、theta_given_phi `(Nphi,Ntheta+1)`、u_edges `(Ntheta+1,)`。
- 最終レコードの共通座標系は粒子集団の輸送参照系。各結晶の body +y と混同しない。
  `shape_polar_axis=[0,-1,0]` は既存 NF の条件付けの軸として定義し、意味を metadata に併記する。
- hg.g は保存セル PDF の一次モーメント。axis_particle_frame は保存入射軸。
- 材質と形状は氷/六角柱として記録する。水や等体積球半径と偽装しない。
- 任意方位角の非canonical入射は既定の rainbow 契約では拒否する。
  旧来の ensemble 契約は compose_atlas.py の `--output-contract ensemble` で明示選択できる。
- NF の既存 reader は `provenance.source_repository` に rainbow をハードコードしている。
  元の metadata.source_repository は正しく ice-crystal-phase と記録する。学習コード自体は変更しない。

**CDF単位のreader/サンプラー契約を検査するもので、条件付き学習全体やGPU学習の実行保証ではない。**
`batch_plan.json` は既存の M3 インデックスを維持し、各レコードに condition を明示追加する。
現在のNFが1レコード指定の学習CLIを使っている場合も、そのレコードは直接指定できる。
複数レコードの列挙は学習側の既存ワークフローに渡す。

## 3. 確認後に波長・入射条件だけを増やす

`configs/spectral_dataset.json` は同じ高価な atlas 設定・形状・配向を使用し、
390..830nm を10nm間隔 (45波長)、-90..90度を5度間隔 (37入射条件) の1665レコードを作る例。
波長/入射間隔は実験計画用の明示値であり、必要精度を満たすとの保証ではない。
対称性を利用した負入射角の自動生成/コピーは行わず、指定条件をそのまま合成する。

```bat
python hex_atlas\tools\run_spectrum.py --exe build\hex_atlas\Release\ice_precompute_atlas.exe --plan hex_atlas\configs\spectral_dataset.json --atlas-root out\m4_atlases --out out\m4_dataset --reuse
python hex_atlas\tools\check_nf_compatibility.py --dataset out\m4_dataset --nf-root C:\Users\sy415\workspace\spherical-hg-flow --out out\m4_dataset_nf_check.json
```

previewの450/550/650nmは再追跡せず利用し、他の波長だけ追加生成する。
既存の一定屈折率 atlas を分散したものとして再利用しない。
`--reuse` は実行引数・実行ファイル/隣接fatbinのバイト一致・metadata一致を照合し、
さらに合成器が全 node の光束を監査する。数値配列を無条件に信用してスキップしない。
同じビルドを別PCへコピーした場合は再利用できるが、再ビルドでバイナリが変わった場合は
保守的に拒否する。設定変更時に既存 atlas を黙って流用しない。
未完了 .part は自動削除・採用しない。新たなCDF出力先を使い、失敗ログを先に確認する。

生 atlas は約0.355GiB/波長、45波長で約16.0GiB。CDF配列本体は全体で約0.82GiB。
metadata/ログ/一時ファイルには別の領域を見込む。圧縮ZIPの大きさで必要ディスクを見積もらない。
合成はCPU/NumPyであり、GPU事前計算だけの時間では全生成時間を見積もれない。

## 保持したモデル上の制限

本パッチでは相対入射間の双線形補間を変更していない。光路の切替わりやコースティックが
格子より細かい場合、ピークの位置/幅/弱い成分に誤差が残る。CDF integrity/光束収支が正常でも
物理的な分布形状まで正しいとは限らない。高解像度の数条件を確認してから全条件へ展開する。
今回の配向分布の幅を shape/size/気象条件に結び付ける作業は、合意した通り保留する。
