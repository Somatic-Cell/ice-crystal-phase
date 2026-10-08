# Phase CDF 解析ツール v1

## 目的と範囲

`metadata.json`, `phi_cdf.npy`, `theta_given_phi_cdf.npy`, `u_edges.npy` を読み，
保存された「セル内で立体角密度が一定」の分布を検査・可視化・比較する．
ソルバ，CDF，reader，ビルド設定は変更しない．

入力は ZIP を展開した1レコードのディレクトリ．`PhaseRecord(validate=True)` を使う．
`.part`，不正 CDF，g 不整合，未知の方向規約は拒否する．
CDF の修復，再正規化，微小値 floor，平滑化，角度の自動縮約は行わない．
本ツールは Maxwell 解・連続配向積分・幾何光学デルタ分布の検証ソルバではない．

対応する座標契約は `ice.phase_cdf.ensemble_coordinates.v1` と
`rainbow.phase_cdf.coordinates.v1`．比較は同じ契約・座標空間・サンプリングフレームに限る．

## インストールとテスト

リポジトリのルートで：

```bat
python -m pip install -r tools\phase_analysis_requirements.txt
python -m unittest discover -s tests -p test_phase_analysis.py -v
```

Python 3.10 以上．CMake/CUDA の再ビルドは不要．既存 CTest は変更しないので，
新しい Python テストは上のコマンドで別に実行する．

## 1レコードの検査

```bat
python tools\phase_analysis.py inspect out\phase_m2 --out analysis\phase_m2
```

`--out` は存在しないディレクトリにする．既存ファイルへの上書きや出力の混在を拒否する．
生成物は `analysis.json`, `theta_profile.json`，5枚の PNG，最後に `analysis_complete.json`．
例外後の出力には完了ファイルが存在しない場合があり，完成した解析結果として扱わない．

`--no-plots` を付ければ NumPy のみで数値検査が可能．入力ファイルの SHA256 を記録する．
JSON は有限数のみ・ASCII エスケープを使用し，コンソール診断も CP932 で出力可能な形式にする．

### 図の意味

- `pdf_linear.png`, `pdf_log.png`：横が phi，縦が入射軸からの散乱極角 theta．
  theta=0 は前方．色はセル内の PDF [sr^-1] で，確率質量ではない．
  セル境界を使う `pcolormesh(shading='flat')`．補間なし．対数表示で厳密なゼロはマスクし，
  最小正値への置換はしない．線形・対数で元の数値は同一．
- `azimuth_average.png`：方位角平均 PDF [sr^-1]．theta 方向の周辺 PDF とは異なる．
- `theta_marginal.png`：theta の周辺確率密度 [degree^-1]．球面ヤコビアンを含む．
- `cone_cdf.png`：前方軸を中心とする円錐内の累積確率．セル内部も u で積分する．

色やカラーマップは Matplotlib のデフォルト．線形・対数とも元のゼロ，集中成分を残す．
極冠の帯は保存規約による有限角度化を含み，回折を意味しない．

## 確率測度

A は phi 周辺 CDF，T は条件付き CDF，U は u 境界とする．

```python
P = np.diff(A)[:, None] * np.diff(T, axis=1)  # [phi, theta]
dOmega = (4*np.pi/Nphi) * np.diff(U)        # [theta]
pdf = P / dOmega[None, :]
g = sum(P[j,i] * (1-U[i]-U[i+1])) / sum(P)
```

質量を密度と混同しない．phi のセル幅は規約どおり 1/Nphi を共通に使い，
丸めた描画座標 v_j=j/Nphi の差分で PDF を再定義しない．

theta バンド i の確率を m_i とすると，方位角平均 PDF は m_i/(4*pi*du_i)．
周辺 PDF は radian 測度で m_i*sin(theta)/(2*du_i)，degree 測度ではさらに pi/180 を掛ける．
後者はセル内でも theta に関して一定ではない．

## 2レコードの比較

```bat
python tools\phase_analysis.py compare out\reference out\candidate --out analysis\comparison
```

異なる角度解像度にも対応．粗い格子へ両方を縮約してから比較するのではなく，
保存 U 境界の和集合，および phi の有理数境界の和集合を用いた共通細分で積分する．
したがって保存分布のセル内情報を意図的に捨てない．
U の1 ULP の違いも区間として残す．丸めで中点が隣接セルへ入らないよう左境界でセルを同定する．
phi の分割は最小公倍数に基づく整数座標を用い，j/Nphi の丸めに由来する偽の細区間を作らない．

出力は TV (= L1/2)，L2 ノルム [sr^-1/2]，Hellinger 距離，Jensen-Shannon divergence [nats]，
片方のみで支持される確率質量，theta 周辺 TV，g の差，入力ハッシュ．
JS はほぼ等しい密度での桁落ちを避ける形で評価する．真の KL が無限大になる場合を
有限の floor で隠さないため，この版では KL 自体を指標にしていない．
L1/TV は密度の二乗誤差ではなく，セル積分された確率質量の絶対差に対応する．

### 条件が違う比較

屈折率・波長・粒子形状・主要な光学モデルが違う場合は既定で拒否する．
物理条件を意図的に比較する場合に限り `--allow-physical-differences` を付ける．
配向 CSV のハッシュが違う／ない場合は `--allow-orientation-plan-difference` が必要．
このフラグは「同じ配向分布である」と認定するものではない．異なる配向数の収束実験では，
同一の連続配向分布を使ったことを別の実験記録で保証する．CSV の改行コードだけでも
ハッシュは変わるので，無条件にハッシュ不一致を無視しない．

フレームは最大絶対差 2e-12 以下を丸めによる等価として扱い，差をレポートする．
それを超える場合は明示的に拒否する．自動回転や補間による救済は実装しない．

## メモリと計算量

入力は最大 `--max-cells`（既定1000万セル）．mmap の reader を使用するが，
解析では質量・密度の dense 配列を作るためメモリ量は O(Ntheta*Nphi) であり，
巨大データを定数メモリで解析する実装ではない．
比較の中間配列は共通細分の1 phi 行ずつ処理する．共通細分の総セル数は
`--max-refined-cells`（既定5000万）以下に制限する．超過時は自動で解像度を下げず拒否する．

## 評価上の注意

1. `analysis.json` では保存配列から再計算した量と，生成器の metadata に報告された量を分ける．
   生のレイや縮約前ヒストグラムがないため，未回収光束や生成時の縮約 TV をこの1ファイルから
   独立に再測定することはできない．
2. `coarsening_tv` の警告閾値0.05は表示用の助言で，物理的許容誤差でも受理条件でもない．
   ソルバの停止基準・実験条件を変更しない．TV が大きくても質量積分が1であることは矛盾しない．
3. 2姿勢の離散混合をレイ数だけ増やしても，連続的な配向分布への積分は実現しない．
4. 1組の比較では収束の認定や信頼区間を出さない．同じ配向リストでレイ数を増やす軸，
   同じ連続分布で配向数を増やす軸，角度格子を変える軸，乱数の独立反復を分ける．
5. 生のデルタ成分と有限セル PDF は別の数学的対象．保存 PDF の TV を幾何光学のデルタ分布
   への TV 収束と読み替えない．角度窓の確率・モーメントも併せて評価する．

## 検証履歴

実データ phase_m2.zip と，前回配布 CPU 参照CDFの読み込み・図生成・比較を実行済み．
19件のユニットテストを追加．一様球面，前方極冠，支持が交わらない分布，異なる非入れ子格子，
共通細分での情報保存，1 ULP の境界，極小 JS，フレーム不一致，メモリ制限，破損 CDF，
g 不整合，入力不変性，既存出力拒否，CSV ハッシュ不一致などを含む．
環境は Linux / Python 3.13.5 / NumPy 2.3.5 / Matplotlib 3.10.8．
Windows 実機，CUDA ソルバの追加実行はこの解析ツールの検証に含めない．

API資料：
- NumPy load: https://numpy.org/doc/stable/reference/generated/numpy.load.html
- Matplotlib flat shading: https://matplotlib.org/stable/gallery/images_contours_and_fields/pcolormesh_grids.html
