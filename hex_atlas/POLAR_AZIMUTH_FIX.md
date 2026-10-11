# M3 polar-azimuth hotfix 1

## 修正対象

付属 `horizontal_gaussian.json` の世界入射傾斜90度では，柱軸方位角 alpha が一様な
配向分布と鉛直入射の組合せなので，スカラー位相関数は出射 phi に依存しない．
M3 は alpha を16点の中点求積で近似していたため，5度幅の出射セルへ疎な単粒子応答を
シフト合成すると，周期模様と40--45度帯のゼロセルが残った．
配布済み CPU example atlas の同条件で，ゼロ24セルを再現した．ゼロは CDF 前の光束にも
存在し，CDF 丸めで正値を失ったものではなかった．この原因は CPU/CUDA 共通の合成器にある．

## 数学的な変更

R=Ry(alpha) Rz(beta) Ry(gamma)，ki=(0,+/-1,0) では R.T ki と投影面積は alpha に
依存しない．alpha は出射 phi の平行移動だけを与える．したがって，

    (1/(2*pi)) integral_0^(2*pi) K(theta,phi-alpha) d alpha

は各 theta における方位角平均であり，区分一定の有限セルモデルに対して厳密に評価できる．
修正は画像や最終 PDF にフィルタをかけず，配向成分ごとの循環畳み込みカーネルを一様にする．
単粒子の生の非正規化セル光束，投影面積の重み，未回収光束は保持する．

- `horizontal_gaussian` / `isotropic` の一様な世界方位角だけを対象とする．
- `fixed` 成分は，混合されている場合も平均しない．
- 入射 ki の x,z が厳密に0の場合だけ適用する．近鉛直を鉛直へ丸めない．
  `inclinations_deg` の +/-90 は既存 CLI でも厳密な軸方向になる．
- その他の入射方向は元のコード経路のまま．横方向の模様が物理的な場合もあるため，
  全入射角の PDF を一律に方位平均することはしない．
- beta/gamma の求積と単粒子の入射格子補間の誤差は，この修正では消えない．
- sigma，個数混合比，角度格子，光路打ち切り，CDF 許容誤差は変更しない．
- 分布に新たな一様散乱成分や PDF floor を加えない．本来ゼロの theta 帯はゼロを保つ．

`coupling()` は元の有限節点積分として残し，合成 CLI は宣言された配向モデルを参照する
`coupling_for_rule()` を使う．metadata の `orientation.azimuth_integration` に，各成分の
積分方法と実際の代表姿勢数を記録する．元の配向配列は再現用の参照求積列として残すが，
鉛直条件では一様 alpha を解析積分したことを，レコードと manifest に明示する．
`source_sha256` / `generator_version` も更新される．

## 適用・実行

前回配布した M3 を適用済みの作業ツリー向け．変更は5ファイル．
C++/CUDA/CDF形式は変更しないので再ビルド・atlas再生成は不要．
GitHub 上に直接 push はしていない．

```bat
git apply --check ice-crystal-phase-m3-polar-azimuth-hotfix1.patch
git apply ice-crystal-phase-m3-polar-azimuth-hotfix1.patch
python -m unittest discover -s hex_atlas\tests -p "test_*.py" -v

python hex_atlas\python\compose_atlas.py --atlas out\atlas_550 --config hex_atlas\configs\horizontal_gaussian.json --out out\deterministic_phase_polar_fixed

python tools\phase_analysis.py inspect out\deterministic_phase_polar_fixed\records\i0006\w0000 --out analysis\m3_polar_fixed --label "M3: vertical incidence, analytic azimuth"
start "" "analysis\m3_polar_fixed\pdf_log.png"
```

出力先と `.part` は新しい名前にする．既存データへの上書きを許可しない．
標準の全7条件 config で生成すると，0--75度の6条件の NPY は旧合成器と同じで，
90度の条件だけが修正対象になる（数値実行環境が同一の場合）．
NumPy依存以外に追加依存はない．新しい CTest `ice_atlas_polar_azimuth` の登録には
`cmake -S . -B build` が必要だが，上記の unittest コマンドには不要．

## 検証の範囲

- 既存11件＋新規8件のPython単体テストを実行．
- 配布済み CPU atlas から sigma=10/3度，各7入射角，計14レコードを再合成し，
  Git blob一致を確認した無変更の `PhaseRecord(validate=True)` で検証．
- 0--75度の6条件のNPYは，同じ環境で実行した旧合成器とバイト一致．
- +/-鉛直，非鉛直，近鉛直，固定姿勢との混合，未回収光束，ゼロtheta帯を回帰試験．
- analytic alpha と，出射セルの全循環シフトを含む72点 alpha 求積との一致を検査．
- 新しいユーザーのCUDA出力NPYは受け取っていない．提示画像の数値をOCRしたものではなく，
  前回同梱したCPU例のatlasを使った再現である．CUDA実行とWindows実機は未検証．
- 正規化・再現性だけでは対称性の正しさを検出できなかったため，新たな回帰試験を追加した．
