# 軸対称 hoop ソースの自由流保持 — 圧力ゲージ整合と離散閉性面積

## メタ

- **area**: `axisymmetric`
- **status**: `in-progress`
- **related_docs**:
  - [`methods/axisymmetric/implementation.md`](../../methods/axisymmetric/implementation.md)
- **related_plans**:
  - [`precision-mixed-axisym.md`](../accepted/precision-mixed-axisym.md) (float32 起因の近軸問題。本 plan は「倍精度化せずに直す」側)
  - [`architecture-node-centroid-value-position.md`](architecture-node-centroid-value-position.md) (case/43 の親調査)
- **created**: `2026-08-16`
- **owner**: `sano` (指示)・Claude (実装)

## 1. 目的

r 重み方式 (`axisymMethod: 0`) の半径運動量は

$$\mathrm{res}_{\rho u_r}=-\sum_f p_f\,r_f S_{f,y} \;+\; p\,A$$

で、一様圧なら $A=\sum_f r_f S_{f,y}$ のときだけ厳密に打ち消える。forge は $A$ に**解析 `A_planar`** を使う一方、
面側は float32 メトリック (`geom_float = float`) なので、**高 AR 壁 CV で両者が最大数十 % ずれ**
偽の半径力が立つ (case/43: 生産 NS y+1.5 で 59%、自由流テストで 3.4e4 m/s²)。倍精度化はコストが高い。
本 plan は**倍精度に頼らず**自由流を回復する 2 手を入れる。

## 2. 方針

1. **圧力ゲージの整合 (バグ修正)**: `space.pRef` は対流流束を $(p-p_{\rm ref})S$ で組むが、**hoop ソースは
   絶対 $p$ のまま**だったため、軸対称 × `pRef>0` では一様圧 $p=p_{\rm ref}$ で source だけが残り
   **偽半径力 $p_{\rm ref}A$** が立っていた (実測 1.25e6 m/s²)。ソースも $(p-p_{\rm ref})A$ にする。
   これで大きな $p$ が悪条件な metric 和に掛からなくなり、誤差は $|p-p_{\rm ref}|$ に比例する。
2. **離散閉性面積 (`mesh.hoopAreaFromClosure: 1`)**: $A$ を解析値でなく**同じ面ベクトルの和**
   $\sum_f r_f S_{f,y}$ (setup 時に計算、`A_closure_y`) にする。構成上、丸めに依らず flux と source が
   同じ量を参照する。`axisRFloor>0` の既存経路と同じ配列を流用し、床なしでも使えるようにする
   (x 成分の閉性補正も同時に有効化)。既定 0 でビット不変。

## 3. 検証 (case/43.node_axis_dof)

自由流テスト (`optest/freestream_test.py`, 一様圧・静止, 直管 AR 0.5–1250) の壁 CV 偽半径加速度 [m/s²]:

| AR(壁) | 解析 A | 離散閉性 A | 解析 A + pRef | 離散閉性 A + pRef |
|---|---|---|---|---|
| 2.5 | 1.33e+02 | 1.73e+01 | 8.0e−06 | 5.9e−07 |
| 62.5 | 1.57e+03 | 1.11e+02 | 1.0e−04 | 7.3e−06 |
| 1250 | 3.44e+04 | 2.14e+03 | 2.1e−03 | **1.5e−04** |

- 生産 Euler (case/41 モデル) で `hoopAreaFromClosure: 1` は場を 9e−5 しか動かさず収束も同等 (run_0021)。
- 既定 OFF はビット不変 (run_0022, 差 1.8e−6 = node の atomicAdd 非決定性レベル)。
- **NS (y+1.5) の `nodeValueAtNode` 発散は pRef でも閉性でも直らない** (step 64→65) = 別原因。

## 4. 残課題

- `pRef` は定数ゲージなので、圧力比の大きいノズルでは全域最適にはならない (誤差は $|p-p_{\rm ref}|$ 比例)。
  区間別ゲージや $p_\infty(x)$ ゲージは将来課題。
- `axisymMethod: 1` (SU2 流 planar+1/y ソース) は hoop 源を持たないので本問題自体が無いが、別の未解決
  (喉部リミットサイクル) がある。
- 既定化は生産系列 (case/41) の再収束確認後に判断する。
- **2026-10-10 の観察** ([architecture-float-state-double-geometry](architecture-float-state-double-geometry.md) §6.5): case/45 の B0 (`hoopAreaFromClosure 0`・`pRef 0`) の一様な静止場で、**FP64 でも**軸の近く (j 2〜4) に `3.3e-3·P·A_planar` の偽の半径方向の力が出る。`hoopAreaFromClosure 1` で丸めの水準まで消える。float32 のメトリックだけの問題ではない。
  - codex の第 1 仮説 (確度 中): 変換器が折れた双対面の区間を 1 本の面ベクトルにまとめ、その重心の半径を掛ける (r̄·Σ S_k ≠ Σ r_k S_k)。§1 の「多角形では厳密」は、区間ごとに半径を掛けた場合に限る。
  - B0 の残差の床への影響は、[time_integration-line-implicit-speed](time_integration-line-implicit-speed.md) §6.23 の短い比較で確かめる。

### 4.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | case/45 の B0 での効き | line-implicit-speed §6.23 (`hoopAreaFromClosure` 0/1 の 2000 step、事前登録済み) の結果で、床への影響と解への影響のフラグを見る **結果 (2026-10-10、line-implicit-speed §6.24)**: 判定は棄却 (Ω の Σ\|res_ro\| の比 0.94〜1.00)。j 2〜4 の局所の床だけ 0.34〜0.41 倍。Q_w に −0.11 % のフラグ (過渡の途中)。同じ状態の ΔR は式どおり (差 2.3e-8)。FP64 の `A_closure` は生の double の幾何からの Σ r_f S_f と一致 (float 化 plan §6.5 の 4b) なので、j 2〜4 の欠損は `A_planar` 側と Σ r_f S_f の差である | O |
| 2 | FP64 でも残る閉性の欠損の原因 | 実メッシュの双対区間から Σ r_k S_k を作り直し、今の r̄·Σ S_k との差が j 2〜4 の欠損の符号・大きさを説明するか確かめる (codex の反証条件)。説明しなければ面積・向き・まとめ方を点検する | F |
| 3 | 既定化の判断 | 1・2 の後。生産系列 (case/41・case/45) の再収束の確認を含む | F |

### 4.2 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan 免除 | 2026-08-16 | — | — | 2026-09-09 以前に起票・実装済み (AGENTS.md の codex レビュー節の経過措置) |

## 5. 変更ログ

- `2026-08-16` — 起案 + 実装 + case/43 で検証 (上表)。`pRef` 整合はバグ修正、`hoopAreaFromClosure` は既定 0 の新フラグ。
- `2026-10-10` — §4 に case/45 の観察と codex の第 1 仮説を追記し、残作業表 (§4.1) とレビュー記録 (§4.2) を足した。
- `2026-10-10` — §4.1 #1 の結果を記入 (line-implicit-speed §6.24: 床は短い時間では下がらない)。
