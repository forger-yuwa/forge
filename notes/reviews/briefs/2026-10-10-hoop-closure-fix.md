# 諮問ブリーフ: 軸対称の hoop の閉性の欠損の直し方 (2026-10-10)

- 依頼者: 主セッション (AGENTS.md のエスカレーション条件 1・6: plan の方針を決める、数値の振る舞いを変える前)
- plan: `plans/active/axisymmetric-freestream-hoop-gauge.md` (§1・§2・§4.1 #2〜#4)。関連: `plans/active/time_integration-line-implicit-speed.md` §6.23・§6.24、`plans/active/architecture-float-state-double-geometry.md` §6.5
- ユーザの決定: 「次に取り組む課題」(2026-10-10)。
- コード: ブランチ `feature/nozzle-wall-fit-and-pipeline` HEAD (`cb6a0d8b`)。

## 0. 問い

1. 直し方を (a) と (b) のどちらにすべきか。ほかの案はあるか。
   - (a) 変換器で、区間ごとに半径を掛けた面ベクトル Σ_k r_k S_k (x・y 両成分) を作り、ソルバの r̄·S (面ベクトル × 面重心の半径) を置き換える。
   - (b) `hoopAreaFromClosure 1` を既定にする (hoop ソースの面積を、流束と同じ r̄·S の和にする)。
2. (a) のとき、流束の r の重みを区間ごとの値に替えると、何が変わるか (保存・面の流束の向き・再構成・ラインや壁の扱い・3D/平面のケースへの影響)。どう事前登録して確かめるべきか。
3. (b) のとき、残る「流束側の近似」(r̄·ΣS_k ≠ Σ r_k S_k) は、解の精度 (軸のマッハ数など) にどれくらい効く見込みか。測る必要があるか。

## 1. 観測事実

- 一様な静止場 (P 1e5 Pa・T 300 K・u 0、`pRef 0`) の case/45 の FP64 で、`hoopAreaFromClosure 0` (既定) だと、軸の近く (j 2〜4) に |res_roUy|/(P·A_planar) 最大 3.3e-3 の偽の半径方向の力が出る。10 step で偽の速度 0.257 m/s (`run_0401〜0403`)。
- `hoopAreaFromClosure 1` で丸めの水準 (8.7e-10) に消え、偽の速度は 6.4e-4 m/s (`run_0408`)。float では 4e-3 m/s (`run_0405〜0407`)。
- **原因の確認** (`case/45.isobutane_m6_d155/dual_segment_closure.py`): 実際の格子 (`run_0252` の `nozzle.h5`、FP64 の変換器) で双対の面を区間ごとに作り直した。
  - 変換器は、内部の双対面を「エッジの中点 → 隣の primal セルの重心 (頂点平均)」の区間の和で作り、重心は区間の中点の長さ加重 (`gmshReader.hpp:1232-1244, 1475-1531`、境界の半割面は `:1641-1699`)。
  - ソルバは r_f = max(pcy, 1e-20) を面ベクトルに掛け、A_planar は双対の面積 (`variables.cpp:751-780`)。
  - 区間ごとの Σ r_k S_k,y は双対の面積と ≤ 1e-12·A で一致する (j 2〜8)。今の r̄·ΣS との差は、j 2〜4 の欠損を全量・符号・位置とも説明する (相関 1.000、残り ≤ 2e-11)。
  - 欠損はほぼ全部、軸方向の辺 (同じ j) の双対面 (法線 ~x) から来る。軸の近くの行は半径方向に長い (x/r_t −3 で Δr ≈ 0.2 r_t、Δx ≈ 0.0056 r_t) ので、上下 2 区間の中点の半径の差が大きい。j 6〜8 にも 4〜7e-4。収縮部では def < 0 (軸へ向かう偽の力)。
- **B0 の床への効き** (line-implicit-speed §6.24、事前登録の A/B): `hoopAreaFromClosure` 1 にしても、Ω = x/r_t ∈ [−5, 0) の Σ|res_ro| は 0.94〜1.00 倍 (棄却)。j 2〜4 の局所の床だけ 0.34〜0.41 倍。2000 step の時点で Q_w −0.11 % (フラグ、過渡の途中)、θ_r +0.002〜0.003 %。同じ状態の ΔR は式どおり。
- 軸の近くの格子は粗い (j 1〜4 で dr 0.10〜0.69 r_t)。別のセッションが格子を見直している (ユーザの指摘)。

## 2. 期待値と出典

- 多角形の CV では、区間ごとの r の重みで Σ r_k S_k,y = A_planar が厳密に成り立つ (§1、確認済み)。今の離散化は面の重心の半径で 1 点近似している。
- 既存の `hoopAreaFromClosure` は 2026-08-16 に入れた既定 0 のフラグ (plan §2・§3)。生産 Euler (case/41) で場を 9e-5 しか動かさなかった。

## 3. 呼び出し側の案 (検証していない)

- (a) が根本。流束の r の重みを区間の積分にすると、一様な圧力の釣り合いが hoop ソースの面積 A_planar のままで厳密になり、流束の幾何の精度も上がる。変換器が r の重みの面ベクトル (x・y) を新しいデータセットとして書き、ソルバはそれがあれば r̄·S の代わりに使う。
- (b) は 1 行で済むが、流束の側の近似は残り、ソースだけを合わせる。

## 4. 読んでよいファイル

- `plans/active/axisymmetric-freestream-hoop-gauge.md`、`plans/active/time_integration-line-implicit-speed.md` (§6.23・§6.24)
- `case/45.isobutane_m6_d155/dual_segment_closure.py`
- `solver_density_cuda/mesh/gmshReader.hpp` (双対の構築、grep で該当箇所)、`solver_density_cuda/variables.cpp` (`:590-680`、r の重み・A_planar・closure)、`solver_density_cuda/cuda_forge/axisymmetricSource_d.cu`
- 巨大なファイル (`*.h5`・ログ) は読まないこと。
