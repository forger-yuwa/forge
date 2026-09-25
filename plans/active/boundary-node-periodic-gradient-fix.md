# node 周期の継ぎ目勾配の修正 (LSQ の 2 重計上、SST k/ω 勾配の未合算)

## メタ

- **area**: `boundary`
- **status**: `draft`
- **related_docs**:
  - [`methods/gradient.md`](../../methods/gradient.md) (「境界寄与」node × 周期の継ぎ目)
  - [`methods/discretization.md`](../../methods/discretization.md) §7.3 (LSQ)、§4.5 (周期の DOF 同一視)
- **related_plans**:
  - [`boundary-node-rotational-periodic.md`](boundary-node-rotational-periodic.md) (回転周期。本 plan の完了を前提にする。本 plan は同 plan の §4.0 / §5.1 #0 を切り出したもの)
  - [`discretization-median-dual-3d.md`](discretization-median-dual-3d.md) §4.5 (node 周期の骨格)
- **created**: `2026-09-26`
- **owner**: `CFD Dev`

## 1. 目的

node の周期境界の継ぎ目にある既存の欠陥 2 件を直す。どちらも回転と無関係に、**今の並進周期の node 計算すべてに効く**。

1. **LSQ 勾配の 2 重計上**: node の勾配は `gradLSQ: 2` 固定で、継ぎ目で割れた各部分 CV が片側の隣接だけで完全な LSQ 勾配を解き、
   `periodicGradientGather` が和を取る。**線形場で継ぎ目の勾配が正確に 2 倍** (2026-09-26 実測、case/09 TGV 32³、`Ux = 10+y` で x 継ぎ目 1458 点の
   $\partial U_x/\partial y$ = 2.000000、内部 1.000000。`case/09.Taylor-Green/_g0_lsq_seam/`)。
2. **SST の $k,\omega$ 勾配が合算されていない**: `periodicGradientGather` (`main.cpp:1459`) の後で `ransGradient` (`main.cpp:1481`) がゼロから作り直し
   (`ransTransport_d.cu:211`)、その後に gather が無い。F1 (`ransSource_d.cu:344`) は片側の勾配を読む。`ransGradient` は周期半割面も除外していない。

**スカラー全体の点検 (2026-09-26、ユーザ質問「その他スカラーは大丈夫?」)** — main で勾配を作る 4 か所 (`calcGradient`・`speciesGradient`・
`passiveGradient`・`ransGradient`) と継ぎ目の gather の順序を確認した:

| 量 | 勾配 | 継ぎ目 | 判定 |
| --- | --- | --- | --- |
| 速度・密度・圧力・温度 | LSQ | gather は和 | **欠陥 1 (2 倍)** |
| $k,\omega$ | Green–Gauss (`ransGradient`) | gather の**後**に作り直され、周期半割面も除外していない | **欠陥 2 (片側)** |
| 化学種 | Green–Gauss (`speciesFaceReconstruction ≥ 1` のみ) | 周期半割面を除外 (`excludePeriodic`) し、合併体積で割った部分寄与の和 | 正しい |
| 受動種・凝縮モーメント | Green–Gauss (`passiveGradient`) | 同上 | 正しい |
| $\gamma,\ Re_{\theta t}$ | 拡散は 2 点差分 (勾配配列なし)。生成項は速度勾配を読む | 状態は root からミラー | 拡散は影響なし。**生成項は欠陥 1 の影響を受け、欠陥 1 の修正で直る** |

## 2. スコープ

- **やる**: 合併 stencil の LSQ 係数の事前計算、$k,\omega$ 勾配の合算位置の修正、並進周期での検証 (G0 拡張、二次場、case/39、case/09)。
- **やらない**: 回転周期 (`boundary-node-rotational-periodic` で本 plan の後に)。Green–Gauss の勾配 (化学種・受動種・凝縮モーメント。現行の和が正しい)。
  陰解法の行縮約 (別 plan)。継ぎ目の部分双対面の合併 (並進の押し出しでは部分面が同一平面で差が無い)。

## 3. 関連 docs と前提

- LSQ の事前計算: `calcGradient_d.cu:577-920` (重み $w=1/|\Delta\mathbf x|^2$ は :668、スペクトル打ち切りは :671)。
- 継ぎ目の gather: `periodicNode_d.cu:166-212`。
- 双対面は primal edge を一意化して生成し、**継ぎ目の接線方向のエッジは両側に存在する** (`gmshReader.hpp:1841`)。

## 4. 設計方針 (2026-09-26 `diagnostician`、codex plan (回転 plan の 2 回目) を採用)

### 4.1 合併 stencil の LSQ

- 各 group (root と member) で、全 member の隣接を**同じ物理隣接ごとに同定**する (隣接の `periodicRoot` と、root 系での $\Delta\mathbf x$ の一致で重複を除く。
  `periodicRoot` が同じでも方向の違う隣接は潰さない)。
- 同じ物理隣接が複数の member に現れるとき、配分係数 $\alpha_{mj}$ (同一隣接で総和 1) を付ける。
- $M_r=\sum_{m,j}\alpha_{mj}w_{mj}\,d_{mj}d_{mj}^{\mathsf T}$ ($d$ = 隣接へのベクトル、並進では root 系と同じ)。**スペクトル打ち切りは合併した $M_r$ に 1 回だけ**。
- 各部分 CV の係数 $c_{mj}=M_{r,\tau}^{+}\,\alpha_{mj}w_{mj}d_{mj}$ を焼き込む。**毎 step の gather は現行の和のまま** (部分和が合併 LSQ になる)。
- 壁∩継ぎ目では、実在する内部隣接だけを合併する (壁の疑似点は足さない。現行仕様と整合)。
- **不採用**: root の値での上書き、一律 0.5 倍 (4・8 member の角、非対称 stencil、部分的な rank 欠損を扱えない。codex の反例)。

### 4.2 SST の $k,\omega$ 勾配

`ransGradient` の直後 (`ransBlendF1` の前) に $k,\omega$ 専用の gather を置く (Green–Gauss、周期半割面を除外して積算、合併体積で割った部分寄与の和)。
早い方の gather (`main.cpp:1459`) からは $k,\omega$ を外す。

### 4.3 影響

継ぎ目の節点の勾配が変わる (意図した修正)。粘性応力・2 次再構成・リミタ・SST 生成・F1 に効く。並進の他の演算経路は不変。

## 5. 実装ステップ

1. codex plan 段。
2. G0 の拡張 (試験を先に: 2/4/8 member、非対称 stencil、壁∩継ぎ目、root 交換) と再現物の保存。
3. 合併 LSQ の実装 (`calcGradient_d.cu`、`mesh.cpp`)。
4. $k,\omega$ 勾配の gather (`main.cpp`、`periodicNode_d.cu`、`ransTransport_d.cu`)。
5. 二次場・case/39・case/09 の検証。
6. docs、codex result。

### 5.1 残作業 (優先順)

**計算資源**: 小規模はローカル、重いものは AWS (ユーザ指示: ローカルで大きな計算をかけない。AWS は他セッションと共有、手動 stop しない)。
**ユーザ指示: 実装前・区切りごとに codex に諮る**。

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | codex plan 段 | Critical/Major は `diagnostician` に諮る | F |
| 2 | G0 拡張と再現物 | 線形場を生成スクリプトで焼く (`initial` ではなく)。入力 h5 生成・集計スクリプト・修正前の勾配出力・revision を `case/09.Taylor-Green/_g0_lsq_seam/` に保存し README に登録。2/4/8 member (面・辺・角)、非対称 stencil、壁∩継ぎ目、root 交換のメッシュ | O |
| 3 | 合併 LSQ | `calcGradient_d.cu` の事前計算、`mesh.cpp` の group 情報。合格: §6 G0 | O (F レビュー、条件 6) |
| 4 | $k,\omega$ 勾配の gather | `main.cpp`、`periodicNode_d.cu`、`ransTransport_d.cu`。合格: §6 G1 | O (F レビュー) |
| 5 | 検証 | §6 G2・R1・R2。**区切りで codex** | O (結論 F) |
| 6 | docs + codex result | `methods/gradient.md` の「修正中」を外す | F |

## 6. 検証 (測る前に固定)

| # | 試験 | 合格 |
| --- | --- | --- |
| G0 | 線形場 (局所作用素試験) — 2/4/8 member (面・辺・角)、非対称 stencil、壁∩継ぎ目、root 交換 | 非退化方向の最大誤差 ≤ 1e-5、ゼロ成分は絶対誤差 ≤ 1e-5 × |真勾配|。退化方向は同じ打ち切りを施した参照解と比較 |
| G1 | $k,\omega$ 勾配: 線形の $k,\omega$ 場で継ぎ目の勾配 | 真値との最大誤差 ≤ 1e-5 (相対) |
| G2 | 二次場: CPU double の一意 stencil の LSQ を参照にした作用素試験と、格子細分 (3 水準) | 作用素試験は参照と float32 丸め内。細分では継ぎ目の勾配誤差が内部と同じ次数 (勾配は一般に O(h)。面再構成の O(h²) と混同しない) |
| R1 | case/39 周期丘: 現行レシピに設定を整備 (`wallTreatmentSST: 1`・旧 `kInf/omegaInf` を落とす) し、修正前後 | 継ぎ目の速度勾配・壁応力・$k,\omega$ 勾配・F1 の不連続が消える (継ぎ目と隣接列の差が内部の列間差と同水準)。定常結果は全残差と対象量の VERDICT を併記 |
| R2 | case/09 TGV | 一意 DOF で質量・運動量・全エネルギーの保存、KE・エントロピー履歴。修正前後の差を記録 (定常 PASS は要求しない) |

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |

### 6.2 結果

- **G0 (修正前、2026-09-26)**: 並進、case/09 TGV 32³、線形場 `Ux = 10+y` で x 継ぎ目 1458 点の $\partial U_x/\partial y$ = 2.000000 [1.999990, 2.000012]、
  内部 19683 点 1.000000。`case/09.Taylor-Green/_g0_lsq_seam/G0_translational.txt`。

## 7. 影響範囲

- `solver_density_cuda/cuda_forge/calcGradient_d.cu`、`periodicNode_d.cu`、`ransTransport_d.cu`、`main.cpp`、`mesh/mesh.cpp`。
- **既存の node 周期 run (case/39、case/09、周期翼列など) の結果は継ぎ目付近で変わる** (修正)。再現には修正前の commit のバイナリが要る。
- docs: `methods/gradient.md`。

## 8. 完了条件

- [ ] `methods/gradient.md` の「修正中」を外す
- [ ] 実装・検証完了 (§6)
- [ ] codex レビュー 2 回 (`plan` / `result`) を §6.1 に記録し、Critical / Major の採否を残作業表に反映済み
- [ ] `status` を `done` に変更し、§9 に変更ログを記載
- [ ] `plans/active/` → `plans/accepted/` へ移動、[`plans/README.md`](../README.md) を同期

## 9. 変更ログ

- `2026-09-26` — 初稿。回転周期 plan の G0 で見つかった LSQ の 2 重計上と、codex (同 plan の 2 回目) が見つけた SST $k,\omega$ 勾配の未合算を、
  回転と独立に先に直すため切り出した (`diagnostician` 判断)。
