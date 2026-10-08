# SERN カウル後縁の下流の格子を後縁の接線から出す (te_wake_blend_H) の生産受入れ

## メタ

- **area**: `tooling` (SERN メッシャ・設計チェーン)
- **status**: `draft`
- **related_docs**:
  - [`methods/convection/implementation.md`](../../methods/convection/implementation.md) (後縁の冷点と $w$ 処置の節)
- **related_plans**:
  - 起点: [`convection-zero-thickness-edge-reconstruction.md`](convection-zero-thickness-edge-reconstruction.md) §4.1・§4.2・§6.0 (格子の判別 A/B で結果 A、2026-10-08)
  - 親: [`tooling-nozzle-sern-chain.md`](tooling-nozzle-sern-chain.md) (R7b の再開条件)、[`tooling-nozzle-sern-3d.md`](tooling-nozzle-sern-3d.md) (旧メッシャ・R5h)
- **created**: `2026-10-08`
- **owner**: `CFD Dev (Claude セッション: SERN 3D)`

## 1. 目的

旧メッシャはカウル後縁の下流の中間線を MOC の自由境界の終端角 θ_b (−40.2°) で引くので、境界層の格子線が後縁で約 36° 折れ、
m10_on で後縁の 1 節点が温度床 50 K に張り付いて評価が成立しない。中間線を後縁の接線から出して 1.0 H で旧線へ戻す局所変形
(`mesh3d.te_wake_blend_H: 1.0`) は、g3 で冷点を消した (起点 plan §6.0 の結果 A)。これを**生産の評価系列に入れてよい状態**にする
(格子感度・来歴・旧評価の扱い・2D)。

## 2. スコープ

- **やる**: g3 → g4 の格子感度、毎更新の床補正の事象カウンタ、評価系列への格子の識別 (曲線版・実格子署名・双対幾何・変換器) の結び付け、
  2D メッシャへの同じ方式の実装と独立した検証、m6_on・m4_off の再評価と旧評価・設計差の持ち越しの判定、生産 YAML。
- **やらない**: $w$ 処置 (起点 plan で opt-in の予備として保持、§6.1 は保留)、遷移長・曲線の再設計 (結果を見て調整しない)、
  有限厚の全ヘキサ模型 ([`tooling-sern-mesh-blocking.md`](tooling-sern-mesh-blocking.md))、ランプ後端の x 間隔の跳び (親 plan R14)。

## 3. 関連 docs と前提

- g3 の判別 A/B: A′ = `case/46.sern_design/run_1078_tewake_A0_m10` (症状を再現、GATES FAIL FLOOR_STUCK)、B = `run_1079_tewake_B10_m10`
  (床 0・R_TE 最低 227.4 K・GATES PASS・4 力量 STEADY)。写し `notes/investigations/2026-10-08-te-wake-ab/`。投入条件と来歴は B の `TE_WAKE_ADMISSION.txt`・`TE_WAKE_PROVENANCE*.json`。
- 道具: メッシャ `mesh_sern3d._te_wake_midline` (3 次 Hermite、曲線版 `TE_WAKE_CURVE_VERSION`)、格子の判定 `case/46.sern_design/diag/te_wake_grid_check.py --admission`、
  初期場 `solver_density_cuda/tools/restart_field_deformed.py`、温度の監視 `case/46.sern_design/diag/te_monitor.py`、壁終端の点検 `solver_density_cuda/tools/check_wall_end_geometry.py`。
- 判断の記録: codex diagnose `notes/reviews/2026-10-08-te-wake-ab-result-diagnose.md` (結果 A を採用・生産採用の完了は認定しない・既定は 0 のまま・受入れは別 plan・力の差を冷点除去の物理効果と認定しない)。

## 4. 設計方針

**判断: 2026-10-08 codex (diagnose) `2026-10-08-te-wake-ab-result-diagnose.md`、Major 4・Minor 2 を全件採用。**

- **既定は 0 のまま**。生産候補の YAML に `mesh3d.te_wake_blend_H: 1.0` を明示する。メッシャの既定を 1.0 に変えない (評価系列の識別が格子の変更を区別できるまで)。
- **評価系列の識別**: 曲線版・$L_b$・実格子署名 (`mark_zero_thickness_edges.mesh_signature`、座標・接続) に加え、**双対幾何 (`PLANES/surfVect`・`CELLS/volume`) のハッシュと変換器の版**を
  `stage_key`・`metrics.json`・設計 DB (driver_sern の台帳) に入れる。座標・接続の署名だけでは変換器の修正 (`ef74e342`) のような双対幾何の変更を検出できない (run_1055 と A′ は署名が同じで双対の閉性が FAIL/PASS)。
  属性の転記でなく実入力格子から再計算する。旧データは削除せず保持し、修正後の学習系列への混入を止める (識別が違えば別系列)。
- **床の判定**: 生産受入れでは**毎更新の床補正の事象カウンタ** (起点 plan §5 の 3 の定義: 更新後の保存量の内部エネルギーが $e(\rho, Y, T_{min})$ を下回り最終状態を床へ補正した事象を、実節点・更新単位で数える。反復中のクランプは数えない) を使う。
  保存標本の「床近傍 0 節点」を「全更新の床補正 0 件」と読み替えない。計測の追加でバイナリが変わったら、比較する g3 側も同じバイナリで回し直す。
- **力の差の扱い**: g3 の B − A′ (末尾 20 標本の平均差: C_T +6.54e-7、C_T_with_shear −7.90e-5、C_L +2.08e-4、C_M −4.53e-3) は離散化の変更への応答であって、冷点除去の物理効果とは認定しない。
  採用判断に原因別の分解は要らない。旧 m10_on は床到達で評価不成立のまま。m6_on・m4_off の旧評価の持ち越しは §6 #6・#7 で判定する (格子感度の許容を持ち越し基準に使わない)。
- **監視領域**: 起点 plan §6.2 の $R_{TE}$・$R_{SE}$ と復帰区間 RET (x ∈ [L_cowl, L_cowl + L_b]) と全域。マスク外周 $R_M$ は $w$ 固有なので外す。
- **2D**: 2D メッシャ (`mesh_sern.py:256`) も同じ作り。同じ方式を実装し、2D で独立に検証する (3D の成功で 2D を受理しない)。

## 5. 実装ステップ

1. g4 の格子感度 (§6 #4): `te_wake_blend_H: 1.0`・$w$ 無効の m10_on を g4 で生成し、投入条件 (`te_wake_grid_check.py --admission`、g4 は A′ に相当する g4 の 0 格子を基準に) を通し、`run_1079` の最終場から cross-mesh restart。
2. 毎更新の床補正の事象カウンタ (§6 #1): `solver_density_cuda/` の EOS・保存量の更新後の補正箇所。無効化不要の出力専用 (数値はビット不変を確認)。
3. 評価系列の識別 (§4): `stage_manifest.py`・`runner_sern3d.py`・`runner_sern.py`・`driver_sern.py`・`te_wake_grid_check.py` (双対幾何のハッシュ)。
4. 2D メッシャへの実装と 2D の判別 A/B (§6 #10)。
5. m6_on・m4_off の修正格子での再評価と旧評価・設計差の判定 (§6 #6〜#8)。
6. 生産 YAML (`problem_3d_prod_*` に `te_wake_blend_H: 1.0`) と親 plan R7b の再開。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | codex plan レビュー | 本 plan の §4・§6 | F (採否) |
| 2 | g4 の格子感度 (§6 #4) | AWS (起動はユーザ確認)。g4 の 0 格子と 1.0 格子を現行コードで生成、`--admission` PASS、`run_1079` の最終場から `interp_field` (cross-mesh)、20000 step (+20000 を 1 回まで)。合格は §6 #4 | O (解釈は F) |
| 3 | 毎更新の床補正のカウンタ | 起点 plan §5 の 3 の定義。出力専用、無効時の残差・保存量がビット一致 (case/66 の手順)。単体試験 (正常・床未満・下限近傍) | O (cuda_forge を触るので実装前に F) |
| 4 | 評価系列の識別 | 双対幾何のハッシュ・変換器の版・曲線版・$L_b$ を stage_key・metrics・設計 DB へ。旧系列と区間・学習が分かれる試験 | O |
| 5 | 2D メッシャ | `mesh_sern.py` に同じオプション (既定 0 でビット一致)、2D の判別 A/B (2D m10_on/m6_on) | O (解釈は F) |
| 6 | m6_on・m4_off の再評価 | 修正格子 (g3、1.0) で 20000 step、§6 #6・#7・#8 | O (解釈は F) |
| 7 | 生産 YAML と R7b 再開 | §6 の受入れが揃ってから | F |

## 6. 検証

**事前登録 (2026-10-08、codex diagnose `2026-10-08-te-wake-ab-result-diagnose.md` と起点 plan §6.2 から格子修正に必要な項目を引き継ぐ)**。各 run は 20000 step、未定常なら +20000 を 1 回。採否は末尾 10000 step (500 step 間隔の 20 標本)。
全残差の `check_convergence.py` と対象量の `check_quasisteady.py` の VERDICT を区間付きで残す (プラトー受理の方針は維持、収束とは呼ばない)。$a$ = 判定区間の $\max|q - \bar q|$、$d$ = 前後 10000 step 窓の平均の差の絶対値。

| # | 量 | 合格 | 判定不能 |
| --- | --- | --- | --- |
| 1 | 床補正の事象 (毎更新のカウンタ) | 判定区間の全実節点・全更新で 0 件 | カウンタの欠落・切り詰め |
| 2 | 端の近くの低温 ($R_{TE}$・$R_{SE}$・RET・全域) | $T < 150$ K の節点なし、各領域の最低温度が STEADY、$T < 200$ K の体積と位置を記録 | 保存場の欠落 |
| 3 | 窓 (4 力量) | C_T・C_T_with_shear・C_L: $\max(a, d) \le 5\times10^{-5}$、C_M: $\le 5\times10^{-4}$、4 量 STEADY | +20000 後も未達 |
| 4 | **g3 → g4 の格子感度** (同じ形状・作動点・1.0 H 曲線、BC・数値設定・物性・介入の有無を固定) | 両側で #1〜#3 が合格し、4 力量の平均差が C_T・C_T_with_shear・C_L ≤ 0.002、C_M ≤ 0.05。**既存の g4 仕様は複数の格子パラメータを変えるので「複合的な格子感度」と呼び、収束次数は推定しない** | g4 の投入条件の不合格・転送の不備・延長後も未定常 |
| 6 | 旧評価の持ち越し (m6_on・m4_off) | $D_q = |\text{平均差}| + \max(a_A, d_A, 10^{-6}) + \max(a_B, d_B, 10^{-6})$ が C_T・C_T_with_shear・C_L ≤ 5e-5、C_M ≤ 5e-4 なら旧評価を持ち越す。超えたら旧評価と学習データを失効させ取り直す | 片側が未定常 |
| 7 | 設計差 ($L_{sw}$ 0.8/1.0) | 起点 plan §6.2 #7 と同じ式・予算 | 同上 |
| 8 | 設計への利得 | 起点 plan §6.2 #8 と同じ | — |
| 10 | 2D | 2D の判別 A/B (0 と 1.0) で #1〜#3 を 2D の領域に読み替えて判定。3D の結果から 2D を受理しない | 2D で症状が再現しない (対照の不成立) |

- **#4 の結果の読み方**: g4 で冷点なし・感度が許容内 →「改善は g3 だけの現象」をこの比較範囲で退ける。g4 で冷点が再発・移転、または感度超過 → 生産候補としての十分性を退ける (g3 で観測した冷点の消失自体は撤回しない)。
- やらないこと: 既定値の一括変更、結果を見て遷移長を調整すること、$w$ との同時変更、3D の成功による 2D の自動受入れ。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| diagnose (起票の判断) | 2026-10-08 | [2026-10-08-te-wake-ab-result-diagnose.md](../../notes/reviews/2026-10-08-te-wake-ab-result-diagnose.md) | 結果 A を採用、次は g4 の格子感度 (C0/M4/m2) | **全件採用** → §4・§6。M1 (生産採用の完了は認定しない・床カウンタの代用を全更新 0 と読み替えない) → §4・§6 #1、M2 (既定は 0・評価系列の識別) → §4・§5.1 #4、M3 (受入れは別 plan・$w$ は予備・§6.1 保留) → 本 plan の起票、M4 (力の差を物理効果と認定しない・持ち越しは #6・#7) → §4・§6、m1 (A′ の rms_roK は 0.5 桁低下中) → 起点 plan の記録を訂正、m2 (力は末尾窓平均で比べる・$a$ = max\|q − 平均\|) → §4・§6 |

## 7. 影響範囲

SERN 3D (と 2D) の生産評価系列。メッシャの既定は 0 のままなので、`te_wake_blend_H` を書かない既存の YAML・run には影響しない。

## 8. 完了条件

§6 の #1〜#4・#6・#7・#10 が合格 (または #6・#7 の判定に従って取り直しの段取りが親 plan に入る) し、codex result レビューを経ること。

## 未確定事項

- g4 の 0 格子 (A′ 相当) を g4 の判定に使うか (投入条件の基準として作るが、g4 の 0 格子で冷点が再現するかは測らない予定)。
- 毎更新のカウンタを入れたバイナリで g3 の B を回し直すか (計測の追加でバイナリが変わる場合)。

## 9. 変更ログ

- `2026-10-08` — 起票。起点 plan §6.0 の格子の判別 A/B で結果 A (g3、run_1078/1079)。codex diagnose で生産受入れを本 plan に分け、次の手を g4 の格子感度に。
