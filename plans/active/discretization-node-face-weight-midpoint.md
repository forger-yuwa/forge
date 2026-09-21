# node 内部双対面の面補間重み `fx` を中点 (0.5) に固定する

## メタ

- **area**: `discretization`
- **status**: `in_progress`
- **related_docs**:
  - [`methods/discretization.md`](../../methods/discretization.md) (「node 内部双対面の面補間係数」の節)
  - [`methods/boundary.md`](../../methods/boundary.md) (壁熱流束 `iface_q_eff` の定義)
- **related_plans**: [`boundary-conjugate-heat-transfer.md`](boundary-conjugate-heat-transfer.md) §5.1 #57–#61 (発見の経緯)、
  [`architecture-node-option-consolidation.md`](architecture-node-option-consolidation.md) (`nodeMidpointFx` を撤去した plan。その前提を本 plan が覆す)
- **created**: `2026-09-21`
- **owner**: Claude (ユーザ指示「熱流束のガタつきに対処する」)

## 1. 目的

冷却翼 (case/53 C3X) の壁熱流束 `iface_q_eff` に、SU2 には無い**約 25 節点周期・空間固定のうねり** (±1.5–2 kW/m²) が乗る。
真因は面補間重み `fx` で、node の内部双対面でも幾何から計算しており、(i) 式が回転不変でなく、(ii) 正しい射影で計算しても
高アスペクト比の曲面壁層では 0.5 から大きく外れる。node では辺の中点で値を取る (`fx = 0.5`) のが本来の離散化であり、これに固定する。

## 2. スコープ

- **やる**: node の内部双対面 (`ip < nNormalPlanes`) で `fx = 0.5` を固定スキームにする。`methods/discretization.md` の該当節を現状に合わせる。
  cell 側の `fx` の式 (回転不変でない) の修正要否を記録する。
- **やらない**: 境界半割面 (`ip >= nNormalPlanes`) の `fx`。cell モードの挙動変更 (cell は使わない方針。式の欠陥は記録のみ)。
  壁熱流束診断 `iface_q_eff` の蓄積項補正 (CHT plan §4.3 / §5.1 #61 で扱う)。

## 3. 関連 docs と前提

- `fx` は `calcStructualVariables_d.cu` で 1 回だけ計算される。**値を実際に読むカーネル** (codex plan レビュー M3 で補完):
  粘性流束 `viscousFlux_d.cu` (面値 $U_f$, $\mu_f$, $k_f$ と面勾配の補間)、SST の $k$/$\omega$ 面補間と GG 勾配 `ransTransport_d.cu`、
  スカラー拡散 `scalarTransport_d.cu` (面係数と陰的対角)、多成分拡散 `speciesTransport_d.cu`、受動スカラー/FCT
  (`passiveKernels_d.cuh`, `passiveFct_d.cuh`)、CFL 評価 `setDT_d.cu`、`interpVelocity_c2p`。
  `limiter_d.cu` と `ducrosSensor_d.cu` は引数に取るだけで**読まない**。
- 2026-06-14 に opt-in `nodeMidpointFx` として同じ固定を入れ、2026-08-16 の node オプション整理で
  「値=ノード座標では幾何 fx が自動的に中点相当」として**撤去**した (カーネルの `nodeMode==1` 分岐は残っている)。
- SU2 は辺の両端の算術平均 (0.5)。

## 4. 設計方針

### 4.1 実測した欠陥

`case/53.c3x_vane_cht/run_0121_resid_split` (C3X、一様壁温 566 K、第一層 2 µm・壁沿い 0.73 mm = AR 約 350) で、
壁半 CV に入るエネルギー残差を対流・熱伝導・粘性仕事に分けた (`interfaceDiag` の `ifaceRconv`/`ifaceRpre`、
`FORGE_WI_FORCE_DIAG=1` の `wi_eheat`/`wi_ework`。和は 0.01 W/m² で閉じる)。負圧面 $s/S$ 0.45–0.95 のうねり rms:

| 成分 | 平均 [W/m²] | うねり rms [W/m²] |
| --- | --- | --- |
| 対流 | −0 | 70 |
| 内部面の熱伝導 − 壁面流束 | 761 | 109 |
| **内部面の粘性仕事 $\tau\cdot U_f$** | 4106 | **571** |

粘性仕事の面速度は $U_f=(1-f)U_1$ (壁節点は $U=0$)。`fx` の式は

$$d_0=\sqrt{\textstyle\sum_i (n_i\,\Delta_{0,i})^2},\qquad f=d_1/(d_0+d_1),\qquad \Delta_0=x_{pc}-x_0$$

で、**法線への射影 $|n\cdot\Delta|$ ではなく成分ごとの積のノルム**である。$\Delta$ の接線成分 $a$ は射影なら消えるが、
この式では $\sqrt2\,|n_xn_y|\,a$ として残る。面重心の接線ずれは壁沿い間隔の不均一から µm 級で、$d_1/2$ = 1 µm と同程度。
この式を後処理で再現すると**カーネルの粘性仕事と 0.0 W/m² で一致**し、射影で計算した重みでは相関 0.14 しかない。
壁側重みは負圧面後半で 0.42–0.75、翼全周で **0.07–0.96**。

**式を射影に直すだけでは足りない**。射影でも翼全周で 0.03–1.00 に散る。曲率 $\kappa$ の壁では面重心が弦のたるみ
$\kappa\Delta s^2/8$ だけ沈み (C3X で 0.1–1 µm)、これが $d_1/2$ と同程度になるためである。「幾何 fx は中点相当」は
等方的なセルでしか成り立たず、壁解像の高 AR 層で破れる。`methods/discretization.md` は式を「法線方向に射影した距離比」と
記述しており、**実装が文書と違っていた**。

### 4.2 対策

node の内部双対面で `fx = 0.5` に固定する (`calcStructualVariables_d_wrapper` で `nodeMode = (discretization=="node")`)。
オプションにはしない (node は固定スキームという整理方針に従う)。SU2 と同じ辺中点評価。
内部面流束は両 CV に同じ値を逆符号で足すので**保存性は `fx` に依らない**。

**精度について主張する範囲** (codex M2 採用): 辺中点の算術平均は「線形場で厳密、滑らかな場の**値補間**として 2 次」であり、
それ以上は主張しない。伸長格子上の 2 点差分の拡散演算子は `fx` と無関係に点値に対する局所打切り誤差が 1 次
(節点 $(-a,0,b)$, $T=x^3$ で $L_hT=2(b-a)$) で、双対面のパッチ重心と辺中点も一般に一致しない。
**解の収束次数は格子系列で測る** (§5.1 #7)。

### 4.3 試行結果 (環境変数 `FORGE_NODE_FX_HALF=1` での A/B)

`run_0108_cf0_su2turb_long` の最終場から同一設定で 40000 step。対照 `run_0123_fxgeom_ctrl` / 試行 `run_0122_fxhalf`。
どちらも `check_convergence.py` は **NOT CONVERGED (stalled/plateau)** — 後縁の渦放出で翼の全 run が同じ。NaN なし。

| 量 | 幾何 fx (対照) | fx = 0.5 (試行) | SU2 (同一メッシュ) |
| --- | --- | --- | --- |
| `q_eff` うねり rms, $s/S$ 0.45–0.95 [kW/m²] | 1.22 | **0.46** | 1.17 (壁面勾配流束。後縁側の大きな構造を含む) |
| 壁面勾配流束のうねり rms [kW/m²] | 0.98 | 0.62 | — |
| $h$ の偏差 PS / 負圧面層流域 / 遷移後 / 全体 [%] | +19.0 / +43.3 / +8.1 / +22.0 | +18.8 / +46.0 / +8.0 / +22.6 | +21.3 / +45.4 / +6.7 / +23.1 |

25 節点周期のうねりは消え、$h$ の偏差は 1 pt 以内 (負圧面層流域のみ +2.7 pt で SU2 に近づく)。

## 5. 実装ステップ

1. `methods/discretization.md` の該当節を更新 (式の実態・高 AR 曲面壁での破綻・固定に戻す判断)。
2. `calcStructualVariables_d.cu`: 環境変数を外し `nodeMode = (cfg.discretization == "node")`。
3. 検証 (§6) → codex result レビュー → `accepted/`。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 |
| --- | --- | --- |
| 1 | ~~codex plan レビュー~~ | 済 (2026-09-21、GO-with-changes C0/M4。§6.1)。4 件とも採用 |
| 2 | 回帰範囲の補完 (codex M3) | `procedures/verification/README.md` の node 一覧に沿って **`case/09` (周期・受動スカラー)** と **`case/44` (軸対称・多成分 TP・凝縮)** を `FORGE_NODE_FX_HALF` の A/B で回す。済: `case/48` 平板 (`run_0020`/`run_0021`: $q$ 最大 0.0031 %・$\tau_w$ 0.0020 % = 不変)、`case/52` スラブ (`run_0006`/`run_0007`: 5e−7) |
| 3 | 合格条件の固定 (codex M4) | §6 の通り。A/B は同じ熱流束定義で比べる / 壁熱流束の領域別偏差・うねり rms の**時系列**を `check_quasisteady.py --series-csv` で判定 / 未収束の run は機構診断として扱う |
| 4 | 恒久実装 | #2・#3 を満たしてから。`nodeMode = (discretization=="node")`、環境変数 `FORGE_NODE_FX_HALF` は撤去 |
| 5 | 翼の生産 run の更新 | 恒久実装と result レビューの**後**に報告の数字を差し替える。それまでは `FORGE_NODE_FX_HALF=1` の試行 run (`case/53` `run_0125`、`case/54` `run_0024`) を「試行」と明記して併記する |
| 6 | cell の `fx` の式 | 回転不変でない式のまま (同じ幾何を 45° 回すと重みが 0.8 → 0.637 に変わる: codex の代数チェック)。cell は使わない方針なので**修正せず記録のみ** |
| 7 | 格子系列での次数測定 (codex M2) | 成長率 1.1/1.2 の伸長格子・曲面高 AR 格子・回転した同一格子で製造解 (または解析解のある層流) を 3 格子。温度・速度・壁熱流束の誤差次数を測り、2 次を主張する範囲は $p\ge1.8$ を要求。細分化時の成長率の扱いを明記する |
| 8 | codex result レビュー | `done` にする前 |

## 6. 検証

- **ビルド**: `tools/build_native_wsl.sh`。
- **検証ケース**: (a) `case/53.c3x_vane_cht` 一様壁温 (§4.3 の A/B を恒久実装で再現)、(b) node の平板 (`procedures/verification/README.md` の
  node 乱流平板。$C_f$ が SU2/理論比で悪化しないこと)、(c) `case/52.conjugate_slab` (1 次元。`fx` に依らず結果不変のはず)。
- **判定基準**: (a) `q_eff` うねり rms が対照の 1/2 以下、$h$ の領域別偏差の変化 3 pt 以内。(b) $C_f$ の変化 1 % 以内。(c) 壁熱流束の変化 0.01 % 以内。
  いずれも NaN なし。
- **比較の作法** (codex M4 採用): A/B は**同じ熱流束定義**で比べる (`iface_q_eff_raw` と `iface_q_eff` を両方保存)。
  §4.3 の表は蓄積項補正を入れる前のバイナリで両 run とも旧定義。定常回帰 (b)(c) は `check_convergence.py` の `PASS`
  (収束場からの継続は `--from-floor`) を要求する。(a) の翼は全 run が `NOT CONVERGED (stalled/plateau)` なので
  **機構診断**として扱い、領域別偏差とうねり rms の時系列を `check_quasisteady.py --series-csv` で `STEADY` と確認した量だけを引用する。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan | `2026-09-21` | [`notes/reviews/2026-09-21-discretization-node-face-weight-midpoint-plan.md`](../../notes/reviews/2026-09-21-discretization-node-face-weight-midpoint-plan.md) | GO-with-changes, C0/M4/m0 | **全件採用**。M1 (蓄積項の係数は $H_w$ でなく $e_w$) → CHT plan §5.1 #61 を訂正し `conjugateWall.cpp` を修正。M2 (値補間の 2 次と演算子の次数の混同) → §4.2 を限定し §5.1 #7。M3 (`fx` の実使用一覧と回帰範囲) → §3・§5.1 #2。M4 (未収束同士の比較・定義変更の混入) → §6・§5.1 #3、生産値の差し替えは result レビュー後 (#5) |

## 7. 影響範囲

- `solver_density_cuda/cuda_forge/calcStructualVariables_d.cu` (1 行)。**node の全 run の数値が変わる** (等方セルでは丸め程度、壁解像の高 AR 曲面層では上記の通り)。
- `methods/discretization.md`、`plans/active/architecture-node-option-consolidation.md` の `nodeMidpointFx` 行 (前提の訂正)。
- 報告 (Cooled Vane CHT Validation) の壁熱流束の図と数字。

## 8. 未確定事項

- 境界半割面の `fx` は現状のまま (壁半割面の粘性流束は `fx` を読まない)。入口・出口面で問題が出るかは未確認。

## 9. 変更ログ

- `2026-09-21` codex plan レビュー (GO-with-changes)。4 件採用。平板・スラブの回帰は不変を確認。
- `2026-09-21` 起票。真因の特定 (`run_0121_resid_split`) と環境変数での A/B (`run_0122`/`run_0123`) まで。
