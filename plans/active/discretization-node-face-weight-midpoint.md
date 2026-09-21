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

で、**法線への射影 $|n\cdot\Delta|$ ではなく成分ごとの積のノルム**である。$\Delta=b\,n+a\,t$ (法線成分 $b$、接線成分 $a$) と分けると、
射影は $b$ だけを返すが、この式は

$$d^2=b^2(n_x^4+n_y^4)+2ab\,n_xn_y(n_y^2-n_x^2)+2a^2n_x^2n_y^2$$

となり $a$ に依存する。**中央の交差項は辺の両端で符号が逆** ($b\to-b$) なので、面が厳密に中点にあっても $d_0\ne d_1$ になる。
$a=b$・面が中点のとき、壁が x 軸と 0°/45°/90° なら重みは 0.500、22.5° で **0.366**、67.5° で **0.634** (射影なら全角度で 0.500)。
(2026-09-21 の初版は「$\sqrt2|n_xn_y|a$ が残り 45° で最大」と書いたが、それは両端に等しく乗る対称項で、重みを 0.5 から外す主因は交差項。訂正。)
面重心の接線ずれ $a$ は壁沿い間隔の不均一から µm 級で、$b=d_1/2$ = 1 µm と同程度。
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

**精度について主張する範囲** (codex M2 採用、2026-09-22 に製造解の結果で更新): 辺中点の算術平均は「線形場で厳密、滑らかな場の
**値補間**として 2 次」。製造解 (§5.1 #7) で測れたのは: **節点値 (温度・速度) は $p$ = 2.00**、**壁の検査体積に入る熱伝導・粘性仕事は
1 次** (双対面の面積重心が弦のたるみの分だけ辺中点より壁側にあるため。誤差は たるみ / 第一層の半分 に比例)、**2 点差分の壁面値は
重みに依らず 1 次**。「`fx=0.5` で壁熱流束が 2 次になる」とは主張しない。主張するのは「旧式は細分化しても壁 CV の粘性仕事が
25–33 % ずれたまま収束しない (不整合) が、`fx=0.5` は 1 次で収束する」こと。

### 4.3 試行結果 (環境変数 `FORGE_NODE_FX_HALF=1` での A/B)

`run_0108_cf0_su2turb_long` の最終場から同一設定で 40000 step。対照 `run_0123_fxgeom_ctrl` / 試行 `run_0122_fxhalf`。
どちらも `check_convergence.py` は **NOT CONVERGED (stalled/plateau)** — 後縁の渦放出で翼の全 run が同じ。NaN なし。

| 量 | 幾何 fx (対照) | fx = 0.5 (試行) | SU2 (同一メッシュ) |
| --- | --- | --- | --- |
| `q_eff` うねり rms, $s/S$ 0.45–0.95 [kW/m²] | 1.22 | **0.46** | 1.17 (壁面勾配流束。後縁側の大きな構造を含む) |
| 壁面勾配流束のうねり rms [kW/m²] | 0.98 | 0.62 | — |
| $h$ の偏差 PS / 負圧面層流域 / 遷移後 / 全体 [%] | +19.0 / +43.3 / +8.1 / +22.0 | +18.8 / +46.0 / +8.0 / +22.6 | +21.3 / +45.4 / +6.7 / +23.1 |

この 2 本の量ごとの準定常判定 (`tools/h_series.py` → `check_quasisteady.py --series-csv`、step 5000–40000 の 8 ダンプ; codex result M4 採用):
$h$ の領域別偏差は両方とも全量 `STEADY`。うねり rms は対照が `STEADY` (1218 W/m² ± 4.8 %)、試行は **`DRIFTING`** (496 W/m² ± 17 % =
0.40–0.58 kW/m² を往復)。この 2 本は蓄積項補正**前**の熱流束定義で、衝撃足の質量残差のダンプ依存がうねり指標に混ざるためである。
同じ場を補正後の定義で 60000 step 継続した `run_0128_cf0_fx05` では **475 W/m² ± 2.1 % で `STEADY`**。
**この節の 2 本は「旧定義での有限区間の低減傾向」を示すもので、合格の根拠ではない** (合格根拠は同一定義の生産 A/B = §5.1 #3)。
主張は「うねり rms が 1.2 → 0.5 kW/m² 前後に下がる (ばらつきの外)。$h$ の偏差の変化は C3X で 2.7 pt 以内 (負圧面層流域)、
Mark II の層流域は +3.8 pt」に限る。いずれの run も残差は `NOT CONVERGED (stalled/plateau)`。

## 5. 実装ステップ

1. `methods/discretization.md` の該当節を更新 (式の実態・高 AR 曲面壁での破綻・固定に戻す判断)。
2. `calcStructualVariables_d.cu`: 環境変数を外し `nodeMode = (cfg.discretization == "node")`。
3. 検証 (§6) → codex result レビュー → `accepted/`。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 |
| --- | --- | --- |
| 1 | ~~codex plan レビュー~~ | 済 (2026-09-21、GO-with-changes C0/M4。§6.1)。4 件とも採用 |
| 2 | 回帰範囲の補完 (codex M3) | **比較は済・合格証拠は不足 (codex result M1)**。対照 2 本 + 試行 2 本を `check_field_regress.py` でノイズ床と比較。**`case/09.Taylor-Green`** (周期・受動スカラー、`run_0169_fx_ctrl_a`/`run_0170_fx_ctrl_b` 対 `run_0171_fx_half_a`/`run_0172_fx_half_b`): `VERDICT: PASS` (最大比 1.39)。**`case/44.vitiated_air_wt`** (軸対称・多成分 TP・凝縮、`run_0505`/`run_0506` 対 `run_0507`/`run_0508`、`--boundary`): `VERDICT: PASS` (最大比 1.88 < 2)。ただし `case/44` は Euler で `fx` を読む粘性・拡散経路をほとんど通らないので**不変の確認にしかならない**。`case/48` 平板 (`run_0020`/`run_0021`: $q$ 最大 0.0031 %・$\tau_w$ 0.0020 %)、`case/52` スラブ (`run_0006`/`run_0007`: 5e−7) も不変。**ただし §6 が要求する `check_convergence.py` の `PASS` は平板・スラブとも取れていない** (4 本とも `NOT CONVERGED (stalled/plateau)`。`--from-floor case/48…/run_0011_Bplain_tw300_y3` は参照 run 自身が継続 run で `PASS` でないため `REFUSED`)。差が小さいことは回帰合格の証明ではないので、**平板を `procedures/verification/48-flat-plate-cooled.md` の全段起動で新バイナリから回し直して `PASS` と基準値 (C_f/VD-II・2St/C_f・エネルギー閉合) を取る**のを #2c として残す |
| 2c | ~~平板の全段起動による回帰合格 (codex result M1)~~ | **済 (2026-09-22)**。`procedures/verification/48-flat-plate-cooled.md` の手順どおり新バイナリで冷間起動: `case/48.flat_plate_cooled_m4/run_0024_A_ad_y3_fx05` (断熱、全段 + 本段 48000) → `run_0025_B_tw300_y3_fx05` (冷却壁 300 K、soft + ランプ + 本段 48000)。`cooled_plate_eval.py --ref … --closure --series --integrals`: $C_f$/VD-II **0.969–0.994** (基準 0.97–0.99)、$2St/C_f$ **1.15–1.16** (1.16)、エネルギー閉合 **1.016** (1.017)、θ の CONTUR 比 0 %、`SERIES VERDICT: STEADY`。2026-09-12 の基準 run `run_0005_B_tw300_y3` を同じツールで評価し直した値との差は $C_f$ +0.1〜+0.2 %、$q_w$ +0.1〜+0.2 %、δ\* −0.9〜+0.8 %、積分 $C_D$・熱量とも +0.18 % = **手順の回帰基準 (3 %) の 1/15 以下**。この差には 2026-09-20 のリミッタ既定変更も含まれる。**`check_convergence.py` の `PASS` はこのケースでは取れない** — 前縁特異点による残差床 (`rms_ro`≈1e−7、`rms_roe`≈0.2) が手順書に既知の床として明記されており、基準 run 自身も `NOT CONVERGED (stalled/plateau)` である。そのため §6 の合格条件を「手順書の物理量基準 + 系列 `STEADY`」に改めた (下)。付随: `gen_runs.py` の廃止キー (`meshFormat`・`isCompressible`・`ro`・`last.control`・`kInf`/`omegaInf`) と、`interp_field.py` より後に config を書いていた順序を修正 |
| 2b | ~~回転不変性の直接試験~~ | **済 (2026-09-21)**。`case/48` の平板のメッシュ・場・入口速度を **z 軸まわりに 30° 回した**同一問題 (`mesh.h5` の座標・面ベクトル・運動量を回転) を 4000 step。回す前の同じスキームの解との差 ($x/L$ 0.3–0.95 の rms): 幾何 `fx` (`run_0022_rot30_fx_ctrl`) は $\tau_w$ **0.085 %**・壁面勾配流束 **0.254 %**・`q_eff` の平均 **+0.28 %**、`fx=0.5` (`run_0023_rot30_fx_half`) は **0.024 %**・**0.062 %**・**+0.01 %**。**4000 step 時点で、回転前後の差が `fx=0.5` で 3.5 倍 ($\tau_w$)・4.1 倍 (壁面勾配流束) 小さい**。両 run とも `NOT CONVERGED` (`rms_roOmega` の低下 0.6–0.7 桁で停滞) なので、残る差を座標の float32 丸め ($x\sim1$ m で 0.06 µm = 第一層 3 µm の 2 %) に帰属させることは**まだできない** (codex result M3 採用。原点移動または倍精度座標の対照が要る)。回したメッシュでは `q_eff` の節点間ノイズが両スキームとも 0.9 % 出る — 壁 CV の質量残差 (float32 の面ベクトル閉性) 由来で `fx` とは無関係。旧定義 `iface_q_eff_raw` なら 3.1 % |
| 3 | ~~合格条件の固定 (codex M4)~~ | **済**。§6 に反映。生産設定の A/B (`case/53` `run_0125`/`run_0126`、`case/54` `run_0024`/`run_0025`) は領域別偏差・うねり rms・2 節点振幅の時系列 (`tools/h_series.py`) を `check_quasisteady.py --series-csv` にかけ 4 本とも `ALL STEADY` |
| 4 | ~~恒久実装~~ | **済 (2026-09-22)**。`calcStructualVariables_d_wrapper` で `nodeMode = (discretization=="node")`、環境変数 `FORGE_NODE_FX_HALF` は撤去。環境変数版との照合 `case/53.c3x_vane_cht/run_0127_fxhalf_permanent_check`: 4000 step 後の最大相対差 ρ 4.6e−5・T 2.0e−5 (後流の非定常の範囲) |
| 5 | ~~翼の生産 run の更新~~ | **済 (2026-09-22)**。報告に載る run を新スキームで 60000 step 継続し直した: `case/53` `run_0125` (生産)・`run_0128` (SU2 対照)・`run_0129` (節点配置)・`run_0130`–`run_0133` (入口乱流)・`run_0134` (層流)・`run_0135` (1 µm)、`case/54` `run_0024` (生産)・`run_0026` (層流)・`run_0027`/`run_0028` (格子)。引用する量は `tools/report_numbers.py` で `check_quasisteady --series-csv` の判定つきで出す。連成は CHT plan §5.1 #62 |
| 6 | cell の `fx` の式 | 回転不変でない式のまま (同じ幾何を 45° 回すと重みが 0.8 → 0.637 に変わる: codex の代数チェック)。cell は使わない方針なので**修正せず記録のみ** |
| 7 | ~~格子系列での次数測定 (codex M2 / result-2 M1・M2)~~ | **済 (2026-09-22 改訂: 熱 + 運動量 + 粘性仕事)**。[`notes/investigations/2026-09-22-mms-node-face-weight.md`](../../notes/investigations/2026-09-22-mms-node-face-weight.md)。内部面の式は共通関数 `solver_density_cuda/tools/node_visc_face.py` で、**実機の場に当てるとカーネルが積んだ熱伝導・粘性仕事を 480 壁節点で 0.28 / 0.00 W/m² の差で再現する** (旧式の重みのとき。射影・0.5 では仕事項が 2.0 kW/m² ずれる = 照合は重みを見分ける)。製造解は半径 0.2 m の 15°–30° 扇形 (旧式の欠陥が最大の向き)、最粗 AR 350、成長率 1.1 / 1.2、壁沿い間隔に変調 ±30 % + 交番 ±0.3 %、4 水準 (両方向 2 倍、第一層 1/2、成長率は平方根)。**結果**: (i) **旧式は不整合** — 壁 CV に入る粘性仕事の誤差が細分化しても 25–33 % から動かず、重みが 0.5 に寄らない (0.34–0.69)。実機のうねりと同じ大きさ。節点値も $p$ = 1.8–1.9。(ii) **`fx=0.5` の節点値 (温度・速度) は $p$ = 2.00**。(iii) **壁 CV に入る流束は `fx=0.5` でも 1 次** (粘性仕事 31 % → 2.9 %)。双対面の面積重心が弦のたるみの分だけ辺中点より壁側にあるため。射影の重みはこの項に限れば誤差が半分だが、節点値の誤差は 1.25 倍で、実機の格子では 0.03–1.00 に振れるので採らない。`q_eff` への影響は +0.8 % 程度 (たるみ / 第一層の半分 = 0.3)。(iv) 2 点差分の壁面熱流束・壁せん断は重みに依らず 1 次 ($p\approx0.95$)。初版の「壁面熱流束も 2 次」は、厳密温度に同じ 2 点差分を当てた値との差を見ていた誤りで撤回 (codex result-2 M1)。**限界**: SST・化学種の拡散は未測定、実カーネルそのものは呼んでいない (上の照合で担保)、節点勾配は解析値、対流・リミッタ・float32 なし |
| 8 | codex result レビュー | **1 回目 2026-09-22: NO-GO (C0/M4/m1)** — 恒久実装は設計どおりだが `accepted` へは移さない。#2c・#7 を終えてから再レビュー |

## 6. 検証

- **ビルド**: `tools/build_native_wsl.sh`。
- **検証ケース**: (a) `case/53.c3x_vane_cht` 一様壁温 (§4.3 の A/B を恒久実装で再現)、(b) node の平板 (`procedures/verification/README.md` の
  node 乱流平板。$C_f$ が SU2/理論比で悪化しないこと)、(c) `case/52.conjugate_slab` (1 次元。`fx` に依らず結果不変のはず)。
- **判定基準**: (a) `q_eff` うねり rms が対照の 1/2 以下、$h$ の領域別偏差の変化 3 pt 以内。(b) $C_f$ の変化 1 % 以内。(c) 壁熱流束の変化 0.01 % 以内。
  いずれも NaN なし。
- **比較の作法** (codex M4 採用): A/B は**同じ熱流束定義**で比べる (`iface_q_eff_raw` と `iface_q_eff` を両方保存)。
  §4.3 の表は蓄積項補正を入れる前のバイナリで両 run とも旧定義。定常回帰の合格条件は **(b) 平板 = 手順書 `48-flat-plate-cooled.md` の物理量基準 (基準 run 比 3 % 以内) と `SERIES VERDICT: STEADY`**、
  **(c) スラブ = 解析解との一致が不変** とする。当初ここに書いた「`check_convergence.py` の `PASS`」は、**今回の基準 run も試行 run も満たさない** (平板は手順書が既知の床として明記する前縁特異点の残差床、
  スラブは純伝導で残差が丸め床に張り付く) ので 2026-09-22 に改めた。**これは物理量の回帰合格であって残差収束の合格ではない**
  (新 run の `rms_roUy` は 1.38e−4 で基準 run の 1.60e−5 の 8.6 倍あり、残差の床が同じとも言えない)。(a) の翼は全 run が `NOT CONVERGED (stalled/plateau)` なので
  **機構診断**として扱い、領域別偏差とうねり rms の時系列を `check_quasisteady.py --series-csv` で `STEADY` と確認した量だけを引用する。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| result (2 回目) | `2026-09-22` | [`notes/reviews/2026-09-22-discretization-node-face-weight-midpoint-result-2.md`](../../notes/reviews/2026-09-22-discretization-node-face-weight-midpoint-result-2.md) | **NO-GO**, C0/M2/m3 | **全件採用、`active` のまま**。平板の合格条件を手順書の物理量基準に変えたことは「妥当、逃げではない」との判断 (残差収束とは区別して書くこと、との但し書き)。M1 (製造解の壁面熱流束の誤差が解析値に対するものでない) → 指標を 3 つに分けて測り直し、「壁面熱流束も 2 次」を撤回 (§5.1 #7)。M2 (主因の粘性仕事を測っていない) → 運動量と粘性仕事を追加し、面の式を実カーネルの出力と照合した共通関数に (同)。m3 (「$p\ge1.8$ は `fx=0.5` のみ」は数値と矛盾) → 撤回。m4 (§4.3 の合格根拠) → 同一定義の生産 A/B (#3) を合格根拠とし、§4.3 は「旧定義の有限区間での低減傾向」に限定。m5 (#2 の文言・「原理的に取れない」) → 下記のとおり訂正 |
| result | `2026-09-22` | [`notes/reviews/2026-09-22-discretization-node-face-weight-midpoint-result.md`](../../notes/reviews/2026-09-22-discretization-node-face-weight-midpoint-result.md) | **NO-GO**, C0/M4/m1 | **全件採用、plan は `active` のまま**。M1 (平板・スラブに `PASS` が無い) → §5.1 #2 を差し戻し #2c。M2 (回帰は次数測定の代わりにならない) → #7 を `accepted` の前提に。M3 (回転試験の残差を丸めに帰属) → #2b の表現を限定。M4 (§4.3 の準定常判定が別 run) → `run_0122`/`run_0123` 自身の判定を §4.3 に追記 (うねり rms は試行側 `DRIFTING`、補正後定義の `run_0128` で `STEADY`)。m1 → `plans/README.md` 同期。実装 (node 内部面のみ 0.5、`roe/ro × Rro` の係数・符号) は「設計どおり」との確認 |
| plan | `2026-09-21` | [`notes/reviews/2026-09-21-discretization-node-face-weight-midpoint-plan.md`](../../notes/reviews/2026-09-21-discretization-node-face-weight-midpoint-plan.md) | GO-with-changes, C0/M4/m0 | **全件採用**。M1 (蓄積項の係数は $H_w$ でなく $e_w$) → CHT plan §5.1 #61 を訂正し `conjugateWall.cpp` を修正。M2 (値補間の 2 次と演算子の次数の混同) → §4.2 を限定し §5.1 #7。M3 (`fx` の実使用一覧と回帰範囲) → §3・§5.1 #2。M4 (未収束同士の比較・定義変更の混入) → §6・§5.1 #3、生産値の差し替えは result レビュー後 (#5) |

## 7. 影響範囲

- `solver_density_cuda/cuda_forge/calcStructualVariables_d.cu` (1 行)。**node の全 run の数値が変わる** (等方セルでは丸め程度、壁解像の高 AR 曲面層では上記の通り)。
- `methods/discretization.md`、`plans/active/architecture-node-option-consolidation.md` の `nodeMidpointFx` 行 (前提の訂正)。
- 報告 (Cooled Vane CHT Validation) の壁熱流束の図と数字。

## 8. 未確定事項

- 境界半割面の `fx` は現状のまま (壁半割面の粘性流束は `fx` を読まない)。入口・出口面で問題が出るかは未確認。

## 9. 変更ログ

- `2026-09-22` codex result 再レビュー (result-2) NO-GO (C0/M2/m3) → 全件採用。製造解を熱 + 運動量 + 粘性仕事に拡張し、面の式を実カーネルと照合した共通関数に載せ替え。「壁面熱流束も 2 次」は撤回。
- `2026-09-22` #2c (平板の冷間起動: 基準 run と 0.2 % 以内) と #7 (製造解: `fx=0.5` は $p$=2.00) を実施。§6 の合格条件を手順書の基準に改めた。codex result 再レビューへ。
- `2026-09-22` codex result レビュー 1 回目 NO-GO → 全件採用。実装は維持、plan は active のまま #2c・#7 を残す。
- `2026-09-22` 恒久実装。報告の run を新スキームで回し直し、報告を改訂。codex result レビューへ。
- `2026-09-21` 回帰 (`case/09`, `case/44`: `check_field_regress` PASS) と回転不変性の試験 (`case/48` 30° 回転) を実施。恒久実装へ。
- `2026-09-21` codex plan レビュー (GO-with-changes)。4 件採用。平板・スラブの回帰は不変を確認。
- `2026-09-21` 起票。真因の特定 (`run_0121_resid_split`) と環境変数での A/B (`run_0122`/`run_0123`) まで。
