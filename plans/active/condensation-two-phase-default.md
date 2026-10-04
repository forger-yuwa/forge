# 凝縮の二相拡散を既定にする (段階既定化)

## メタ

- **area**: `condensation`
- **status**: `draft`
- **related_docs**:
  - `methods/condensation.md` §7c (二相拡散の現在仕様)
  - `procedures/solver-settings.md`・`procedures/recommended-settings.md` (既定変更の記載先)
- **related_plans**:
  - `plans/active/condensation-two-phase-transport.md` (親。初版の実装・検証、§5.1 #9 のユーザ決定)
  - `plans/accepted/species-passive-scalar-unification.md` (dual-time の受動種 FCT、§5.1 #26 の緩和の教訓)
  - `plans/accepted/condensation-source-limiter-steady.md` (`condLimiterMode 1`)
- **created**: `2026-10-04`
- **owner**: `CFD Dev`

## 1. 目的

**ユーザ決定 (2026-10-03)**: 「物理的には二相拡散は ON にしないといけない。既定にしていってほしい」。
既定経路 (`condTwoPhaseDiffusion 0`) は総水分 $Y_w$ だけを拡散し、液 $\rho g$・モーメントは拡散しない (液が渦で混ざらない)。
二相拡散 (ON) は液・Q を他の輸送量と同じ渦拡散で混ぜ、蒸気は気相組成 $z_v$ の勾配で分子拡散する。
これを、検証済みの範囲から段階的に既定にする。完了時: 検証済みの定常包絡で既定 ON、包絡外は黙って OFF にならず、実効の作用素が記録に残る状態。

## 2. スコープ

- **やる**: 既定 ON の段 (S0・S1)、包絡外の扱いの 3 分類、実効作用素の記録、既定変更の手順。
- **やる (後段 S2)**: dual-time への拡張 (緩和・θ を sub-iter の不動点から外す設計と検証)。
- **やらない**: CPG carrier (空気凝縮) の乱流部 (親 plan #7 → 別 plan)、RK・`speciesImplicitCoupling 2`・平衡凝縮・凝縮種 2 以上 (需要が出たとき)、有限液シュミット数 (親 #8)。

## 3. 関連 docs と前提

- 判断: 2026-10-04 diagnostician (Fable、ユーザ指定)。ブリーフ `notes/reviews/briefs/2026-10-04-twophase-default-plan.md`。§4・§6 は応答の骨子を採用して書いた。
- 検証済みの ON の実体は **`condTwoPhaseSolver 1` (緩和整合 scalar-DPLUR) + `condTwoPhaseNonnegLimit 0`** (親 #4g・#4h・#4j、case/16 `run_0516`/`0519`/`0521`)。
  ところがキーの既定は Solver 0 (点対角)・NonnegLimit 1 (`input/solverConfig.cpp:1159-1162`) で、その組は親 #4g の A 側 (θ=0 セル 359/更新・残差停滞) と #1b の `run_0506` (残差上昇停滞) で悪い。
  **`condTwoPhaseDiffusion` の既定だけを 1 にすると、未検証かつ既知の悪い組が既定になる。**
- 検証済みの数値設定は case/16 `run_0482` 由来 (cfl 2.0、implicitRelax 1.0、nStepInner 5、coupling 1、SFR 2)。生産レシピ (cfl_pseudo 4〜6 + implicitRelax 0.7、`procedures/recommended-settings.md:86`) での ON は未検証。DPLUR の sweep の ω は `cfg.implicitRelax` (`condensationTransport_d.cu:852,932`)。
- 影響範囲 (既存 run の棚卸し, diagnostician): 凝縮 run は case/16 (43 本、TP carrier・定常 NS)、case/34 (CPG 1 本)、case/44 (3 本、`unsteady 1` だが `viscMethod 0`)。
  `condTwoPhaseDiffusionValidate` は viscMethod 0・CPG・pure を不活性で通す (`condensationTransport_d.cu:729-736`) ので、**既定 ON で結果が変わるのは case/16 型だけ**の見込み (G0 で確認)。
- 拒否条件の既定値: `passiveScalarScheme` 1・`condLimiterMode` 1・`condEquilibrium` 0・`speciesImplicitCoupling` 0 → 普通に書いた凝縮 NS config は拒否条件に当たらない。
- **dual-time の教訓 (ユーザ注意)**: species-passive-scalar-unification §5.1 #26 — dual-time の sub-iter 数依存の原因は `implicitRelax 0.7` 単独 (倍精度 nSub 40 vs 80 の終了場 L2 差が 9.0e-6 → 5.65e-4)。緩和は不動点を変えないが nSub 内で収束しきらず、残差ノルムには出ない。
  二相更新には緩和の入口が 2 つ (DPLUR の ω = implicitRelax と `condTwoPhaseRelax`、`twoPhaseDiffusion_d.cuh:188`) あり、θ_vg (非負、:201-204) と commit の床 (:223-227) は不動点そのものを動かしうる。設計メモ §14.3 の 3 セル 1000 物理更新の総液量 3.14e-6 (>1e-6) も同型の警告。
- **温度への効果の書き方 (diagnostician 指摘)**: −L·J_l は EOS の −L δ(ρg) と相殺し、温度には −R_wT·r_g (L の約 4 %) しか効かない (`methods/condensation.md` §7c)。
  −L·J_l は「液を運ぶと偽の加熱が出ない」ための項であり、ON−OFF の壁温差 −2.4 K を「潜熱輸送の効果」と書かない。既定化の根拠は「液・Q が他の輸送量と同じ渦拡散で混ざり、蒸気は気相組成の勾配で拡散する」混合作用素の整合に置く。

## 4. 設計方針

1. **既定の意味**: `condTwoPhaseDiffusion` 省略時 = 1。キー ON の実効条件は `condTwoPhaseDiffusionValidate` の包絡: TP carrier (`condGasSpecies ≥ 0`、`thermalMethod 2`、`nSpecies ≥ 2`)・`viscMethod ≠ 0`・`unsteady 0, timeIntegration 11`・`speciesImplicitCoupling ≠ 2`・`passiveScalarScheme 1`・`condEquilibrium 0`・`condLimiterMode 1`・`nCondSpecies 1`。
2. **包絡外の 3 分類** (「同じ物理が構成で違う式になる」危険を、黙って OFF に落ちる構成を作らないことで防ぐ):
   - (a) **物理が同一 → 不活性 (INFO)**: `condensation 0`、`viscMethod 0` (拡散自体が無い)。ON/OFF でビット一致。
   - (b) **モデルが構造的に適用できない → 不活性 + 毎回 WARNING + 出力に記録**: CPG carrier・pure 凝縮。「液は拡散しない (旧近似)」を `[twophase]` 行と `res_*.h5` 属性 (`twophase_diffusion_effective = 0`) に書く。
   - (c) **実装が未対応 → エラー終了** (黙って OFF にしない): dual-time・RK・coupling 2・平衡凝縮 (g が輸送変数なら; 代数量なら (b) — 要確認)・凝縮種 2 以上。メッセージに「旧作用素で回すなら `condTwoPhaseDiffusion: 0` を明示」。明示 0 は起動時 WARNING「legacy operator: liquid not diffused」。
3. **ON の既定の組を検証済みに揃える**: ON のとき `condTwoPhaseSolver` 既定 1・`condTwoPhaseNonnegLimit` 既定 0。点対角・非負 θ は診断用 opt-in として残す (削除しない)。
4. **実効作用素の記録**: `[twophase]` 起動行・`res_*.h5` 属性・`stage_manifest.json` の方程式署名に実効の `condTwoPhaseDiffusion` を入れる (今の `tools/stage_manifest.py:73` は二相キーを見ないので、OFF の場から ON で継続すると同一区間に連結されうる)。`check_solver_config.py` は 0 明示を WARN、dual-time + ON を FAIL。
5. **dual-time 拡張 (後段 S2)** — sub-iter が収束したとき BDF 解が nSub・ω・θ に依存しない契約:
   - 更新順は定常と同じ。R_v = R_w − R_g と D_v は BDF 対角込みで線形整合 (`speciesTransport_d.cu:1645-1653`、`twoPhaseDiffusion_d.cuh:149-154`)。
   - θ_thr (dg_max/dT_max) は**各物理 step の初回 sub-iter だけ** (既存の受動種と同じ規則、`condensationTransport_d.cu:466`)。
   - **緩和は dual-time では使わない**: `condTwoPhaseRelax < 1` + dual-time は起動拒否、`implicitRelax < 1` は既存どおり FAIL (#26)。
   - θ_vg と commit の床は sub-iter 内では「記録付きの最後の砦」。作動量を `[twophase-corr-gate]` に積算し、総量比 ≤1e-6 を要求 (作動 = 不動点が動いた証拠)。有界化の本体は物理 step 末尾の FCT。
   - FCT: 液 g と Q2/Q1/Q0 を既存の受動種 FCT に乗せ、二相カーネルの面の乱流係数を **F_H と F_L の両方**に入れる (今の低次作用素の拡散はトレーサ 1 本だけ `speciesTransport_d.cu:1808`)。総水分 ρY_w は化学種経路のまま FCT を通さない。蒸気の非負は step 末尾の実現可能性射影 (記録) に任せる。
6. **既定変更の手順**: `solver-settings.md` §9 旧設定に「0 = 既定変更日以前の作用素」を日付付きで、`recommended-settings.md` の凝縮 NS レシピに ON 既定と検証済みの組を明記。case README の該当 run 行に「既定 ON 後の再現には 0 を明示」。OFF の場から ON での restart は同一メッシュ restart として許し、manifest では新しい区間に切る。

## 5. 実装ステップ

段 S0 → S1 → (S2) の順。各段の合格は §6。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | codex plan 段レビュー | §4・§6 を `codex_review.py --stage plan` で点検し、採否を §6.1 と本表へ | F |
| 2 | S0: 既定の組・記録・ツール | ON 時の `condTwoPhaseSolver` 既定 1・`condTwoPhaseNonnegLimit` 既定 0、包絡外 3 分類のメッセージと属性、`stage_manifest` の方程式署名、`check_solver_config.py`。**`condTwoPhaseDiffusion` の既定はまだ 0**。合格: G0 (旧 config のビット不変 = ノイズ床内) と負例試験 | O |
| 3 | S1-a: G1 (0 step 作用素 A/B) | 親 #4k(1)。`run_0520` res_16000 を共通入力に、拡散作用素だけ A=OFF / B=ON を forge の診断出力で評価。合格は §6 G1 | O (解釈 F) |
| 4 | S1-b: G2 (生産レシピ ON) | case/16、run_0482 由来 IC、cfl_pseudo 4〜6 + implicitRelax 0.7 + nStepInner 5 + Solver 1 + nn0。合格は §6 G2・G3 | O (解釈 F) |
| 5 | S1-c: 既定を ON に | G0〜G3 全 PASS のときだけ。docs (§4-6)・変化量の転記 | O |
| 6 | S2: dual-time 拡張 | §4-5 の設計。μt > 0 の NS 凝縮 dual-time ケースを用意 (case/44 は Euler なので不可)。合格は §6 G4 | F (設計) / O |
| 7 | 平衡凝縮の分類 | `condEquilibrium 1/2` で g が輸送変数か代数量かを確認し、§4-2 の (b)/(c) に振り分ける | O |

## 6. 検証

事前登録 (結果を見てから変えない)。

- **G0 不活性クラスの不変**: case/44 dual-time Euler (`run_0487` 系 config)・case/34 CPG・case/16 乾き (`run_0524` config) を変更前後のバイナリで、親 #5a と同じプロトコル (旧 ×3・新 ×2、短 step、`check_field_regress` ノイズ床 × 2) → PASS。S0 の後と S1-c の後の 2 回。
- **G1 (0 step 作用素 A/B、既定化のゲート)**: (i) A の総水分流束が乱流域で |J_w^A| ≤ 0.05·max|J_l^B| (旧経路は ∇Y_w ≈ 0 で流束ほぼ 0)、(ii) B の J_l = −(μt/Sc_t)∇g が格納場から float64 で組んだ値と 8ε₃₂ で一致し、L·|J_l| / 熱流束の p50 が `analyze_liquid_diffusion_error.py` の 0.77 の 0.5〜2 倍、(iii) B の分子蒸気流束が ∇z_v で評価した値と一致。**どれかが外れたら既定化を止める** (符号・大きさの説明を後から変えない)。
- **G2 生産レシピでの ON**: 親 #4j の (A)〜(D) (7 報告量 STEADY `--tail 0.5 --drift 0.0001 --osc 0.0001 --min-snaps 21`、非有限/負値 0、RISING 0、corr-gate ≤ κ) + NaN 0。加えて**レシピ感度 < 既定化の効果**: 7 量の |ON@生産 − ON@cfl2 (`run_0521`)| ≤ 0.1 × |ON − OFF| (親 #4j の差: g_exit_mw 3.24 %、pw42 0.78 %、Tw −2.44 K など)。超えたら「レシピ依存」として既定化を保留。
- **G3 収支**: ON run の `[twophase-corr-gate]` commit・floor・clamp ≤ κ、射影は記録。再正規化 max|f−1| は OFF と同程度 (記録)。
- **G4 dual-time (S2)**: species-passive-scalar-unification #28 のゲートを流用 — nSub 倍増比 ≤0.1、3 水準次数 (流れ・化学種 BDF2 1.7–2.3、受動種 ≥1.3)、`rms_roYv` を含む全列の sub-iter 低下 ≥2 桁、`check_passive_budget.py` PASS (floor+lim+射影+二相 withheld/vround の総量比 ≤1e-6)、**ω = 1 固定**。緩和の A/B は診断として 1 組だけ記録: |q(ω=0.7, nSub) − q(ω=1, 2nSub)| が |q(ω=1, nSub) − q(ω=1, 2nSub)| の 2 倍を超えたら「dual-time で緩和は使えない」を本 plan に書く (#26 の再現)。
- **G5 記録・ツール**: 実効作用素が `[twophase]`・h5 属性・manifest に出る。`check_solver_config.py` の負例試験 (dual-time + ON → FAIL、0 明示 → WARN)。
- **外部参照の限定**: CFD レベルの外部参照 (解析解・文献・SU2) は無い。物理の向きの根拠は (a) 単体のエネルギー接線試験 (親 plan、済)、(b) G1、(c) 作用素の構造 (乱流拡散が全輸送量で共通、分子拡散は気相内で Σj=0) の 3 つまでと明記する。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| diagnostician (§4・§6 の方針) | `2026-10-04` | ブリーフ [`notes/reviews/briefs/2026-10-04-twophase-default-plan.md`](../../notes/reviews/briefs/2026-10-04-twophase-default-plan.md)、応答の要旨は §3・§4・§6 | 段階既定化 (定常包絡で ON、包絡外は 3 分類)、既定の組を検証済みに揃える、G1・G2 を既定化のゲートに、dual-time は後段 | 全件採用 (本 plan の §3〜§6) |

## 7. 影響範囲

- `input/solverConfig.{hpp,cpp}` (既定値・包絡外の扱い)、`cuda_forge/condensationTransport_d.cu` (`condTwoPhaseDiffusionValidate`・記録)、`output/` (h5 属性)、`tools/stage_manifest.py`・`tools/check_solver_config.py`。S2 で `speciesTransport_d.cu` (FCT の拡散)・`main.cpp` (dual-time)。
- 既存ケース: 既定 ON で結果が変わるのは case/16 型 (TP carrier NS 定常) の見込み。再現には `condTwoPhaseDiffusion: 0` の明示が要る。
- docs: `methods/condensation.md` §7c、`procedures/solver-settings.md` (§9 旧設定)、`procedures/recommended-settings.md` (凝縮 NS レシピ)。

## 8. 完了条件

- [ ] S0・S1 を実施し、G0〜G3・G5 を満たして既定を ON にした (または G1/G2 の FAIL で既定化を止めた記録)
- [ ] S2 (dual-time) は別途 G4 で判定 (本 plan で行うか後継に送るかを S1 完了時に決める)
- [ ] 関連 docs を更新済み
- [ ] codex レビュー 2 回 (`plan` / `result`) を §6.1 に記録し、Critical / Major の採否を残作業表に反映済み
- [ ] `status` を `done` にし §9 に変更ログ、`plans/accepted/` へ移動、`plans/README.md` を同期

## 9. 変更ログ

- `2026-10-04` — 起票。ユーザ決定 (二相拡散を既定にしていく) と diagnostician の判断 (段階既定化、既定の組を検証済みに揃える、dual-time は後段で緩和を使わない設計) を §3〜§6 に反映。
