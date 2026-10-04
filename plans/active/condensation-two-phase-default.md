# 凝縮の二相拡散を既定にする (段階既定化)

## メタ

- **area**: `condensation`
- **status**: `in_progress`
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

1. **既定の意味と包絡**: `condTwoPhaseDiffusion` 省略時 = 1。キー ON が実際に作動する包絡 (検証済み・包絡に含めるもの):

   | 項目 | 包絡 | 検証の状態 |
   | --- | --- | --- |
   | 熱物性・凝縮形 | TP carrier (`condGasSpecies ≥ 0`、`thermalMethod 2`、`nSpecies ≥ 2`)、非平衡 (`condEquilibrium 0`)、`nCondSpecies 1`、`condLimiterMode 1` | case/16 で検証 |
   | 方程式 | NS (`viscMethod ≠ 0`)、定常 (`unsteady 0`、`timeIntegration 11`) | 同上 |
   | 離散化・幾何 | node、平面 2D / 3D、**非周期** | 平面 2D のみ実行検証 (3D は未検証と明記) |
   | 化学種の更新 | `speciesImplicitCoupling` 1 (検証済み)・0 (包絡に含めるが実行検証なし; 二相更新は非水種の更新方式に依存しない)、`speciesFaceReconstruction` 2 (検証済み)・0 (同上) | 表のとおり |
   | 受動種 | `passiveScalarScheme 1` | 検証済み |

   **軸対称・周期は、二相拡散が実際に作動する ON の NS 試験が通るまで (c) (エラー)** とする (既存の入力が無い: case/16 は平面、case/44 の軸対称は Euler で不活性、周期の凝縮ケースは無い)。cell は「未検証」を起動 WARNING で明記。
2. **判定表 (省略 / 明示 0 / 明示 1 × 実効状態)** をソルバ・`check_solver_config.py`・`res_*.h5` 属性・`stage_manifest.json` で共有する。実効状態は次の 4 つ:
   - **active**: 包絡内。省略・明示 1 で作動。明示 0 は旧作用素 + 起動 WARNING「legacy operator: liquid not diffused」。
   - **inactive-(a) 物理が同一**: `condensation 0`、`viscMethod 0` (拡散自体が無い)。どの指定でもビット一致、INFO。dual-time + Euler もここ (FAIL にしない)。
   - **inactive-(b) モデルが構造的に適用できない**: CPG carrier・pure 凝縮・**平衡凝縮の EOS 拘束形 (`condEquilibrium 2`、g は EOS の状態量で `res_rog = 0`)**。不活性 + 毎回 WARNING + `twophase_diffusion_effective = 0` を記録。**「物理が同一」とは書かない**: 既定経路のエネルギー流束は液を蒸気として数え (h_v J_w、−L·J_t(g) が無い)、分子拡散は ∇Y_w 駆動のままという近似が残る。
   - **unsupported-(c) 実装が未対応 → エラー終了**: dual-time・RK・`speciesImplicitCoupling 2`・**平衡凝縮の緩和形 (`condEquilibrium 1`、g を輸送 + 緩和ソース)**・凝縮種 2 以上・軸対称・周期 (実効判定で。省略・明示 1 のとき)。メッセージに「旧作用素で回すなら `condTwoPhaseDiffusion: 0` を明示」。
3. **ON の既定の組を検証済みに揃える**: ON のとき `condTwoPhaseSolver` 既定 1・`condTwoPhaseNonnegLimit` 既定 0。点対角・非負 θ は診断用 opt-in として残す (削除しない)。
4. **実効作用素の記録**: `[twophase]` 起動行・`res_*.h5` 属性・`stage_manifest.json` の方程式署名に実効状態を入れる (今の `tools/stage_manifest.py:73` は二相キーを見ないので、OFF の場から ON で継続すると同一区間に連結されうる)。
5. **dual-time は本 plan の範囲外 (後継 plan に分離、codex plan M3)**: 本 plan では起動拒否 (c) を維持する。後継 plan の設計要件:
   - 蒸気・液・他の気相種を**連成して制限**し、総水分は蒸気 + 液から戻す (「液だけ FCT・蒸気は射影」は撤回: 総水分 [0.10, 0.20]・液 [0.09, 0.11] に保存的補正 [+0.02, −0.02] で蒸気が負になり、clamp `condensationRealizability_d.cuh:125-127` が液を削って総液量 −5 %。FCT の制限係数は成分別で ρY_w を見ない `passiveFct_d.cuh:290-300`)。
   - エネルギー流束と BDF 履歴 (G/H) の定義、FCT の作動条件 (SFR ≥ 2・SLAU、`speciesTransport_d.cu:1781`) を包絡に明記。
   - #26 の教訓を持ち越す: 緩和 (DPLUR の ω と `condTwoPhaseRelax`) は dual-time で使わない、θ_thr は各物理 step の初回 sub-iter だけ、θ_vg・床は記録して収支で判定。検証は nSub 倍増差・3 水準次数・sub-iter 低下・収支 (旧 G4 の骨子) + 緩和の A/B を 1 組 (|q(ω=0.7,nSub) − q(ω=1,2nSub)| が |q(ω=1,nSub) − q(ω=1,2nSub)| の 2 倍超なら「dual-time で緩和は使えない」と書く)。
6. **既定変更の手順**: `solver-settings.md` §9 旧設定に「0 = 既定変更日以前の作用素」を日付付きで、`recommended-settings.md` の凝縮 NS レシピに ON 既定と検証済みの組を明記。case README の該当 run 行に「既定 ON 後の再現には 0 を明示」。OFF の場から ON での restart は同一メッシュ restart として許し、manifest では新しい区間に切る。

## 5. 実装ステップ

段 S0 → S1 の順。各段の合格は §6。dual-time は後継 plan。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | ~~codex plan 段レビュー~~ | **完了 2026-10-04**: GO-with-changes (C0/M5/m2)、採否は diagnostician に諮り全件採用 (M1・m6・m7 は内容を修正して採用)。§4・§5.1・§6 を改訂 (§6.1) | F (完了) |
| 2 | S0-前提: 平衡凝縮の分類 | **確定 (2026-10-04、コード確認 + diagnostician)**: `condEquilibrium 1` → (c)、`condEquilibrium 2` → (b) (§4-2)。判定表を S0 で実装 | O (完了) |
| 3 | ~~S0: 既定の組・判定表・記録・ツール~~ | **完了 2026-10-04** (commit c910f308、implementer 実装、CUDA は AWS でクリーンビルド)。ON 時の既定 Solver 1・NonnegLimit 0 (明示値が優先)、判定表 `condTwoPhaseDiffusionClassify` (`cuda_forge/condensationTransport_d.cuh`) と Python 版 `tools/twophase_state.py` (定数 `kCondTwoPhaseDiffusionDefault`/`DEFAULT` = 0 を一致試験)、`[twophase]` 起動行・`res_*.h5` 属性 (`twophase_diffusion_effective/state/requested`)・`stage_manifest` の hard キー、`check_solver_config.py` の G5 負例、単体試験 `tests/unit/test_twophase_state.py` ALL PASS。補足: 周期は bcond に `periodic` があれば cell/node とも (c)、既存の拒否 (`passiveScalarScheme 0`・`condLimiterMode 0`) は (c) に残した、3D の未検証は docs のみ (cfg が次元を知らない)。**G0 PASS (S0 後)**: 旧 649f2d77 ×3・新 c910f308 ×2、200 step、`check_field_regress` (液・モーメント込み): case/44 dual-time 非粘性凝縮 `case/44.vitiated_air_wt/run_0524`–`0528` (入力は run_0482_prune_regress_float、species_db の旧 H2O エントリが #10 で両バイナリとも拒否されたので外して移行)、case/34 CPG `case/34.arthur_n2_nozzle/run_0110`–`0114` (run_0108 の入力)、case/16 乾き `case/16.nozzle_wys/run_0555`–`0559`。新バイナリの記録: case/44 は `[twophase] ... state inactive-a, effective 0`・h5 属性あり、凝縮なしの case/16 乾きは行なし (仕様どおり) | O (完了) |
| 4 | S1-a: G1 (0 step 作用素 A/B) — 診断 D1 | 判断: 2026-10-04 codex diagnose (`notes/reviews/2026-10-04-twophase-g1-g3-diag-design-diagnose.md`) — D1 は条件付き採用。**設計 (確定)**: 環境変数 `FORGE_DIAG_TP_FACES=<出力 h5>` (既定 off)。組立前処理 (射影・EOS・BC・勾配) を一度だけ通した同一状態から、全通常面について OFF / ON の面作用素を評価して書き、**更新へ進まず明示終了** (nStepOuter の値に頼らない)。前処理のクランプ・境界上書きは前後差を記録。**本番カーネル (`species_diffusion_d`・`twophase_diffusion_d`) は書き換えない**: OFF の面計算は診断専用の写しで評価し、その面流束から倍精度で組み直した節点の拡散寄与が、本番の組立の拡散段の残差寄与と atomicAdd の雑音内で一致することを確かめる (写しの正しさの検査)。ON は本番と同じ `tp_build_face_in`/`tp_face_flux<float>` と、格納値 (ρ、ρY、ρg、ρQ) から double で差・除算・正規化をやり直した `tp_face_flux<double>` (監査 #4f と同じ経路)。出力: 面 id・両端節点・向き・面重心・幾何 (geo, geo_abs, 面積)・f・ct・D_k (OFF と ON の実係数を別々に; 定数 Schmidt では分母が ρ_f と ρ_g,f で違う)・h_k・L・vis_lam・vis_turb・両端の状態・skip 理由、OFF の J_s・q、ON の J・Jv・Jl・JQ・q・Sm・up0 (float と double)、**分子蒸気流束 j_v = j_v⁰ − z_up,v·Σj⁰ を補正前・補正項と分けて**、誤差尺度用の差し引き前の項。`TpFaceIn` 全体を保存して参照を再現可能に。境界半割面は拡散 skip (記録)、周期は包絡外で拒否。**G1 の事前登録の詳細 (実行前に固定)**: 符号 = J はセル 0 へ入る向きが正の面積分値 (外向きへの変換は解析側)。乱流域 = 面の μt,f/μ_f ≥ 1。(i) 乱流域の全面で abs(J_w^A) ≤ 0.05·max abs(J_l^B)。(ii) abs(Jl_float − Jl_double) ≤ 8ε₃₂·A_l、A_l = abs(ct·geo)·(abs(ρg₀/ρ₀) + abs(ρg₁/ρ₁))、分子蒸気は abs(j_v,float − j_v,double) ≤ 8ε₃₂·A_v、A_v = abs(ρ_g,f D_w geo)·(abs(z_v0)+abs(z_v1)) + abs(z_up,v)·Σ_k abs(ρ_g,f D_k geo)·(abs(z_k0)+abs(z_k1)) (結果を見て広げない)。double 照合が外れたら診断未成立。(iii) **第 1 仮説の判別 (記録、既定化のゲートにはしない)**: x 35〜70 mm の壁距離帯 0〜0.1・0.1〜0.4・0.4〜1.6 mm (上下壁を壁距離差と面の向きで判定し、y の符号でまとめない) で、壁法線方向の面積分和が 分子蒸気 = 壁から外向き・液乱流 = 壁向き、かつ誤差尺度の和を超える → 支持、そうでなければ棄却。入力 = `case/16.nozzle_wys/run_0520_twophase_off_16k` res_16000 (主判定)。ON の最終場 (`run_0521`) でも出してよいが主判定と混ぜない **実装 (2026-10-04, implementer, 本番カーネル不変・追記のみ)**: `cuda_forge/twoPhaseFaceDiag_d.cuh`・`speciesTransport_d.cu` 診断節・`main.cpp runTpFacesDiag`・判定 `notes/investigations/2026-10-04-twophase-g1/g1_judge.py`。**実行前に固定した実装上の選択 (結果を見て変えない)**: (1) 面の拡散係数が読む μ_t は組立後半 (`turbulent_viscosity`) で作られるので、前半 → 後半 → 面評価の順に通す (後半が面の入力 ρ・T・P・μ・ρY・液・Q を変えていないことをバイト比較で記録し、変わっていれば (i) は UNDETERMINED)。前処理の差の基準は初期化後の状態。(2) OFF の写しの一致検査 P0: 本番 `species_diffusion_d` を同じ状態のまま 0 初期化した別配列へ流した節点寄与と、写しの面流束から倍精度で組んだ節点寄与の差 ≤ n_face·ε₃₂·Σabs(J) (P0 不成立なら G1 は UNDETERMINED)。(3) (i) の max abs(J_l^B) は乱流域の面での最大 (全面の最大は参考)。(4) (iii) の面の壁距離は両端節点の wall_dist の算術平均、壁法線成分は J·Δwd/abs(Δx_cc)、誤差尺度の和も同じ abs(cos) で重み、上下壁は壁距離が増える向きの y 成分の符号で分け、面のある全ての (帯 × 壁) が支持のときだけ総合で支持。(5) judge の終了コード 0 = (i)(ii) PASS、1 = FAIL、2 = UNDETERMINED。run ディレクトリは「診断 run (更新なし)」として収束ゲートの対象外 | O (解釈 F)
| 4r | #4 結果 (2026-10-04, AWS, build 3d6e853f) | 診断 run `case/16.nozzle_wys/run_0560_g1_faces_off0520` (更新なし; 成果物の写し `notes/investigations/2026-10-04-twophase-g1/evidence/`)。前提: P0 (OFF の写し vs 本番 `species_diffusion_d`) は全成分で許容の 0.28〜0.30 倍 → ok。面 83745 (評価 82813、境界半割面 932 skip、OFF/ON の skip 不一致 0、Sm/up0 の自己検査不一致 0)。前処理の追加変化は roe 1921 節点 (最大 3.9e-3)・ghost の roUx/roe・roQ1 13 節点 (8e-39) — 記録のみ。
**(ii) float vs double: 事前登録の判定は FAIL** — 液: 574 面で abs(Jl_f − Jl_d) > 8ε₃₂·A_l (最大比 1.05e6)、分子蒸気: 0 面 (最大比 0.068)、非有限 0 (NaN は skip 面 932 のみ)。**観測**: 574 面はすべて真値 abs(Jl_d) < FLT_MIN (1.18e-38) の非正規化数の範囲で、絶対誤差は最大 5.4e-45 (全面の最大 abs(Jl) 1.44e-5 に対し 40 桁下)。正規化数の範囲の面では最大比 0.19・超過 0。→ 事前登録の誤差尺度は非正規化数の範囲 (相対精度が失われる) を考慮していなかった。**事前登録の規則では「double 照合が外れたら診断未成立」なので、(i)(iii) の判定も未成立のまま**。(i) の生の比較: 乱流域 50446 面で max abs(J_w^A) 1.42e-10 vs 0.05·max abs(J_l^B) 7.18e-7 (比 0.000)。(iii) 記録: 帯 0.1〜0.4・0.4〜1.6 mm は上下とも支持 (分子蒸気は壁から外向き、液は壁向き、誤差尺度の 10³〜10⁵ 倍)、帯 0〜0.1 mm は棄却 (分子蒸気が壁向き −8.9e-9、液は ~2e-23 で液がほぼ無い) → 事前登録の「全ての帯で支持」を満たさず総合は棄却。**扱いは上位に諮る (条件 3)** | O (完了; 解釈 F) |
| 4s | G1(ii) の判別 A/B (指数スケール) | 判断: 2026-10-04 codex diagnose (`notes/reviews/2026-10-04-twophase-g1-result-diagnose.md`)。採否 (全件採用): 誤差尺度の改訂は**条件付き** (実際の演算と FTZ 条件から伝播した絶対誤差項を導く、観測最大値に合わせて床を選ばない、旧基準 FAIL と改訂を併記)、「FLT_MIN 未満は表現不能で対象外」は却下 (非正規化数は FLT_TRUE_MIN 刻みで表現できる、除外しない)、壁直近帯の解釈を訂正 (分子蒸気は誤差尺度の 0.02 倍で**向きは確定できない**、液は微小だが壁向きの流束が尺度の 10⁶ 倍で検出されている)、G1 成立前に G2 へ進まない、本番の global `atomicAdd` (f32) は非正規化の入力・結果をゼロ化する仕様なので面評価の誤差と本番組立で失われる寄与は別物 (後者は G3 の組立誤差として計上)。**事前登録 (実行前に固定)**: ビルド条件 = nvcc 既定 (`-ftz` 指定なし = 非正規化数を保持、fast-math なし; CMakeLists に FTZ 系の指定が無いことを確認)。保存した `tp_faces.h5` (run_0560) の面入力から、**本番の `tp_face_flux<float>` と double 参照を GPU の試験プログラム** (`tests/unit/tp_jl_scale_ab.cu`、既定のコンパイル条件) で評価する。A = 保存した `rg0, rg1` のまま、B = `rg0, rg1` を正確に 2⁶⁴ 倍 (ρ・ct・geo ほかは固定; 2⁶⁴ で最小の非零液入力 ~1e-45 の流束が正規化数に入り、最大の液 ~1e-2 は溢れない) して評価し、B の Jl は double 上で 2⁶⁴ で割り戻す。対象 = FAIL の 574 面 + 正規化数の範囲で abs(Jl_d) 上位の対照 574 面 + 全評価面。前提: A の Jl_f は保存値とビット一致 (試験プログラムが本番の評価を再現している)。判定: **A が FAIL・B が全件旧尺度 (8ε₃₂·A_l、A_l は B の入力で計算) 内 → 第 1 仮説 (非正規化数の範囲だけが原因) を支持し、アンダーフロー込みの誤差尺度の導出へ進む**。B でも FAIL → 第 1 仮説を棄却し尺度の緩和を止める | O (解釈 F) |
| 4g3 | G3 の再設計 (D2 は却下) | codex 2026-10-04: 残差差分だけを正本にしない (小さい加算項が float で消え、望遠鏡和しか確かめない)、時間当たり (流束・S·V) と更新当たり (補正 Δq·V) を同じ収支に足さない、0 step では更新補正が生じない。→ **作用素の収支** (Σ_Ω R_final = −F_adv − F_diff + S⁺ + S⁻ + B_pin + R_other + E_assembly; 面流束は atomic 前、ソースは実際に加算する S·V を保存、記録点は wrapper 内部: speciesTransport の移流後・拡散後、化学ソース、ピン前後、condensationTransport の移流後・二相拡散後、condensationSource 前後; 蒸気は double で R_v = R_w − R_l 等) と **更新写像の収支** (Σ_Ω V(q_after − q_before) = Σ V·δq_update + Σ_a Σ V·C_q,a + E_update; 補正は操作ごとの符号付き格納差、液滴消滅・境界上書き込み、twoPhaseHoldWater の二重計上なし; 一更新の診断) に分ける。合否の尺度 (成分・単位・比較する液流束差の定義) を実装前に固定する。G1 の後に設計を書いて諮る | F |
| 5 | S1-b: G2 (生産レシピ ON) | case/16、run_0482 由来 IC、cfl_pseudo 4 + implicitRelax 0.7 + **nStepInner 4** (生産レシピ) + Solver 1 + nn0。先に同条件反復でノイズを測る。合格は §6 G2 | O (解釈 F) |
| 6 | S1-c: 既定を ON に | G0〜G3・G5 全 PASS のときだけ。docs (§4-6)・変化量の転記。G0 をもう一度 | O |
| 7 | 軸対称・周期の ON NS 試験 | 包絡に入れるための試験ケースを用意 (ケース未定; 軸対称 NS の TP carrier 凝縮、node 周期の凝縮)。通るまで (c) | O (ケース選定 F) |
| 8 | 平衡形 (EOS 拘束形) の二相輸送 | `condEquilibrium 2` のエネルギー流束 (−L·J_t(g)) と ∇z_v 駆動の扱いを別途設計 (g が EOS で再決定されるので非平衡形の非分割更新は使えない。親 §4.2 の「エネルギーだけ直すのは禁止」は g を輸送する非平衡形の話) | F |
| 9 | dual-time の後継 plan | §4-5 の要件で起票 | F |

## 6. 検証

事前登録 (結果を見てから変えない)。

- **G0 不活性クラスの不変**: case/44 dual-time Euler 凝縮 `case/44.vitiated_air_wt/run_0482_prune_regress_float` (同条件反復 `run_0483_prune_regress_float_rep` をノイズ床に)・case/34 CPG・case/16 乾き (`case/16.nozzle_wys/run_0524_floor_dry_L1` の config) を変更前後のバイナリで、親 #5a と同じプロトコル (旧 ×3・新 ×2、短 step、`check_field_regress` ノイズ床 × 2) → PASS。S0 の後と S1-c の後の 2 回。入力 (config・IC・メッシュ) は完全パスで固定し、番号だけで参照しない。
- **G1 (0 step 作用素 A/B、既定化のゲート)**: 同一面・同一格納入力・同一係数で独立な double 参照を組む。絶対誤差尺度 = 8ε₃₂ × (差し引き前の項の大きさ: ρ_g,f D abs(z)・c_t abs(g/ρ) など) × 面の幾何量 (流束自身を分母にしない — 小さい勾配では差し引きの丸めが支配し、正しい演算でも相対 0.3 ずれる)。評価領域・符号規約・面積重みを事前に固定。
  合格: (i) A の総水分流束が乱流域で abs(J_w^A) ≤ 0.05·max abs(J_l^B)、(ii) B の J_l = −(μt/Sc_t)∇g と分子蒸気流束 (∇z_v) が double 参照と上の尺度内で一致。**外れたら既定化を止める**。0.77 (`analyze_liquid_diffusion_error.py` の節点回復勾配の統計) は参考値で合否に使わない。
- **G2 生産レシピでの ON**: 親 #4j の (A)〜(D) (7 報告量 STEADY `--tail 0.5 --drift 0.0001 --osc 0.0001 --min-snaps 21`、非有限/負値 0、RISING 0、corr-gate ≤ κ) + NaN 0。
  **レシピ感度**: 先に同条件反復 (2 本以上) で量ごとのノイズ σ を測る。abs(ON−OFF) > 3σ の量だけに「abs(ON@生産 − ON@cfl2 `run_0521`) ≤ 0.1 × abs(ON − OFF)」を課す。それ以外の量は絶対許容 (onset ±0.01 mm、その他はノイズ × 3)。分解能が足りなければ「判定不能」(既定化は保留、「レシピ依存」とは書かない)。7 量の抽出マスク・断面の数値照合を前提に入れる。
- **G3 閉じた収支 (親 #4k(2))**: ON・OFF それぞれ同じ制御体積で、移流 + 拡散の境界流束・相変化ソース (正負別)・残差・成分別 (ρY_w, ρv, ρg, ρQn) の数値補正 (再正規化・射影・床・clamp) を積算し、収支の不整合 < 観測された液流束差の 10 %。加えて `[twophase-corr-gate]` の commit・floor・clamp ≤ κ。準定常性だけで輸送の収支を代替しない。
- **G5 記録・ツール (既定変更前の必須ゲート)**: §4-2 の判定表どおりに実効状態が `[twophase]`・h5 属性・manifest に出る。`check_solver_config.py` の負例試験 (dual-time + 凝縮 NS + 省略/明示 1 → FAIL、dual-time + Euler → 通る、明示 0 → WARN)。
- **外部参照の限定**: CFD レベルの外部参照 (解析解・文献・SU2) は無い。物理の向きの根拠は (a) 単体のエネルギー接線試験 (親 plan、済)、(b) G1、(c) 作用素の構造 (乱流拡散が全輸送量で共通、分子拡散は気相内で Σj=0) の 3 つまでと明記する。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| diagnose (G1 結果の扱い) | `2026-10-04` | [`notes/reviews/2026-10-04-twophase-g1-result-diagnose.md`](../../notes/reviews/2026-10-04-twophase-g1-result-diagnose.md) | G2 保留、指数スケール A/B で判別, M4 | 全件採用 → #4s |
| diagnose (G1/G3 診断の設計) | `2026-10-04` | [`notes/reviews/2026-10-04-twophase-g1-g3-diag-design-diagnose.md`](../../notes/reviews/2026-10-04-twophase-g1-g3-diag-design-diagnose.md) | D1 条件付き採用・D2 却下, M5/m1 | 全件採用: D1 の設計と G1 の詳細 → #4、本番カーネルは書き換えず写し + 一致検査、G3 は作用素収支と更新写像収支に分けて再設計 → #4g3 |
| diagnostician (§4・§6 の方針) | `2026-10-04` | ブリーフ [`notes/reviews/briefs/2026-10-04-twophase-default-plan.md`](../../notes/reviews/briefs/2026-10-04-twophase-default-plan.md)、応答の要旨は §3・§4・§6 | 段階既定化 (定常包絡で ON、包絡外は分類)、既定の組を検証済みに揃える、G1・G2 を既定化のゲートに、dual-time は後段 | 全件採用 |
| plan | `2026-10-04` | [`notes/reviews/2026-10-04-condensation-two-phase-default-plan.md`](../../notes/reviews/2026-10-04-condensation-two-phase-default-plan.md) | GO-with-changes, C0/M5/m2 | 採否は diagnostician に諮った (2026-10-04、全件採用・M1/m6/m7 は修正して採用)。M1 包絡を検証範囲に絞る (軸対称・周期は (c)、cell は未検証明記、coupling/SFR を表で区別、G2 は nStepInner 4) → §4-1・#7。M2 G3 に閉じた収支 → §6 G3・#4。M3 dual-time を後継 plan に分離・「液だけ FCT」撤回 → §4-5・#9。M4 G1 を double 参照と絶対誤差尺度に → §6 G1。M5 G2 はノイズ σ を先に測り有意な量だけに 10 % 条件 → §6 G2・#5。m6 分類を S0 前提に・判定表・G5 を必須に → §4-2・#2・#3。m7 G0 の参照 run を完全パスに (case/44 `run_0482_prune_regress_float` と反復 `run_0483`) → §6 G0 |

## 7. 影響範囲

- `input/solverConfig.{hpp,cpp}` (既定値・包絡外の扱い)、`cuda_forge/condensationTransport_d.cu` (`condTwoPhaseDiffusionValidate`・記録)、`output/` (h5 属性)、`tools/stage_manifest.py`・`tools/check_solver_config.py`。S2 で `speciesTransport_d.cu` (FCT の拡散)・`main.cpp` (dual-time)。
- 既存ケース: 既定 ON で結果が変わるのは case/16 型 (TP carrier NS 定常) の見込み。再現には `condTwoPhaseDiffusion: 0` の明示が要る。
- docs: `methods/condensation.md` §7c、`procedures/solver-settings.md` (§9 旧設定)、`procedures/recommended-settings.md` (凝縮 NS レシピ)。

## 8. 完了条件

- [ ] S0・S1 を実施し、G0〜G3・G5 を満たして既定を ON にした (または G1/G2 の FAIL・判定不能で既定化を止めた記録)
- [ ] dual-time の後継 plan を起票 (#9)
- [ ] 関連 docs を更新済み
- [ ] codex レビュー 2 回 (`plan` / `result`) を §6.1 に記録し、Critical / Major の採否を残作業表に反映済み
- [ ] `status` を `done` にし §9 に変更ログ、`plans/accepted/` へ移動、`plans/README.md` を同期

## 9. 変更ログ

- `2026-10-04` — G1 結果の扱いを codex に諮問: G2 保留、液入力を 2⁶⁴ 倍する A/B (#4s) を事前登録。壁直近帯の解釈を訂正。
- `2026-10-04` — #4 G1 診断を実行: P0 ok、(ii) は非正規化数の範囲の液流束 574 面で事前登録 FAIL (正規化数の範囲は最大比 0.19)、(i) の生比較は満たす、(iii) は壁直近帯で棄却。扱いを諮問へ。
- `2026-10-04` — #4 D1 実装 (本番カーネル不変)。実装上の選択 5 点と P0 の許容を実行前に固定。
- `2026-10-04` — G1/G3 診断の設計を codex に諮問: D1 (面流束ダンプ、更新なし) を採用し G1 の誤差尺度・マスクを固定、D2 は却下して G3 を作用素収支と更新写像収支に分けて再設計 (#4g3)。
- `2026-10-04` — S0 (#3) 完了・G0 PASS (3 ケース)。既定は 0 のまま。
- `2026-10-04` — codex plan 段 (GO-with-changes, M5/m2) を diagnostician に諮って全件採用: 包絡を検証範囲に絞る (軸対称・周期は (c))、判定表、G1 の誤差尺度、G2 のノイズ基準、G3 の閉じた収支、dual-time を後継 plan に分離、平衡凝縮の分類確定 (1 → (c)、2 → (b))。
- `2026-10-04` — 起票。ユーザ決定 (二相拡散を既定にしていく) と diagnostician の判断 (段階既定化、既定の組を検証済みに揃える、dual-time は後段で緩和を使わない設計) を §3〜§6 に反映。
