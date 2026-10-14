---
name: nozzle-design
description: 軸対称風洞ノズル (wind_tunnel_axisym_axismach) の設計・計算の作業手順。(1) 軸長・出口径・入口径・入口条件を指定して非粘性 (MOC) ノズルを新規設計する、(2) 出口・入口・条件を固定して軸長を MOC でパラスタし、有望形状を Euler → NS (境界層補正) → 形状調整 → 凝縮計算する、(3) 作成済みのノズルで Euler/NS を回す、(4) 実寸の要求 (全圧・全温・出口マッハ・スロート径・入口半径・入口長さ・出口長さ・壁温) から問題を作り、STEP・gmsh・Fluent メッシュ・Excel 表の成果物一式を出す — のどれかをユーザが頼んだとき (「ノズルを設計したい」「軸長を振りたい」「このノズルで Euler/NS を回して」「凝縮を見たい」「この寸法で風洞ノズルを作って」「STEP/Fluent/gmsh/Excel で欲しい」) に使う。
---

# /nozzle-design — 軸対称風洞ノズルの設計・計算

**棲み分け**: 要件の確定 (目的・Pt/Tt/組成・M_design・寸法・凝縮の扱い・コスト承認) と problem YAML・方針メモの作成は
skill [`design-intake`](../design-intake/SKILL.md) が担う (ここでは繰り返さない)。本 skill はその YAML ができた後の
**axismach ノズルの作り方・回し方**だけを扱う。YAML が無い新規の依頼は、まず `design-intake` を通す。

**手順の正本は [`procedures/nozzle-design-workflow.md`](../../../procedures/nozzle-design-workflow.md)**。
ここには AI が判断を誤りやすい点だけを置く。作業前に正本の該当節 (§1 新規設計 / §2 パラスタ〜凝縮 / §3 既存形状) を読む。

## 実寸の要求から作る・成果物を出す (2026-10-13; plan [`tooling-nozzle-wt-deliverables`](../../../plans/accepted/tooling-nozzle-wt-deliverables.md))

ユーザが寸法を実寸で言ったら、**要求カード** (`forge.wt_card/v1`、例 `design/forge_design/intake/examples/m5_air_demo.yaml`) に写す。
ユーザの言葉とキーの対応 (2026-10-12 ユーザ回答で確定):

| ユーザの言葉 | カードのキー | 注意 |
| --- | --- | --- |
| 全圧・全温 | `conditions.Pt` [Pa]・`conditions.Tt` [K] | 絶対圧 |
| 出口マッハ数 | `M_exit` | = `spec.M_design` |
| スロート径 | `throat.diameter_physical_m` または `throat.diameter_design_m` | **物理 (境界層込みの実際の穴) か設計 (非粘性) かを必ず聞く**。物理なら r_t は CFD 前の δ で解く |
| スロート曲率 | `throat.curvature_radius_m` (設計壁のスロートの曲率半径 [m]) か `throat.curvature_ratio` (= 曲率半径 / r_t = `geometry.R`) | インテークで必ず聞く (無ければ雛形の R = 2 になることを伝えて了承を取る)。実寸で与えると r_t を解き直すたびに割り直す。R < 0.5 は遷音速解の精度が落ちるので拒否。物理壁の曲率は境界層の分だけ違う (DIMENSIONS.json の info に STEP から測った値) |
| 入口半径 | `inlet.radius_m` | 縮流部の上流の直管の半径 |
| 入口長さ | `inlet.pipe_length_m` と `inlet.contraction_length_m` | **両方を指定** (直管長と縮流部長)。縮流部が長すぎると μ = L_U²/(R Δr) > 20 で単調に絞れず例外 |
| 出口長さ | `divergent_length_m` | **設計スロート (x=0) → 物理出口**の長さ。出口の後ろの直管ではない |
| 壁温 (一様) | `wall: {thermal: isothermal, Tw}` | 断熱なら `adiabatic` |
| 壁温 (分布) | `wall: {thermal: isothermal, Tw_profile: {s: [...], Tw_K: [...]}}` (s: 入口直管 −2→−1・縮流部 −1→0・拡大部 0→1) か `{x_m: [...], Tw_K: [...]}` (実寸、原点は設計スロート)、または `Tw_profile_csv: <カードからの相対パス>` | **長さを振るなら s** (分布が長さに追従)。表は全長を覆うこと (足りなければ止まる)。s の表は s = 0 を点として含み、直管長 > 0 (直管なしは x_m)。NS は `wall_profile_3.csv` (`wallProfile: 1`)、積分法は同じ表を実効の r_t で読む。run の `wall_temperature_profile.json` と CSV が正本 (書き換えると段の起動・成果物で止まる)。Fluent は `fluent/wall_temperature.prof` (壁面中心の値; Fluent での読み込みは未確認) |
| ガス | `gas: {preset: air}` か `gas: {species: {...}}` | 空気の M5 以上は出口の N₂ 液化を生成器が警告する |
| 寸法の許容差 | `tolerances: {length_m, diameter_m}` | 既定値は仮。ユーザの実用公差を聞いて書く (出荷ゲートの基準) |
| knot のマッハ数 M_knot | 書かない (既定 = 自動) か `overrides: {geometry.M_knot: auto}`。固定するなら `overrides: {geometry.M_knot: 2.1}` | **既定は自動** (2026-10-14): 長さを固定して最大壁角が最小の M_knot を、網の検査・余裕 ≥ 1°・成立の境目から 0.05 以上の範囲で選ぶ (plan `tooling-nozzle-mknot-auto`)。カードの生成に数十分かかる (M5 で 45 分以上) — **AWS で回す**。記録は YAML の `_card.M_knot_rule` (要約) と card_log |

流れ (正本は手順書 §1・§2。ここでは新しい入口と出口だけ):

1. `cd design && /home/sano/work/forge/design/.venv-opt/bin/python -m forge_design.intake.wt_card CARD --out ../case/NN/problem_<tag>.yaml`
   — カード → problem YAML。物理スロート径なら「r_t を解く → 長さを r_t で割り直す」を自動で繰り返す。
   標準出力の対応表・警告 (液化・格子の見積もり) をユーザに見せる。**格子の見積もりは初期値**で、合格の根拠は NS 後の `check_wall_resolution.py`。
2. 以降は手順書 §1 ③〜§2 (Euler → NS/δ* → C2 → 最終 NS)。**Euler は生成器が一緒に書く `<stem>_euler.yaml` で回す** (最終の雛形の
   `physical_wall_repr`・`pw_upstream` は NS 専用で、Euler の prepare が拒否する。写しは設計壁がビット同一)。CFD 前の物理壁とメッシュだけなら
   `deltastar_loop --problem YAML --euler-ref none_precfd --init-integral --prepare-only --run-dir RUN` (`FORGE_BIN`・`FORGE_CONVERTER` を新しいビルドに)。NS 後に r_t を解き直したら、同じカードに `--r-throat R` を付けて YAML を作り直す
   (長さの割り直し。r_t 単位の値を手で直さない)。
3. 最終 NS の run で `-m forge_design.export.deliverables RUN --card CARD --problem YAML` → `RUN/deliverables/` (README・MANIFEST・cad/・mesh/・fluent/・excel/)。
   **応答には `DIMENSIONS.json` の判定 (寸法の出荷ゲート) と、熱量の `thermal_use` の状態を必ず書く**。ゲートが FAIL なら成果物を渡さず原因を報告する。
4. Fluent での読み込みはこちらでは確認できない — 初回は「Fluent で `Mesh > Check` を見てほしい」とユーザに頼み、結果を plan §9 に記録する。

罠:
- **2D の Fluent メッシュの面の向きは「n0→n1 の左側が c0」**。3D の規則からの類推は逆 (2026-10-13 に一度誤った)。書き出し器を触ったら実物の 2D メッシュで測り直す。
- **軸対称の `iface_q_eff` は r·q** (熱流束ではない)。熱量は `iface_Qf_eff` × 2π から作る。`qwall` は診断値。
- **新しい設計の最初の NS は Euler の場から始められない** (2026-10-13): `--ic-from <Euler run>` は interp_field が化学種の照合で拒否する
  (Euler は輸送の記録が無い)。`--ic-from` を付けずに 1 次元等エントロピー場 + 段階起動 (`--stages full`) で始める。
  `FORGE_ALLOW_UNVERIFIED_SPECIES` で黙らせない。
- 熱量の列の状態が `valid` 以外、または `thermal_use` が `PASS` でない表を、1 次元熱伝導の入力として「使える」と書かない。
- STEP を CAD で編集した後のメッシュの作り直しは `deliverables/mesh/rebuild_from_step.py` (x(u) を解いて station に合わせる)。
  編集後の形状の格子品質は `check_mesh_quality.py` で確かめる。

## 計算の長さの既定 (2026-10-14)

- NS の本段は**停止規則つきで既定 20 万 step を上限**に回す (2026-10-14 ユーザ決定; `run_staged_ns`、`prepare_info.euler_ref` が要る)。2 万 step の窓で、
  入口−出口の質量の収支 r (|平均| ≤ 0.05 %・振幅 ≤ 0.5 %) と θ_r (設計壁の出口 x_w に対する x/x_w 0.4・0.7・0.95、平均の変化 2 万 step あたり ≤ 0.05 %・
  振幅 ≤ 0.5 %) を満たした後、**さらに 4 万 step その条件を保った出力で止める** (`evaluate.stop_rule_hold_steps`、既定 40,000; 途中で外れたら数え直す。
  上限 20 万は到達 + 保持を含む)。上限に届けば「未到達」(合格でない)。記録は run の `stop_rule.json`・`STOP_RULE_STATE.json`。
  キー無しで使えない (euler_ref が無い等) ときは理由を記録して `ns_steps` の長さで回し、`evaluate.stop_rule: true` を明示したときだけ準備の前に例外。
  使わないときは `--no-stop-rule` か `evaluate.stop_rule: false` で、本段は `--steps` > `evaluate.ns_steps` > 60,000。
  `evaluate.nStepOuter` は Euler の本段だけ。冷却壁の case/70 は入口の直管の質量の過渡で未到達になる (plan `tooling-nozzle-cooled-wall-ns-settling` §6 S-stop)。
  未到達の run から次の壁を作るのは `deltastar_loop --allow-not-reached "<理由>"` の明示の上書きだけ。窓の質量の収支・直近 2 窓の δ_E・δ_target の
  実測が閾値を超える run (case/70 の r 0.7〜1.4 % など) はさらに `--allow-unsettled "<理由>"` が要る (plan `tooling-nozzle-isothermal-wall-chain`
  §5.1 #36h)。理由は必須で、`HANDOFF_GATE.json` (version 3) に理由・実測の値と `settled: false` が残る。判定不能 (記録の欠け等)・被覆・held・次の壁・
  局所 y₁⁺ は上書きでは通らない。
- **計算は手元の PC でなく AWS で** (2026-10-14 ユーザ「ローカルで計算かけない方がいいかも」、手元の PC がメモリ不足で落ちた)。skill `forge-aws-run`。

## メッシュの規模 (作る前に示す; AGENTS.md「メッシュの規模は作る前に示す」)

- 要求カードの生成・`prepare`・`prepare_ns` の前に、**節点数** (Euler・NS。生成器は YAML のコメントと標準出力に `# 節点数` を出す) と時間の見込みを応答に書く。
  目安 (手元の GPU、2026-10-14): Euler 17.5 万点で 1 万 step 約 36 s (3.6 ms/step)、NS 15〜17 万点の本段で 1 万 step 約 80 s。
- 上限は 2D で **30 万点** (Euler・NS 共通)。超えると生成器・prepare が格子を作る前に止まる。超える格子が要るなら**ユーザの了承を取ってから** `overrides: {mesh.max_nodes: N,
  mesh.max_nodes_reason: "…", mesh_euler.max_nodes: N, mesh_euler.max_nodes_reason: "…"}` を書く。NS の ni は壁際の AR (≤ 0.97 × 5000) から決まり、第 1 層が薄い (Re が高い)
  と上限 6000 まで増える — まず第 1 層・`ar_max`・`throat_refine` を見直す。
- Euler の ni は NS と同じ (減らすと Euler の解が動く: case/70 で ni 600 の出口 M +3.8e-4、plan `tooling-nozzle-mesh-size-guard` §4.1)。明示の `mesh_euler.ni` は保持される。

## まず使うツール (新しく書かない; AGENTS.md「既存ツールの再利用を最優先」)

検証・判定・後処理のスクリプトを新しく書く前に、この表から探す。足りないときだけ書き、応答に「既存の何が足りなかったか」を 1 行書く。

| やりたいこと | 使うもの |
| --- | --- |
| 実寸の要求 → 問題 YAML (r_t の固定点・Euler 用の写し) | `python -m forge_design.intake.wt_card CARD --out P.yaml` (NS 後の δ_E で作り直すなら `--ns-correction RUN --delta-next CSV`) |
| Euler を回す | `python -m forge_design.evaluate.runner_axismach P_euler.yaml RUN --steps N --stages soft` |
| NS (CFD 前の積分法の壁) | `python -m forge_design.feedback.deltastar_loop --problem P.yaml --euler-ref EULER --init-integral --run-dir RUN [--ic-from RUN] [--stages full\|none] [--prepare-only]` (`interfaceDiag` を足す台本の例: `case/70.m5_air_wt_demo/run_ns_ifacediag.sh`) |
| NS (補正の表の壁) | `deltastar_loop --problem P.yaml --euler-ref EULER --delta-r-csv P_delta_r.csv --ic-from RUN --run-dir RUN` |
| 同じメッシュで続きを回す | `solver_density_cuda/tools/restart_field.py RES.h5 NEW/nozzle.h5` → `runner_axismach.run_staged_ns(RUN, "none")` (Euler は `run_staged`) |
| メッシュ品質・収束・準定常・壁解像 | `check_mesh_quality.py`・`check_convergence.py RUN --segment`・`check_quasisteady.py`・`check_wall_resolution.py RUN --mesh RUN/nozzle.h5 --groups wall` |
| 出口コア M・流量の系列と準定常 | `case/70.m5_air_wt_demo/qs_series.py RUN [--tail 0.4 --drift 5e-5 --osc 5e-5]` (出口 M の合否は専用の許容) |
| 残差の図 | `solver_density_cuda/tools/plot_residual.py RUN/residual_history.csv -o RUN/residual_history.png --all` |
| 標準の報告 (図・評価量・pptx) | `python -m forge_design.report.nozzle_report RUN --euler EULER` → `build_pptx.py RUN/report` (`.venv-pptx`) |
| 成果物一式 (STEP・gmsh・Fluent・Excel) | `python -m forge_design.export.deliverables RUN --card CARD --problem P.yaml` |
| δ* の抽出 (1 枚・指定のスナップショット) | `deltastar_loop.extract_and_merge(NS, EULER, band_select="edge", out_dir=...)`; スナップショットを指定する例 `case/70.m5_air_wt_demo/band_c_ab.py`・`band_persist_ab.py` |
| M_knot を自動で選ぶ (要求カードの既定) / 長さ固定で振る | 自動はカードの生成器 (`forge_design.intake.wt_card`、`evaluate/mknot_select.py`)。振るだけなら `case/70.m5_air_wt_demo/mknot_sweep.py PROBLEM_euler.yaml OUT.json [--n-axis N] M_K ...` |
| 遅い量 (熱量・排除厚さ・流量・出口 M) の時系列 | `case/70.m5_air_wt_demo/slow_mode_probe.py OUT.json RUN[:OFFSET] ...` |
| 寸法の逆算の A/B・壁温の δ_r → 物理壁の A/B | `case/70.m5_air_wt_demo/sizing_ab.py`・`walltemp_wall_ab.py` (手本) |

## 最初に決めること

1. **どのユースケースか** を確定する (§1 / §2 / §3)。曖昧ならユーザに 1 問だけ聞く。
2. **新規設計 (§1・§2) で problem YAML が未確定なら `/design-intake` に渡す** (要件インタビューはそちらの担当)。
3. **最新の実例を手本にする**: `case/45.isobutane_m6_d155` (M6、2026-10)。生産レシピは
   [`plans/accepted/tooling-nozzle-cfd-pinned-initial-line.md`](../../../plans/accepted/tooling-nozzle-cfd-pinned-initial-line.md) §4・§5.1 #11〜#11f と §9 末尾。
   古い run の YAML をそのまま複製しない (試作の `M_design` 較正値・廃止キーが残っている)。

## 判断の分かれ目

- **既存ノズルで条件 (Pt・Tt・H2O 分率・背圧) だけ変えたい** → 正本 §3a: `solver_density_cuda/tools/rerun_conditions.py REF NEW ...` (形状固定、forge は起動しない)。Pt を変えるなら `--Ps`/`--keep-Ps` 必須、回し方は `RERUN_CONDITIONS.json` の `recommended_stages` (Pt 変更は full + 本段 cfl 1 + `--scale-ic pt` を推奨; scale なしは入口壁際に逆流域が残り準定常に達しなかった)。
  **problem YAML の spec/gas を変えて runner/deltastar_loop に通すと壁が作り直される** (それは §1・§2 の設計作業)。
  組成の lump を変えると化学種の定義が変わり restart_field が拒否する → convert_species_field か等エントロピー初期場。

- **探索段階** (軸長・R・M_knot を振る) は Hall 初期線 + 補間壁 (YAML の既定) と粗い格子で十分。
  **形状を詰める段階**で `wall_repr: joint`・`initial_line: cfd`・`Md_moc_offset` を入れる (ユーザ決定 2026-10-05)。
- 壁の評価観点 (ユーザの方針): **壁がなめらか (r″ が連続) / 試験部に圧力波が無い / オーバーシュートが小さい** (凝縮を強めるので)。
  軸 M は形状の生成器にすぎない。Pa 単位の波の追い込みはしない (Euler の数値床と同じ桁)。
- 試験部・膨張部 (MOC の領域) に Euler のずれを返さない。CFD に合わせるのは遷音速のスタート部 (初期線) と出口の 1 係数だけ。
- **出口較正 (`Md_moc_offset`) は生産 NS と同じ格子パラメータの Euler で決める** (粗い Euler 格子だと NS で出口 M が約 −0.0008 ずれた)。
- NS の格子は**壁解像 PASS** (y1+ ≤ 1、`check_wall_resolution.py --over-frac 5`) と**格子ゲート** (δ_E(x_F) の細分前後差 ≤ 1 %) を通したものを使う。
  細分格子は本段 CFL 5 で発散する → 段階起動 + cfl 1・60000 step。

## 落とし穴 (実際に踏んだもの)

- `delta_r_next.csv` の δ は **`delta_E` 列 (未緩和)** を使う。第 2 列 `delta_target` は ω=0.5 緩和後 (取り違えて k_f・r_t を誤って解いた)。
- 出口径は spec に書けない。`solve_rt` で r_t を解いて `spec.r_throat` を手で書き換える。r_t を変えても `r_inlet`・`L_total` は r_t 単位のまま。
- `design_chain` は fold・特性線トポロジ・壁マージンを検査しない → パラスタでは `geometry/moc_diagnostics.py` を併用する。
- Euler の `prepare` と NS の `prepare_ns` は同じ mesh キーを読む (2026-10-06 に `mesh_params()` で統一。以前の Euler は一部のキーを黙って無視していた)。
- 同じ格子の続きは `restart_field.py`、格子が違えば `interp_field.py`。凝縮 ON への引き継ぎは `convert_species_field.py --mode conserve` (restart_field は拒否する)。
- case スクリプトに `/home/sano/work/forge/...` の絶対パスが書かれているものがある (worktree で走らせると別ツリーを読む)。使う前に確認する。
- AWS の共有バイナリは他セッションが作り直すことがある。自分の worktree でビルドし、`RUN_PROVENANCE.txt` の sha256 を記録する。
- 判定は判定ツールの VERDICT で行う。残差 plateau (`NOT CONVERGED`) はこの系の常態 — そのまま書き、**報告する量の時系列**に
  `check_quasisteady.py --series-csv` をかける (`nozzle_report` の末尾幅は VERDICT の代わりにならない)。
- 判定基準は回す前に plan に書く (事前登録)。外れたら閾値を動かさず、延長 1 回 → なお未達なら上位の判断役に諮る (AGENTS.md「モデル分担とエスカレーション」)。

## 未実装 (ユーザに聞かれたら正直に言う)

汎用の MOC パラスタ CLI・凝縮引き継ぎの自動化・外部壁 CSV の取り込みは無い (case スクリプト頼み)。壁の形状表・STEP・メッシュの書き出しは
`forge_design.export.deliverables` (上の節)、k_f は YAML の `deltastar_initializer.cf_scale` に書ける。
統一は [`plans/active/tooling-design-problem-campaign-recipe.md`](../../../plans/active/tooling-design-problem-campaign-recipe.md) の残作業。

## 終わったら

`nozzle_report` → pptx ([`procedures/nozzle-design-outputs.md`](../../../procedures/nozzle-design-outputs.md))、case README の run 一覧、plan §9 に記録。
