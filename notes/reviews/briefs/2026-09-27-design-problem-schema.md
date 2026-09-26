# 諮問ブリーフ: 設計ツール (forge_design) の問題定義・設計変数・詳細解析指定をどう構造化するか (2026-09-27)

## 問い (1 つ)

forge_design を「ユーザと対話しながら任意の機種 (風洞 axismach / ベル+スラスタ / デュアルベル+スラスタ /
SERN / 将来は熱制約・FEM 連成) の最適化フローを組める」基盤にしたい。
**problem YAML・設計変数 (dv)・最適化後の詳細解析の指定を、どういう構造に寄せるべきか。推奨を 1 つに絞って示してほしい。**

具体的にユーザが迷っている論点 (これに縛られなくてよい):

- (a) problem を簡素化して axismach に絞る / (b) problem を機種ごとに色々な形に対応させる /
  (c) problem は薄い「契約」だけにして、機種ごとの中身はプラグイン (runner) 側が持つ / (d) かっちり決めない (スクリプト駆動) 方がよい
- dv (設計変数・範囲) をどこで指定するか (problem YAML 内 / キャンペーン定義の別ファイル / 最適化スクリプト内)
- 最適化後の詳細解析 (NS δ* 補正・凝縮・入口分布・3D・FEM 等) をどう指定するか
- 熱的制約・FEM 連成など「評価の段」が増えるときに破綻しない形は何か

回答に欲しいもの: 推奨する構造 1 つ (スキーマの骨格を YAML 例で)、その理由 (下の観測事実の `ファイル:行` に紐づけて)、
**最初の一歩として何を 1 つだけ変えるべきか** (大改修の前に効果を確かめられる最小の変更)、捨てるべきもの。

## 観測事実 (2026-09-27 のセッションで確認)

1. **problem YAML は 1 ファイルに「物理条件・形状パラメータ・メッシュ・評価設定・dv」を全部持つ**。
   例 `case/44.vitiated_air_wt/problem_va3_M4.19_Lc8_noneq_lumpX.yaml` (55 行): `gas` / `spec` / `dv` / `geometry` / `mesh` / `evaluate`。
   `type:` で runner を選ぶ (`probdef.py` の `KNOWN_TYPES`、`design/CAPABILITIES.md` §1)。
2. **semiperfect でも `gas.gamma` / `gas.cp` を書かせている**。熱力学の本体は NASA-9 だが値は次で読まれる:
   forge `physProp` への素通し (TP では無効)、`collect` の出口一様性診断 `exit_uniformity(..., gamma=p.gamma)`
   (`runner_axismach.py:581-582`、定数 γ で特性線追跡)、NS 経路の ω 床の R 概算 (`runner_axismach.py:821-822`)、
   `evaluate.cfd_gas: cpg` のときだけ実際の CPG 定数。**省略すると黙って γ 1.4 / cp 1004.5** (`probdef.py:167-168`)。
3. **`dv.L_c.min/max` は axismach では探索に使われない**。読むのは `probdef.py:221-226` の「value が [min,max] 内か」検査と、
   MOO ドライバ (`opt/driver.py:77-78`, `opt/driver_sern.py:69`) だけ。`runner_axismach` は設計チェーン自身が計算した許容窓
   (`runner_axismach.py` の `design_chain` 内 `lo <= L_c <= hi`、va3 では `Lc_window [0.01, 10.94]`) で検査する。
   YAML の `max: 14.0` は設計上ありえない値のまま放置されていた。
4. **MOO ドライバは機種ごとに別実装**: `opt/driver.py` は `from ..evaluate import runner` (ベル専用, `driver.py:40`) で
   DV_ORDER を固定、`opt/driver_sern.py` は SERN 専用 (`driver_sern.py:41` DV_ORDER = M_c, f, θr0, θc0, L_cowl)。
   axismach (風洞) は MOO に未接続。length-dv plan (`plans/accepted/tooling-nozzle-axismach-length-dv.md:30-33`) は
   「ドライバは dv dict に対して汎用だから YAML を書けばそのまま乗る」と書くが、コード上は未接続 (driver.py が runner を固定 import)。
5. **風洞の L_c 最適化は problem YAML を読まない使い捨てスクリプトで行われた**: `case/42.isobutane_wt/optimize_axislaw_A_shortest.py`
   → `sweep_axislaw_A_MK.py` → `compare_axislaw_D.py` と import 連鎖し、条件 (R=3, M_d=6, Tt=1550 K, 組成) も探索範囲
   (L_c 35–40, M_K 2.2–3.2) もスクリプト内定数。設計段のみ (CFD なし) の 600/1200 点探索。
   ほかに設計チェーン内の自動決定 `Lc_mode: max` (窓上限×0.98) / `from_length` (全長から逆算) がある。
   case/44 の R×L_c×L_U 選定は problem YAML を並べて `run_study_va.py` で一括評価するスイープだった。
6. **Euler か NS (粘性 δ* 補正) かは YAML で表現されず、エントリポイントで決まる**: `python -m forge_design.evaluate.runner_axismach`
   の CLI (`runner_axismach.py:608-634`) は常に Euler (`prepare`: `_config_euler_node` + `_bcond(p, euler=True)` = slip 壁)。
   粘性経路 `prepare_ns` / `run_staged_ns` は CLI に出ず `forge_design.feedback.deltastar_loop` からだけ呼ばれる
   (`deltastar_loop.py:138-192`)。δ* 生産レシピは CLI 引数の組 (`--init-integral --euler-ref --ic-from --stages --cfl --implicit-relax`)。
7. **runner の段階起動が規約とずれていた**: `run_staged` は同一メッシュの段間引き継ぎに `interp_field.py` (cross-mesh 用) を使う
   (`runner_axismach.py:538-540`)。AGENTS.md は同一メッシュは `restart_field.py` と定めている。今回の再計算
   (case/44 run_0509–0511) では runner を `--prepare-only` で止め、段階起動を case 内の自作スクリプト `run_lumpX_staged.py` に置き換えた。
   入口組成を bcond の lump モル分率 `X0/X1` に置換するのもこのスクリプトで行った (runner は常に質量分率 `Y{s}` を書く, `runner_axismach.py:175-180`)。
   → **「決められた runner では足りず、case ごとに段や前後処理を差し替えたい」需要が既に発生している**。
8. **推奨数値設定の正本は別文書** (`procedures/recommended-settings.md`) で、runner は config テンプレートを内蔵しており
   現行レシピ (cfl 6 + implicitRelax 0.7, nStepInner 4) とずれる (runner 既定 cfl 4 / nStepInner 5; CLI で上書き)。
9. 機種の現状 (`design/CAPABILITIES.md` §1): `thruster_bell` (MOO 検証済), `wind_tunnel_axisym_axismach` (生産),
   `sern_2d` (多作動点 MOO 検証済), 旧 `wind_tunnel_axisym` / `walldriven` (回帰対照)。デュアルベル・熱構造 FEM・CHT は 📋 未実装。
   インテーク手順は `.claude/skills/design-intake/SKILL.md` (対話で problem YAML を「契約書」として確定し、キャンペーンは決定的に回す)。
10. ユーザの希望: 最適化フローを「ユーザの意見を聞きながら好きに組んでいく」。対話型で都度組み替えたいが、再現性 (ledger・VERDICT) は保ちたい。

## 期待値と出典

- 収束・準定常・メッシュ品質・発散は判定ツールの VERDICT をファイルで残す契約 (`design/CAPABILITIES.md` §5, AGENTS.md)。
  どの構造にしてもこの契約は維持する前提。
- ノズル設計の仕様・チェーンの説明は `methods/design/overview.md` (§問題定義 YAML は 48-59 行)。

## 実施済みの操作と結果

- 2026-09-27: case/44 va3 M4.19 L_c8 を lump モル分率入力で dry / noneq / eq 再計算 (`run_0509`–`0511`, commit 070d9003)。
  runner は prepare のみ使用し段階起動は自作スクリプト。全 NaN 0・残差 plateau (NOT CONVERGED)・報告量 ALL STEADY。
  この作業でユーザが上の 2〜6 に気づき、設計ツールの構造を見直したいとなった。
- まだ plan は起票していない (本諮問の結論を受けて起票する)。

## 仮説 (当方の見立て。棄却してよい)

- H1: problem YAML を「何を作るか (機種・仕様・物理)」と「どう探すか (dv・目的・制約・探索法)」と
  「どう評価するか (段の並び: 設計段 → Euler → NS δ* → 凝縮 → FEM …、各段の VERDICT ゲート)」の 3 層に分け、
  評価は段 (stage) の DAG/リストとして宣言し、段の実体は機種非依存のプラグインにするのがよい。
- H2: 逆に、いま汎用スキーマを先に作ると機種ごとの差 (SERN の多作動点、風洞の軸 M 則、ベルの TOP) を吸収しきれず
  過剰設計になるので、「機種ごとの runner + 共通の小さな契約 (dv 宣言・評価段の入出力ファイル・VERDICT)」に留めるべきかもしれない。
- どちらが良いか当方では決めきれていない。

## 禁止事項 (厳守)

- ファイルを変更しない (read-only サンドボックスで動いている)。
- **`*.log`, `residual_history.csv`, `res_*.h5`, `*.vtu`, `*.json` (大きいもの), `plans/README.md` を読まない** (巨大)。
- 下に列挙した `sed -n 'A,Bp' <file>` 以外のファイル読みをしない。grep は可。
- 両論併記で逃げず、**推奨は 1 つに絞る**。根拠は `ファイル:行` か本ブリーフの番号で示す。

## 読んでよいもの

- `sed -n '1,237p' design/forge_design/probdef.py`
- `sed -n '1,85p' design/CAPABILITIES.md`
- `sed -n '1,60p' methods/design/overview.md`
- `sed -n '1,126p' .claude/skills/design-intake/SKILL.md`
- `sed -n '1,55p' case/44.vitiated_air_wt/problem_va3_M4.19_Lc8_noneq_lumpX.yaml`
- `sed -n '1,44p' case/40.nozzle_design_tool/problem_bell_rao_rd15.yaml`
- `sed -n '1,117p' case/46.sern_design/problem_3d_prod_m6on_g1.yaml`
- `sed -n '86,183p' design/forge_design/evaluate/runner_axismach.py`
- `sed -n '435,560p' design/forge_design/evaluate/runner_axismach.py`
- `sed -n '600,660p' design/forge_design/evaluate/runner_axismach.py`
- `sed -n '1,120p' design/forge_design/opt/driver.py`
- `sed -n '1,110p' design/forge_design/opt/driver_sern.py`
- `sed -n '1,60p' design/forge_design/feedback/deltastar_loop.py`
- `sed -n '130,200p' design/forge_design/feedback/deltastar_loop.py`
- `sed -n '1,40p' case/42.isobutane_wt/optimize_axislaw_A_shortest.py`
- `sed -n '1,104p' case/44.vitiated_air_wt/run_lumpX_staged.py`
- `sed -n '15,45p' plans/accepted/tooling-nozzle-axismach-length-dv.md`
