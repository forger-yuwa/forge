# 諮問ブリーフ: 設計チェーンへの壁温分布 T_w(x) の配管 (plan tooling-nozzle-wall-temperature-profile §4・§6)

AGENTS.md 条件 1 (新規 plan の §4・§6)。判断役のスイッチが codex。ユーザ指示 (2026-10-13): 壁温分布を入れられるように、設計開始時に尋ねるように、デモする。

## 観測事実

- ソルバの `wall_isothermal` + `ints: {wallProfile: 1}` + `wall_profile_<physID>.csv` は実装・検証済み (case/48 `run_0015_wallprofile` で場の壁ノードの T が入力と最大差 0.068 K)。
  CSV は先頭の座標列 + `Ts`。座標列 1 つなら 1D 線形補間・範囲外は端値 (`procedures/inlet-profile.md` 25 行・`solver_density_cuda/boundaryCond.cpp:471`)。
- 積分法 `integral_bl` は `thermal_bc.Tw_table [[x_rt, Tw]]` に対応 (`design/forge_design/feedback/deltastar_integral.py:132-140`)。
- 設計チェーンの単一ソース `probdef.Problem.wall_thermal` は一様のみ (`design/forge_design/probdef.py:60-95`)。
- CONTUR の熱閉包は温度形で冷却壁の δ_r が最大 2.8 % 変わりうる既知の危険 (isothermal plan §5.1 #10〜#12)。
- case/70 (M5 空気、Tt 1250 K、Pt 3 MPa、等温 350 K) は C2 方式で NS まで回っている (case/70 README)。

## 計画の要点

plan §4.1〜4.4 を読むこと。表の x は実寸 [m]・原点は設計スロート (要求カード・STEP と同じ)、壁の全長を覆わなければ例外、
bcond の `Ts` の既定値は x=0 の値、一様・断熱はバイト同一。

## 問い

**この配管のうち、実装後に「分布を与えたつもりが NS と積分法・成果物で別の壁温になっている」形で失敗する可能性が最も高い点を 1 つ挙げ、防ぐ変更を 1 つ推奨してほしい。**
候補 (縛られなくてよい): (a) x の原点・単位 (実寸 + 設計スロート原点) と、物理壁 (δ_r を足した壁) の x の対応、(b) 壁の端 (入口・出口の角の節点) の扱いと
ソルバの端値クランプ、(c) `run_staged_ns` の段の切り替え・restart・`rerun_conditions` で CSV が落ちる/古いまま残る、(d) V2 の合格条件 (場の壁ノードの T と表の線形補間の最大差 ≤ 0.1 K) が検証として十分か。

## 禁止事項 (厳守)

- ファイルを変更しない (read-only サンドボックスで動いている)。
- **`*.log`, `residual_history*.csv`, `res_*.h5`, `*.vtu`, `*.msh`, `plans/README.md` を読まない** (巨大)。
- 下に列挙した `sed -n 'A,Bp' <file>` 以外のファイル読みをしない。grep は可。
- 両論併記で逃げず、**推奨は 1 つに絞る**。根拠は `ファイル:行` か提示した数値で示す。

## 読んでよいもの

- `sed -n '1,200p' plans/active/tooling-nozzle-wall-temperature-profile.md`
- `sed -n '55,100p' design/forge_design/probdef.py`
- `sed -n '120,145p;190,235p' design/forge_design/feedback/deltastar_integral.py`
- `sed -n '260,330p;440,480p' solver_density_cuda/boundaryCond.cpp`
- `grep -n "wall_bcond_line\|wall_thermal\|bcondConfig" design/forge_design/evaluate/runner_axismach.py | head -30`
- `grep -n "def run_staged_ns" -A60 design/forge_design/evaluate/runner_axismach.py | head -90`
- `sed -n '575,612p' methods/boundary.md`
