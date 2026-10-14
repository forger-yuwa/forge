# 諮問ブリーフ: 風洞ノズルの実寸カードと成果物一式 (plan tooling-nozzle-wt-deliverables §4・§6)

AGENTS.md「モデル分担とエスカレーション」条件 1 (新規 plan の §4 設計方針・§6 検証計画)。判断役のスイッチが codex。

## 背景 (ユーザの要求, 2026-10-12)

- 入力: 全圧・全温を決め、入口半径・スロート径・出口マッハ 5・入口長さ・出口長さ・壁温 350 K で風洞ノズルを作る。
- 欲しいもの: 解析メッシュで Fluent 計算をするので gmsh の入力、ノズル形状 (STEP)、**STEP から gmsh で全く同じ構造のメッシュ**、
  Excel でノズル形状を扱う方法 (1 次元の壁熱伝導を評価するかも)。
- 「いろいろなマッハ数・形状で実施するので準備を整えてほしい。Claude に作成を指示させる可能性が高い」。
- ユーザ回答: 出口長さ = スロート→出口の長さ (既存 `Lc_mode: from_length` + `dv.L_total`)、入口長さ = 直管長 `L_pipe` と縮流部長 `L_U` の両方。

## 観測事実 (出典つき)

1. 設計チェーンの入力は $r_t$ 単位 (`procedures/nozzle-design-workflow.md` §0 の表、落とし穴「r_t を変えても r_inlet・L_total は r_t 単位のまま (実寸がずれる)」)。
   物理スロート径から $r_t$ を解く `solve_rt_throat` (`design/forge_design/feedback/deltastar_loop.py:288`)。
   `Lc_mode: from_length` は x_F(L_c) = L_total を固定点反復で解き、既定許容差 0.05 r_t (`design/forge_design/evaluate/runner_axismach.py:562-620`)。
2. メッシュ `design/forge_design/meshing/mesh2d.py:143-258`: 格子線は $x=x_i$ の縦線、$r_{ij}=r_w(x_i) s_{ij}$、msh4.1 を `.10g` で書く。
   別ブランチ (未マージ) には壁法線の近壁層 `wall_normal_layer` がある (縦線でない格子)。
3. STEP: `design/forge_design/export/wall_step.py` (物理壁 `single_bspline` の 1 本の 5 次 B-spline、助変数 u = x [mm]、FreeCAD で書き読み直し、転送誤差 位置 1e-6 mm)。
   3D 回転面は読み直し検査 (`revolve`) だけで STEP には書いていない。
4. Fluent の面の向き: 実物の Fluent メッシュ 2 本 (3D) で、節点順の右手系の法線 (Newell) が c0 側を向く割合を測った:
   `case/32.fan_duct_fluent/mesh/fan.msh.h5` 内部面 24,542/24,542 が c0 向き・境界面 1,484/1,484 が内向き、
   `case/30.mixing_elbow_fluent/mesh/mixing_elbow.msh.h5` 86,782/86,782・4,799/4,799。
   読み込みは `solver_density_cuda/tools/fluent_h5_to_forge.py` の `FluentMesh` (節点順は変えない、232-270 行)。
   2D は「辺 n0→n1 を +z に押し出した面の法線 (n1−n0)×ẑ が c0 を向く = 進行方向の右側が c0」と推定した (2D の実物は手元に無い)。
5. 壁面量: `res_wall_<physID>_*.h5` に Ps・Ts・qwall (流体向き正、診断)・twall_{x,y,z}・utau・ypls。`output.interfaceDiag: 1` で iface_q_eff (正本、固体向き正) ほか
   (`methods/boundary.md` 400-410 行付近と「界面熱量の 4 つの定義」節)。現行の設計チェーンの NS は interfaceDiag を付けていない。
6. 手元に Fluent は無い。ParaView 6.1.1 (pvpython, Fluent 読み込み器あり)、gmsh 4.15.2 (msh・CGNS・UNV・BDF・STEP を書ける)、FreeCAD 1.1.1、
   `.venv-pptx` に xlsxwriter がある。

## 計画の要点 (plan §4)

- §4.1 実寸カード → YAML。物理スロート径のとき「solve_rt_throat → 長さを r_t で割り直す」を |Δr_t|/r_t ≤ 1e-6 まで反復。y1+ ≤ 0.8 の第 1 層を 1 次元等エントロピー + Eckert 参照温度 + 平板 Cf で見積もる (合格は NS 後の check_wall_resolution)。
- §4.3 構造ファイルを run の nozzle.msh から作り (メッシャを再実装しない)、STEP の曲線を u=x_i で評価して r_ij = s_ij r_w(x_i) で作り直す。縦線でなければ拒否。
- §4.4 Fluent ASCII .msh の自前書き出し (c0 = 右側) + gmsh 経由の CGNS。検証は自前の面積・閉路検査と ParaView の読み直し。Fluent での確認はユーザ。
- §4.5 Excel: ASCII 列名の CSV。壁への熱流束は iface_q_eff、無ければ −qwall を「診断」と明記。

## 問い

**この計画のうち、実装後に「成果物は出たが使うと誤る」形で失敗する可能性が最も高い設計判断を 1 つ挙げ、それを防ぐ変更を 1 つ推奨してほしい。**
候補として当方が気にしているのは (a) r_t の固定点と L_total (from_length の許容差 0.05 r_t、NS 後の C2/k_f の解き直しで長さが動く) の連成、
(b) 2D の面の向きを 3D の実測から推定していること (V4 の検証で十分か)、(c) Fluent に渡す熱流束・境界条件表が forge と別物になる危険、
(d) 壁解像の事前見積もり — だが、これに縛られなくてよい。あわせて、§6 の合格条件で**検証になっていないもの** (自己一致しか見ていない等) があれば指摘してほしい。

## 禁止事項 (厳守)

- ファイルを変更しない (read-only サンドボックスで動いている)。
- **`*.log`, `residual_history.csv`, `res_*.h5`, `*.vtu`, `*.msh.h5`, `plans/README.md` を読まない** (巨大)。
- 下に列挙した `sed -n 'A,Bp' <file>` 以外のファイル読みをしない。grep は可。
- 両論併記で逃げず、**推奨は 1 つに絞る**。根拠は `ファイル:行` か提示した数値で示す。

## 読んでよいもの

- `sed -n '1,400p' plans/active/tooling-nozzle-wt-deliverables.md`
- `sed -n '1,120p' procedures/nozzle-design-workflow.md`
- `sed -n '143,258p' design/forge_design/meshing/mesh2d.py`
- `sed -n '1,60p;138,152p;290,354p' design/forge_design/export/wall_step.py`
- `sed -n '200,360p' design/forge_design/feedback/deltastar_loop.py`
- `sed -n '455,625p' design/forge_design/evaluate/runner_axismach.py`
- `sed -n '395,450p;600,680p' methods/boundary.md`
- `sed -n '1,80p' solver_density_cuda/tools/check_wall_resolution.py`
- `sed -n '200,270p' solver_density_cuda/tools/fluent_h5_to_forge.py`
