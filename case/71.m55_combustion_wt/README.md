# case/71 — M5.5 の燃焼ガスの風洞ノズル (要求カードから NS まで)

2026-10-14 のユーザ依頼:「mach 5.5、入口配管長 0.5 m、収縮部 0.5 m、スロート曲率 5 mm、スロート直径 Φ200、スロート〜出口 8 m のノズルを作って、
凝縮計算はせずとも、NS 計算までやりきる」。うまくいったら、何も知らない低級なモデル (Haiku) に同じことをさせ、つまずいた所を流れの側で直す。

## 要求 (カード `card_m55.yaml`)

| 項目 | 値 | 出所 |
| --- | --- | --- |
| 出口マッハ数 | 5.5 | ユーザ指示 |
| ガス | イソブタンの燃焼ガス (当量比 0.9、case/45 と同じ組成: 質量分率 CO₂ 0.1677・H₂O 0.0858・O₂ 0.0220・N₂ 0.7245) | ユーザ選択 (提案値) |
| 全圧・全温 | 3.5 MPa・1,250 K | ユーザ指示 |
| スロート | 物理の直径 0.200 m、曲率半径 0.500 m (R ≈ 5) | 直径はユーザ指示 (物理)、曲率は**暫定** (ユーザ「500 mm かも、自信なし」; 5 mm は R = 0.05 で拒否される) |
| 入口 | 半径 0.300 m (提案値、ユーザ了承)、直管 0.50 m、縮流部 0.50 m | 長さはユーザ指示 |
| スロート〜出口 | 8.0 m | ユーザ指示 |
| 壁 | 等温 350 K | ユーザ指示 |
| 許容差 | 雛形の仮の値 (長さ 1 mm・径 0.02 mm) | ユーザ「よくわからん」 |
| 凝縮 | 計算しない (dry の NS まで) | ユーザ指示 |
| 格子の節点数 | NS・Euler とも ni 4,800 × nj 97 = 465,600 (手組みの格子 `problem_m55_ns.yaml`)。#3b (最終の雛形) から NS は ni 4,800 × nj 92 = 441,600 (径方向 `axis_cap_frac 0.03`・比 ≤ 1.2、ユーザ 2026-10-15「上限3%」; `problem_m55_final.yaml`)、Euler は 465,600 のまま (`problem_m55_final_euler.yaml`) | ユーザ決定「a」(2026-10-15): 見積もりの偏りを M6 の実績で補正した約 45 万点で Euler → NS (上限 30 万を超えることは了承済み)。それ以前に ni 10,200・約 99 万点を一度了承 → 同日取り下げ (「あまりにメッシュ多すぎる、古いやり方で作ってない？」) |

## 計算 run 一覧

| run | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_euler` | Euler (slip) の参照 (AWS A)。`problem_m55_final_euler.yaml` (最終の雛形の Euler の写し)、格子 ni 4,800 × nj 97 = 465,600 (全域 0.005)、`runner_axismach` の `--stages full` (soft 3,000 → mid 3,000 → 本段 cfl 1・12,000)、FP64 のビルド (`~/forge-integ-dflow2-fp64`、2fd36b0c)・`FORGE_CUDA_BLOCKSIZE=128` | 完走 (227 s): 出口の軸の M 5.5002、dM 最大 0.06 % Md、NaN なし。`check_convergence --segment` NOT CONVERGED (残差 2.1〜2.5 桁で停滞; M6 の Euler の参照と同じ型 → δ* の参照としては未適格と記録、plan §5.1 #3b)。δ* の参照は `res_12000.h5` (sha256 93cb88b0…) に固定。記録 [`_runs/run_0001_euler/`](_runs/README.md) | ref (run_0002_ns・run_0003_ns_cont の Euler の参照) |
| `run_0002_ns` | NS pass 0 (AWS A): `deltastar_loop --init-integral --stages full` (CONTUR の積分法の初期壁、停止規則・上限 200,000)、`problem_m55_final.yaml` (手組みの格子 ni 4,800 × nj 92 = 441,600、`axis_cap_frac 0.03`)、Euler の参照 `run_0001_euler`、FP64 のビルド | 本段 0〜183,998 step でディスクが満杯 (OSError) になり終了、停止規則は未到達。本段の `check_convergence --segment` (manifest の写しに本段を足した事後復元) NOT CONVERGED (各列 158k〜172k で最小の後に上昇)、停止規則の系列 DRIFTING (θ_r 0.4・0.7 TRANSIENT-UNSETTLED)、壁解像 PASS (最大 1.004)、res_183000 で P0 > Pt が 50,441 点・最大 +0.24 %。全圧の A/B (#3c) と診断は plan §5.1 #3b・#3c。記録 [`_runs/run_0002_ns/`](_runs/README.md) | ref (run_0003_ns_cont の初期場 res_183000) |
| `run_0003_ns_cont` | NS の延長を全部 double で (AWS A、plan §5.1 #3d): `run_0002_ns/res_183000.h5` を同じ格子に index restart (`restart_field.py --keep-src-dtype`、9 量ビット一致)、前段なし (`run_staged_ns` stages none)。本段の設定から `physProp.thermoFloat: 0`・`time.deltaT.implicitSolvePrecision: 1`・`outStepInterval` 2,000 だけを変え、環境変数 `FORGE_TRANSPORT_TABLE=0`・`FORGE_DIAG_FACE_H_DOUBLE=1` (`RUN_PROVENANCE.txt` に記録)。停止規則は明示 (窓 2 万・保持 4 万・上限 = 追加 200,000 step)、Euler の参照は `run_0001_euler/res_12000.h5` (sha256 93cb88b0…) に固定、帯はこの run の最初の出力 res_2000 から (θ_r は run_0002_ns と別の系列)。残る float: 化学種の拡散の D_s・h_s・X の float の写し (`speciesTransport_d.cu:240-356`)・リミッタの長さの `sqrtf`・XMF の型の宣言 | 投入 2026-10-15 00:24 UTC、26.9 ms/step。**判定 (本段 0〜199,999)**: 停止規則 CENSORED (200,000; θ_r 0.4・0.7 の D 6.3e−4・5.3e−4 > 5e−4、単調に減少)、`check_convergence` NOT CONVERGED (ρ・k・Y falling、運動量・ρE flat; float の末尾より 1.6〜3.2 桁下)、θ_r 3 断面・ṁ・出口コア M (5.52045) STEADY、壁解像 PASS (y1+ 最大 1.002、超過 0 %)、P0 > Pt 46,936 点・最大 +0.462 %、δ の 2 窓 FAIL (measured、δ_target の振れ 0.031 mm = δ の 0.3〜0.5 %、x 1.5〜2.2 m)。ユーザ決定「a」で関門の明示の上書きで δ* の補正 (plan §5.1 #3e) へ。判定の原文 `_runs/run_0003_ns_cont/judge/`、記録 [`_runs/run_0003_ns_cont/`](_runs/run_0003_ns_cont/README.md) | ref (δ* の補正の 1 回目の入力) |
| `run_0004_euler` | δ* の補正の 1 回目の Euler の参照 (AWS A、plan §5.1 #3e (4)): `problem_m55_p1_euler.yaml` (`wt_card --ns-correction run_0003_ns_cont` で物理スロート径 200 mm に r_t を解き直したカード `problem_m55_p1_card.yaml` → `_mesh/hand_mesh.py build`、r_t 0.10024473616136853 m = 前の +0.34 %、M_knot 2.80)、格子 466,570 点 (品質 PASS、AR 最大 ≤ 1000)、`runner_axismach` の `--stages full` (本段 cfl 1・12,000)、FP64 のビルド。台本 `_runs/run_0004_euler/p1_euler.py`・`_runs/run_0005_ns/p1_chain.sh` | 完走 (228 s、非有限なし): 軸の出口 M の最大 5.5013、dM 最大 0.058 % Md。`check_convergence` NOT CONVERGED (1.2〜1.8 桁で横ばい; run_0001_euler と同じ型)。参照は `res_12000.h5` (sha256 78294389…) | ref (run_0005_ns の Euler の参照) |
| `run_0005_ns` | δ* の補正の 1 回目の NS (AWS A、plan §5.1 #3e (5)): run_0003_ns_cont の δ_E から作った補正の表 `problem_m55_p1_card_delta_r.csv` の物理壁、`problem_m55_p1.yaml` (手組みの格子 442,520 点、品質 PASS・AR 緩和 違反 0)、初期場は run_0003_ns_cont の最後の場から補間、`--stages full`、**精度は既定 (float; ユーザ決定 2026-10-15)**・FP64 のビルド、出力 2,000、停止規則 (窓 2 万・保持 4 万・上限 200,000)、Euler の参照 `run_0004_euler/res_12000.h5`。台本 `_runs/run_0005_ns/p1_ns.py`・`p1_follow.sh` (終わったら成果物を作る) | 起動 2026-10-15 05:29 UTC。S1 (3,000 step) は非有限なし | 走行中 |

- **AWS A の途中の全場の削除 (2026-10-15、判定済み)**: `run_0003_ns_cont` の 96 枚 (res_0・140000・160000・180000・200000 は #3f の再生用に残す)、`run_0002_ns` の 21 枚 (res_0・res_183000 = run_0003 の初期場を残す)、`run_0001_euler`・`run_0004_euler` の各 11 枚 (res_0・res_12000 = Euler の参照を残す)。壁・出口の出力・残差履歴・判定の記録は残した。

## 経緯

- 2026-10-14〜15: 旧の M_knot の規則 (`min_theta_max_v1`) でカードを生成 (AWS A `~/forge-m55`、約 98 分) → 最終の組の 4800 の網の長さ合わせが全点未収束で失敗 (plan `tooling-nozzle-mknot-fast` §3)。
- 2026-10-15: 新の規則 `coarse_min_theta_v2` (feature/mknot-fast) で M_knot 2.80 (144 s) → NS の格子の節点数の上限 (30 万) で停止 → ユーザ了承で ni・上限を明示。
- 2026-10-15: ユーザが 99 万点の了承を取り下げ。今のカードの格子の見積もりは第一層を 2 値で持つ古いやり方で、拡大部全体が最も厳しい場所の薄い第一層になる。場所ごとの第一層 (`wall_first_frac_table`) と流れ方向の密度 (`x_density_table`) で見積もり直してから格子を決める。新の規則での生成は M_knot の最終の組の途中で止めた (A `~/m55_card/`、結果なし)。
- 2026-10-15 (plan `tooling-nozzle-card-local-mesh` §5.1 #3a): カードに `mesh.ni: 3000` (生成のためだけ) を足して AWS A `~/forge-cardmesh` (479723ab) で生成 → `problem_m55_card.yaml` (M_knot 自動 2.80、r_t 0.09990877093936801 m、177 s)。その `mesh`・`mesh_euler` を `_mesh/hand_mesh.py build` で手組みの格子 (場所ごとの第一層の表・流れ方向の密度の表・壁法線の近壁層・17 桁) に置き換えて `problem_m55_ns.yaml` (Euler・NS の両方に使う)。CFD 前の物理壁で作った格子の検査 (forge は起動しない、A `~/m55_mesh_hand/check/`): NS `check_mesh_quality --ar-max 5000` PASS (AR 最大 4849.4)・AR > 1000 の 47,607 セルのスキュー最大 0.0015 (違反 0; 出口面は非平面 0.62 mm で記録だけ)・変換後の第一層厚 ≤ c·y₁,req (全 4,800 断面)・見積もりの y₁⁺ 最大 0.797。Euler PASS (AR 最大 6.1)。再現の入力と記録は `_mesh/` ([`_mesh/README.md`](_mesh/README.md))。
- 2026-10-15 (plan §5.1 #3b): カードの雛形を `axismach_final_v1` に替えて (探索の雛形の物理壁では NS の停止規則が使えない) AWS A `~/forge-cardmesh` (94790833) で生成 → `problem_m55_final_card.yaml`・`_euler.yaml` (M_knot 自動 2.80 で不変、r_t 0.0999084914252905 m、寸法 PASS、194 s)。`_mesh/hand_mesh.py build` で手組みの格子に置き換えて `problem_m55_final.yaml` (NS)・`problem_m55_final_euler.yaml` (Euler)。NS の径方向はユーザ決定 (「上限3%」) で `axis_cap_frac 0.03`・nj 92 (隣接比の最大 1.1974)。検査 (FP64 の変換器、forge は起動しない): NS 品質 PASS (AR 最大 4849.0)・AR 緩和 違反 0 (AR > 1000 の 36,115 セル、スキュー最大 0.0015)・第一層 実寸 / (c·y₁,req) 最大 0.9999997 (超過 0)・見積もりの y₁⁺ 最大 0.797・軸側のセル / r_w 0.030 (全断面); Euler PASS (AR 最大 6.1)。`run_0001_euler` を準備 (未起動)。記録は `_mesh/hand/final/` ([`_mesh/README.md`](_mesh/README.md))。
