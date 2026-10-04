# codex レビュー: tooling-nozzle-sern-chain (plan)

- **plan**: [`plans/active/tooling-nozzle-sern-chain.md`](../../plans/active/tooling-nozzle-sern-chain.md)
- **stage**: `plan`
- **date**: 2026-09-19
- **commit**: `879b05bf` (feature/sern-design)
- **codex**: effort `high`, 5.4 min, rc=0
- **判定**: **NO-GO**, 指摘 C1/M4/m0
- **focus**: §4.14 の 3D 領域トポロジに集中。(1) 機体幅内の下面は A/B/C/D のどれか (2) 幅外を fluid にする際の見落とし (3) 床 20Pa と DEPVAR_TMAX 6000K が m10_on で物理解を切る件は幾何修正だけで足りるか

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**

§4.14 の「幅外に架空の固体がある」という診断は正しい。ただし、温度上限を発散原因とする結論は実測で裏付けられていません。  
推奨は **D：有限幅ノズルの基準形**です。体積トポロジ、壁熱条件、EOS の許容範囲を修正対象として明文化してから着手してください。  
ファイル変更・新規 CFD 実行は行っていません。以下はコード、既存保存場、読み取り専用の再計算に基づきます。

1. **Critical — D 案にも、有限厚機体の側面と幅外の流体接続が欠けている**

   **根拠:** [`mesh_sern3d.py:228`](/home/sano/work/forge/design/forge_design/meshing/mesh_sern3d.py:228) の `W_vehicle` は境界タグを選ぶだけです。上面生成も [`同:288`](/home/sano/work/forge/design/forge_design/meshing/mesh_sern3d.py:288) で全幅に押し出されています。メモリ上で `W_vehicle=4→2` を比較した結果、**座標とセル接続は完全同一**でした。したがって、パラメータ変更だけでは D になりません。

   また、D でも `z=W/2` に、ランプ下面から機体上面までを閉じる**機体側面**が必要です。この面は `L_sw` より下流にも存在します。`sidewall_in/out` が表すダクト側壁とは別物です。「横方向レリーフは `L_sw` が決める」だけでは形状仕様が足りません。

   **対案:** D の固体体積を明示し、次を実装前に固定してください。

   - `z>W/2` では旧下面と旧上面の間を流体セルで埋め、上下の外部流と内部面で接続する。
   - 機体側面を専用タグにし、ダクト側壁と力の帳簿を分ける。
   - ランプ端・機体側面・上面・側壁後縁が交わる各線で、共有ノード／壁両側の重複ノードを定義する。
   - 新規流体領域の組成・初期値を領域情報から与える。現行の index 算術による分類 [`runner_sern3d.py:72`](/home/sano/work/forge/design/forge_design/evaluate/runner_sern3d.py:72) は、新トポロジへ無条件に流用しない。

   既存のメッシュテストは再実行で **30 項目 ALL PASS**でした。つまり現在のテストは、架空の固体を含む領域も合格させます。

2. **Major — 「1 次だからスキームではない」「6000 K クランプが原因」は断定できない**

   **根拠:** `case/46.sern_design/run_0118_r5_3d_sst_frozen_m6on/` を再確認しました。

   | 確認対象 | 結果 |
   |---|---|
   | `check_convergence.py` | **DIVERGED (NaN/Inf)** |
   | CSV の最初の非有限値 | step **179**, `outer_begin`：`ro`, `roe`, `roY0`, `roY1` の残差 |
   | `res_nan_180.h5` | `ro` NaN **116 ノード**。位置は計画記載と整合 |
   | `T=6000 K` | **1 ノードのみ**。その保存量・組成は既に NaN |
   | 密度が有限かつ `P=20 Pa` | **140 ノード**、うち **139 ノードが `T=50 K`** |

   保存場は初期場と発散ダンプだけで、6000 K 到達が NaN に先行した証拠はありません。1 次でも SLAU、node 境界処理、block-DPLUR、EOS 射影は働くので、除外できるのは「SST が必要条件」「2 次再構成が必要条件」という仮説までです。

   さらに [`dependentVariables_d.cu:209`](/home/sano/work/forge/solver_density_cuda/cuda_forge/dependentVariables_d.cu:209) はクランプ後の温度から **`roe` を再構成**します。問題は単に「保存エネルギーと温度が不整合」ではなく、**保存エネルギーを非保存的に変更すること**です。

   **対案:** §4.14 の根因断定を撤回し、最初の床・上限到達時点で、更新前後の保存量、反転残差、組成、補正量を記録する検証を先に置いてください。`DEPVAR_TMAX` の引き上げを発散対策として先行させる根拠はありません。

3. **Major — m10 の物理解を切るのは圧力床だけではない。幾何変更だけでは保証できない**

   **根拠:** 生産 YAML の外気条件と既存 `FrozenGas.air()` を使い、凍結組成・等エントロピーの平面膨張を独立に積分しました。以下は **CFD の定常値ではなく解析参照値**です。膨張で静温・静圧が下がる関係は [NASA の説明](https://www.grc.nasa.gov/WWW/BGH/expans.html)とも整合します。

   | 外気条件・転向角 | 静圧 | 静温 | 密度 |
   |---|---:|---:|---:|
   | m6・20.2° | 50.4 Pa | 69.5 K | `2.53e-3 kg/m³` |
   | m10・20.2° | **0.147 Pa** | **18.0 K** | **`2.83e-5 kg/m³`** |

   m10 は `pMin=20 Pa` に加えて、[`DEPVAR_TMIN=50 K`](/home/sano/work/forge/solver_density_cuda/cuda_forge/dependentVariables_d.cu:8)、[`roMin` の既定値 `1e-4`](/home/sano/work/forge/solver_density_cuda/input/solverConfig.hpp:417) にも抵触します。`run_0118` の config は `roMin` を上書きしていません。

   一方、この膨張は高温化を要求しません。**6000 K 上限が m10 の正常な膨張解を切る、という説明は誤り**です。NASA-9 の外挿が実装されていることも、18 K の物性モデルが物理的に検証済みであることを意味しません。

   **対案:** D で不要な幅外膨張を除いた後も、残る側縁・後縁膨張について `pMin/roMin/Tmin/Tmax` の到達と感度を確認してください。現行ゲート [`sern_gates.py:39`](/home/sano/work/forge/design/forge_design/metrics/sern_gates.py:39) は有限・正値を検査しますが、**床に張り付いた有限解を拒否しません**。EOS 整合性、床による保存量変更、床感度を受理条件へ追加する必要があります。

4. **Major — 3D runner は生産指定の等温壁を無視している**

   **根拠:** [`problem_r5_3d_sst_frozen_doe001.yaml:9`](/home/sano/work/forge/case/46.sern_design/problem_r5_3d_sst_frozen_doe001.yaml:9) は `wall_thermal: {mode: isothermal, Tw: 1000.0}` を指定しています。しかし [`runner_sern3d.py:44`](/home/sano/work/forge/design/forge_design/evaluate/runner_sern3d.py:44) は `kind: wall` を直接生成し、実際の `run_0118` もランプ・カウル・側壁が断熱壁です。

   2D は [`runner_sern.py:269`](/home/sano/work/forge/design/forge_design/evaluate/runner_sern.py:269) で `p.wall_bcond_line()` を使っており、既に共通実装があります。R4b の処方は 3D に届いていません。

   **対案:** 3D の物理壁も共通生成関数へ接続し、生成済み `bcondConfig.yaml` の壁種・`Ts` を検査してください。機体上面・新設側面の slip／粘性壁の選択は別途明示します。この配管修正は R5 の「生産構成再現」の前提です。ただし、これだけで `run_0118` の発散が直るとは断定しません。

5. **Major — 幅外を開いた後の遠方境界・領域独立性を、現在の検証計画では保証できない**

   **根拠:** [`runner_sern3d.py:49`](/home/sano/work/forge/design/forge_design/evaluate/runner_sern3d.py:49) は `side_far` を slip に固定しています。これは横方向流出を遮り、有限幅機体から来た擾乱を反射します。`top_out` も `run_0118` では slip で、面形状は機体上面に追従します。

   [`r4_domain_study.py:22`](/home/sano/work/forge/case/46.sern_design/r4_domain_study.py:22) は加速点の Euler を対象とし、明示された許容値は `|ΔC_T|<0.002` だけです。生産 TP・SST と `C_L/C_M` の独立性を保証しません。引用された `run_0102_r4_domain/` の成果物も、この checkout にはありません。

   **対案:** R4→R5 を次の順に具体化してください。

   - 固体領域と流体領域の断面検査、内部面の共有、境界面の過不足、正体積、双対閉性を検査する。AR/skew の PASS とは別ゲートにする。
   - まず生産条件 m6/m10 の外気膨張を node で検証し、床・反射・境界市松を切り分ける。
   - 固体形状と近傍格子を固定して外部領域だけを拡大し、遠方境界の影響を測る。
   - 全生産作動点で全残差の VERDICT と `C_T_with_shear/C_L/C_M` の **STEADY** を要求し、各係数の格子・領域許容誤差を数値で固定する。

   なお現行手順は [`procedures/verification/README.md:40`](/home/sano/work/forge/procedures/verification/README.md:40) の **node のみ**です。cell 回帰を今回の必須条件にはしません。

**推奨は D に一本化します。** 未知の機体下面を A/B/C で作り、計算が通る角度に調整する根拠はありません。D を「有限幅ノズル単体の基準形」と位置づけ、実機全体の再現とは区別してください。D でも側面・側縁流れは残るため、幾何修正だけで数値健全性が成立したとは扱えません。

実装前の優先順は、**①D の固体体積・接続・境界・帳簿を確定 → ②等温壁配管を修正対象へ追加 → ③温度原因説を訂正して EOS 診断・受理条件を追加 → ④生産条件の検証基準を確定**です。一般 EOS Jacobian、床の config 化、3D median-dual 幾何には既存の accepted plan があり、再実装せず利用します。今回必要なのは、それらの上で成立する有限幅領域と評価器です。

読み取り専用レビューのため、これらの提案は **plan 未反映**です。反映先は対象 plan の §4.14・§5.1 R4/R5・§6 です。

指摘数: Critical 1 / Major 4 / Minor 0
