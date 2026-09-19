# codex レビュー: convection-node-wall-reconstruction (plan)

- **plan**: [`plans/active/convection-node-wall-reconstruction.md`](../../plans/active/convection-node-wall-reconstruction.md)
- **stage**: `plan`
- **date**: 2026-09-19
- **commit**: `3d4b0d25` (feature/sern-design)
- **codex**: effort `high`, 7.2 min, rc=0
- **判定**: **NO-GO**, 指摘 C1/M7/m1
- **focus**: §4.11 W1b: Venkatakrishnan (生産既定) が密度を一切制限しない件。eps2 = K^3*volume が変数と無関係な次元で K=1 固定。案 S (無次元化: 既存の受動スカラー用 limiter_r1_scaled_d と同じ phi_ref スケーリングを流れ5変数へ) / 案 K (K を config キーに, SU2 既定 0.05) / 両方。論点 (i) 基準値の取り方 (速度は 0 を跨ぐ) (ii) eps2 の局所 volume 依存を残すか SU2 のような大域参照長か (iii) limiterMatchRecon と同キーか別キーか (iv) 既定を変えるか opt-in か。あわせて §4.10 の A/B 実測と既定パス非退行の扱いも見てほしい

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**

W1b の単位依存という問題同定は正しいですが、案 S の単純移植では幾何・速度基準の問題が残ります。  
さらに、診断の誤合格、出力の時相混在、非退行判定の未達を実測で確認しました。実装前に診断と設計契約を修正すべきです。

1. **Critical — 面診断が NaN を「非物理な再構成なし」と判定する**

   根拠は [check_face_reconstruction.py:63](/home/sano/work/forge/solver_density_cuda/tools/check_face_reconstruction.py:63)。判定が `<= 0` だけなので NaN を見逃します。実際に次を再現しました。

   ```text
   case/46.sern_design/run_0218_w1_venkat_mr1_long/
   --steps 31
   非物理な面 0、再構成 P の最小 nan
   VERDICT: OK (非物理な再構成なし)
   ```

   同じ run の `check_convergence.py` は **`DIVERGED (NaN/Inf)`** です。また、指定 snapshot が存在しなくてもスキップされ、検査件数ゼロで合格できます（同ファイル:90）。

   **対案:** 入力・増分・面状態の `isfinite` を必須にし、非有限値、欠落、検査対象ゼロを非ゼロ終了コードで拒否する。この診断を修正するまで G3 合格を受理しない。

2. **Major — 診断の密度と勾配が別の反復時点に属し、G1 も検査していない**

   定常陰解法は残差構築後に保存量を更新し、そのまま出力します。[main.cpp:1543](/home/sano/work/forge/solver_density_cuda/main.cpp:1543)、[main.cpp:1701](/home/sano/work/forge/solver_density_cuda/main.cpp:1701)。`updateVariablesOuter` は勾配・リミッタを再評価しません。

   したがって `res_n.h5` の `ro` は更新後、勾配・リミッタは更新前です。現在の診断はそれらを組み合わせています。`case/46.sern_design/run_0219_w1_barth_mr1_long/` の step 20 で、密度の近傍有界性を丸め許容付きで再計算すると、

   | 密度の参照 | 範囲外の面側数 |
   |---|---:|
   | 同じ snapshot の `ro` | 5,259 |
   | 前 snapshot の `ro` | 0 |

   これは W1 の失敗ではなく、**診断が流束計算時の状態を再現していない証拠**です。また現行ツールは `P,ro > 0` しか見ず、近傍 min/max も速度も検査しません。§4.10 の「G1/G3 合格」は成立していません。

   **対案:** 流束評価直前の状態で、共有再構成関数を使った診断を行う。G1 は流れ5変数の範囲逸脱、G3 は有限性・物理性として分離する。保存 snapshot による診断は評価時点を明示し、`convMethod: 2` にも対応させる。

3. **Major — 案 S＋K でも、寸法付き `volume` を残せば単位依存は解消しない**

   [limiterFunctions_d.cuh:13](/home/sano/work/forge/solver_density_cuda/cuda_forge/limiterFunctions_d.cuh:13) の問題は確認できました。ただし `delta` だけを無次元化しても、右辺の `K³·volume` は無次元になりません。

   加えて forge の `volume` は共通の長さ尺度ではありません。

   - 平面2Dでは面積×単位厚み：[gmshReader.hpp:44](/home/sano/work/forge/solver_density_cuda/mesh/gmshReader.hpp:44)。
   - 軸対称ではさらに半径を掛ける：[variables.cpp:543](/home/sano/work/forge/solver_density_cuda/variables.cpp:543)。
   - 周期 node では合併体積になる：[mesh.cpp:710](/home/sano/work/forge/solver_density_cuda/mesh/mesh.cpp:710)。

   特に軸対称でこれを使うと、同じ再構成増分でも半径によって平滑化が変わります。

   SU2 の係数も直接移植できません。ローカルソースの `CConfig.cpp:5021` は `RefElemLength=1.0`、Venkat はその参照長から閾値を作ります。既定 `K=0.05` の閾値は `1.25e-4`。対象ノードで forge の K だけを同じ値にすると約 `1.28e-10` で、約100万倍違います。[SU2 v8.5.0 の実装](https://raw.githubusercontent.com/su2code/SU2/v8.5.0/SU2_CFD/include/limiters/CLimiterDetails.hpp)

   **対案:** 局所依存は残すが、`ε̂²=(K h_i/L_ref)³` と定義する。`h_i` は平面・軸対称では半径重み前の面積の平方根、3Dでは体積の立方根。周期では合併 CV に対応させる。単位変換・相似拡大・格子細分に対する試験を必須にする。

4. **Major — 速度への局所絶対値スケーリングは、ゼロ交差以外にも問題がある**

   移植元は [passiveLimiter_d.cuh:43](/home/sano/work/forge/solver_density_cuda/cuda_forge/passiveLimiter_d.cuh:43) の `max(|Qc|,|Qmax|,|Qmin|,floor)` です。速度に適用すると、一様な速度加算だけで平滑化強度が変わります。

   例えば `δp=0, |δm|=0.01, ε̂²=1e-6` の同じ局所差分でも、参照値が `0.02 → 1000.02` になると、ψ は **約 `2e-6 → 0.9998`**。成分ごとの絶対速度を基準にする案は採用すべきではありません。

   また §4.11 の「密度はどこでも制限されない」は過大です。`run_0218` の step 20 には、`limiter_ro < 0.9` が **4,661 / 62,317 節点**あり、最小値は **0.002031** でした。問題は局所状態と尺度による不適切な効き方です。

   **対案:** run 中固定の正の物理参照値 `ρ_ref, p_ref, a_ref` を使い、速度3成分は共通の `a_ref` で無次元化する。`p_ref` はゲージ差圧ではなく絶対圧力。係数を公開するなら、現行式の余分な `delta_m` 因子も約分する。現行の3次積は、float32・`K=0, δp=0, δm=1e-19` で **`0/0 → NaN`** を再現しました。

5. **Major — 「密度が頭打ち」「排出が止まる」は実測に反する**

   [plan:273](/home/sano/work/forge/plans/active/convection-node-wall-reconstruction.md:273) の主張を、`case/46.sern_design/run_0219_w1_barth_mr1_long/` の節点62104で確認しました。

   | step | `ro` |
   |---:|---:|
   | 100 | 0.00115554 |
   | 240 | 0.00101269 |
   | 300 | 0.000906855 |
   | 400 | 0.00103937 |

   全400 snapshot の系列を `check_quasisteady.py --series-csv` に渡した結果は、

   ```text
   rho_62104: drift=5.4%/tail、fluct=40.4%
   VERDICT相当: DRIFTING
   OVERALL: NOT ALL STEADY
   ```

   収束判定も **`NOT CONVERGED`**。低下は `rms_roe` が1.2桁、`rms_roK` が0.4桁です。

   **対案:** 結論を「400 step まで有限状態を維持した」に限定する。排出停止の判定には壁 CV の密度時系列と符号付き質量残差を用い、正式な定常性判定まで継続する。run 索引の該当記述も修正する：[README.md:270](/home/sano/work/forge/case/46.sern_design/README.md:270)。

6. **Major — 既定パス非退行は、指定された判定ツールで未達**

   `run_0214_w1_base_mr0`・`run_0222_w1_repeat_mr0` を反復基準、`run_0221_w1_bitcheck_baseline` を候補として、`check_field_regress.py --boundary` を実行しました。すべて `case/46.sern_design/` 配下です。

   | step | 既定のノイズ床×2判定 |
   |---:|---|
   | 1 | **FAIL：4量**。`roUy`、`roe` の L∞差はノイズ床の約2.92倍 |
   | 10 | **FAIL：3量**。`h0` は約3.62倍 |
   | 27 | **PASS** |

   差は小さく、これだけで実装退行とは断定できません。しかし「同オーダー」と「宣言したゲートに合格」は別です。コンパイラのレジスタ割り当てが原因という説明も未検証です。

   **対案:** §4.10 を「短期差は小さいが非退行未確定」に修正。旧・新双方の反復でノイズ床を測り、許容倍率を事前固定する。発散直前の run に加え、安定した標準ケースでも全保存量・原始量・壁出力を判定する。

7. **Major — W1 は設計どおりに共通化されておらず、適用範囲も揃っていない**

   増分式は現在も [limiter_d.cu:309](/home/sano/work/forge/solver_density_cuda/cuda_forge/limiter_d.cu:309)、[limiterPeriodic_d.cuh:103](/home/sano/work/forge/solver_density_cuda/cuda_forge/limiterPeriodic_d.cuh:103)、[convectiveFlux_common_d.cuh:99](/home/sano/work/forge/solver_density_cuda/cuda_forge/convection/convectiveFlux_common_d.cuh:99) に別々に存在します。

   さらに周期側の共通 wrapper は `limiterMatchRecon` を化学種・受動スカラーにも渡しますが、非周期側の `limiter_r1_d` / `limiter_r1_scaled_d` は旧経路のままです。cell も `mr:1` なら増分の算術が変わります。「流れ5変数の node のみ」とは一致しません。

   **対案:** W1b の前に共有増分関数と対象変数を確定する。初回は流れ5変数・node・`convMethod: 0/1/2` に限定し、その他の経路への変更を明示的に遮断する。周期試験は `mr:1` でも実施する。既存 [test_periodic_limiter.cu:142](/home/sano/work/forge/solver_density_cuda/tests/unit/test_periodic_limiter.cu:142) は `mr:0` 固定です。

8. **Major — 検証計画に、変更を実際に作動させる試験と数値基準が不足している**

   [plan:373](/home/sano/work/forge/plans/active/convection-node-wall-reconstruction.md:373) は許容差を先送りしています。G2 の「逸脱が減る」には、対象量・正規化・比較時点がありません。各方式の発散直前を比較すると状態差が交絡します。

   また、標準ケース一覧は [verification/README.md:9](/home/sano/work/forge/procedures/verification/README.md:9) と異なります。Taylor–Green の標準手順は純粋 KEEP であり、これだけでは Venkat の周期経路を検証できません。

   **対案:** 実装前に以下を固定する。

   - 同一の凍結状態で G1/G2 を評価し、その後に起動 A/B を行う。
   - 平滑な極値を含む格子細分試験を入れ、単に強く制限して一次精度へ落ちる案を排除する。
   - `case/20` または `case/08` の標準回帰、SLAU＋Venkat の周期移流、`convMethod: 2`、軸対称、非一様組成の TP を明示する。
   - TP の組成再構成は実際に `limiter_ro` を読むため、圧力・全エンタルピー・種収支も検証する：[convectiveFlux_slau_d.inc.cuh:275](/home/sano/work/forge/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:275)。
   - W2 は有限性検査、EOS再評価、同一面流束の共有、フォールバック保持期間を設計する。保持期間を設けるなら「新規検出数」と「一次化中の面数」を分ける。面状態の正値性だけで block-DPLUR 更新後の正値性まで保証したことにしない。

9. **Minor — 現行方針と撤回済み説明が同じ文書に残っている**

   [plan:6](/home/sano/work/forge/plans/active/convection-node-wall-reconstruction.md:6) は実装着手前の `draft`、§4.10 は実装済み。§2 は更新制限・ロールバックを対象に含め、W3 は別 plan へ分離しています。§4.1 の一方弁説明、§6 の削除済み測定器、§9 の撤回済み SU2 比較も残っています。

   **対案:** 現在の結論を本文に一本化し、旧説は撤回済みと明示する。`plans/README.md` と accepted plan を照合した限り、W1/W1b は未解決の独立課題です。一方、更新ガードには既存の [accepted plan](/home/sano/work/forge/plans/accepted/time_integration-update-positivity-guard.md) があり、W3 を分離する判断は妥当です。

推奨は **「案 S＋案 K。ただし固定物理参照値と無次元の局所長さを使う方式」** に絞ります。速度は共通音速尺度、幾何は上記 `h_i/L_ref`、平滑化設定は `limiterMatchRecon` と別キーにします。初期導入は **opt-in** とし、SU2 の `0.05` を同等設定として移植しません。

優先順は、**診断の誤合格・時相修正 → W1 の共通化と範囲確定 → W1b の尺度・係数・ゲート確定 → 実装・回帰 → W2 の独立 A/B** です。読み取り専用レビューのため、これらは **plan 未反映**、ファイル変更なしです。

指摘数: Critical 1 / Major 7 / Minor 1
