# codex レビュー: performance-3d-node-sst-speedup (plan)

- **plan**: [`plans/active/performance-3d-node-sst-speedup.md`](../../plans/active/performance-3d-node-sst-speedup.md)
- **stage**: `plan`
- **date**: 2026-09-12
- **commit**: `0512823d` (feature/perf-3d-speedup)
- **codex**: effort `high`, 5.6 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M7/m2

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
高速化の目的は妥当ですが、現行案には離散式の変更と「ビット同一」の誤分類があります。  
計測・検証基準を修正し、まず離散式を保つ最適化に絞ることを推奨します。

レビュー対象は貼付 plan と開始時の `HEAD=0512823d` です。AWS の `run_0400_perf_baseline`、`run_0234_user_node3d_sst_dry_half_ext_ps_newdef` の実データ・生ログ・`CONVERGENCE_VERDICT.txt` はこの作業ツリーに無く、**82.85 ms/step と基準場の収束は独立検証できていません**。ファイルは変更していません。

1. **Major — 化学種物性のセル前計算＋補間は、精度変更ではなく離散式変更です。**

   根拠: 現行コードは面の `T/P/Y` を作り、その状態で `D_s` と `h_s` を評価します。plan §4.2-2 はこれをセル評価値の補間に置き換えます。一般に  
   `h(f*T0+(1-f)*T1) ≠ f*h(T0)+(1-f)*h(T1)`  
   であり、`D(T,P,X)` も同様です。[speciesTransport_d.cu:214](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/speciesTransport_d.cu:214)、[同:240](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/speciesTransport_d.cu:240)、[同:257](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/speciesTransport_d.cu:257)

   内蔵 `N2` 係数をメモリ上で倍精度評価すると、`T0=300 K, T1=900 K, f=0.5` で両者の差は **8,716.88 J/kg** でした。これは CFD run の結果ではなく、**丸め誤差がなくても変更が生じる反例**です。[係数:62](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/thermo_d.cu:62)

   **対案:** 本 plan では面状態での評価を維持し、重複演算削減と演算精度変更を個別に評価してください。セル物性の補間は別のスキーム変更として切り出します。一様組成に近い dry ノズルだけでは、変更した拡散流束を十分検証できません。

2. **Major — float Newton の安全性を、出力値の解像度だけでは保証できません。**

   根拠: plan §4.2-3 の `e≈2e5` に対する float 間隔は、実際には **0.015625 J/kg** です。さらに問題は間隔そのものより、NASA 多項式、混合和、`h-R*T`、Newton 残差の演算誤差です。既存仕様は生成エンタルピーを含む大きな値の桁落ちを理由に double を保持しています。[thermo_d.cuh:414](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/thermo_d.cuh:414)、[thermophysics.md:341](/home/sano/work/forge-perf/methods/thermophysics.md:341)

   また、温度反転後に `roe` を再構成して上書きするため、反転誤差は表示温度だけに留まりません。凝縮有効時でも液相が少ない分岐は通常の `thermo_T_from_e` を使うので、`cond_*` を double のまま残すだけでは影響を隔離できません。[dependentVariables_d.cu:138](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/dependentVariables_d.cu:138)、[同:167](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/dependentVariables_d.cu:167)

   **対案:** Newton の float 化は後段へ送り、まず double 参照に対する単体検証を追加してください。対象は使用 DB、組成端点・微量成分、50–6000 K、温度区間境界・外挿、datum の有無。温度誤差に加えてエネルギー残差、反復上限到達、反転→再構成のドリフトを測り、未検証条件は double に戻す設計にします。

3. **Major — 「100 step 後の全場ビット一致」は、現行ソルバの非決定性と整合しません。**

   根拠: 面流束の蓄積には浮動小数 `atomicAdd` が使われ、既存 M6 検証にも同一バイナリ間の非決定性が記録されています。[speciesTransport_d.cu:255](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/speciesTransport_d.cu:255)、[既存 plan:122](/home/sano/work/forge-perf/plans/accepted/thermophysics-multicomponent-tpgas.md:122)

   加えて、`pow` の double 評価を float の `x*sqrt(x)` に変える処理はビット同一ではありません。現在の 5×5 解法はピボット付き消去と RHS 更新を交互に行うため、LU 保存も演算順序の同一性を別途確認する必要があります。[gasProperties_d.cu:62](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/gasProperties_d.cu:62)、[timeIntegration_d.cu:177](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/timeIntegration_d.cu:177)

   **対案:** ビット一致は固定入力に対する局所カーネル出力へ限定します。全計算は旧×旧反復でノイズ幅を測ってから A/B 判定してください。float 化の比較対象には、現行リストから漏れている `Uy/Uz/roe/roY*/Y*/vis_turb` も含め、ゼロ近傍を扱える量別の絶対・相対許容差を定義します。

4. **Major — 現行ベンチは、本番速度を測る基準として不十分です。**

   根拠: `bench_steps.sh` は既定で `FORGE_PROFILE=1` を設定し、`RuntimeProfiler::measureCuda` はセクションごとに `cudaEventSynchronize` を挿入します。本番の非同期実行とは条件が違います。[bench_steps.sh:27](/home/sano/work/forge-perf/solver_density_cuda/tools/bench_steps.sh:27)、[main.cpp:839](/home/sano/work/forge-perf/solver_density_cuda/main.cpp:839)

   また、スクリプトは `Time` を「起動込み」と表示しますが、この時計は初期出力後に開始されています。プロセス全体の `/usr/bin/time` と混同しています。[main.cpp:1717](/home/sano/work/forge-perf/solver_density_cuda/main.cpp:1717)、[bench_steps.sh:35](/home/sano/work/forge-perf/solver_density_cuda/tools/bench_steps.sh:35)

   **対案:** 性能合否は `FORGE_PROFILE=0`、GPU warm-up 後、複数回の A/B 交互実行で判定してください。初期化・I/O・時間ループを分離し、バイナリ・config・IC・メッシュのハッシュ、ビルド条件、GPU 状態、生ログを保存します。nsys/ncu は原因分析用に別測定し、FP64 稼働率から「42 ms がそのまま削減可能」とは解釈しません。

5. **Major — 「同桁の残差＋12000 step＋壁圧差」では収束解不変を判定できません。**

   根拠: plan §4.3・§6 には、比較する両 run の `PASS` を必須とする条件と、壁圧自体の準定常判定がありません。`check_convergence.py` は各保存量のピークからの低下とトレンドを判定するため、単なる最終残差の桁比較とは異なります。[plan:95](/home/sano/work/forge-perf/plans/active/performance-3d-node-sst-speedup.md:95)、[check_convergence.py:91](/home/sano/work/forge-perf/solver_density_cuda/tools/check_convergence.py:91)

   **対案:** 同一 IC から基準・変更版を走らせ、両者の全保存量の `PASS` と、報告する壁圧分布・流量等の `STEADY` を必須にしてください。壁圧に対応する時系列判定を整備し、`check_quasisteady.py` の閾値も壁圧差 **0.1%** と整合させます。12000 step は初回投入長とし、未収束なら比較を確定しません。CFL/sweep 比較も、共通の到達条件までの壁時計で評価します。

   §5.1 の先頭には、基準 run の VERDICT 確認、メッシュ品質、旧キー監査、IC 確認、初期 NaN 検査、スナップショット・残差図保存・case README 索引更新を追加してください。

6. **Major — 対角キャッシュの適用範囲と回帰試験が不足しています。**

   根拠: 状態凍結中の対角再利用自体は妥当です。しかし現行には `double` solve、line-implicit、軸・壁・等温壁の拘束行、sweep ごとの周期補正ミラーがあります。既存 `diag_block_*` は line 経路で保存する配列です。[timeIntegration_d.cu:962](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/timeIntegration_d.cu:962)、[同:989](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/timeIntegration_d.cu:989)、[同:1419](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1419)、[main.cpp:1282](/home/sano/work/forge-perf/solver_density_cuda/main.cpp:1282)

   plan の float キャッシュを double 経路にも使えば精度が変わります。また、対角を再利用しても **RHS の拘束処理は各 sweep に必要**です。現行の検証一覧ではこれらを網羅しません。共有コード変更時の Taylor–Green 検証も手順に明記されています。[verification/README.md:27](/home/sano/work/forge-perf/procedures/verification/README.md:27)

   **対案:** 最初は float・point-DPLUR 経路に限定し、キャッシュは各 `blockDPLURSolve` 呼び出しで更新する、と明記してください。拘束行、周期ミラー、LU ピボット・失敗処理を保持し、対象経路の検証を追加します。node/cell の標準ケースに加え、周期保存、軸近傍、等温壁、非一様組成拡散を検証してください。凝縮へ共有変更が届く場合はその回帰も必要です。

7. **Major — 計測スクリプトが検証証拠と restart 入力を削除し得ます。**

   根拠: 同じ run を label 違いで再利用し、最後に無条件で `rm -f res_*.h5 res_*.xmf` を実行します。既存の restart ファイルや `res_nan_*` も対象です。さらに `set -e` がなく、config 変更用 Python の失敗を明示確認しないため、意図しない設定で計測を続行できます。[bench_steps.sh:7](/home/sano/work/forge-perf/solver_density_cuda/tools/bench_steps.sh:7)、[同:9](/home/sano/work/forge-perf/solver_density_cuda/tools/bench_steps.sh:9)、[同:37](/home/sano/work/forge-perf/solver_density_cuda/tools/bench_steps.sh:37)

   **対案:** A/B・反復ごとに専用 run を作り、入力を不変にしてください。設定生成失敗時は即終了し、実行 step 数も確認します。ワイルドカード削除は撤去し、NaN ダンプと各試行の残差 CSV を保存します。

8. **Minor — 既存の TP 高速化との関係が整理されていません。**

   根拠: accepted の M6 ですでに輸送係数の FP32 化、`cp+h` 融合、`Rmix` キャッシュが完了しています。現在の `thermo_Dbinary` と `thermo_Dmix_species` も float 実装です。「化学種拡散は全演算 double」は正確ではありません。[既存 plan:118](/home/sano/work/forge-perf/plans/accepted/thermophysics-multicomponent-tpgas.md:118)、[thermo_d.cuh:370](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/thermo_d.cuh:370)

   **対案:** M6 と既存 profiling plan を関連計画に追加し、「実装済み」「今回残る double 演算」「別の近似を導入する案」を整理してください。今回の目的全体が解決済みという意味ではありません。

9. **Minor — GPU 性能の前提値と現行宣言に誤記があります。**

   根拠: plan §3 の CC 8.6 の FP64/FP32 比 `1/32` は、通常の加算・乗算・FMA のスループット表では **2/128＝1/64** です。[plan:36](/home/sano/work/forge-perf/plans/active/performance-3d-node-sst-speedup.md:36)、[NVIDIA CUDA Programming Guide](https://docs.nvidia.com/cuda/archive/12.6.3/cuda-c-programming-guide/index.html#arithmetic-instructions)。また、現行カーネルは `__launch_bounds__(128)` 相当であり、`(128,2)` ではありません。[timeIntegration_d.cu:649](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/timeIntegration_d.cu:649)

   **対案:** plan と `performance.md` を訂正し、占有率変更は現行ビルドのレジスタ数・spill・実測時間から評価してください。

**推奨は、離散式を維持する最適化を先行させる一案です。** 実装前の修正順は次のとおりです。

1. ベンチの削除・再利用を修正し、基準バイナリ・入力・生ログ・VERDICT を固定する。
2. ビット一致、許容誤差、収束・準定常・速度の合否条件を定義する。
3. セル物性補間を本 plan から外し、Newton float 化を検証付きの後段へ移す。
4. リミッタ特殊化、確認済みのリテラル昇格除去、適用範囲を限定した対角再利用から実装し、効果が測定誤差内または悪化する施策は採用しない。

反映先は対象 plan の §4.2・§4.3・§5.1・§6 です。read-only 指示に従い、**plan 未反映**です。

指摘数: Critical 0 / Major 7 / Minor 2
