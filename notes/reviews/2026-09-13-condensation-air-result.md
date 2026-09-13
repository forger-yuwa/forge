# codex レビュー: condensation-air (result)

- **plan**: [`plans/active/condensation-air.md`](../../plans/active/condensation-air.md)
- **stage**: `result` (diff base `feature/sern-design`)
- **date**: 2026-09-13
- **commit**: `9167bf35` (feature/condensation-air)
- **codex**: effort `high`, 8.1 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M5/m2
- **focus**: 実装 diff (condensationProperties の低温整合と CondPropOpts、括弧付き cond_T_from_e_cpg、CPG carrier 形の EOS/source/SLAU 面状態/実現可能性、slip 状態保持、config 検査) が plan §4 (v3) と codex plan レビュー 2 回の採用事項に一致するか。§9 の結果 (N2 新 +2.1 K / 空気 +1.7 K で理論線 ±3 K 以内、空気−N2 0.7 K、node/cell 一致、slip 回帰、Arthur 実験との ~9 % 過大は物性前後で不変) が run (case/34 run_0014–0026, onset_analysis.py --series) の実測で裏付けられるか。accepted にしてよいか
- **extra**: `case/34.arthur_n2_nozzle/README.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**NO-GO**

主要な cell 結果（新 N2＋2.14 K、空気＋1.69 K、onset 差0.73 K）は再現できました。
ただし、温度反転失敗時の処理、蒸発時の実現可能性判定、受付範囲、node 後処理に不備が残り、現状の「実装・検証完了」は承認できません。

1. **Major — 温度反転失敗時に、失敗した温度をそのまま採用している。**

   **根拠:** [dependentVariables_d.cu:230](/home/sano/work/forge-cond/solver_density_cuda/cuda_forge/dependentVariables_d.cu:230) は `ok=false` でも `T/P/sonic` を更新します。条件付きなのは `roe` の代入だけで、失敗フラグも外部へ残しません。[plan §4.1:80](/home/sano/work/forge-cond/plans/active/condensation-air.md:80) の「前ステップの T を保持しフラグを診断に出す」が未実装です。保存エネルギーと不整合な温度・圧力が、次の流束・核生成評価へ渡ります。

   **対案:** `ok=false` を原始量更新より前で処理し、計画どおり温度を保持して診断を残す。反転不能・非有限入力を含む試験で、保存量保持と失敗検出を確認してください。

2. **Major — CPG carrier の液滴消滅判定だけ、N2 分圧ではなく全圧を使っている。**

   **根拠:** [condensationTransport_d.cu:91](/home/sano/work/forge-cond/solver_density_cuda/cuda_forge/condensationTransport_d.cu:91) は `roY_w==nullptr` なら `pv=P` とします。CPG carrier はこの分岐へ入り、`opts.Yw` は上限処理にしか反映されません。source 側の分圧式と不一致です。

   コードの物性・EOSで、`T=50 K`、`Yw=0.7671`、`g=0`、`ρ=0.0462648 kg/m³` とすると：

   - 全圧：666.68 Pa
   - N2 分圧：526.67 Pa
   - 飽和圧：606.07 Pa

   N2 は未飽和ですが、この処理は過飽和として戻ります。`g=0` で `Q0/Q1/Q2` が残った「モーメント塵」を消せません。

   **対案:** source と実現可能性処理で蒸気状態の評価を共通化し、CPG carrier では `ρ max(Yw−g,0) Rw T` を使う。上記の全圧と分圧で判定が分かれる状態を回帰試験に追加してください。

3. **Major — 計画で限定した境界構成を、実装は制限していない。**

   **根拠:** [solverConfig.cpp:712](/home/sano/work/forge-cond/solver_density_cuda/input/solverConfig.cpp:712) は solver・物性モデル等を検査しますが、[plan §4.1:88](/home/sano/work/forge-cond/plans/active/condensation-air.md:88) の境界制限を実装していません。例えば CPG carrier と `wall` を組み合わせられますが、[boundaryCond_d.cu:235](/home/sano/work/forge-cond/solver_density_cuda/cuda_forge/boundaryCond_d.cu:235) の CPG 壁は単相 EOS で圧力を再構成します。

   `ρ=1`、`T=40 K`、`g=0.1`、静止した空気 carrier 状態では、二相 EOS の圧力 **10,340.8 Pa** に対し、この ghost 式は **2,474.8 Pa** になります。

   また `condVaporMassFraction` の検査を `>0` の内側に置いているため、`.nan` は検査を通らず pure 経路へ落ちます。

   **対案:** 境界読込み後に受付範囲を検査し、未対応境界を拒否する。超音速出口という動的条件には実行時検査を設ける。`Yw` は分岐前に有限性と許容値を検査してください。

4. **Major — node の「中心線」に中心線外の点が混入し、膨張率と旧 N2 onset を誤算している。**

   **根拠:** [onset_analysis.py:31](/home/sano/work/forge-cond/case/34.arthur_n2_nozzle/onset_analysis.py:31) は全点を等間隔の x ビンへ分け、その中の最小 `|y|` を選びます。中心線点のないビンでも内部点を採用します。

   `case/34.arthur_n2_nozzle/` の保存場から再計算すると：

   | run | 現行 `--series` の Ṗ | 実際の `y=0` 列の Ṗ |
   |---|---:|---:|
   | `run_0020_n2_ref_node` | 7.17×10⁴ /s | 1.75×10⁴ /s |
   | `run_0022_n2_new_node` | 5.16×10⁴ /s | 1.95×10⁴ /s |

   `run_0020` の報告 onset 点は **y=2.229 mm** です。中心線だけで評価すると onset は **2.3566 in／37.768 K** となり、記載の2.336 in／37.89 Kから変わります。現行ツールは誤った抽出結果にも `VERDICT(series): STEADY` を返します。

   **対案:** node は中心線境界のノード列、壁・出口は `BCONDS` の面と隣接 CV から抽出する。修正後に全時系列を再判定し、Ṗを含む表を保存してください。空気用の `--R` と各 run の `--dry` も再現コマンドに明記し、報告量を `check_quasisteady.py` の判定へ接続してください。

5. **Major — 回帰合格の根拠と「全件採用・検証完了」の主張が不足している。**

   **根拠:** [plan §9:206](/home/sano/work/forge-cond/plans/active/condensation-air.md:206) の場差を、`max|新−旧| / max|旧|` で再計算しました。

   - `case/34.arthur_n2_nozzle/run_0015_dry_slip_cell/` 対 `run_0015o_dry_slip_cell_oldbin/`：ρ **3.003×10⁻⁴**、Uy **1.996×10⁻³**。
   - `case/34.arthur_n2_nozzle/run_0014_n2_ref_cell/` 対 `run_0008_ref_n2/`：Uy **2.264×10⁻³**。「場差≤1e−3」は成立しません。
   - node の新旧差 **8.636×10⁻⁶** は再現できました。

   cell の差を、ρについて示された約1e−3のノイズから全変数へ拡張して「ノイズ以内」とは判定できません。

   また [test_cond_air.cpp:72](/home/sano/work/forge-cond/solver_density_cuda/tests/unit/test_cond_air.cpp:72) の SLAU・slip 試験は、テスト内に書いた double の代数式を比較するだけです。実カーネルや float32 保存量からの温度復元を通しておらず、前回レビューの要求を満たしていません。

   **対案:** 同じ初期場・同じ反復数の同一バイナリ反復で、変数別・同一ノルムのノイズを測る。実装を通す float32 試験と、採用事項だった node/cell の質量・エネルギー流束収支を補い、R0/R1 の合否を再判定してください。

6. **Minor — 「約9%過大の原因はレート側」とする判断は、今回の比較からは導けない。**

   **根拠:** [plan §5.1:149](/home/sano/work/forge-cond/plans/active/condensation-air.md:149) は「物性修正前後で差が変わらないので物性ではなくレート側」と断定しています。確認できるのは、この物性変更で壁圧偏差が解消しなかったことです。これだけでは、メッシュ依存性・形状近似・境界条件・他のモデル近似を排除できません。

   **対案:** 後続項目を「レート較正」から「壁圧過大の原因切り分け」へ変更する。併せて今回の未解決事項を §5.1 に戻し、未掲載の準定常ツール統合も追加してください。採用済みの Iland 使用方針が [§10:232](/home/sano/work/forge-cond/plans/active/condensation-air.md:232) で未確定のままなのも同期が必要です。

7. **Minor — 現在仕様・検証条件の文書に、古い値と実装と異なる説明が残っている。**

   **根拠:**

   - [methods/condensation.md:496](/home/sano/work/forge-cond/methods/condensation.md:496)：実装済み物性を「計画中」、飽和圧変化を旧記述の0.59倍と説明。
   - [同:646](/home/sano/work/forge-cond/methods/condensation.md:646)：CPG 音速を `√(γp/ρ)` と記載。実装は `√(γRgas T)` で、二相では異なります。
   - [case README:9](/home/sano/work/forge-cond/case/34.arthur_n2_nozzle/README.md:9)：「凝縮モデルは無い」が残存。
   - [plan:164](/home/sano/work/forge-cond/plans/active/condensation-air.md:164)：空気入口324.48 m/s・6.161 kg/m³と記載。実 config は **325.12 m/s・6.137 kg/m³**。実設定からは `M=1.05004`、`T0=290.002 K`、`P0=844.128 kPa` となり、設定側は妥当です。
   - [plan:221](/home/sano/work/forge-cond/plans/active/condensation-air.md:221)：`cl ±500` に対して onset `∓0.7 K` とありますが、実測38.95→39.67→40.25 Kは同符号です。

   **対案:** 実装と実 config を正本に修正し、新キー・受付範囲・旧物性への戻し方を `procedures/` にも記載する。`methods/index.md` と `plans/README.md` の項目追加は確認できました。

再実行で裏付けられた部分も明確です。対象15 runの全保存場で NaN/Inf はなく、`check_convergence.py` は `run_0021_dry_slip_node` のみ **PASS**、他は **NOT CONVERGED**。現行 `--series` は全対象 **STEADY**、`check_quasisteady.py --quantity machmax,pmax` も **STEADY**でした。ただし後者は onset・壁圧比の判定を代替しません。

`case/34.arthur_n2_nozzle/run_0017_n2_new/` と `run_0024_air_cpgcarrier/` の理論線との差 **＋2.141／＋1.692 K**、壁圧比 **1.2421/1.3715/1.4758／1.2138/1.3514/1.4547** は再現できています。主要6 runの実出口面でも全保存時刻で `u_n/c>5.44`。メッシュ品質は cell/nodeとも **VERDICT: PASS** でした。これらは[case の run 一覧](/home/sano/work/forge-cond/case/34.arthur_n2_nozzle/README.md:35)に紐づく**未収束の準定常比較**として扱えます。

**推奨は、`active` に留め、1・2の状態整合修正 → 3の受付制限 → 4の再解析 → 5の回帰検証 → 文書同期の順で完了させ、result レビューを再実施することです。** ファイルは変更しておらず、指摘は **plan 未反映**です。既存単体実行ファイル4本の PASS は確認しましたが、再ビルドは実施していません。

指摘数: Critical 0 / Major 5 / Minor 2
