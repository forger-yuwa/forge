# codex レビュー: condensation-air (result)

- **plan**: [`plans/active/condensation-air.md`](../../plans/active/condensation-air.md)
- **stage**: `result` (diff base `feature/sern-design`)
- **date**: 2026-09-13
- **commit**: `4390f17c` (feature/condensation-air)
- **codex**: effort `high`, 7.1 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M3/m2
- **focus**: 2026-09-13 の result レビュー 1 回目 (NO-GO, M1–M5/m1–m2) の採用・反映 (§6.1 / §9 2026-09-13) が実装 diff・単体試験・再取得 run (run_0027–0032) と整合しているか。特に M1 (T 反転失敗時の凍結)、M2 (分圧判定)、M5 (同一バイナリ反復ノイズ床と --tolfile 判定) の裏付け。accepted にしてよいか
- **extra**: `case/34.arthur_n2_nozzle/README.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**NO-GO**

M2 の分圧化と通常状態での主要結果は確認できました。  
ただし、M1 は `Inf` 入力を反転成功と誤判定し、M5 の「全変数 ok」は比較ツールの **FAIL** と矛盾します。  
現状の「前回指摘を全件反映済み」は承認できません。

1. **Major — M1：`e_in=±Inf` が成功扱いになり、凍結・警告を迂回する。**

   **根拠:** [condensationEOS_d.cuh:101](/home/sano/work/forge-cond/solver_density_cuda/cuda_forge/condensationEOS_d.cuh:101) は `tol=1e-9*abs(e_in)+0.05`、[同:114](/home/sano/work/forge-cond/solver_density_cuda/cuda_forge/condensationEOS_d.cuh:114) は `abs(G)<=tol` で成功を返します。有限の `T_guess=40 K`、`g=0.1` で、この成功条件を実評価すると：

   | `e_in` | `tol` | 初回 `G` | 成功判定 |
   |---|---:|---:|---|
   | `+Inf` | `Inf` | `−Inf` | `true` |
   | `−Inf` | `Inf` | `+Inf` | `true` |

   したがって [dependentVariables_d.cu:226](/home/sano/work/forge-cond/solver_density_cuda/cuda_forge/dependentVariables_d.cu:226) の凍結分岐に入らず、成功側で `roe` を再構成します。**異常な保存エネルギーを有限値で上書きし、`g_condTinvFail` にも数えない経路**が残っています。

   `ok=false` 時に `T/P/sonic/Ht` を保持する修正自体は確認しました。しかし、再取得 run の警告0件は、この異常入力の検出能力を保証しません。

   **対案:** 反転入口で非有限入力を拒否し、成功条件にも残差・温度の有限性を要求する。`±Inf`、NaN、有限だが反転不能な入力を実装経由で試験し、呼出し側の保存量保持と診断発火まで確認してください。

2. **Major — M5：「2倍ノイズ床内・全変数 ok」は実測と違う。**

   **根拠:** [plan §9:231](/home/sano/work/forge-cond/plans/active/condensation-air.md:231) と [case README:67](/home/sano/work/forge-cond/case/34.arthur_n2_nozzle/README.md:67) の主張を、次の保存場で再判定しました。

   - 参照：`case/34.arthur_n2_nozzle/run_0024_air_cpgcarrier/res_12000.h5`
   - 比較：`case/34.arthur_n2_nozzle/run_0027_air_cpgcarrier_v2/res_12000.h5`
   - 許容値：`case/34.arthur_n2_nozzle/noise_cell_air_12000.json`、`--factor 2`

   `Ux` の相対場差は **1.683983837120607e−4**、許容値は **1.665054617830001e−4**。許容値を **1.14%超過**し、`diff_res.py --tolfile` は **FAIL、exit 1** を返します。「1.7e−4ちょうど」は表示丸めを合否判定に持ち込んだ説明です。

   ノイズ床 JSON は `run_0027/0030` から完全に再現できました。旧 dry の床も再現でき、`run_0015o_dry_slip_cell_oldbin` 対 `run_0015_dry_slip_cell` は **PASS** です。

   **対案:** E1 の判定を未達へ訂正する。独立した反復を追加し、変数別ノイズの評価方法を明示して再判定してください。この小超過だけで物理的な退行とは断定できませんが、現在の基準で合格とは書けません。

3. **Major — 新しい面エンタルピー関数の絶対床が、受付可能な状態で EOS を変える。**

   **根拠:** [condensationEOS_d.cuh:75](/home/sano/work/forge-cond/solver_density_cuda/cuda_forge/condensationEOS_d.cuh:75) は `R_eff<1` を一律に1へ変更します。一方、[solverConfig.cpp:723](/home/sano/work/forge-cond/solver_density_cuda/input/solverConfig.cpp:723) の受付条件は正値であることだけです。

   コードの式で、`gamma=1.4`、`cp=1038.67062845`、`Yw=g=0.999`、`ρ=1`、`T=40 K` を評価すると、受付検査を通り、`R_eff=0.2598367`、`P=10.39347 Pa` です。しかし面側では：

   - EOS 温度：**40 K**
   - 面復元温度：**10.39347 K**
   - 静エンタルピー：EOS **−196.417 kJ/kg**、面側 **−255.598 kJ/kg**

   標準空気 `Yw=0.7671` の今回の run では発火しませんが、§4.1 が受付対象とする状態範囲で「同じ面状態の `h=e+p/ρ`」が成立しません。

   **対案:** 正で有限な `R_eff` をそのまま使用し、非正値・非有限値を異常として扱う。枯渇上限付近と `R_eff<1` を単体掃引へ追加してください。

4. **Minor — m2 の文書同期が未完了で、実装と異なる説明が残る。**

   **根拠:**

   - [methods/condensation.md:646](/home/sano/work/forge-cond/methods/condensation.md:646)：音速を `sqrt(γp/ρ)` と記載。実装は [dependentVariables_d.cu:244](/home/sano/work/forge-cond/solver_density_cuda/cuda_forge/dependentVariables_d.cu:244) の `sqrt(γRgas T)`。二相では異なります。
   - [case README:61](/home/sano/work/forge-cond/case/34.arthur_n2_nozzle/README.md:61)：入口が依然 `324.48 m/s・6.161 kg/m³`。実 config は `325.12 m/s・6.137 kg/m³`。
   - [procedures/recommended-settings.md:128](/home/sano/work/forge-cond/procedures/recommended-settings.md:128)：`cl ±500` に対して onset `∓0.7 K`。実測は **38.95→39.67→40.25 K** で同符号です。
   - [同:123](/home/sano/work/forge-cond/procedures/recommended-settings.md:123)：`condKantrowitz` の既定を1と説明しますが、コードの省略時既定は0です。推奨設定と区別が必要です。

   **対案:** 実装・実 config・再解析値に同期してください。`methods/index.md` と `plans/README.md` の対象項目追加は確認できました。

5. **Minor — §5.1 に、採用した検証要求と後続作業が残っていない。**

   **根拠:** [前回レビュー:78](/home/sano/work/forge-cond/notes/reviews/2026-09-13-condensation-air-result.md:78) は実装を通す境界試験と node/cell の質量・エネルギー流束収支を要求しています。しかし [test_cond_air.cpp:79](/home/sano/work/forge-cond/solver_density_cuda/tests/unit/test_cond_air.cpp:79) の slip 試験は、依然テスト内の double 代数式の比較です。

   また、[plan:132](/home/sano/work/forge-cond/plans/active/condensation-air.md:132) で後続とした `check_quasisteady.py` 統合が [§5.1:136](/home/sano/work/forge-cond/plans/active/condensation-air.md:136) にありません。#10 にも「物性ではない」という未立証の除外判断が残っています。

   **対案:** 未実施の採用事項を残作業表へ戻し、境界試験・流束収支を完了させる。後続統合も表へ追加し、原因については「今回の物性変更では偏差を解消しなかった」に限定してください。

確認できた成果もあります。M2 の分圧式、M3 の境界種別・`Yw` 有限性検査、M4 の node 抽出修正、SLAU 面エンタルピーの共通関数化は確認しました。既存単体実行ファイル4本は再実行で PASS、float 入力試験は **90/90 PASS、最大誤差0.0406 J/kg** でした。ただし再ビルドは実施していません。

[case の run 一覧](/home/sano/work/forge-cond/case/34.arthur_n2_nozzle/README.md:35) に対応する `run_0014–0032` と追加対照の計21 run・273保存場で、`VALUE/*` の NaN/Inf は0でした。収束判定は `run_0021_dry_slip_node` のみ **PASS**、他は **NOT CONVERGED**。専用 `onset_analysis.py --series` は全対象 **STEADY**、再取得6 run の `check_quasisteady.py --quantity machmax,pmax` も **STEADY** です。後者は onset・壁圧比の判定を代替しません。

最終保存場では、`case/34.arthur_n2_nozzle/run_0027_air_cpgcarrier_v2/` の onset **38.94 K**・壁圧比 **1.2134/1.3514/1.4547**、`case/34.arthur_n2_nozzle/run_0029_n2_new_v2/` の **39.52 K**・**1.2422/1.3715/1.4758** を再現しました。これは**未収束の準定常比較**です。使用メッシュの品質判定は cell/node とも **VERDICT: PASS** でした。

**推奨は `active` に留め、1・3のコード修正と異常系試験 → 2の回帰再判定 → 5の検証補完 → 4の文書同期の順に完了し、result レビューを再実施することです。** ファイルは変更しておらず、今回の指摘は **plan 未反映**です。

指摘数: Critical 0 / Major 3 / Minor 2
