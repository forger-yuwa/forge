# codex レビュー: tooling-nozzle-cfd-pinned-initial-line (plan)

- **plan**: [`plans/active/tooling-nozzle-cfd-pinned-initial-line.md`](../../plans/active/tooling-nozzle-cfd-pinned-initial-line.md)
- **stage**: `plan`
- **date**: 2026-10-05
- **commit**: `c0f37689` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 8.5 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M6/m2
- **extra**: `case/45.isobutane_m6_d155/cfd_initial_line.py`, `case/45.isobutane_m6_d155/test_cfd_initial_line.py`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**GO-with-changes**

CFD 初期線の抽出は実現可能ですが、採用判定と入力場の適格性確認が不足しています。  
再計算では V1 の点上ゲートも不合格でした。以下を修正し、形状ゲートを通すまで Euler 投入には進めません。

目的と差分は妥当です。[既存 A8](/home/sano/work/forge-integ-1005/plans/accepted/tooling-nozzle-axismach-throat-characteristic.md) は Hall 場からの特性線生成であり、今回の CFD 初期線取得は未解決部分です。変更を `design/forge_design` に閉じる構成も自然で、CUDA・block-DPLUR・周期境界の改修は不要です。node 限定検証は現行の検証手順に合っています。

なお、作業ツリーには貼付版より新しい試作がありました。保存処理を除いてメモリ内で再実行し、合成 Hall 試験の全項目合格と、以下の V1 不合格を確認しました。ファイルは変更していません。

1. **Major — V3 は「劣らない」を証明できない採用判定になっている。**

   **根拠:** [plan:85](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-cfd-pinned-initial-line.md:85) は、改善・悪化を検出できなかった場合を「差なし」とし、形状で採用できます。しかし、差 `+0.02 %pt`、不確かさ `U=0.04 %pt` でもこの分岐に入り、許容差 `0.03 %pt` 内の非劣化は保証できません。

   また、指定評価器の [eval_wallfit_euler.py:152](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/eval_wallfit_euler.py:152) は `P_slope` を判定対象へ割り当てておらず、傾きの `U`・合否を生成しません。出口 Mach で規格化したオーバーシュートも未実装です。

   **対案:** 各量に許容悪化幅 `Δq` を登録し、全量で `差＋U ≤ Δq` を採用条件にする。不確かさが大きい場合は「判定保留」とする。傾きの絶対値化、規格化オーバーシュート、準定常の前提条件を含め、評価器改修を §5・§7 に追加する。

2. **Major — 凍結する入力の時間安定性と、V2 の自己整合判定が不足している。**

   **根拠:** 入力元  
   `case/45.isobutane_m6_d155/run_0062_euler_wallfit_fit_r1_ext6k/CONVERGENCE_VERDICT.txt`  
   の保存判定は **`NOT CONVERGED (stalled/plateau — needs scheme change, not more steps)`**。全残差の低下は約 0.3〜0.6 桁です。手元に `residual_history.csv` がなく、再実行した判定器も `NO residual_history.csv` を返しました。

   一方、対照 `run_0056_euler_wallfit_fit_r1/wallfit_series.csv` の再判定は **`OVERALL: ALL STEADY`**。これは試験部指標の判定であり、初期線・軸微分の定常性は保証しません。[plan:86](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-cfd-pinned-initial-line.md:86) は初期線の時間変動を別項へ送り、[V2:84](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-cfd-pinned-initial-line.md:84) には線の座標 `x(r)`・`x₀`・`M′_A` の判定がありません。

   **対案:** 線・アンカー・m* の時系列確認を凍結前へ移す。V2 は凍結した全出力を比較し、`D₁>D₀` は各許容差で無次元化した同一ノルムで定義する。未収束場による実験と、採用資格の検証を分け、後者には残差判定の成立を要求する。

3. **Major — 固定した壁表現では、登録済み V1 が既に不合格。**

   **根拠:** `run_0062/.../res_6000.h5` を入力に [cfdpin_pass1.py:37](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/cfdpin_pass1.py:37) の計算を再実行しました。

   | 指標 | 再計算値 | 登録基準 |
   |---|---:|---:|
   | 最大位置誤差 | `7.7746e-6 r_t` | `≤5e-6 r_t` |
   | 最大角度誤差 | `0.0124427°` | `≤0.005°` |
   | r″ の最大値 | `0.510160` | Hall の `0.534085` 未満 |

   **`V1_point_gate: false`、`all_pass: false`** です。曲率の改善では位置・角度ゲートの不合格を相殺できません。

   **対案:** §5.1 に「V1 不合格時の MOC 点群と壁当てはめの整合確認」を追加し、これを Euler 前の必須作業にする。閾値の事後緩和や曲率だけによる通過は認めない。壁表現の変更が必要なら、除外したスコープを先に改訂する。

4. **Major — `_flux_along` の m* を CFD の実流量と同一視できず、出口 Mach の符号による棄却にも根拠がない。**

   **根拠:** [moc_inverse.py:335](/home/sano/work/forge-integ-1005/design/forge_design/geometry/moc_inverse.py:335) は M と θ から、指定ガスの**等エントロピー流束**を再構成します。CFD の保存量 `roUx`・`roUy` を積分していません。

   単一 snapshot の変換誤差診断として、`run_0062/.../res_6000.h5` の同じ線を積分すると、保存量補間からの流束とこの再構成流束には相対 **`1.88〜1.90e-4`** の差がありました。これは V2 の m* 許容差 `2e-4` と同程度です。定常値の主張ではありません。

   また、再計算では m* が `+4.20594e-4`、出口面積が `+4.20645e-4` とほぼ同率で増えます。面積も変える逆設計なので、[plan:85](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-cfd-pinned-initial-line.md:85) の「出口 M が下がれば m* の符号が逆」は導けません。

   **対案:** m* を「MOC ガスモデルへ写像した流束」と明記し、保存量からの独立積分との差、組成・h0 の整合、格子感度を誤差予算へ入れる。出口 M は目標値からの偏差で判定し、変化の符号だけによる棄却条件を削除する。

5. **Major — 抽出不能な場を拒否する契約がない。**

   **根拠:** [試作:58](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/cfd_initial_line.py:58) は `max(M,1+1e-9)` により亜音速を隠します。同じ格子へ一様 `M=0.9` を与えた実行でも、例外なく `x₀=4.47213e-5` の「特性線」を返しました。既存 [HallThroat:297](/home/sano/work/forge-integ-1005/design/forge_design/geometry/transonic.py:297) は `M≤1` を明示的に拒否しています。

   読込側も [試作:19](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/cfd_initial_line.py:19) では node/Euler/軸対称、ガス条件、snapshot とメッシュの対応を確認しません。

   **対案:** 対応範囲を構造化 node・軸対称 Euler に限定し、不適合入力を入口で拒否する。RK の全評価点で有限値・超音速・領域内を検査し、軸未到達も失敗にする。`start_line: vertical` や旧アンカー指定との併用規則、明示 snapshot・入力ハッシュも定義する。

6. **Major — 固定量と比較座標の宣言が実装と一致していない。**

   **根拠:** [plan:53](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-cfd-pinned-initial-line.md:53) は `x_F 不変` としていますが、既存実装は特性線と壁の交点から出口を決めます。実際の再計算では壁末尾が `95.22667765 → 95.22725644 r_t` に動き、特性線交点はさらに別の `95.29650227 → 95.29708106 r_t` です。

   [評価器:18](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/eval_wallfit_euler.py:18) は旧 `X_E/X_F` を固定しています。また「M″_A 据え置き」は現試作では**同じ Hall 式を新座標で評価**する意味で、値は `0.0213539 → 0.0245527`、約 **15%** 変わりました。「変数は初期線・m*・M_A だけ」という V3 の説明は不正確です。

   **対案:** 今回は `L_c` 固定・出口位置可変と明記する。特性線交点とメッシュ出口を別記録し、両腕に共通する物理区間で比較する。M″ の評価座標と、変更するアンカー全成分を明記する。

7. **Minor — 動機の因果帰属を強く書きすぎている。**

   **根拠:** [verification plan:405](/home/sano/work/forge-integ-1005/plans/active/verification-m6-axis-wave-mesh-su2.md:405) 自身、V0→V4b を「スロート直上流の接続＋始点処理」の変更としています。上流壁だけの寄与を分離した比較ではありません。また、CFD ピンで傾きの差が検出されないことから、V4b の効果を上流 blend に帰属させる [plan:85](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-cfd-pinned-initial-line.md:85) の解釈も成立しません。

   **対案:** 動機を「スロート近傍の複合変更への感度」と記述する。差が不確かさ以下なら「今回の比較では効果を分離できない」と報告する。

8. **Minor — スコープ・工程・仕様更新の順序が矛盾している。**

   **根拠:** §2 は最大3回の固定点反復、§4 は凍結＋1回検証、§5 は再び固定点 pass です。[plans/README.md:66](/home/sano/work/forge-integ-1005/plans/README.md:66) も反復方式のままです。さらに §5 は methods 更新を最後に置き、実装前に仕様を更新する開発規則と逆順です。

   **対案:** 全節を凍結＋1回検証へ統一する。methods 更新、入力契約・評価器仕様、provider、形状判定、Euler の順へ直し、失敗時の停止条件も残作業表へ入れる。

推奨は、**Hall を既定に保ち、CFD ピンを「凍結＋1回検証」の限定実験として進める**ことです。実装前の優先順位は、①採用判定の修正、②入力場と自己整合の資格条件、③m*・固定量・入力拒否規則、④工程の同期です。その後も、上記 V1 不合格の解消を Euler 投入の条件としてください。生産採用と格子誤差の帰属は、その検証後に判断すべきです。

本レビューは **plan 未反映**です。参照 run の恒久索引は [case/45 README](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/README.md:65) です。

指摘数: Critical 0 / Major 6 / Minor 2
