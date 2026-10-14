# codex 諮問 (diagnose): vgc-result

- **brief**: [`notes/reviews/briefs/2026-10-13-vgc-result.md`](../../notes/reviews/briefs/2026-10-13-vgc-result.md)
- **plan**: [`plans/active/tooling-nozzle-isothermal-wall-chain.md`](../../plans/active/tooling-nozzle-isothermal-wall-chain.md)
- **date**: 2026-10-13
- **commit**: `90c5c58e` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 4.0 min, rc=0
- **結論**: **次は壁温と物理壁を一つに固定し、新しい G0・上限 1 の参照を伴う G0c・上限 1/50 の A/B を事前登録して、壁更新を始める前に適用範囲を確かめる。**
- **extra**: `case/45.isobutane_m6_d155/vgc_judge.py`, `case/45.isobutane_m6_d155/_band_ab/cold_pair/vgc_judge.json`, `design/forge_design/feedback/deltastar_loop.py`, `plans/active/architecture-float-state-double-geometry.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：

**V-g0c-chain の結果 A は採用する。#29〜#31 は「登録した 300 K・固定形状・開始手順での実装と検証完了」として閉じてよい。壁温変更後と δ* 反復への適用は別の残作業として残す。** Critical に相当する根拠は確認できなかった。

| 重大度 | 論点・採否 | 根拠と対案 |
|---|---|---|
| Major | **結果 A：採用** | 保存 JSON の数値を `sp_judge.py` の閾値で再判定し、T2a・T2b の双方について全比較が許容内と確認した。θ_r の最大差は **0.075843 %**、代理 Q_w は約 **−0.040 %**、出口の最大 \|ΔM\| は **0.000503**。「登録した開始手順から停止規則に到達し、停止出力の比較要件を満たした」と記す。[判定実装](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/sp_judge.py:44) |
| Major | **#29〜#31 の完了：範囲を限定して採用** | 登録自体が壁温変更・全 pass・float への拡大を除外している。[plan §6](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-isothermal-wall-chain.md:309)。3 項目の完了と、plan 全体・設計チェーン全体の検証完了を混同しない。停止監視のチェーン組み込みも未完のまま残す。 |
| Major | **「準定常の設計入力が得られた」：却下** | 本段の保存判定は `OVERALL: ALL PASS` だが、100,000〜120,000 step の θ_r 3 断面は **DRIFTING**、代理 Q_w は **TRANSIENT-UNSETTLED**。[保存判定](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/cold_pair/vgc_judge.json:201)。残差の合格は保持し、派生量の定常化を別に判定する。起点残差が大きいことだけで PASS を無効化する必要もない。 |
| Major | **新壁温で G0c の上限 1/50 だけを比較すれば十分：却下** | この対では、両腕に共通する粗格子の影響を検出できない。従来の採用条件も G0・上限 1 との直接比較を要求している。[float plan §6.58](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:1507)。新壁温の G0・上限 1 を共通参照として追加する。旧 T2 は新壁温の合否判定には使わない。 |
| Major | **既存判定器の壁温・形状変更への流用：要再検証** | `cold_xcheck.py` は旧 run の `scale` と旧断熱／300 K の抽出帯を読み、`sp_judge.py` は出口を同じ行番号で比較する。[抽出条件](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/cold_xcheck.py:46)、[出口比較](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/sp_judge.py:55)。新しい r_t・節点配置・境界層に対して抽出帯と座標対応を確認する。代理 Q_w・τ_w は実物性による熱負荷・壁応力の保証には使わない。 |
| Major | **`run_pass` をそのまま自動反復させる：却下** | `run_pass` は solver の戻り値を表示した後に `collect` を行い、`check_convergence` も `check=False`、その後に次壁を抽出する。準定常判定や今回の停止監視は通過条件になっていない。[deltastar_loop.py:381](/home/sano/work/forge-integ-1005/design/forge_design/feedback/deltastar_loop.py:381)。当面は pass ごとに監視・判定し、壁更新に使う **δ_E・δ_target 自体**の時系列を確認してから次へ進む。 |
| Major | **FP64 採用と `thermoFloat` 基準の変更を同時に行う：却下** | `thermoFloat=1` の hybrid 経路は実際に float Newton を呼ぶ。[thermo_d.cuh:728](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/thermo_d.cuh:728)。したがって FP64 ビルドを「全演算 double」と呼べない。一方、これだけで今回の差の原因とも言えない。新壁温試験では検証済みバイナリと実効設定を固定し、float の置換・main 統合・`thermoFloat=0` の基準化は保留する。 |
| Minor | **復元・時間・AR の但し書き：採用** | `run_0062` の復元は今回生成した .msh と変換後格子の照合で裏付ける範囲に限る。本段 **8,884 s** は共有 GPU 上の到達時間で、速度倍率に使わない。品質の保存記録は **`--ar-max 10000`: PASS**、AR>5000 は壁から 0〜5 層の **1,325 セル**、スキュー最大 **0.002435**。一般の AR 基準変更ではない。[結果記録](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-isothermal-wall-chain.md:317) |

根拠 run は `case/45.isobutane_m6_d155/run_0810_vgc_chain/`、比較対象は同 case の `run_0752_t2_d1/`・`run_0754_t2_d2/`。保存記録では停止比較に `res_120000.h5` を使う一方、残差履歴は step **121389** まである。したがって `ALL PASS` は「本段の保存履歴に対する判定」と記し、120,000 step で切った判定と同一視しない。

結論: **次は壁温と物理壁を一つに固定し、新しい G0・上限 1 の参照を伴う G0c・上限 1/50 の A/B を事前登録して、壁更新を始める前に適用範囲を確かめる。**

第 1 仮説: 新壁温でも、G0c・B0 は G0・上限 1 に対する停止時の比較要件を満たす。確度: **低**。  
根拠: 300 K では `run_0810_vgc_chain/` が異なる開始手順から比較要件を満たした。ただし壁温を変えた証拠はなく、初期場に対する頑健性と壁温に対する頑健性は別である。  
反証条件: 新壁温で G0c・上限 50 が予算内に到達しない、発散する、または G0・上限 1 との直接比較で登録許容を超える。これは運用要件の反証であり、原因を上限や格子へ即断しない。

第 2・第 3 仮説: 追加しない。低マッハ前処理や float の原因調査は今回へ混ぜない。

判別 A/B: **主 A/B で変える設定は `lineDtDirectionalCap: 1 → 50` の一点。両腕とも G0c。新壁温の G0・上限 1 を共通参照として 1 本作る。計 3 本で、初回から反復実行は足さない。**

- **固定するもの**：新壁温 T_w₁、物理壁、r_t、入口・出口条件、物性、FP64 バイナリ、リミッタ基準、段階起動、出力間隔。δ* 壁更新は試験中に行わない。共通ドナーを SHA で固定し、G0c の対は同じ S2 後入力から開始する。
- **長さ**：各本段は最大 **200,000 step**、出力 **2,500 step** ごと、既存 §6.11 の停止規則。短い固定窓の近さで早期合格にしない。入力・抽出不成立と数値的失敗を分ける。
- **見る量**：既存 M1〜M7 を同じ許容で比較する。θ_r・δ_loc・代理 Q_w・流量差は **0.1 %**、壁分布は既存の混合許容、出口 M は最大 **0.002**・平均 **0.001**。G0c の両腕を G0 参照へ**それぞれ直接比較**し、比較の連鎖で代用しない。ゼロ近傍の δ_loc・Q_w では相対誤差が破綻するため、該当時の絶対許容を結果を見る前に定める。
- **A**：全腕到達、主 A/B と G0 への直接比較がすべて許容内 → **この壁温・固定壁の停止時運用**への拡大を支持する。
- **B**：到達腕の許容超過、発散、予算内未到達 → 拡大を見送る。G0c の上限間だけ合格して G0 参照との比較が落ちた場合も、G0c 採用は合格にしない。欠損・設定不一致・抽出不能は判定不能。

この試験の A と、δ* 壁更新に使えることは別判定にする。全残差の `check_convergence --segment`、θ_r・Q_w に加え **δ_E(x)・δ_target(x)・出口 M の `check_quasisteady`** を残す。壁更新の許容誤差に対して時間変動が十分小さいことを、絶対寸法で事前登録する。現在の M1〜M7 には、次壁として使う δ_target 全分布の保証がない。

各 pass では、最低限次を確認する。

1. **格子と形状**：新しい壁に対して生成し、実際の区間端・列数・r_t・第一層距離を記録する。`x_respace` は壁曲率や流れを見て区間を選び直す機構ではない。[mesh2d.py:176](/home/sano/work/forge-integ-1005/design/forge_design/meshing/mesh2d.py:176)。区間外の不変性は、その pass の G0 と G0c の間で確認する。
2. **品質と壁解像**：品質 VERDICT、高 AR セルの位置・直交性、局所 y₁⁺ の最大・超過面積を確認する。AR≤10,000 だけで旧 G0c と同じ資格とはしない。
3. **設定と移植**：S1・S2／本段の実効設定、壁温、物性、移植元出力と SHA を保存する。形状変更時は新格子の壁距離を使う。
4. **次壁への受け渡し**：δ_E・δ_target の定常性、抽出の不合格断面、平滑化、物理壁の単調性を確認する。形状変化の許容範囲は未検証なので、正規化区間が同じという理由だけで全 pass を認定しない。

やらない方がよいこと: **新壁温を旧 T2 と比較して採否を決めること、G0c 内の上限比較だけで粗格子の影響を除外すること、DRIFTING の場から自動で壁更新を重ねること、今回同時に float・`thermoFloat`・既定値を変えること。**

呼び出し側の前提への異議: **「出口 M が高いのは壁未調整だから」は今回の検証からは確定できない。** 壁温変更後の固定壁評価と、壁更新後の評価を分けて確認する。また、M7 は腕間の流量差であって収支閉合ではなく、代理熱流束は実熱負荷の保証ではない。

不足情報: **新しい壁温の値、壁更新に許す絶対寸法誤差、各 pass の適用範囲**。対象 run の実ディレクトリは当地になく、HDF5・残差履歴からの独立再判定は未実施。今回確認できたのはコード、plan、保存 JSON とその数値の閾値照合である。

**plan 未反映。** 呼び出し側の反映先は `tooling-nozzle-isothermal-wall-chain.md` §4.8・§5.1 #29〜#31・§6・§8、および float plan §5.1 #28。ファイル変更・forge 起動は行っていない。
