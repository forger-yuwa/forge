# codex 諮問 (diagnose): mknot-min-theta

- **brief**: [`notes/reviews/briefs/2026-10-13-mknot-min-theta.md`](../../notes/reviews/briefs/2026-10-13-mknot-min-theta.md)
- **plan**: [`plans/active/tooling-nozzle-mknot-auto.md`](../../plans/active/tooling-nozzle-mknot-auto.md)
- **date**: 2026-10-13
- **commit**: `0e9ec170` (feature/nozzle-cooled-wall-settling)
- **codex**: effort `high`, 1.7 min, rc=0
- **結論**: **境目を「不成立の標本点」から「成立／不成立を挟んだ区間」に変更し、M5 の低端で N/2N の境目と最小点を同時に測る A/B を先に行う。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **境目を「不成立の標本点」から「成立／不成立を挟んだ区間」に変更し、M5 の低端で N/2N の境目と最小点を同時に測る A/B を先に行う。**

第 1 仮説: **［Major・要再検証］標本点から測る ΔM_safe = 0.05 が安全距離を過大評価し、解像度変更後の境目への接近を見逃す。** 確度: 中（距離の定義の欠陥は明確、実際の選択失敗は未確認）。

  根拠:
  
  - [plan:62](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-mknot-auto.md:62) は「最も近い不成立点」を境目としている。しかし [low_end_N2400.json:28](/home/sano/work/forge-deliv/case/70.m5_air_wt_demo/_mknot_sweep/low_end_N2400.json:28) では 1.97 が失敗、[low_end2_N2400.json:3](/home/sano/work/forge-deliv/case/70.m5_air_wt_demo/_mknot_sweep/low_end2_N2400.json:3) では 1.98 が成功。局所的に成立性が単調だとしても、分かるのは境目がこの間にあることだけ。1.97 基準で選ぶ 2.02 は、実際の境目から 0.05 離れているとは限らない。
  - 失敗条件は [wall_axismach.py:320](/home/sano/work/forge-deliv/design/forge_design/geometry/wall_axismach.py:320) の `abs(k0_ratio - 1.0) > 0.15`。これは QA の閾値であり、物理的な成立限界そのものではない。
  - 観測された 2.05 の最小点については、境目まで約 0.07〜0.08 あり、**現時点で危険とは判定できない**。N/2N の既存検証は 2.40〜2.60 の別領域だけである（[plan:119](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-mknot-auto.md:119)）。

  対案: **境目の区間同定を解像度ゲートに組み込む。** 成立／不成立の隣接区間を幅 ≤ 0.001 まで詰め、N/2N のうち成立領域を狭める側の「成立側端点」から距離を測る。0.05 は当面の設計要求として保持してよいが、測定済みの安全保証とは扱わない。

  反証条件: 固定した M5 入力で、N/2N の選択差 ≤ 0.01、両方で全制約成立、選択点から上記の保守的な境目まで ≥ 0.05 を満たすこと。この場合、この入力で境目への接近が選択を不安定化する仮説は棄却する。

第 2・第 3 仮説: 追加しない。M6 の同点規則による失敗を最優先にする証拠は弱い。N=600 の L_c=45 では M_K=2.6 と 2.7 の θ_max 差は約 0.00724°で、同点許容 0.001°より大きい（[summary.csv:91](/home/sano/work/forge-deliv/case/42.isobutane_wt/sweep_axislaw_A_MK_summary.csv:91)、同ファイル:92）。

判別 A/B: `r_t` と他の geometry を固定し、**変える設定は `geometry.n_axis_inv` の 2400 → 4800 だけ**。M_K=1.90〜2.20 を 0.01 刻みで設計評価し、成立性が切り替わる区間だけ幅 ≤ 0.001 に詰める。CFD は不要。全制約、境目区間、θ_max、同点処理後の選択値を比較する。  
→ **A：選択差 ≤ 0.01・全制約成立・保守的距離 ≥ 0.05**なら、この固定入力の境目起因の不安定性は除外。  
→ **B：どれか不成立**なら、現行規則の安定性は未成立。境目が動いたのか、θ_max の最小点だけが動いたのかを同じ評価表で区別する。

やらない方がよいこと: 根拠なく 0.05 を大きくすること。最大壁角の最小化から選択を遠ざけても、境目の測定誤差は解決しない。

呼び出し側の前提への異議: 「境目 ≈ 1.975」は測定値ではなく区間内の推定。「解像度 A/B 合格」は旧目的の領域に限る。また M6 の実装は制約内の単純な `min(theta_max_deg)`（[sweep_axislaw_A_MK.py:118](/home/sano/work/forge-deliv/case/42.isobutane_wt/sweep_axislaw_A_MK.py:118)）であり、安全距離と許容同点処理まで含めて「M6 と同じ規則」とは言えない。

不足情報: 低端領域の N=4800 の評価、両解像度の曲率 QA の未丸め値、M″ 制約の判定。今回は固定 `r_t` の診断であり、`r_t` 変化への頑健性は未検証。**plan 未反映（依頼どおり read-only。反映先は §4.1・§4.2・§6）。**
