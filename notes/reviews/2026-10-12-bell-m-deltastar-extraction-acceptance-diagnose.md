# codex 諮問 (diagnose): bell-m-deltastar-extraction-acceptance

- **brief**: [`notes/reviews/briefs/2026-10-12-bell-m-deltastar-extraction-acceptance.md`](../../notes/reviews/briefs/2026-10-12-bell-m-deltastar-extraction-acceptance.md)
- **plan**: [`plans/active/tooling-nozzle-bell-two-mode.md`](../../plans/active/tooling-nozzle-bell-two-mode.md)
- **date**: 2026-10-12
- **commit**: `0f9e2753` (feature/thruster-nozzle-design)
- **codex**: effort `high`, 5.8 min, rc=0
- **結論**: ?
- **extra**: `notes/reviews/2026-10-12-bell-m-deltastar-calibration-diagnose.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0、Major 4、Minor 1）。**現行抽出は不受理を維持し、較正と L70 の開封は保留する。**

| 論点 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| 0.75倍を除き、外向きだけで受理する | **却下／Major M1** | 再抽出で789断面すべて不受理、最大61.96%を再現した。0.75倍帯は718断面で既存の勾配条件1%を超えるため、境界層内への侵入という説明は支持される。ただし、固定帯では再試行が停止し、この条件違反が不合格理由に残らない。[`deltastar.py:267`](/home/sano/work/forge-thruster/design/forge_design/metrics/deltastar.py:267)、[同:300](/home/sano/work/forge-thruster/design/forge_design/metrics/deltastar.py:300)。**無効な帯の値を通常の不確かさと混ぜる問題と、主系列の帯選択の問題を分け、帯選択を検証する。** |
| x < 0.9を較正区間から外す | **却下／Major M2** | x≈0.50024では帯を2倍にすると抽出値が9.29%増える。これは除外の物理的根拠ではなく、抽出不足を疑う証拠。登録区間はスロートから出口である。[`bell_deltastar.py:63`](/home/sano/work/forge-thruster/design/forge_design/metrics/bell_deltastar.py:63)。**まず全区間で抽出を直す。** |
| 定数1個の較正方式／現時点でkδを確定する | **方式は採用、確定は要再検証／Major M3** | 未受理の抽出と時系列未判定の場への当てはめでは、抽出偏りを係数へ吸収する。さらに現在の参考最良定数は断面数による等重みで、格子の集中に依存する。[`bell_deltastar.py:282`](/home/sano/work/forge-thruster/design/forge_design/metrics/bell_deltastar.py:282)。**受理後、固定区間・x方向の積分重みで較正し、許容される抽出系列による係数の幅を併記する。** |
| 代表点と現在のL2だけでM3を判定する | **要再検証／Major M4** | 現在の`dr_l2`は断面値の単純RMSで、分布の局所変化を隠せる。[`bell_deltastar.py:169`](/home/sano/work/forge-thruster/design/forge_design/metrics/bell_deltastar.py:169)。**代表点・出口・x積分によるL2に加え、全断面の時系列を判定する。** |
| Pr=0.72を燃焼ガスと完全に整合した物性と扱う | **却下／Minor m1** | 回復係数は既定Pr=0.72の立方根。[`deltastar_integral.py:74`](/home/sano/work/forge-thruster/design/forge_design/feedback/deltastar_integral.py:74)、[同:204](/home/sano/work/forge-thruster/design/forge_design/feedback/deltastar_integral.py:204)。閉包の近似として明記し、較正前に扱いを固定する。燃焼室のPr≈0.709へ置換するだけで全温度域の整合が取れるとはしない。 |

**結論:** L100の固定場で、帯選択だけを`adaptive`から既存の`edge`へ変更し、元の倍率・5%線・全補正区間を維持した比較を行う。

**第1仮説:** 現行の`adaptive`は抽出値が十分に安定する前に停止しており、内向き感度の大半は境界層内で参照分布を当てはめ直した影響である。 **確度: 高**

根拠となる場は次の2本。以下はすべて**最終スナップショットの診断値であり、定常値とは認定しない**。

- `case/68.bell_two_mode/run_0019_tp_m_ns_l100/res_60000.h5`
- `case/68.bell_two_mode/run_0018_tp_m_euler_l100_pin/res_30000.h5`

x/r_t=10.887305では、選択帯の`band_slope`は0.00973で既存の停止条件0.01を満たす。一方、0.75倍帯では0.17802となり、δ_rが61.96%減る。1.5倍・2倍帯ではδ_rがそれぞれ3.148%、3.160%増える。停止条件は帯内の比の変化を見ており、δ_r自身の帯位置依存性を保証しない。[`deltastar.py:304`](/home/sano/work/forge-thruster/design/forge_design/metrics/deltastar.py:304)。

反証条件: 境界層外の有効な帯であることを確認しても同程度の大差が残り、連続した帯位置の範囲で抽出値が安定しないこと。その場合、単なる帯の早期停止では説明できない。

**第2仮説:** 帯を遠ざけると境界層外の分布変化も参照外挿へ取り込み、「外向きなら有効」という扱いも破綻する。 **確度: 中**

x/r_t≈0.751669で、元の抽出値に対する増加は2倍帯で6.25%、3倍で5.77%、4倍で21.46%。この観測は外向き感度の非単調性を示す。原因を波・曲率・補間のどれかに確定する証拠はまだない。

**第3仮説:** δ_rに時間変化が残り、較正係数がまだ動く。 **確度: 低・未確認**

手元の両runには最終場しかなく、δ_rの時系列VERDICTは作れない。

**判別A/B:** 同じ2枚の場を使い、変更は`band_select="adaptive" → "edge"`の1点。**追加CFDは0 step**。倍率`(0.75, 1.5, 2.0)`、受理線5%、区間x≥0、補間・積分・平滑化設定は固定する。CONTURとの一致は帯選択に使わない。

- **A：** 全789断面で感度≤5%、有効な局所縁が存在し、主系列と感度帯の診断も適合するなら、現行`adaptive`の選択位置が不受理の主因という説明を支持する。新しい抽出をM3判定へ進める。
- **B：** 1断面でも満たさなければ、「既存`edge`への切替だけで解消できる」を棄却する。較正は保留し、不成立位置を記録する。

`edge`には欠測縁の平滑補間とx<3の端値保持があるため、**有限な出力が出ただけではAとしない**。`edge_y_star`の欠測、平滑化前後の位置、最終的な帯の診断を確認する。[`deltastar.py:415`](/home/sano/work/forge-thruster/design/forge_design/metrics/deltastar.py:415)。このA/Bは今回未実施で、呼び出し側で実行前に登録する提案である。

**やらない方がよいこと:** 0.75倍系列を削除して旧結果をPASSへ変更すること、x≥1への切詰めで受理すること、CONTURに近い帯を選ぶこと、現在の参考最良定数1.0279を確定値としてL70へ進むこと。

**呼び出し側の前提への異議:**

1. **旧FAILは保持する。** 61.96%を「有効な抽出の誤差幅」と解釈するのは不適切だが、それは現行抽出を受理する根拠にはならない。今回は受理線を緩めず、抽出方法の改訂として扱う。倍率や区間を将来変更するなら、旧登録の不受理と新しい定義を別記する。

2. **kδの不確かさは空間的なばらつき1.3%ではない。** 抽出・M3受理後、固定区間で  
   `kδ = Σ(wᵢ Cᵢ dᵢ / sᵢ²) / Σ(wᵢ Cᵢ² / sᵢ²)`  
   とする。CはCONTUR、dは受理したNS抽出、wはΔxの台形積分重み、sは事前に固定した相対誤差の分母と床。許容帯系列・定常末尾系列ごとに同じ式で求めたkδの最小～最大を、系統的な感度幅として報告する。信頼区間とは呼ばない。

3. **M3の新規登録案は`--tail 0.4 --drift 0.005 --osc 0.01 --min-snaps 20`。** 同じ実効設定の区間で20点以上、末尾に8点以上を確保する。対象は既存の代表x・出口、x積分によるL2、および全補正断面。全対象で`STEADY`を要求する。単調系列は漸近値との差も1%以内とし、減衰外挿が成立しないものを「残差0」と扱わない。これは抽出の5%予算に対し、時間変化を1%以内に抑えるための**未検証の提案値**である。

4. **L70の開封条件を固定する。** L100の抽出受理、M3受理、kδと感度幅、CONTUR設定、補正区間、誤差の重み・分母を凍結してから開封する。L70自身にも抽出・M3ゲートを適用し、その後に既登録の出口10%・RMS10%・局所20%で判定する。抽出ゲートに落ちた場合はモデル移植性を判定不能とし、L70への再較正で成功へ変えない。

**不足情報:** L100のNS中間場、Euler参照の局所分布の時間安定性、改訂後の帯診断、δ_rの準定常VERDICT。

確認済みの判定は、両runの保存残差VERDICTが判定区間`S3_main`で **`NOT CONVERGED (stalled/plateau…)`**。残差CSVは手元になく再実行していない。推力履歴は今回ツールを再実行し、両runとも **`OVERALL: ALL STEADY`**。これをδ_rへ流用できない。

ファイル変更・forge起動・L70のδ_r閲覧は行っていない。**plan未反映**。呼び出し側の反映先は[`tooling-nozzle-bell-two-mode.md`](/home/sano/work/forge-thruster/plans/active/tooling-nozzle-bell-two-mode.md:147)の§4.4・§5.1 #7・§6。
