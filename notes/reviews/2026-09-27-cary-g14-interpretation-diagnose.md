# codex 諮問 (diagnose): cary-g14-interpretation

- **brief**: [`notes/reviews/briefs/2026-09-27-cary-g14-interpretation.md`](../../notes/reviews/briefs/2026-09-27-cary-g14-interpretation.md)
- **date**: 2026-09-27
- **commit**: `192d0448` (feature/gap-heating-precision)
- **codex**: effort `high`, 4.1 min, rc=0
- **結論**: **追加の CFD 計算より先に、G14 の固定 `Re_H` 点で前縁入熱の復元 A/B と補間範囲の監査を完了する。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **追加の CFD 計算より先に、G14 の固定 `Re_H` 点で前縁入熱の復元 A/B と補間範囲の監査を完了する。**

採否表:

| 対象 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| R を条件付き予測偏差として報告 | **採用** | 記録値は R = 0.998–1.363。H の端点比は 0.896–1.028、L は 0.754–0.840。この記述まではよい。ただし **H に壁温依存がない、L では SST が冷壁を過大評価する、とは言えない**。原因帰属を禁じる [acceptance.json:406](/home/sano/work/forge/case/59.flatplate_cary_m6/acceptance.json:406) を維持する。 |
| H/L の効果量差を Reynolds 数依存と読む | **却下・Major** | 比較域が H: 4000–5600、L: 3050–3350 と異なり、壁温比の上端も 0.7／0.6 と異なる。H 内の単位 Re も一定ではない（[acceptance.json:134](/home/sano/work/forge/case/59.flatplate_cary_m6/acceptance.json:134)）。対案は**別の条件群として併記**し、Re と壁温の交互作用とは呼ばない。 |
| 非単調性を「実験の反復性10%以内」で説明 | **要再検証・Major** | H の `Re_H=5200` では R(0.4)/R(0.2) = 1.363/1.128 = **1.208**（[README.md:58](/home/sano/work/forge/case/59.flatplate_cary_m6/README.md:58)）。反復性10%は、系列間比の信頼区間でも総不確かさでもない。対案は、実測 St が縦軸と積分座標の双方に入る相関を含めて扱い、現状は「説明未確定」とする。 |
| 格子差 +5% を全系列の誤差幅とする | **却下・Major** | H・Tw/Tt=0.2 の**旧17点平均**だけの差であり、G14各点や壁温間比の格子感度ではない。さらに `run_0005_tw02_re027_fine` は壁解像のツール判定が **FAIL** のまま（[README.md:32](/home/sano/work/forge/case/59.flatplate_cary_m6/README.md:32)）。対案は暫定表にこの制約を明記する。格子に依存しない壁温傾向を主張するには、少なくとも L の両端条件で同じ格子変更を比較する必要がある。3段目は格子収束・誤差推定に進む場合に必要で、現状の記述的な表のために全9系列へ一律追加する必要はない。 |
| D-7275 を位置 x で比較 | **採用。ただし受け渡しの過大解釈は Major** | **実測位置 x を一次座標**とし、位置 II の局所 `B=q_CFD/q_exp`、パネル平均は別量とする。前縁／トリップ起点は感度であり、物理的上下界ではない（[plan:902](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:902)）。欠測上流入熱を CFD で補った座標を「実験 Re_H」と呼ばない。Cary の R を補正倍率として持ち込まず、すきまでは実測分母と同条件の forge 平板分母を併記する（[plan:909](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:909)）。 |

第 1 仮説: **Tw/Tt=0.2 の旧17点で成立した復元感度の確認を、G14 の固定 `Re_H` 点へ一般化できない。**  
確度: **高。比較方法の問題として確認できたもので、物理的偏差の原因という主張ではない。**

  根拠: **Major — 比較する量が変わっている。** `cary_energy_compare.py` は実験位置を固定し、復元によって CFD の対応位置を動かす。一方、G14 は `Re_H` を固定し、実験側 St の対応位置が動く（[cary_energy_compare.py:103](/home/sano/work/forge/case/59.flatplate_cary_m6/tools/cary_energy_compare.py:103)、[cary_g14.py:62](/home/sano/work/forge/case/59.flatplate_cary_m6/tools/cary_g14.py:62)）。

  `conditions.json` と生成器の公称条件から、前縁外挿係数を最小→最大へ変えて独立再計算すると：

  - H・Tw/Tt=0.2：実験の有効下限が **3885.7 → 約4026.0**。固定点 **4000 は範囲外**。
  - L・Tw/Tt=0.2：同下限が **2997.4 → 約3157.4**。固定点 **3050 は範囲外**。
  - H・Tw/Tt=0.3：範囲内の5600でも **R_B/R_A ≈ 0.9484、約−5.16%**。

  現行の `np.interp` 呼び出しには範囲検査がなく、範囲外でも端点値が返る（[cary_g14.py:67](/home/sano/work/forge/case/59.flatplate_cary_m6/tools/cary_g14.py:67)）。**現行の腕Aの表が誤っていると断定する証拠ではないが、「外挿感度0.7%以下で確認済み」とする根拠は失われる。**

  反証条件: AWS の実効条件による同じ再計算で、両腕とも全固定点が有効範囲内にあり、点別感度も小さいことが示されること。遷移前1点の系列は、この方法では感度を評価できない。

第 2・第 3 仮説: **L 群の下流2点の減少は、今回の最小／最大係数による復元差だけでは消えない。確度: 中。** 両系列を最大係数へ変更した計算では、端点比は3200／3350で約 **0.822／0.859**（元は0.840／0.835）。ただし格子感度・実験不確かさ・遷移履歴は未分離。第3仮説は置かない。

判別 A/B: **変更するのは実験側の未計測前縁入熱 Q₀ だけ。追加 forge 計算は0 step。**

- 腕A＝事前登録どおりの最小係数、腕B＝同じ遷移前点の最大係数。他の積分法・物性・比較点は固定する。
- AWS の既存壁出力で、CFD側 `Re_H` の単調性と固定点の包含を確認する。実験側も両腕で包含を確認し、範囲外は数値を返さず**判定不能**とする。
- **全対象点で範囲内かつ差≤5%なら**、この復元範囲に限って第1仮説を退ける。**範囲外または差>5%なら**、G14全体へ方法確認を一般化する判断を退ける。上の独立計算は既に後者を示している。
- 5%は既存 T4-0a-E の方法感度基準を使う**今回の事後診断案**で、G14の事前登録済み合否基準とは呼ばない。1点しかない系列の A=B を「感度ゼロ」と扱わない。

やらない方がよいこと: **Major — 比較点の事後削除、実験に近い外挿腕の選択、Prₜ の調整、Cary の倍率による TP-1187 の補正。** 根拠は [acceptance.json:111](/home/sano/work/forge/case/59.flatplate_cary_m6/acceptance.json:111) と [plan:910](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:910)。対案は元の表を維持し、感度・未評価・範囲外を追加表示すること。

呼び出し側の前提への異議: **「全9系列の収束・定常性が確認された」ことと「比較方法・格子・壁解像まで成立した」ことを分ける。** また、Cary 原報の「冷壁で Van Driest が過大」という結論には比較規約の条件がある。原報本文p.12自身が、エネルギー厚さ座標と Reynolds analogy factor=1.0を用いた別整理では予測評価が変わると記している。現在の `Re_H` 比較の期待値へ、そのまま移植しない。

不足情報: AWS の実効設定・`stage_manifest.json`・残差・壁時系列・判定ファイルがローカルにない。**`PASS (converged)`／`ALL STEADY` は [run索引](/home/sano/work/forge/case/59.flatplate_cary_m6/README.md:25) の記録の引用であり、独立再判定していない。** 特に、判定した St 系列が G14 の固定点に対応することと、CFD側の補間範囲を確認する必要がある。

ファイル変更・forge 起動なし。**plan 未反映。**
