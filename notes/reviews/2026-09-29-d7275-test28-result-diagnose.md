# codex 諮問 (diagnose): d7275-test28-result

- **brief**: [`notes/reviews/briefs/2026-09-29-d7275-test28-result.md`](../../notes/reviews/briefs/2026-09-29-d7275-test28-result.md)
- **plan**: [`plans/active/case-hypersonic-gap-heating-validation.md`](../../plans/active/case-hypersonic-gap-heating-validation.md)
- **date**: 2026-09-29
- **commit**: `9f776451` (feature/gap-heating-precision)
- **codex**: effort `high`, 3.8 min, rc=0
- **結論**: **追加 CFD より先に、`case/60.flatplate_d7275_m7/run_0025_t28_y2/` の壁解像を全点で再集計し、「前縁1点だけ」という完了判断の根拠を監査する。**
- **extra**: `case/60.flatplate_d7275_m7/README.md`, `case/60.flatplate_d7275_m7/acceptance.json`, `case/60.flatplate_d7275_m7/digitize_d7275_fig19c_test28.json`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical なし）

| 判断・重大度 | 採否 | 根拠と対案 |
|---|---|---|
| 平均比から「関係は変わらない」と確定する — **Major** | **算術的な読みだけ採用** | 記録値の 1.411/1.321＝1.06813 は登録条件内。ただし、1.10 までの余裕は相対 **2.98%**。格子変更による +2.2〜2.4% は誤差上限ではなく、入口モデルの影響も未定量。「今回のモデル・格子の点推定では登録した10%以内」と限定する。[acceptance.json:274](/home/sano/work/forge/case/60.flatplate_d7275_m7/acceptance.json:274) |
| 壁解像の超過は前縁1点だけ — **Major** | **要再検証** | ツールの `index` は**最大値の位置**であり、超過点の一覧ではない。さらに「面積割合」は実装上、評価点数の割合。実メッシュの壁 DOF は2,149点なので、1点なら約0.047%で、記録の0.3%と整合しない。全超過点の座標・個数・長さ重みを再集計する。[check_wall_resolution.py:369](/home/sano/work/forge/solver_density_cuda/tools/check_wall_resolution.py:369) |
| G15 を条件付き達成として締める — **Major** | **現状では却下** | G15 は位置 II の局所比を一次とし、壁解像も要求する。平均比の判定では代替できない。試験26のユーザ決定は収束ゲートの条件付き扱いであり、試験28の壁解像免除ではない。[plan:1030](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:1030)、[acceptance.json:302](/home/sano/work/forge/case/60.flatplate_d7275_m7/acceptance.json:302) |
| 入口モデルと Tₗ −4%を注記する — **Major** | **採用。ただし制約を明確化** | 計算対象は「衝撃波背後の一様状態を与えた平板」。前縁衝撃波・トリップ・フェンスを含む試験全体の再現ではない。Tₗ差をそのまま熱流束の±4%誤差に換算せず、局所状態の不一致と履歴誤差を未定量と記す。[gen_runs.py:59](/home/sano/work/forge/case/60.flatplate_d7275_m7/gen_runs.py:59)、[plan:904](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:904) |
| Eckert 比で偏差増加を説明する — **Minor** | **定性的整理として採用、厳密な分解は要修正** | 実験/Eckert は図から読んだ R*、forge/Eckert は熱電対位置から換算した R*を使用し、分母が異なる。同じ位置基準にすると実験/Eckert は試験26で **0.74955**、試験28で **0.70361**。傾向は残るが、原因の識別にはならない。[eckert_compare.py:75](/home/sano/work/forge/case/60.flatplate_d7275_m7/tools/eckert_compare.py:75) |

結論: **追加 CFD より先に、`case/60.flatplate_d7275_m7/run_0025_t28_y2/` の壁解像を全点で再集計し、「前縁1点だけ」という完了判断の根拠を監査する。**

第 1 仮説: **壁解像の「最大位置 index 0」を「超過は1点だけ」と読み違えている。** 確度: **高**

  根拠: [acceptance.json:291](/home/sano/work/forge/case/60.flatplate_d7275_m7/acceptance.json:291) の「0.3%＝前縁 index 0 の特異点のみ」に対し、ツールは `argmax` の位置を出すだけ。入力メッシュ `fp_d7275_y2_wall_top8.h5` を独立に確認すると壁 DOF は2,149点あり、全点評価なら0.3%は概ね6〜7点に相当する。

  反証条件: 元の壁出力から、未評価点ゼロ・重複なしで再計算して超過が本当に1点と確認され、0.3%との不整合も説明できること。

第 2・第 3 仮説: **追加しない。** 格子依存と入口近似は未評価の誤差候補だが、現時点で原因の順位を付ける証拠はない。

判別 A/B: **同じ壁出力に対し、集計対象だけを A＝全点、B＝最大値の1点を除いた点集合に変える後処理比較。追加計算は0 step。**

- 最終場と保存済み末尾系列で、超過点数・座標・最大 y₁⁺・未評価点数・壁長さで重み付けした超過割合を出す。
- **Bでも超過が残る →「前縁1点だけ」を棄却。**
- **Bで超過ゼロ → 複数点超過という仮説を棄却し、0.3%という記録の算出根拠を訂正。**
- Bは診断用であり、G15の判定から最大点を除外する許可にはしない。

やらない方がよいこと: **試験26をまず2 µmで回すこと、壁解像 FAIL を黙って免除すること、1.068を物理的な不変性の証明にすること。** 試験26の y₁⁺≈0.27は試験28より小さいが、それだけで格子誤差ゼロとは言えない。現段階では、追加格子より既存判定の監査を優先する。G15達成を根拠とする #60 への復帰は保留する。

呼び出し側の前提への異議: **平均比の比較と G15 の局所比較が混同されている。** 登録済みの熱電対データを位置 II＝1.88 mで線形補間し、記録された CFD 熱流束を用いると、局所比は試験26で約 **1.298**、試験28で約 **1.419**、その比は約 **1.093**となる。これは丸め済みの報告値からの算術確認であり、定常値としての独立検証ではない。位置 II の系列を直接評価し、元の平均比による事前登録判定とは別に報告すべきである。

不足情報: 対象 `run_*` の実効設定、壁・全場出力、残差履歴、判定原本がこの環境にない。台帳上の試験28は **収束 `PASS`／St `ALL STEADY`／壁解像 `FAIL`** だが、今回は再判定できていない。一方、入力メッシュ品質は独立に実行し、**`VERDICT: PASS (AR<=5000, skew<=0.90)`**、最大 AR＝1280.8を確認した。恒久台帳は [case README の計算 run 一覧](/home/sano/work/forge/case/60.flatplate_d7275_m7/README.md:22)。ファイル変更・forge起動なし。**plan 未反映。**
