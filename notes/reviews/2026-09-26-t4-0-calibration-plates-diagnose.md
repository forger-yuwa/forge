# codex 諮問 (diagnose): t4-0-calibration-plates

- **brief**: [`notes/reviews/briefs/2026-09-26-t4-0-calibration-plates.md`](../../notes/reviews/briefs/2026-09-26-t4-0-calibration-plates.md)
- **date**: 2026-09-26
- **commit**: `0abbb1e3` (feature/gap-heating-precision)
- **codex**: effort `high`, 5.2 min, rc=0
- **結論**: **Cary Table II の Tw/Tt=0.2 の条件を固定し、冷壁平板 1 run を使った「比較原点だけの A/B」を最初に行ってください。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

**現状のまま設計確定は不可です。** 校正平板を先に調べる方針は採用しますが、G14/G15 の原因帰属は却下します。Critical はありません。以下を Major とします。

| 重大度・対象 | 採否 | 根拠と対案 |
|---|---|---|
| **Major 1：G14/G15 の原因判定** | **却下** | [plan:911](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:911) は「温度依存→モデル側」「施設間の比の差→施設側」としているが、この対応は一意でない。測定・遷移履歴にも壁温依存があり、モデル誤差にも Mach 数・Reynolds 数・組成依存がある。**R は条件付きの予測偏差、B/R はその施設間比較として報告し、原因ラベルを付けない。** |
| **Major 2：仮想原点と自然遷移の扱い** | **要再検証** | [plan:897](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:897) の「前縁から完全乱流」と「実測ピークを仮想原点」が、比較座標としてどう接続するか未定義。Cary 原報本文 p.13 は、原点の変更だけで予測が **7–18 %** 変わると記載している。**x_exp、x_peak、x_CFD の対応を明記する。横軸を Re に書き替えるだけでは発達履歴は変わらない。** |
| **Major 3：Cary の条件・10 % 閾値** | **要再検証** | [plan:891](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:891) の単位 Re 固定は Table II と一致しない。同表は Tw/Tt=0.7 で **3.1×10⁵/cm**、0.5 で **2.6×10⁵/cm**、0.2 で **2.7×10⁵/cm**。原報 p.10 の反復性は熱伝達 **10 %以内**、遷移位置 **2 cm以内**。**列ごとの条件を再現し、共通 Re_v 範囲・比較点・集約重みを固定する。10 % は効果量の基準に限定し、不確かさ区間が閾値をまたぐ場合は判定不能とする。** |
| **Major 4：トリップ無し 2D を上下界と扱うこと** | **却下** | [plan:904](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:904) の起点 2 通りは有効長の感度であり、トリップ後の厚さ・速度／温度分布・遷移履歴を挟む保証はない。「厚化の影響は数 %」にも根拠がない。**2D は中心線の探索用に採用し、上下界とは呼ばない。** 定量比較には原報の中心線分布と局所外縁条件を使い、残る履歴誤差を未定量とする。 |
| **Major 5：G15 の比較対象と温度比** | **要再検証** | [plan:892](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:892)、[plan:906](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:906)。TP-1187 本文は q_FP を**タイル中心位置**の値と説明しており、carpet plot のパネル平均説は未確認。さらに **300/1867=0.161** で、指定壁温は Tw/Tt≈0.2 ではない。**局所値を一次比較、パネル平均を別の診断量とする。B/R(0.2) の直接比較を止め、壁温・局所 Mach 数・物性・St の還元定義を揃える。** |
| **Major 6：前提ゲートの不足** | **要再検証** | [plan:1032](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:1032) の壁解像コマンドは不十分。[check_wall_resolution.py:207](/home/sano/work/forge/solver_density_cuda/tools/check_wall_resolution.py:207) は既定で超過 **2 %** を許し、同ファイル:368 は面積ではなく点数割合を計算する。**`--over-frac 0`、未評価点ゼロ、面積重みの確認を明記する。** 格子 2 水準の差≤3 %は感度確認であり誤差上限ではない。また準定常判定は一般的な流れ量でなく、比較する各点の St・q の系列に適用する。 |
| **Major 7：「forge は理論側」を確立済みとすること** | **却下** | [plan:883](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:883) の根拠 run を再判定したところ、下記のとおり両方未収束。**既存値は予備観測へ降格し、モデル誤差の帰属や G14 の精度保証に使わない。** |

再実行した `check_convergence.py` の結果は、各 run の現在の本段 `residual_history.csv` 単独で次のとおりです。

- `case/56.gap_tp1187/run_0002_fp_t8_long/`、最終 step 149999：**`NOT CONVERGED (stalled/plateau)`**
- `case/48.flat_plate_cooled_m4/run_0005_B_tw300_y3/`、最終 step 47999：**`NOT CONVERGED (stalled/plateau)`**

Cary の追加の注意点は次のとおりです。

- 原報 p.6 の壁温は平均の周りに概ね **±15 K**。等温近似は出発点として使えますが、冷壁端では壁温そのものに対して大きい変動です。局所壁温と上流の熱履歴の不確かさを残してください。
- 熱量的完全気体の採用と、輸送物性の妥当性は別です。[gasProperties_d.cu:58](/home/sano/work/forge/solver_density_cuda/cuda_forge/gasProperties_d.cu:58) の Sutherland 定数は直書きで、`visc` を変更しても変わりません。65 K から境界層内の高温域までの μ・λ と実効 Re を確認する必要があります。
- 「SST の誤差」には乱流熱流束閉鎖も含まれます。[viscousFlux_d.cu:264](/home/sano/work/forge/solver_density_cuda/cuda_forge/viscousFlux_d.cu:264) は λ_eff に cₚμₜ/Prₜ を加えています。Prₜ・層流 Pr・SST 補正を事前固定し、温度水準ごとに調整しないでください。
- 原報は水分 **1–2 ppm**、霜の無いデータを選別しています。低温というだけで霜を真因候補へ戻す根拠はありません。

**Cary → D-7275 → TP-1187 層流の順序は採用**します。ただし目的は予測偏差の測定であり、この順序だけで H-model と H-facility を識別できるとはしません。

結論: **Cary Table II の Tw/Tt=0.2 の条件を固定し、冷壁平板 1 run を使った「比較原点だけの A/B」を最初に行ってください。**

第 1 仮説: **仮想原点の対応による差だけで、G14 の 10 % 判定を越える。**　確度: **高**  
　根拠: [plan:897](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:897) の対応が未定義。Cary Table II の Tw/Tt=0.2 の加熱ピークは x=22.22 cm。計画自身の St∝x⁻⁰·² を使った概算では、前縁起点とピーク起点の差は x=29.85 cm で **31.4 %**、49.20 cm でも **12.8 %**。これは CFD の実測ではなく、A/B を優先するための感度見積りです。  
　反証条件: 事前固定した下流比較点すべてで、原点変更による R の変化が、不確かさ込みで **10 %未満**となる。

第 2 仮説: SST と一定 Prₜ の熱流束閉鎖に冷却依存の偏差がある。確度: **中、Cary 条件では未確認**。  
第 3 仮説: 実験還元・施設条件に偏差がある。確度: **低、異なる施設に共通する系統誤差は未確認**。

判別 A/B: **変更するのは抽出器の比較原点 x₀ だけ**。同一の完全乱流 CFD 場に対し、A は x_CFD=x_exp、B は x_CFD=x_exp−22.22 cm で評価する。実験側の St、物性、比較点は固定する。これは自然遷移の再現ではなく、比較方法の感度試験です。  
　本段は初回 **20,000 step、1,000 step ごと保存**を試行予算とし、同一設定区間の `PASS` と対象 St 系列の `STEADY` が得られなければ判別を保留する。  
　→ 差が不確かさ込みで **10 %超なら「原点感度は小さい」を棄却**。**10 %未満なら第 1 仮説を棄却**。閾値をまたぐ場合は判定不能。どちらでも施設／モデルへの帰属はまだ行わない。

やらない方がよいこと: **G14/G15 の結果だけで SST または実験を原因と断定すること、実測に合うよう原点・Prₜ を選ぶこと、2D 起点感度をトリップ効果の上下界と呼ぶこと。**

呼び出し側の前提への異議: **文献に記載された偏差、未収束 CFD の予備値、carpet plot 内挿という推測が「観測事実」に混在しています。これらを分離する必要があります。**

不足情報: Cary 各列の条件台帳と壁温分布、比較座標の対応、低温輸送物性の検証、D-7275 の実測時壁温・平均化手順、誤差の相関を含む不確かさ予算。

**plan 未反映。ファイル変更・`forge` 起動は行っていません。**
