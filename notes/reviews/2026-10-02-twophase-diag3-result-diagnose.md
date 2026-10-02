# codex 諮問 (diagnose): twophase-diag3-result

- **brief**: [`notes/reviews/briefs/2026-10-02-twophase-diag3-result.md`](../../notes/reviews/briefs/2026-10-02-twophase-diag3-result.md)
- **plan**: [`plans/active/condensation-two-phase-transport.md`](../../plans/active/condensation-two-phase-transport.md)
- **date**: 2026-10-02
- **commit**: `c7ef3c87` (feature/species-transport)
- **codex**: effort `high`, 4.2 min, rc=0
- **結論**: **次は `run_0511` の既存診断CSVを同じ更新番号・セルIDで独立再集計し、残差内訳と最大セルの矛盾を解消する。**
- **extra**: `notes/investigations/2026-10-02-twophase-diffusion-kernel-design.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表:

| 論点 | 判定 | 根拠・対案 |
|---|---|---|
| 乾燥停止が Q 残差を直接支配する | **棄却を採用（提示集計の範囲）** | `case/16.nozzle_wys/run_0511_twophase_relax05_diag3/` の更新1800–1999で、停止集合内の Q 残差/許容 ≤3e−25、集合外790〜1000、最大セルは200/200更新で集合外。ただし「全域未収束と無関係」への拡張は **Major・却下**。蒸気更新停止の間接影響は未検証。 |
| 残る原因は下流の移流・拡散・ソースの拮抗 | **Major・要再検証** | 提示された残差合計・内訳・最大セルが相互に矛盾する。下記の再集計を先に行う。 |
| 乾燥セルは非負だから放置してよい | **Major・却下** | 非負性と原方程式の残差合格は別条件。さらに [plan:126](/home/sano/work/forge-species/plans/active/condensation-two-phase-transport.md:126) は末尾10%で θ<1 件数0を要求する。実装修正は保留してよいが、未解決事項と受入未達の記録は残す。 |
| #1bを区切り、opt-inのまま凍結する | **採用** | **影響評価は判別不能・受入未達**として区切る。追加の長時間計算は提案せず、既存診断出力の整合確認だけを次の一手とする。完了・問題なしとは扱わない。 |

結論: **次は `run_0511` の既存診断CSVを同じ更新番号・セルIDで独立再集計し、残差内訳と最大セルの矛盾を解消する。**

第 1 仮説: **ブリーフ作成時の抽出・集計・転記で、異なる更新・セル・成分の値が混ざっている。** 確度: **中**

根拠: **［Major］** [依頼ファイル:63](/home/sano/work/forge-species/notes/reviews/2026-10-02-twophase-diag3-result-diagnose.prompt.md:63) の数値を検算すると、以下になる。単位はすべて ×10⁻⁹。

| セル | 記載された R_g,B | 記載された移流＋拡散＋ソース | 差 |
|---|---:|---:|---:|
| 33877 (`max_g`) | 4.49 | 2.17 | 2.32 |
| 41578 (`max_Q0`) | −2.40 | 2.45 | −4.85 |
| 29485 (`max_v`) | 8.58 | 11.84 | −3.26 |

表示丸めでは説明できず、41578では符号も逆になる。また、同じ最終更新・同じ対象集合なら、`max_v` 行の |R_g,B|＝8.58e−9 が `max_g` の4.49e−9を上回ることはない。[設計メモ:880](/home/sano/work/forge-species/notes/investigations/2026-10-02-twophase-diffusion-kernel-design.md:880) は `max_g` を B の |残差| 最大セルと定義している。

これは**報告の不整合が存在する証拠**であり、まだソルバ本体の不具合の証拠ではない。

反証条件: 元の `twophase_diag3_cells.csv` の同一行でも上記の不一致が再現すること。その場合、「要約段階だけの誤り」は棄却する。

第 2・第 3 仮説:
- **第2：診断CSV自体で、総残差・項別残差・最大セル選択の評価時点または定義が揃っていない。確度：低、未確認。** 第1仮説と以下の再集計で切り分ける。
- 第3は置かない。二相経路固有か共通の問題かを順位付けする証拠が不足している。

判別 A/B: **変更するのは抽出・集計方法だけ。追加計算0 step、対象は既存の1800–1999の全200更新。**

- **A＝現在の集計結果、B＝列名と `(update, cell, group)` を明示して行う独立再集計。** 列ごとの窓内最大値を足し合わせず、必ず同じ行から値を取る。
- Bでは全記録について、`Rg_B_double − (Rg_B_advection + Rg_B_diffusion + Rg_source)` を検査する。列整合の許容はCSVの表示丸め誤差＋既存の `Rg_local_error_bound` とし、超過件数0を要求する。これは**列の整合確認**であり、本番float演算の誤差保証には使わない。
- 各更新の `max_g` が、同じ更新の他の記録行の |R_g,B| 以上であることを確認する。Qについても最大セル・停止集合所属を `qcompare` と照合し、不整合件数0を要求する。

**Bで整合し、現在の要約だけが不一致なら第1仮説を支持し、第2仮説を今回の矛盾の原因から除外する。Bでも不一致なら第1仮説を棄却し、診断出力側を調べる。** 必要列や対象集合の定義が欠ける場合は確認未達とし、物理原因の判定に進まない。

やらない方がよいこと:

- 緩和0.5を8000 stepへ延長し、Q1だけの低下で方針を決めること。全Q・流れ・再正規化・制限解除の受入条件が残る。
- 旧作用素と新作用素の残差最大位置が近いだけで、共通原因と断定すること。[設計メモ:612](/home/sano/work/forge-species/notes/investigations/2026-10-02-twophase-diffusion-kernel-design.md:612) にあるとおり、作用素に加えて前処理・更新経路も変わる。
- θの正の下限、負の液増分の無条件破棄、受入閾値の緩和で停止を隠すこと。

呼び出し側の前提への異議:

**［Major］「float丸めではない」は広すぎる。** double組立でも負という記録が正しければ、除外できるのは「同じ入力からのfloat組立だけが負を作った」という説明である。診断Bも本番の面値・係数・ソース値を共有するため、その生成段階の丸めや再構成は除外されない。[設計メモ:875](/home/sano/work/forge-species/notes/investigations/2026-10-02-twophase-diffusion-kernel-design.md:875)

また、`max_Q0` セルに載せた **R_gの内訳はR_Q0の内訳ではない**。その値からQ0未収束の支配項は決められない。対案は、まず上記の列整合を確認し、成分を区別して解釈すること。

不足情報: 元診断CSVの該当行、200更新全体の整合検査結果、最大セル選択の対象集合。禁止されたrunファイルとソースコードは読まず、実測の判定は提示記録に依拠した。同runは記録上 **`renorm-gate: FAIL`、独立監査 `NOT CONVERGED`** であり、受入未達という判断は維持する。

ファイル変更・forge起動なし。**plan未反映**。呼び出し側の反映先は [plan §5.1 #1b-r2](/home/sano/work/forge-species/plans/active/condensation-two-phase-transport.md:129) と #1b の処置欄。
