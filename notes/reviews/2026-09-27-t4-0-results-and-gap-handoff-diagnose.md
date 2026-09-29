# codex 諮問 (diagnose): t4-0-results-and-gap-handoff

- **brief**: [`notes/reviews/briefs/2026-09-27-t4-0-results-and-gap-handoff.md`](../../notes/reviews/briefs/2026-09-27-t4-0-results-and-gap-handoff.md)
- **date**: 2026-09-27
- **commit**: `a8d1bf03` (feature/gap-heating-precision)
- **codex**: effort `high`, 4.1 min, rc=0
- **結論**: **T4-0 の完了と #60 の再開は保留し、D-7275 試験26で「下流 slip 区間を壁距離計算に含める／含めない」だけの A/B を先に行う。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表:

| 対象 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| D-7275 を「収束した校正結果」として扱う | **却下・Major** | [case/60/README.md:27](/home/sano/work/forge/case/60.flatplate_d7275_m7/README.md:27) の記録は **`NOT CONVERGED`／St `ALL STEADY`**。報告できるのは「未収束の探索計算で R=1.266–1.442、平均1.313」。G15 の収束条件を満たした結果には昇格させない。 |
| 比較域の変化が小さいので残差ゲートを免除する | **却下・Major** | [README.md:28](/home/sano/work/forge/case/60.flatplate_d7275_m7/README.md:28) が測ったのは **2000 step 間の状態差**であり、空間残差の分布ではない。陰解法の更新は残差と局所刻み・ヤコビアンを介するため、両者は同義でない（[implementation.md:375](/home/sano/work/forge/methods/time_integration/implementation.md:375)）。対案は下記の一因子 A/B。 |
| Cary と D-7275 の偏差を並べる | **条件付き採用。共通原因への帰属は却下・Major** | 条件・比較座標・判定状態を別々に示す参考表ならよい。「共通の SST 冷壁誤差約30%」とは言えない。[plan:904](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:904) のとおりトリップ履歴は未定量で、Tw/Tt、物性、St の還元も異なる。St 化だけで壁温・熱履歴の影響が消えたとは扱わない。 |
| T4-0 を完了して #60 に戻る | **却下・Major** | [G15:1031](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:1031) は試験26・28と G14 同等の前提ゲートを要求する。試験26のゲート未達に加え、試験28も未提示。なお登録済みの次条件は **試験28・α=9.8°**であり、7.5°系列ではない（[plan:892](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:892)）。まず試験26の成立確認を優先する。 |
| 平板の過大だけでは D1 不足を説明できない | **採用。ただし根拠を修正・Major** | 現行 D1 は既に**実測分母**で判定している（[skin_forward3d.py:241](/home/sano/work/forge/case/56.gap_tp1187/tools/skin_forward3d.py:241)）。forge 分母の高さが、その FAIL を作ったわけではない。復帰時は有次元の観測演算後熱流束、実測分母による比、同条件 forge 平板分母による比を併記する。D-7275 の1.313を補正倍率にはしない。 |
| x 対応 ±3% の感度が確認済み | **要再検証・Minor** | 最下流点は計算すると **2.54596 m、+3%で2.62233 m**。平板出力は2.60 mまでだが、[d7275_compare.py:55](/home/sano/work/forge/case/60.flatplate_d7275_m7/tools/d7275_compare.py:55) は範囲外を端点値で返す。対案は範囲検査を入れ、この点を判定不能とすること。公称位置の R 全体が誤りという指摘ではない。 |

結論: **T4-0 の完了と #60 の再開は保留し、D-7275 試験26で「下流 slip 区間を壁距離計算に含める／含めない」だけの A/B を先に行う。**

第 1 仮説: **下流の人工的な壁終端に伴う `wall_dist` の変化が、試験26の残差プラトーを支配している。**　確度: **中**

  根拠: ローカルの `case/60.flatplate_d7275_m7/mesh/fp_d7275_y3.h5` を読み取ると、下流底面149節点の `wall_dist` は **0.587 mm〜200 mm**。生成元の設定には `wallDistExtraPhysIDs` がない（[make_case.py:29](/home/sano/work/forge/case/56.gap_tp1187/tools/make_case.py:29)）。この接続で SST が不安定になる既往と、壁距離を延長する処置は [calcWallDistance_kdtree.cpp:126](/home/sano/work/forge/solver_density_cuda/input/calcWallDistance_kdtree.cpp:126) に明記されている。`run_0003_t26_probe` の下流で大きい状態変化という台帳記録とも整合する。ただし、**今回の残差の発生場所は未確認**。

  反証条件: 延長した `wall_dist` が実際に読み込まれたことを確認しても、十分な末尾窓で未達残差の水準・傾向が変わらなければ、「壁距離の変化が主因」を退ける。

第 2 仮説: **比較熱流束も下流端処理の影響を受けており、短い状態差では検出できていない。**　確度: **低・未確認**。下記 A/B の比較点別熱流束で判別する。第3仮説は置かない。

判別 A/B:

- **変更因子は下流 slip の壁距離への算入だけ**。A＝現行、B＝変換時に下流 slip を `mesh.wallDistExtraPhysIDs` に指定する。現状の `physID: 5` は**前縁上流と下流を兼用**するため、両腕共通で下流だけ別 ID に分ける必要がある（[mesh:36](/home/sano/work/forge/case/60.flatplate_d7275_m7/mesh/fp_d7275_y3.geo:36)）。`[5]` をそのまま指定して上流まで変えない。
- 同一の FP64 バイナリ・節点配置・保存量初期場・BC・CFLを使い、新規2 runで **まず20,000 step、1,000 stepごと保存**。保存量は `restart_field.py --keep-src-dtype` でコピーし、壁距離は各入力側に残す。
- 見る量は**全残差列、固定10点の St、位置IIとパネル平均の qwall**。比較前に上記の範囲検査を直す。
- **Bで未達列の末尾残差がAの1/10以下となり、熱流束差が全比較点で1%以内なら**、第1仮説を支持し、第2仮説の「1%を超える影響」を退ける。**残差水準が変わらなければ**壁距離単独主因説を退ける。**熱流束差が1%を超えれば**「比較域には効かない」を退ける。1%・1桁は今回提案する事前判別基準であり、既存ゲートの代替ではない。過渡が残れば判別保留。
- 各腕の `check_convergence` と `check_quasisteady` の VERDICT を残す。Bは作用素が変わるので、変更前の履歴と連結して収束判定しない。

やらない方がよいこと: **Major — ゲート閾値を緩めて完了扱いにすること、Prₜ を合わせ込むこと、Cary/D-7275 の倍率で D1 を補正すること。** 対案は上記 A/B。層流 T4-0c、迎角系列、H群格子追加は、この成立確認より先に広げない。

呼び出し側の前提への異議: **「x>2.55 m＝slip バッファ」は誤り。** 実メッシュの壁終端は **2.60 m**（[mesh:3](/home/sano/work/forge/case/60.flatplate_d7275_m7/mesh/fp_d7275_y3.geo:3)）で、診断領域には等温壁末尾も含まれる。領域を分け直し、「状態変化の局在」と「残差の局在」を区別する必要がある。

不足情報: AWS の実効設定、バイナリ識別情報、残差・壁時系列、判定区間と VERDICT 原本がローカルにない。**本回答の run 数値と VERDICT は台帳の引用で、独立再判定ではない。** 独立確認したのはコード・登録条件・入力メッシュ・比較座標の算術。ファイル変更・forge 起動なし。**plan 未反映。**
