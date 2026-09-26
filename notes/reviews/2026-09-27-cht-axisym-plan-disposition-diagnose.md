# codex 諮問 (diagnose): cht-axisym-plan-disposition

- **brief**: [`notes/reviews/briefs/2026-09-27-cht-axisym-plan-disposition.md`](../../notes/reviews/briefs/2026-09-27-cht-axisym-plan-disposition.md)
- **plan**: [`plans/active/boundary-cht-axisymmetric-fem2d.md`](../../plans/active/boundary-cht-axisymmetric-fem2d.md)
- **date**: 2026-09-27
- **commit**: `392effba` (feature/cht-axisym-fem2d)
- **codex**: effort `high`, 1.8 min, rc=0
- **結論**: **先に plan の負例式と拒否解除条件を訂正し、その反映後に `status: in_progress` として §5.1 #3 の固体 `r` 重み実装へ進む。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（全件採用の方向は妥当。ただし「現 plan に正しく反映済み」とは認めない）：

| 指摘 | 判定 | 根拠・必要な修正 |
|---|---|---|
| M1 / Major | **採用** | 一様 IC・同一設定の単一区間で `PASS` を要求するのは正当な検証条件。[plan:178](/home/sano/work/forge-cht/plans/active/boundary-cht-axisymmetric-fem2d.md:178)。ただし `run_0002` の README は残差 `PASS` と壁温「まだ緩和中」を併記しているため、連成全体の収束実績には使えない。[README:43](/home/sano/work/forge-cht/case/52.conjugate_slab/README.md:43)、[README:56](/home/sano/work/forge-cht/case/52.conjugate_slab/README.md:56)。解析誤差・全節点の準定常・G-if との**同時合格**を維持する。 |
| M2 / Major | **採用・反映修正必須** | [plan:188](/home/sano/work/forge-cht/plans/active/boundary-cht-axisymmetric-fem2d.md:188) の負例の面積比が逆。正しい旧実装の模擬は **`q_bad = q_eff × A_fluid^r / A_planar`**。登録円環では正値 91.85535 に対し旧値 0.9185535、現記述の模擬値は 9185.535。現記述で `FAIL` しても、旧実装への検出力を試したことにならない。 |
| M3 / Major | **採用** | 幾何方式・単位・契約版の照合と、平面旧状態の互換性維持は妥当。[plan:102](/home/sano/work/forge-cht/plans/active/boundary-cht-axisymmetric-fem2d.md:102)。ただし公開時の拒否解除条件にも V-ax4 を必須として追記する。 |
| M4 / Major | **採用** | 評価器の重み修正・状態照合と、未対応の外部経路の拒否は必要。[plan:106](/home/sano/work/forge-cht/plans/active/boundary-cht-axisymmetric-fem2d.md:106)。現解除条件は V-ax0〜3 しか列挙していないため、V-ax4 を含める。[plan:141](/home/sano/work/forge-cht/plans/active/boundary-cht-axisymmetric-fem2d.md:141)。 |
| M5 / Major | **採用** | 予測列を再計算すると誤差比は **0.31663・0.29487**、いずれも 0.4 以下。要求次数は約 1.322、予測次数は 1.659・1.762 なので整合する。ただし元の数値は**固体への荷重移送試験**であり、連成計算や非一様格子の合格予測ではない。[レビュー:55](/home/sano/work/forge-cht/notes/reviews/2026-09-27-boundary-cht-axisymmetric-fem2d-plan.md:55)。非一様格子は間隔が増える半径方向も固定する。 |
| m6 / Minor | **採用** | FP32 の受理と精度保証を分ける方針は妥当。ただし未実装の現段階では「FP64 で検証済み」を「FP64 を検証対象とする」に直し、合格後に更新する。[plan:84](/home/sano/work/forge-cht/plans/active/boundary-cht-axisymmetric-fem2d.md:84)。 |

結論: **先に plan の負例式と拒否解除条件を訂正し、その反映後に `status: in_progress` として §5.1 #3 の固体 `r` 重み実装へ進む。**

第 1 仮説: 数学的な中核を変更する必要はなく、現在の障害は採用指摘の検証仕様への反映漏れである。確度: 中  
  根拠: [plan:62](/home/sano/work/forge-cht/plans/active/boundary-cht-axisymmetric-fem2d.md:62) の重み付き剛性・Robin 式と、[plan:188](/home/sano/work/forge-cht/plans/active/boundary-cht-axisymmetric-fem2d.md:188) の逆向きの負例式。後者は代数だけで誤りを確認できる。  
  反証条件: 独立積分を通過した実装でも、十分に定常化した V-ax2b が登録した局所誤差・細分化条件を満たさなければ、「検証仕様を補えば足りる」という仮説を棄却する。

第 2・第 3 仮説: 無し。

判別 A/B: **熱流束診断の分母だけ**を変える。固定した合成荷重 `Q_i=q_* A_i^r` に対し、A は `A_i^r`、B は `A_planar,i` で除算し、同じ評価器を各1回実行する。forge の時間積分は不要。→ **A が PASS・B が FAIL なら**旧分母への検出力あり／**それ以外なら**評価器または面積定義に不備があり、検証仕様が完成したという仮説を棄却する。

やらない方がよいこと: 現記述の逆比の負例が落ちただけで M2 を完了扱いすること、`run_0002` の残差 `PASS` だけで連成定常を主張すること、V-ax4 未合格で拒否を解除すること。

呼び出し側の前提への異議: 「一様 IC なら PASS する」は未確認。正当なのは**一様 IC から PASS を要求する設計**であり、軸対称での達成実績ではない。また、継続 run の判定不能を一般化しない。同一設定区間の連結や `--from-floor` は [AGENTS.md:83](/home/sano/work/forge-cht/AGENTS.md:83) に規定されている。

不足情報: 指定された閲覧範囲には `run_0002` の実効設定・残差原本・`CONVERGENCE_VERDICT.txt` が含まれないため、一様 IC・判定区間・現行ツールでの PASS は独立確認していない。今回確認したのは文書間の整合性と数式・数値比である。ファイル変更なし。提案は **plan 未反映**。
