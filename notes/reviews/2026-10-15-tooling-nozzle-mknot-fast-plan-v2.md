# codex レビュー: tooling-nozzle-mknot-fast (plan)

- **plan**: [`plans/active/tooling-nozzle-mknot-fast.md`](../../plans/active/tooling-nozzle-mknot-fast.md)
- **stage**: `plan`
- **date**: 2026-10-15
- **commit**: `8283a114` (feature/mknot-fast)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 別 (依頼側 `claude-opus-5-5`、実効 `gpt-6-astra`)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 3.7 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-15-tooling-nozzle-mknot-fast-plan-v2.codex.log` (codex: ok)
- **判定**: **GO-with-changes**, 指摘 C0/M4/m2
- **focus**: ユーザ決定で規則を単純化した (§4.0)。4800 にだけ現れる不成立の検出を失うのはユーザが受け入れた限界なので、その点の再指摘は不要。§4.1 の規則の手順に穴が無いか (候補の定義・谷の比較・確認とずらし・r_t の固定点)、§6 の試験と合格条件で実装の誤りを捕まえられるかを見てほしい

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

判定: GO-with-changes
単純化の方針は妥当ですが、候補の比較と確認・ずらしの手順に未定義部分があります。
特に、ずらした結果が別の評価済みの谷より θ_max で 0.98° 悪くなる反例を再現しました。
以下を実装前に修正してください。4800 でのみ現れる不成立は指摘対象にしていません。

1. **Major — 生産網で不成立になった中心点と、ずらせない場合の処理が未定義。**

   **根拠:** [plan:63](/home/sano/work/forge-mknot/plans/active/tooling-nozzle-mknot-fast.md:63) は、生産網の θ_max で谷を比較しますが、生産網での成立による絞り込みを明記していません。現行の [mknot_select.py:192](/home/sano/work/forge-mknot/design/forge_design/evaluate/mknot_select.py:192) は、有限の θ_max を記録した後に制約を判定するため、**θ_max が小さくても `ok=False` の行は存在します**。合成評価でも θ_max=14°、壁余裕=0.5°、`ok=False` を確認しました。

   また、§4.1 の4には「中心だけ不成立」「両隣とも不成立」「移動先が候補外」の分岐がありません。「不成立の側から離れる」だけでは方向が決まりません。V1(iii) の片側不成立だけでは捕まえられません。

   **対案:** 生産網で中心が不成立の谷は比較から除外し、全滅なら表付き例外。ずらしは「中心が成立し、片側だけ不成立で、移動先が同じ候補区間にある場合」に限定する。それ以外はその谷を棄却する。初回確認＋最大3移動と回数を定義し、中心のみ不成立・両側不成立・候補端・移動先の再評価を V1 に追加してください。

2. **Major — ずらした後に谷を比較しないため、既に分かっている良い谷を捨てる。**

   **根拠:** [plan:65](/home/sano/work/forge-mknot/plans/active/tooling-nozzle-mknot-fast.md:65) は「3に戻らずに4をやり直す」と規定しています。次の合成評価をメモリ上で実行しました。

   - 粗い網: `min(15 + 2(M−2)², 15.02 + 400(M−3.45)²)`
   - 生産網: `min(15 + 400(M−2)², 15.02 + 400(M−3.45)²)`
   - 生産網で M=1.95 だけ不成立。他の必要点は成立。

   谷の比較では 2.00／15.00° が勝ちます。しかし確認で左隣が不成立なので 2.05 へ移動し、3点成立で **16.00°を返します**。既に評価した 3.45 は15.02°であり、出力は **0.98°悪化**します。現行 V1(i) と片側不成立の試験を別々に通しても、この組み合わせは検出できません。

   **対案:** **谷ごとに確認・上限付きのずらしを完了し、成立した最終代表点の θ_max を比較する**順序に変更してください。この反例では3.45を返すことを V1 の合格条件にします。全域の再掃引や旧ゲートの復活は不要です。

3. **Major — 候補集合の分断と同点処理を、一意に実装できる定義が足りない。**

   **根拠:** [plan:62](/home/sano/work/forge-mknot/plans/active/tooling-nozzle-mknot-fast.md:62) の候補集合は、不成立点によって複数区間に分かれます。一方、既存の [coarse_local_minima:344](/home/sano/work/forge-mknot/design/forge_design/evaluate/mknot_select.py:344) は配列の隣接要素を比較します。

   実際に、候補だけを詰めた配列  
   `M=[2.00, 2.05, 3.00, 3.05]`、`θ=[15, 16, 17, 18]`  
   を渡すと2.00だけが返り、**別区間の端の谷3.00を落としました**。

   また「谷どうしの差≤0.01°」は推移的ではありません。θ が15.000／15.009／15.018°、境目までの距離が0.1／0.2／0.3の場合、逐次比較では評価順によって選択が変わりました。

   **対案:** 格子を整数 `k` で保持し、候補の連続区間ごとに局所極小を抽出する。同点集合は現行 [_pick:401](/home/sano/work/forge-mknot/design/forge_design/evaluate/mknot_select.py:401) と同様に、**全体の最小値との差≤0.01°**で一度だけ作る。境目の距離は粗い成立区間の端からの推定値と明記し、距離も同じなら θ、M_K の順で決着する。分断区間・平坦列・同点3点・入力順を反転した試験を追加してください。

4. **Major — V2が `r_t` の変化に偏り、実効設定の変更による再選択を保護していない。**

   **根拠:** [plan:68](/home/sano/work/forge-mknot/plans/active/tooling-nozzle-mknot-fast.md:68) と [V2:111](/home/sano/work/forge-mknot/plans/active/tooling-nozzle-mknot-fast.md:111) は、最初／最終の `r_t` と M_K の変更を扱っています。しかし現行コードには別の再試行原因があります。

   - [wt_card.py:1706](/home/sano/work/forge-mknot/design/forge_design/intake/wt_card.py:1706): `n_axis_inv` 倍化、`L_total_tol` 緩和。
   - [wt_card.py:2024](/home/sano/work/forge-mknot/design/forge_design/intake/wt_card.py:2024): 最終壁の判定後にも `x_end_margin` を変更。
   - [wt_card.py:2012](/home/sano/work/forge-mknot/design/forge_design/intake/wt_card.py:2012): 変更後は最初の選択からやり直す。

   したがって「M_K が同じなら出す」だけでは不十分です。同じ `r_t`・M_K でも評価対象が変わります。

   **対案:** 終了条件を「固定点の残差が許容内で、**出力と同じ実効入力に対する選択結果が使用中の M_K と一致すること**」とする。上記3設定の変更で評価を無効化し、`N_s`・`N_p`・許容差・並列数を再設定することを明記してください。V2には各自動調整、M_K の往復、2回目の再固定点で成功／次の変更で例外、最終3点と YAML の入力対応を追加します。

5. **Minor — V4の基準は正しいが、既存の計時欄を使うと再固定点を計上し損なう。**

   **根拠:** [plan:115](/home/sano/work/forge-mknot/plans/active/tooling-nozzle-mknot-fast.md:115) は再固定点込みの600秒を要求しています。一方、現行の [wt_card.py:1972](/home/sano/work/forge-mknot/design/forge_design/intake/wt_card.py:1972) は `mk["t"]` を加算した**後**で `solve_loop` を実行しており、その時間はこの欄に入りません。

   **対案:** §5.1に計時の改修を明記する。合成時計で再固定点に既知の時間を与え、報告時間にその全量が入る試験を V2 に追加してください。V4は既存の `elapsed_s` を無条件に信用せず、内訳と壁時計時間を照合します。

6. **Minor — 残作業表とレビュー採否に、廃止した方針が現役の指示として残っている。**

   **根拠:** [plan:97](/home/sano/work/forge-mknot/plans/active/tooling-nozzle-mknot-fast.md:97) は依然として `N_p/2N_p` の既存ゲートを判断として掲げ、存在しない§4.4を参照しています。[plan:123](/home/sano/work/forge-mknot/plans/active/tooling-nozzle-mknot-fast.md:123) の採用欄も同様です。

   **対案:** 過去のレビュー本文は保存し、採否欄と残作業表には「ユーザ決定により失効」と現行節への参照を記載する。`methods/design/overview.md` の更新も、実装開始前の作業として分離してください。

**推奨は、`coarse_min_theta_v2` を維持し、「候補区間の抽出 → 谷ごとの確認・最大3移動 → 成立した最終代表点の比較 → 実効入力を固定した終了判定」に修正することです。** 実装前の優先順は、指摘1・2、3、4、5、6です。

目的は未解決です。既存の [case/70 の記録:26](/home/sano/work/forge-mknot/case/70.m5_air_wt_demo/_mknot_card_v4_3d/check_final.txt:26) でも M_knot 部分5339秒、カード全体5404秒を確認しました。accepted planの残作業を引き継ぐ位置づけも適切です。捕捉を2件に制限する案は、正常終了時に許容差内の評価を採る [from_length:650](/home/sano/work/forge-mknot/design/forge_design/evaluate/runner_axismach.py:650) と整合します。

本変更は設計前処理であり、CUDA／float32、node／cell、block-DPLUR、周期境界の数値処理を変更しません。今回はCFDの収束・準定常性を主張しておらず、それらのVERDICTは判定根拠に使用していません。コード読解、既存成果物の再集計、ファイルを書かない合成評価で計画上の反例を確認しました。新方式の実データ性能は、実装後のV3・V4で別途判定する必要があります。

証拠: 足りる

指摘数: Critical 0 / Major 4 / Minor 2
