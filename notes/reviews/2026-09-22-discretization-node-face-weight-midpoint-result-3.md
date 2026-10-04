# codex レビュー: discretization-node-face-weight-midpoint (result)

- **plan**: [`plans/active/discretization-node-face-weight-midpoint.md`](../../plans/active/discretization-node-face-weight-midpoint.md)
- **stage**: `result` (diff base `18916ed9`)
- **date**: 2026-09-22
- **commit**: `4a070144` (feature/sern-design)
- **codex**: effort `high`, 4.0 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M1/m3
- **focus**: 3 回目の result レビュー。2 回目 (notes/reviews/2026-09-22-discretization-node-face-weight-midpoint-result-2.md、NO-GO C0/M2/m3) への対応を見てほしい: M1 → 製造解の壁面量の指標を 3 つに分けて測り直し (壁 CV に入る熱伝導・粘性仕事を同じ面の厳密流束と比較 / 2 点差分を解析的な壁面値と比較 / 2 点差分を厳密温度の 2 点差分と比較)、「壁面熱流束も 2 次」を撤回。M2 → 運動量 (2 成分、可変粘性、転置項・発散項つき) と粘性仕事を追加。面の式は node_visc_face.visc_face() で、実機 C3X の場に当てるとカーネルが積んだ wi_eheat / wi_ework を 480 壁節点で 0.28 / 0.00 W/m2 の差で再現する (旧式の重みのとき。射影・0.5 では仕事項が 2.0 kW/m2 ずれる)。結果: fx=0.5 は節点値 p=2.00 だが壁 CV の流束は 1 次 (粘性仕事 31 % → 2.9 %)、旧式は 25–33 % のまま収束しない。問い: (1) この測定と主張の範囲 (§4.2 の『主張する範囲』、§5.1 #7) に過大・誤りは無いか。特に『壁 CV の流束が 1 次なのは双対面の面積重心が弦のたるみの分だけ辺中点より壁側にあるため』という説明と、q_eff への +0.8 % という見積りの妥当性。(2) これで accepted に移してよいか、まだ active に残すべきか、推奨を 1 つに。読んでよいのは本 plan、同梱の 2 ファイル、solver_density_cuda/tools/mms_face_weight.py (全文可、約 230 行)、notes/investigations/2026-09-22-mms-node-face-weight.out.txt。禁止: *.log, residual_history.csv, res_*.h5, *.vtu, plans/README.md, case/*/README.md の全文 (grep は可)。
- **extra**: `notes/investigations/2026-09-22-mms-node-face-weight.md`, `solver_density_cuda/tools/node_visc_face.py`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
今回の製造解は、`fx=0.5` を採用する根拠として有効です。ただし「実機の `q_eff` を +0.8 % 押し上げる」という結論は裏付け不足です。  
以下の記述修正後に `accepted` へ移すことを推奨します。追加の大規模計算を移行条件にはしません。

後段の閲覧制限を優先しました。許可されたコードを読み、成長率 1.1 の全 4 水準・3 重みを再実行して保存表を再現しました。実装 diff、実 run の VERDICT、他の仕様文書は今回再監査しておらず、480 壁節点でのカーネル照合値も独立確認済みとは扱いません。ファイル変更はありません。

1. **Major — `q_eff` の「+0.8 %」は、誤差の定義と分母が揃っていない。**

   根拠: [調査メモ:53](/home/sano/work/forge/notes/investigations/2026-09-22-mms-node-face-weight.md:53) は「粘性仕事の 30 % の過大評価 × 寄与率 2.8 %」で換算しています。しかし [測定コード:179](/home/sano/work/forge/solver_density_cuda/tools/mms_face_weight.py:179) の `wallW` は、**符号付き平均ではなく `rms(wh-we)/mean(abs(we))`** です。

   再実行した最粗格子では、厳密仕事の平均が 767.95 W/m²、`half` が 980.58 W/m²。相対 RMS は **31.21 %**ですが、平均の過大評価は **27.69 %**です。また、2.8 % が数値計算の仕事項の寄与率なら、仮に符号付き誤差を 30 % としても、同じ分母への換算は `2.8% × 0.30/1.30 ≈ 0.65%` になります。これも実機の補正値を意味しません。

   さらに、この製造解は熱と運動量を別々に解き、粘性仕事を後評価しています（[コード:156](/home/sano/work/forge/solver_density_cuda/tools/mms_face_weight.py:156)）。実機のエネルギー解が変化した際の熱伝導との相殺・増幅は測っていません。

   **対案:** 「+0.8 %」を実機の誤差見積りとしては撤回し、「仕事項への相対誤差が数十 % なら、数 % の寄与率を介して熱流束に影響し得る」という条件付きの規模評価に限定してください。実機の符号付き誤差を残すなら、同じ領域・重み付け・分母で別途評価が必要です。

2. **Minor — 「粘性仕事が 1 次」「旧式は収束しない」は、相対誤差に限定する必要がある。**

   根拠: [plan:81](/home/sano/work/forge/plans/active/discretization-node-face-weight-midpoint.md:81) と [調査メモ:45](/home/sano/work/forge/notes/investigations/2026-09-22-mms-node-face-weight.md:45)。再実行で得た、壁面積当たりの仕事の絶対 RMS 誤差は次のとおりです。

   | 重み | 水準 0 → 1 → 2 → 3 [W/m²] |
   |---|---|
   | `half` | 239.68 → 60.03 → 15.11 → 3.79 |
   | `code` | 253.05 → 136.94 → 67.44 → 32.56 |

   **`half` の絶対誤差は約 2 次、旧式も約 1 次で減っています。** 厳密仕事自体が第一層厚とともに小さくなるため、相対誤差ではそれぞれ約 1 次、停滞になります。旧式の熱伝導誤差が停滞する証拠は別にありますが、仕事の相対誤差だけから「仕事が絶対量として収束しない」とは言えません。

   **対案:** §4.2・#7 を「壁 CV の粘性仕事の**相対 RMS 誤差**は `half` で約 1 次、旧式では測定範囲で約 25 % に停滞」と統一してください。節点値の約 2 次も、解析勾配と解析温度由来の物性を使う補助問題での結果と明示すると正確です。

3. **Minor — 幾何による説明は主要なスケーリングとして妥当だが、面重心の沈みを面積分の誤差と同一視している。**

   根拠: [調査メモ:49](/home/sano/work/forge/notes/investigations/2026-09-22-mms-node-face-weight.md:49)。コードの厳密値は、面重心での評価ではなく、折れた双対面の各線分上で積分しています（[コード:70](/home/sano/work/forge/solver_density_cuda/tools/mms_face_weight.py:70)、[コード:117](/home/sano/work/forge/solver_density_cuda/tools/mms_face_weight.py:117)）。

   同じ構成の等間隔円環で確認すると、最粗条件の面重心の沈みは **0.29666 µm**ですが、面上の半径方向距離を積分した平均の沈みは **0.19777 µm**です。近壁速度を壁距離に比例すると近似すれば、後者は前者の約 **2/3**になります。

   したがって、`たるみ / 第一層厚 = O(h)` という相対誤差の説明は支持できますが、その係数は面重心だけでは決まりません。実際、製造解の最粗格子で `proj` の仕事の平均誤差は **−14.87 %**。`half` の +27.69 %から、過小評価側へ移っています。

   **対案:** 原因説明を「曲面上の場と、辺中点評価・双対面積分の位置の違い」に修正し、「射影はこの試験で誤差の絶対値を約半減するが、符号は反転する」としてください。この説明を熱伝導の観測次数にもそのまま適用するのは避けるべきです。

4. **Minor — 残作業表が現在の合格条件・レビュー段階と一致していない。**

   根拠: [plan:116](/home/sano/work/forge/plans/active/discretization-node-face-weight-midpoint.md:116) の #2 は、現在も「§6 が `PASS` を要求」「#2c として残す」と記載しています。一方、#2c は完了し、§6 は物理量による回帰判定へ変更済みです。#8 も第1回レビューの状態のままです（[plan:124](/home/sano/work/forge/plans/active/discretization-node-face-weight-midpoint.md:124)）。

   **対案:** #2 の旧条件を履歴として区別し、#8 を今回の結果に更新してください。#7 の測定完了と、本レビューによる主張修正の未完了も分けて記録してください。

**推奨は、`fx=0.5` を維持し、上記を修正してから `accepted` へ移すことです。** 優先順は、①実機の +0.8 %という断定の撤回、②相対・絶対誤差の区別、③幾何説明の修正、④残作業表の同期です。前回問題だった壁面指標の混同と粘性仕事の未測定は、今回のコードで解消しています。本レビューは read-only のため **plan 未反映**です。

指摘数: Critical 0 / Major 1 / Minor 3
