# codex 諮問 (diagnose): twophase-1b-final

- **brief**: [`notes/reviews/briefs/2026-10-02-twophase-1b-final.md`](../../notes/reviews/briefs/2026-10-02-twophase-1b-final.md)
- **plan**: [`plans/active/condensation-two-phase-transport.md`](../../plans/active/condensation-two-phase-transport.md)
- **date**: 2026-10-02
- **commit**: `25332078` (feature/species-transport)
- **codex**: effort `high`, 4.9 min, rc=0
- **結論**: **再正規化の成分別判定を「各更新の相対補正の末尾窓内最大値」に直し、集計間隔によらず合否が一致することを確認してから、#1b の事前登録を確定する。**
- **extra**: `notes/investigations/2026-10-02-twophase-diffusion-kernel-design.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（HEAD `25332078`。禁止された実データは未読。2 step の数値はブリーフ記載値として扱い、コード確認と独立した算術検算を行った。）

| 論点 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| 局所係数偏差 ≤ 2n_sε₃₂ | **採用** | 今回の比較の事前受入上限として採用する。case/16 では **4.7683716e−7**。ただし「固定点なら必ず達成する」という一般的な保証とはしない。 |
| 成分別相対補正に同じ許容を適用 | **現状の区間値への適用は却下・Major** | 分子は全呼出しの絶対補正の累積、分母はログ出力時点の総量である（[計上](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/speciesTransport_d.cu:2141)、[相対値の出力](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/condensationTransport_d.cu:680)）。**各更新の相対補正を計算し、その末尾窓内の最大値**に許容を適用する。累積値は記録用に残す。 |
| 「Δq = (f−1)q なので相対補正 ≤ max｜f−1｜」 | **厳密な上界としては却下・Major** | 実装は `q⁺ = float(double(q)·f)` であり、格納丸めが入る（[該当コード](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/speciesTransport_d.cu:2139)）。局所係数と**実際の格納値の補正**を別々の受入条件として要求する。 |
| 末尾10 %・全更新 | **採用** | 実更新数 N に対して最後の ceil(0.1N) 更新を対象とし、ログの区切りによらず同じ判定にする。窓境界と最終端の集計を欠かさない。 |
| #1b-pre の不足解消・事前登録の確定 | **部分採用、確定は保留・Major** | 旧作用素監査の分岐・式は確認できる（[旧作用素](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/speciesTransport_d.cu:2237)、[選択](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/speciesTransport_d.cu:2333)）。θの件数・呼出し数の集計もある（[集計](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/condensationTransport_d.cu:590)）。成分別計測は追加されたが、受入式との時間方向の不一致が残る。(4)(5) はメモ§16.4の実装報告として確認した。**追加で必要なのは、再正規化ゲートの時間方向の定義と実装を一致させること1つ。** |

結論: **再正規化の成分別判定を「各更新の相対補正の末尾窓内最大値」に直し、集計間隔によらず合否が一致することを確認してから、#1b の事前登録を確定する。**

第 1 仮説: 現案は1更新の補正と区間累積の補正を混同しており、丸め程度の補正でも窓を長くすると不合格になる。確度: 高  
根拠: [設計メモ§16.5](/home/sano/work/forge-species/notes/investigations/2026-10-02-twophase-diffusion-kernel-design.md:792)の上界には更新回数が入っていない。一方、コードは累積量の差分を現在総量で割っている（`condensationTransport_d.cu:644,680`）。

同じ補正標本を繰り返して集計するCPU検算では、f = 1＋ε₃₂、q = 1、各更新の実補正 = ε₃₂ としたとき、現行形式の相対値は次になった。

| 集計する更新数 | 現行形式の相対補正 | 提案の4ε₃₂との比較 |
|---:|---:|---|
| 1 | 1.1920928e−7 | 合格 |
| 10 | 1.1920928e−6 | 不合格 |
| 400 | 4.7683710e−5 | 不合格 |

これは集計式の検算であり、CFDの収束結果ではない。  
反証条件: 現行出力が実際には更新ごとに正規化した値の最大値であり、同じ補正列を再区分しても判定が変わらないこと。

第 2・第 3 仮説: 追加しない。今回の保留理由は計測定義だけで説明できる。

判別 A/B: **同じ100更新分の再正規化補正列を、A＝1更新ごと、B＝10更新ごとに集計する。変更は集計間隔だけ。** 上記の丸め規模の合成入力で実施でき、forge の長時間計算は不要。

採用する判定は次に固定する。

- κ = 2n_sε₃₂。
- 各更新 n で、全対象セルの max｜f−1｜≤ κ。
- 各成分 q について、C_q,n = Σ_i｜q⁺_i,n−q⁻_i,n｜V_i / Σ_i q⁻_i,n V_i を計算し、**末尾窓内の max_n C_q,n ≤ κ**。q⁻・q⁺は再正規化直前・直後の格納値をdoubleに上げた値とする。
- 分母0なら補正絶対量0を要求する。非有限値は不合格。局所条件と成分別条件は**両方要求する独立条件**であり、一方が他方を数学的に保証するとは書かない。

→ **A・Bで一致すれば集計間隔依存の懸念を当該試験で除外。不一致なら計測修正未完了。** 単に累積補正を更新数で割る平均では、途中の大きな補正を隠せるため代用しない。

やらない方がよいこと: 累積相対補正をそのまま4.8e−7と比較して計算を延長すること、観測された累積値に合わせて許容を広げること。

呼び出し側の前提への異議: 格納丸めを除いた上界を実補正へ直接適用できない。例えば f = 1＋ε₃₂、q = 1.5 の算術検算では、実補正比は **1.5894572e−7 > ｜f−1｜=1.1920929e−7**。また、局所偏差が許容内で非ゼロに停滞すること自体は不合格理由ではない。

不足情報: 修正した再正規化ゲートの集計間隔不変性の確認結果。前回不足していた旧作用素監査を再度追加する必要はない。

**plan 未反映。** 読み取り専用の指示に従い変更していない。呼び出し側の反映先は [plan §5.1 #1b・#1b-pre](/home/sano/work/forge-species/plans/active/condensation-two-phase-transport.md:126) と§6。
