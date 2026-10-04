# codex 諮問 (diagnose): twophase-g3-design

- **brief**: [`notes/reviews/briefs/2026-10-04-twophase-g3-design.md`](../../notes/reviews/briefs/2026-10-04-twophase-g3-design.md)
- **plan**: [`plans/active/condensation-two-phase-default.md`](../../plans/active/condensation-two-phase-default.md)
- **date**: 2026-10-04
- **commit**: `4f92f201` (feature/species-transport)
- **codex**: effort `high`, 3.8 min, rc=0
- **結論**: **格納丸め `C_round` を明示した本番一更新の診断を作り、同一入力の float 格納／double 参照 A/B で G3-b の尺度を確定する。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論は、**G3 の二分は採用、提案された丸め誤差ゲートは却下**です。G2 が示すのは「調べた二つの ON レシピの差が、基準 ON−OFF 差に比べて小さい」までです。

採否表:

| 重大度 | 対象・採否 | 根拠と対案 |
|---|---|---|
| **Major** | G2 の「同じ値」「レシピに依存しない」：**却下** | 保存 CSV を再集計すると、生産−cfl2 は `g_exit_mw = −1.5077e−8 ± 4.7571e−10`、壁温 `+3.5950e−4 ± 3.5926e−5 K`。± は腕平均差の標準誤差で、差自体は反復ノイズを超える。対案は「レシピ感度／ON−OFF 差は 7 量で **1.66e−5〜2.58e−2**、事前登録の 0.1 以下」。根拠：`notes/investigations/2026-10-04-twophase-g2/evidence/run_0561_g2_onprod_r1/`〜`run_0569_g2_off_r3/` の `twophase_series.csv`、24000〜48000 step、各 61 点。 |
| **Major** | G2 の完全な PASS 証明：**要再検証** | [g2_judge.py:35](/home/sano/work/forge-species/notes/investigations/2026-10-04-twophase-g2/g2_judge.py:35) は `commit/passive_floor/clamp` の最後の行だけを拾い、段欠落を必須エラーにせず、`projection`・段合計・`nonfinite`・正式な FINAL VERDICT を確認しない。本体の正式ゲートはそれらを判定する：[condensationTransport_d.cu:1013](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/condensationTransport_d.cu:1013)。対案：数値比較の PASS は維持し、正式ゲートの構造化した記録を補って判定器を修正する。現在の run が補正ゲートに失敗したとまでは言えない。 |
| **Major** | G3-a の `64ε₃₂·Σabs(各項)`：**却下** | 境界流束やソースの**領域積分後**の絶対値では、内部面で相殺される大きな加算を数えられない。本番は内部面の両端へ別々に float atomic 加算する：[passiveKernels_d.cuh:32](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/passiveKernels_d.cuh:32)。対案：各節点に入る全加算項の絶対値と加算回数から上界を作る。固定の 64 には現時点で根拠がない。 |
| **Major** | G3-b の `16ε₃₂·Σabs(VΔq)`：**却下** | 本番には `fl(q_N + δq)` の格納丸めがある：[speciesTransport_d.cu:485](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/speciesTransport_d.cu:485)。誤差は微小な δq ではなく基底 q の大きさに依存する。対案：格納丸めを独立項として記録する。下記の反例で提案ゲートは正常な加算を FAIL にする。 |
| **Major** | κ を「既存と同じ基準」として全補正・全 Ω に適用：**却下、分類を修正** | 既存ゲートは全域の成分別総量を分母とし、再正規化は別ゲート、液滴消滅は対象外：[condensationTransport_d.cu:958](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/condensationTransport_d.cu:958)。局所 Ω・再正規化込み・境界上書き込みは別基準になる。対案：既存ゲートはそのまま再現し、局所の数値補正ゲートを追加する。境界からの物理的供給は数値補正と分ける。 |
| — | 作用素／更新写像の分離、wrapper 内の記録点：**採用、以下の補足が必要** | [main.cpp:1905](/home/sano/work/forge-species/solver_density_cuda/main.cpp:1905) の組立順と整合する。前処理・格納丸め・境界拘束の扱いを明文化する。 |
| — | 実装方式：**本番カーネル内の診断出力を採用** | G3 の目的は、本番が実際に加算・格納した値の監査。面・節点専用スロットへ値を書き、積算は後で double で行う。診断専用の写しを収支の正本にする案は採らない。 |
| — | `run_0561`・`run_0567` の最終場：**診断入力として採用** | 各設定での作用素・更新の監査には適切。ただし、新しく測る出口液・Q 流束の定常性まで、既存の 7 量の VERDICT で保証してはいけない。 |

G2 の保存記録は、9 本とも **`CONVERGENCE: NOT CONVERGED`、`QUASISTEADY: ALL STEADY`** です。残差評価は保存記録上 47999 step まで、報告量の比較窓は 24000〜48000 step。したがって「7 報告量の準定常比較」として扱い、残差収束や場全体の一致とは表現しません。

G3 の仕様には、次を入れてください。

- **制御体積と境界**：Ω は実節点 CV の集合とし、ghost は含めない。面ごとの所属指示値から境界面を決め、領域内面は相殺する。node 境界半割面の移流は境界節点自身の組成を使うため残り、拡散は skip する。[passiveKernels_d.cuh:25](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/passiveKernels_d.cuh:25)、[speciesTransport_d.cu:253](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/speciesTransport_d.cu:253)。x 帯の境界は「CV 集合を切る双対面」であり、補間した平面断面とは別です。入口ピン節点も Ω に含め、`B_pin` を別表示します。

- **作用素の恒等式**：提案式を採用します。`R_final` に体積を再乗算しない。ソースは実際に使った制限後・型変換後の加算値を保存し、早期退出では明示的にゼロを書く。`R_other` は命名された操作の和とし、閉じなかった差をそこへ入れてはいけません。

- **組立誤差尺度**：節点 i の加算項を aᵢⱼ、加算回数を mᵢ として、基本形を  
  `B_assembly = Σᵢ γₘᵢ Σⱼ|aᵢⱼ| + B_underflow + B_double`  
  とします。`γₘ = mu/(1−mu)`、`u = 2⁻²⁴`。内部面の両端への加算も数え、float atomic の非正規化数処理も別途含める。ピンなどの上書きは加算列と分離します。

- **更新写像の恒等式**：  
  `ΣV(q_after−q_before) = ΣVδq_limited + ΣVC_round + ΣₐΣVCₐ + E_bookkeeping`  
  を推奨します。`C_round` は、制限後候補を保存量へ格納する際の丸め差。床・再正規化などと分けます。`δq_update` を格納後の差として定義する方法でも恒等式は閉じますが、その場合も、**格納で失われた増分を別記録**しなければ停止の原因を隠します。前処理の射影・境界上書きも測定区間に含めるか、前処理収支として別に閉じる必要があります。[main.cpp:1836](/home/sano/work/forge-species/solver_density_cuda/main.cpp:1836)

- **10% の尺度**：`ΔF_l = |平均 F_l,exit^ON − 平均 F_l,exit^OFF|` とし、出口マーカーの**本番の外向き面流束**を使います。この node 経路では出口半割面の拡散はゼロです。`wflux.py:17` の節点台形積分や `g_exit_mw` は代用不可。Qₙ は各自の `ΔF_Qn` を使い、平面単位奥行きでは液が `kg/(m·s)`、Qₙ が `mⁿ⁻¹/s`。新しい流束系列にも `check_quasisteady.py` を適用し、差が測定誤差以下なら判定不能とします。物理収支は `|ΣR_final| + B_assembly < 0.1ΔF_c` として、丸めの不確かさで合否が逆転しないことを要求します。`0561−0567` を使うなら「生産 ON−基準 OFF」の差であり、純粋な作用素差とは呼びません。

結論: **格納丸め `C_round` を明示した本番一更新の診断を作り、同一入力の float 格納／double 参照 A/B で G3-b の尺度を確定する。**

第 1 仮説: 提案された G3-b ゲートは、正常な保存量格納の丸めを収支欠落として偽 FAIL にする。確度: **高**。実際の G2 場での寄与の大きさは未確認。  
　根拠: `speciesTransport_d.cu:485` の加算形式について、Python の float32 丸めで `q=1、δq=2⁻²⁶、V=1` を確認した。格納増分は 0、`|E_update|=1.4901161e−8` に対し提案許容は `2.8421709e−14`、**524288 倍の超過**。床や射影を一切必要としない反例。  
　反証条件: 実場の計測で未閉鎖が格納丸めでは説明できず、`C_round` を計上しても導出した演算誤差上界を超えるなら、その未閉鎖の原因としては棄却する。

第 2・第 3 仮説: 現時点では追加しない。G2 の残差プラトーの原因を、この反例だけで格納丸めに断定することもできない。

判別 A/B: `case/16.nozzle_wys/run_0561_g2_onprod_r1/res_48000.h5` から、本番の **1 更新**で得た同一の基底・制限後増分を使う。変えるのは診断側の格納精度だけ：A は本番 float 格納、B は double 参照。B は本番状態へ戻さず、補正操作直前の候補を比較する。見る量は成分別・Ω 別の格納丸め、消失した増分、丸め計上後の未閉鎖。  
　→ **A の未閉鎖が A−B で説明され、計上後に上界内なら第 1 仮説を支持**。  
　→ **計上後も上界外なら「格納丸めだけ」を棄却し、補正項・記録点の欠落を追う**。

やらない方がよいこと: 64ε・16ε の係数を実測超過に合わせて増やすこと。G2 の準定常性を閉じた収支の代替にすること。診断専用の写しが本番残差と大まかに合うだけで G3 を合格にすること。

呼び出し側の前提への異議: 「G2 の比較条件を満たした」と「同じ解」「レシピ非依存」は分ける。また、局所・全補正込みの κ は既存ゲートの再現ではなく追加要件である。**plan 未反映**。反映先は `plans/active/condensation-two-phase-default.md` §5.1 #4g3・#5、§6 G3。

不足情報: 指定された二つの `res_48000.h5` はこの作業環境には存在しない。保存 CSV と VERDICT は確認したが、元 HDF5、正式な補正ゲート全項の記録、抽出照合の実データ、新規の液・Q 出口流束系列は未確認。G3 の実測合否はまだ判断できない。
