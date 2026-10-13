# codex 諮問 (diagnose): tw370-registration

- **brief**: [`notes/reviews/briefs/2026-10-13-tw370-registration.md`](../../notes/reviews/briefs/2026-10-13-tw370-registration.md)
- **plan**: [`plans/active/tooling-nozzle-isothermal-wall-chain.md`](../../plans/active/tooling-nozzle-isothermal-wall-chain.md)
- **date**: 2026-10-13
- **commit**: `b2acab60` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 2.7 min, rc=0
- **結論**: **§6 の合格条件・移植／分岐・抽出検査を修正し、入力の事前検査を通してから、登録した 3 本だけを実行する。**
- **extra**: `case/45.isobutane_m6_d155/cold_xcheck.py`, `case/45.isobutane_m6_d155/sp_judge.py`, `case/45.isobutane_m6_d155/vgc_run.py`, `notes/reviews/2026-10-13-vgc-result-diagnose.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：**3 本構成は採用する。ただし、現登録のままの投入は推奨しない。** 主 A/B の比較漏れと、ゼロ近傍の許容の不連続は、結果を見る前に修正が必要。Critical に相当する根拠は確認できなかった。

| 重大度 | 論点・採否 | 根拠と対案 |
|---|---|---|
| **Major** | **現在の結果 A 条件：却下** | [plan:345](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-isothermal-wall-chain.md:345) は C1/R0・C50/R0 だけを要求し、**C50/C1 が抜けている**。前回の諮問は「主 A/B と G0 への直接比較がすべて許容内」を要求している。[前回:45](/home/sano/work/forge-integ-1005/notes/reviews/2026-10-13-vgc-result-diagnose.md:45)。`sp_judge.compare` を合成入力で実行すると、θ_r が R0 比 **+0.09 %／−0.09 %** の両腕はそれぞれ合格するが、C50/C1 は **−0.179838 % で不合格**になる。**3 組すべての合格を結果 A の条件に加える。** |
| **Major** | **ゼロ近傍だけ絶対許容へ切り替える規則：却下** | [plan:344](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-isothermal-wall-chain.md:344)。登録の Q₃₀₀ = 9.2 MW を使った算術では、参照値が 9,190.8 W なら許容 9,200 W、9,200 W なら許容 9.2 Wとなり、境界で約 **1,000 倍**変わる。また [sp_judge.py:47](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/sp_judge.py:47) は相対比較だけで、ゼロ除算も残る。**連続した混合許容 `abs(a−b) ≤ max(0.001·abs(b), a_q)` にする。** 登録した絶対許容を維持するなら `a_q = 0.001·abs(q₃₀₀)` とし、各断面・各量の正確な q₃₀₀ と出典を固定する。この式ではゼロ近傍以外にも絶対許容が効くことを明記する。 |
| **Major** | **ドナーの選択：採用。同一格子の R0 に `interp_field` を使う説明：却下** | 前の壁温の停止場を開始場に選ぶことは、今回の運用試験の意図に合う。ただし [plan:338](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-isothermal-wall-chain.md:338) の「S1 が入るのでビット一致にこだわらない」は、移植誤差が停止時までに消える証拠ではない。[AGENTS.md:135](/home/sano/work/forge-integ-1005/AGENTS.md:135) も同一メッシュへの `interp_field.py` を禁止している。**R0 は `restart_field.py`、G0c は格子間移植、両者とも登録の S1→S2 を通す。** ドナーのファイル名・SHA、移植後の補正内容と入力 SHA を固定する。 |
| **Major** | **C50 を C1 の S2 後から分岐：採用。ただし引き渡しの実装は要再検証** | [plan:342](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-isothermal-wall-chain.md:342) の方針は正しい。しかし既存の通知点では前段の `res_*` が既に消えていることを要求する。[vgc_run.py:358](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/vgc_run.py:358)。**C1 が本段を始める前に、実際に読む入力を固定保存し、その場を C50 に渡す。** 両腕の本段開始時の保存量・乱流量・組成のビット一致と格子一致を検査する。C50 の記録では、S1・S2 は C1 から継承した履歴とし、既存判定器の「各 run に S1・S2 必須」をそのまま適用しない。[vgc_judge.py:205](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/vgc_judge.py:205)。 |
| **Major** | **抽出帯・出口対応を確認済みと扱うこと：要再検証** | [plan:348](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-isothermal-wall-chain.md:348) の R0 単独の感度では、**腕間差が帯に依存しないことを確認できない**。抽出器は旧 run の scale と帯を読む。[cold_xcheck.py:46](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/cold_xcheck.py:46)。出口比較も座標照合なしで配列を引く。[sp_judge.py:55](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/sp_judge.py:55)。**固定帯は維持し、3 腕すべてで帯 1／1.25 倍の腕間差を記録する。帯で採否が変われば抽出依存として判定保留にする。** 生成後の r_t、断面位置、出口全節点の `(x,r)` 対応を検査し、不一致・抽出値の非有限は数値的失敗 B と区別して判定不能にする。 |
| **Major** | **GPU 共有での同時実行：条件付き採用** | step 基準の判定なので、共有自体は否定しない。監視の停止対象も run の cwd で限定されている。[m9_watch.py:94](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9_watch.py:94)。ただし既存の 4.5 GiB 検査は **`/proc/meminfo` のホスト RAM** で、GPU メモリではない。[vgc_run.py:528](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/vgc_run.py:528)。**既存ジョブを含む GPU メモリ・ホスト RAM・出力容量を確認し、各腕は別プロセスで起動する。** `Observer` はモジュールの関数を差し替えるため、同一プロセス内のスレッド並列には流用しない。[vgc_run.py:312](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/vgc_run.py:312)。資源不足による終了は判定不能とする。 |
| **Minor** | **δ_E・δ_target の許容決定を #33 に分離：採用** | 今回を固定壁の停止時運用に限定するなら妥当。[plan:349](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-isothermal-wall-chain.md:349)。ただし **#32 の合格を壁更新の許可にしない**。今回も抽出方法・評価位置・時系列窓・`check_quasisteady` の引数を固定し、#33 では壁更新を始める前に絶対寸法の許容を登録する。 |

結論: **§6 の合格条件・移植／分岐・抽出検査を修正し、入力の事前検査を通してから、登録した 3 本だけを実行する。**

第 1 仮説: **370 K・固定壁でも、G0c の上限 1／50 は相互に、かつ G0・上限 1 に対して停止時の比較要件を満たす。** 確度: **低**。  
根拠: [plan:316](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-isothermal-wall-chain.md:316) に 300 K の別開始手順での結果 A が記録されている。ただし、370 K の実測根拠はまだなく、今回その原データも独立再判定していない。  
反証条件: 入力・抽出・監視が成立した条件で、いずれかの腕が発散／200,000 step まで未到達、または **C50/C1・C1/R0・C50/R0 のいずれか**が登録許容を超える。

第 2・第 3 仮説: 追加しない。今回の登録点検から発散機構や壁温応答の真因までは推定しない。

判別 A/B: **G0c の同一本段入力から `lineDtDirectionalCap: 1 → 50` だけを変える。** 共通参照として同じ壁温・物理壁の G0・上限 1 を置く。対象は次の 3 本。

- R0：`case/45.isobutane_m6_d155/run_0811_tw370_g0_c1/`
- C1：`case/45.isobutane_m6_d155/run_0812_tw370_g0c_c1/`
- C50：`case/45.isobutane_m6_d155/run_0813_tw370_g0c_c50/`

各本段は最大 **200,000 step**、出力 **2,500 step**、既存の停止規則を維持する。M1〜M7 を上記 3 組で比較し、全残差の `check_convergence --segment` と派生量の `check_quasisteady` は別に記録する。  
→ **A**：全腕到達・3 組すべて許容内・抽出検査成立なら、この壁温・固定壁の停止時運用への拡大を支持する。  
→ **B**：有効な試験で発散・未到達・許容超過なら拡大を見送る。上限間だけが外れれば上限変更の互換性を否定できる。両 G0c 腕が R0 から外れる場合は、粗格子と登録した移植手順を含む運用の不成立であり、格子単独の原因とは断定しない。

やらない方がよいこと: **R0 への 2 比較だけで主 A/B を合格にすること、S1 が移植誤差を消すと仮定すること、抽出不成立や資源不足を数値手法の失敗に数えること、今回の合格から自動で δ* 壁更新へ進むこと。**

呼び出し側の前提への異議: 「370 K の境界層は 300 K と断熱の間に入る」は未確認の予測。「300 K の 3 本で帯感度が小数 3 桁まで同じ」はブリーフの観測報告であり、370 K の帯の妥当性を保証しない。「時間を判定に使わない」は共有実行を認める理由になるが、資源確認を省略する理由にはならない。

不足情報: 当地にはドナー・対象 run の実ディレクトリ、および `tw370_run.py`・`tw370_judge.py` がない。ドナー SHA、帯感度の原データ、生成後の座標対応、共有 GPU の余力は未確認。今回実行したのは判定関数の合成入力検査と許容式の算術確認で、forge は起動していない。

**plan 未反映。** 呼び出し側の反映先は `tooling-nozzle-isothermal-wall-chain.md` §6「V-tw370 の事前登録」・§5.1 #32／#33。ファイルは変更していない。
