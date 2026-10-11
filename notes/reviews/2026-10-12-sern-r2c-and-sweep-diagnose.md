# codex 諮問 (diagnose): sern-r2c-and-sweep

- **brief**: [`notes/reviews/briefs/2026-10-12-sern-r2c-and-sweep.md`](../../notes/reviews/briefs/2026-10-12-sern-r2c-and-sweep.md)
- **plan**: [`plans/accepted/tooling-sern-te-wake-grid.md`](../../plans/accepted/tooling-sern-te-wake-grid.md)
- **date**: 2026-10-12
- **commit**: `dfac3c32` (feature/sern-design)
- **codex**: effort `high`, 3.8 min, rc=0
- **結論**: ?

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表

**Q1 は採用、Q3 は現時点では保留です。R2c の合格と、生産経路全体の受入れは分けて判断します。** Critical 相当の根拠は確認できませんでした。

| 対象・重大度 | 採否 | 根拠と対案・行き先 |
|---|---|---|
| **Q1：R2c** | **採用** | §6 #2 は最低温度の `STEADY` を要求し、#11 は「未定常だけ」を延長対象にしている。[`sern_gates.py:582`](/home/sano/work/forge-sern-design/design/forge_design/metrics/sern_gates.py:582) も温度の未定常を `UNSETTLED`、低温・欠測を別扱いにする。**ALL の `TRANSIENT-UNSETTLED` は正当な延長理由**。下記の保存系列の再判定でも結果を再現した。 |
| **R2h(1)：Major** | **採用・最優先、本 plan の完了前** | [`sern_gates.py:149`](/home/sano/work/forge-sern-design/design/forge_design/metrics/sern_gates.py:149) が、入力の判定不能を無視して `ok` を再構成する。メモリ内試験でも SST の `rms_roOmega` 欠落を `ok=True` とした。**入力の妥当性と収束状態を別の構造化結果として渡し、判定不能は受理・延長とも拒否**する。単に `analyze.ok` を必須にすると、登録済みのプラトー受理まで変わるため不可。 |
| **R2h(2)：Major** | **採用、本 plan の完了前** | [`sern_gates.py:80`](/home/sano/work/forge-sern-design/design/forge_design/metrics/sern_gates.py:80) と同ファイル263行は存在する化学種だけを見る。全 `roY*` 欠落なら組成検査自体が消える。**実効 config から必須集合を決め、存在・配列寸法を先に検査**する。ΣY の検査だけでは、欠落種の質量分率がゼロの場合を検出できない。 |
| **R2h(3)：Major** | **採用、3D 生産経路は本 plan の完了前** | [`runner_sern3d.py:398`](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern3d.py:398) は欠落面を足さず、416行で通常の係数を返す。**形状・境界定義から必要面を固定し、各予定 step で完全性を要求**する。全形状に8部品を一律要求するのではなく、存在しない面と出力欠落を区別する。2D 側は移管済みの2D plan に同じ契約を置く。 |
| **R2h(4)：Major** | **採用、3D 生産経路は本 plan の完了前** | [`runner_sern3d.py:432`](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern3d.py:432) は1部品でも摩擦があれば、他部品をゼロ補完して `C_T_with_shear` を作る。R2a の「列がある」検査では止まらない。**目的量に寄与する全粘性壁の必要 traction 成分を要求し、欠落は判定不能**にする。2D 側は2D plan へ。 |
| **R2h(5)：Major** | **採用。ただし一律のコピー義務化は却下** | [`runner_sern.py:1364`](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern.py:1364) は欠落・寸法不一致を飛ばし、標準起動の1617行でも呼ばれる。R2d の1928行の検査は判定段への移行に限られる。**標準生産経路の初期化と SST→SST 継承は完了前に検査**する。一方、[`runner_sern3d.py:308`](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern3d.py:308) の `keep_turb=True` は Euler からの暖機で初期乱流量を保つ明示仕様なので維持する。一般の変形 restart・暖機経路は状態転送の別 plan に集約する。 |
| **R2h(6)：Major** | **採用、共有ツールの別 plan へ** | [`check_quasisteady.py:418`](/home/sano/work/forge-sern-design/solver_density_cuda/tools/check_quasisteady.py:418) は適用不能量を除き、空集合でも成功し得る。**要求量ごとの適用可否を返し、未評価を `ALL STEADY` にしない**。ただし R2c は明示4列の CSV 経路で、同ファイル361行に列欠落の拒否がある。この欠陥を R2c の不合格理由にはしない。 |
| **R2h(7)：Minor〈本 plan への影響〉** | **採用、共有ツールの別 plan へ** | [`check_floor_events.py:149`](/home/sano/work/forge-sern-design/solver_density_cuda/tools/check_floor_events.py:149) は最終場がない場合に照合を省く。ただし生産ゲートは [`sern_gates.py:74`](/home/sano/work/forge-sern-design/design/forge_design/metrics/sern_gates.py:74) と246行で保存場なしを拒否する。**照合を要求したモードでは場なしを判定不能にし、記録単体の監査は照合未実施を明示**する。本 plan を止める独立の穴とは認めない。 |
| **R2h(8)：Minor** | **採用、親 plan R7b の監視一般化へ** | [`te_monitor.py:74`](/home/sano/work/forge-sern-design/case/46.sern_design/diag/te_monitor.py:74)・86行の可変列は修正すべき。ただし固定形状との不一致は `runner_sern3d.py:70`、空領域は `sern_gates.py:551` で拒否する。**固定 CSV 列・空領域の明示・形状由来の領域生成**を移管する。現状の保護を残した固定形状の受入れまで止める根拠はない。 |
| **Q3：result 再レビュー** | **要再検証** | **(1)〜(4) と (5) の標準生産経路を修正し、正常入力の判定不変・欠落入力の拒否を確認した後なら回してよい**。移管項目は追跡先・適用範囲を本 plan に残す。再レビュー前に `done` と結論づけない。 |

R2c の根拠は `case/46.sern_design/run_1100_r2c_lsw10_m6_judgment/`、写しは [`notes/investigations/2026-10-12-r2c/run_1100_r2c_lsw10_m6_judgment/`](/home/sano/work/forge-sern-design/notes/investigations/2026-10-12-r2c/run_1100_r2c_lsw10_m6_judgment/) です。

- 各 session の40標本を再判定し、保存された `judgment` と一致：1回目 `UNSETTLED`、2回目 `PASS`。
- 2回目の判定窓は **(10000, 20000]**。`C_T_with_shear` の a/d は **2.67e−6 / 3.70e−7**、`C_M` は **2.55e−4 / 2.38e−5**で登録許容内。
- `check_quasisteady.py` の再実行：1回目の温度は **`TRANSIENT-UNSETTLED`**、2回目は温度・力とも **`OVERALL: ALL STEADY`**。
- 床の保存判定は **session 2、10000更新、`VERDICT: PASS`**。残差の保存判定は両 session とも **`NOT CONVERGED (stalled/plateau)`**。`rms_roK` は両回とも `falling` で、全列プラトーではない。
- 保存された力の内訳は全40標本で8部品、ノズル5面の `Fx_tau` 欠落は両回とも0。R2h(3)(4) がこの結果で発現した証拠はない。

**結論:** 次は R2h(1) の「残差の判定不能を受理・延長へ通さない」修正を、下記の欠落列 A/B を受入れ条件として行う。

**第 1 仮説:** 残る誤受理の主因は、入力の不完全性を収束未達と区別せず、上位ゲートで捨てていること。 **確度: 高**  
　**根拠:** `check_convergence.py:183` の判定不能が `sern_gates.py:149` で失われる。実施した100行のメモリ内試験では、完全な SST 残差と `rms_roOmega` だけ欠落した残差の両方が `residual_health.ok=True`。欠落側では `analyze` が明確に「必須の保存量残差列が無い」と返している。  
　**反証条件:** 欠落側が上位ゲートまで `JUDGMENT_UNVERIFIABLE` として伝わり、受理・自動延長の両方を拒否すること。

**第 2 仮説:** 部品・摩擦の欠落による部分和は、列の存在と定常性の検査を通過し得る。確度: 高。ただし R2c での発現は支持されない。  
**第 3 仮説:** 設けない。R2c の物理異常や GPU 経路の不具合を追加で疑う根拠はない。

**判別 A/B:** **forge は0 step**。呼び出し側で保存した R2c 1回目の同一 session を読み取り再生し、A は残差 CSV 完全、B は **`rms_roOmega` 列だけ除く**。config・力・温度・床判定は固定し、`residual_health` → `evaluate_gates` → `judgment_extendable` を確認する。  
→ **A のみ延長可、B は判定不能・延長不可なら修正成立**。**B も延長可なら第1仮説の欠陥が残存**。A まで拒否するなら、登録済みの未収束受理方針を誤って変更している。

**やらない方がよいこと:** R2c を取り直してソフトウェアの欠測処理を検証すること、`require_residual_pass=True` で穴を隠すこと、Euler→SST の明示的初期化まで SST→SST restart と同じ必須コピーにすること。

**呼び出し側の前提への異議:** 「正常な保存 run の判定が前後一致」は欠落入力を拒否する証拠ではない。また、run_1099 との少数点の近さを R2c の合格根拠には採らない。合格根拠は登録された窓・床・温度・延長経路の充足である。

**不足情報:** ローカルの R2c 写しには生の残差 CSV・HDF5 がなく、残差と場の再計算はできない。残差は保存 VERDICT、力・温度は保存系列の再判定に基づく。変形 restart 全経路の安全性は今回未確認。**plan 未反映**であり、反映先は本 plan §5.1 R2h・§6・§6.1、移管先の残作業表とする。
