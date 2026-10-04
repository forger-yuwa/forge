# codex 諮問 (diagnose): twophase-g1-interpretation

- **brief**: [`notes/reviews/briefs/2026-10-04-twophase-g1-interpretation.md`](../../notes/reviews/briefs/2026-10-04-twophase-g1-interpretation.md)
- **plan**: [`plans/active/condensation-two-phase-default.md`](../../plans/active/condensation-two-phase-default.md)
- **date**: 2026-10-04
- **commit**: `32312c79` (feature/species-transport)
- **codex**: effort `high`, 5.5 min, rc=0
- **結論**: **G2 はまだ投入せず、G1 判定器の不整合を修正し、保存面入力の 0 step 演算監査で丸め上界の適用条件を閉じる。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表:

| 重大度 | 判断 | 根拠と対案 |
|---|---|---|
| **Major** | G1 の成立確定は**要再検証** | **判定器に不整合がある。** [`g1_judge.py:166`](/home/sano/work/forge-species/notes/investigations/2026-10-04-twophase-g1/g1_judge.py:166) は改訂判定を `v2_used` に保存するが、`:205`・`:257` は旧判定 `v2` を `finish` に渡す。現行関数を読み取り実行し、`(i)=PASS、旧(ii)=FAIL、改訂(ii')=PASS` でも **exit 1** になることを確認した。**終了判定と `(iii)` の注記を選択した基準に統一し、旧基準 FAIL は別記録として保持する。** |
| **Major** | 非正規化域が旧基準 FAIL の原因という解釈は**採用**。「0.995 が導出の完全性を証明する」は**却下** | 保存された [`SCALE_AB_VERDICT.txt:1`](/home/sano/work/forge-species/notes/investigations/2026-10-04-twophase-g1/evidence/SCALE_AB_VERDICT.txt:1) は、A のビット一致、574 面の超過、B の超過 0・最大比 0.1889 を記録する。ただし、上界への接近は省略項の不存在を証明しない。**演算ごとの誤差伝播と適用範囲を明記する。安全率を観測値に合わせて増やさない。** |
| **Major** | G1(i) を作用素変更の根拠に使うことは**採用**。物理精度・壁温変化の説明への拡張は**却下** | [`G1_VERDICT_revised.txt:26`](/home/sano/work/forge-species/notes/investigations/2026-10-04-twophase-g1/evidence/G1_VERDICT_revised.txt:26) の数値から、乱流域の OFF 総水分流束最大値は ON 液流束最大値の **9.92×10⁻⁶ 倍**。言えるのは、この入力で旧作用素がほぼ運ばない状態でも新作用素は液を運ぶこと。**既定化の根拠は混合作用素の整合に限定し、CFD の物理精度は未検証と書く。** |
| **Major** | G2 の現行事前登録だけで投入することは**却下** | [`condensation-two-phase-default.md:104`](/home/sano/work/forge-species/plans/active/condensation-two-phase-default.md:104) は σ の定義・比較値の集約方法・判定区間を固定していない。親 [`condensation-two-phase-transport.md:130`](/home/sano/work/forge-species/plans/active/condensation-two-phase-transport.md:130) も固定マスクの数値照合を未実施としている。**下記の条件を先に固定する。** |
| **Minor** | G3 の設計・診断実装を G2 と並行することは**採用** | [`condensation-two-phase-default.md:89`](/home/sano/work/forge-species/plans/active/condensation-two-phase-default.md:89) の作用素収支と更新写像収支の分離は妥当。**G2 用バイナリを固定して並行開発し、最終 G3 判定には G2 の受入状態を使う。** 診断追加が場を変えるなら同じ検証系列として扱わない。 |

G1 の導出はもう一段だけ厳密にする必要がある。実装は [`twoPhaseDiffusion_d.cuh:92`](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/twoPhaseDiffusion_d.cuh:92) で逆数と `ctg` を計算し、`:114` で液流束を作る。一方、判定器の `ctg` は入力を double に上げた積である（[`g1_judge.py:155`](/home/sano/work/forge-species/notes/investigations/2026-10-04-twophase-g1/g1_judge.py:155)）。

必要なのは、①逆数・`ctg` が正規化数という前提の全対象面での確認、②それらの相対丸めを含む上界、③絶対誤差項の後続演算による増幅、④実際の FMA 縮約への対応である。`8ε₃₂ = 16u` には相対誤差側の余裕があるが、plan の「6u 以下」という記述だけでは厳密な上界になっていない。**0.995 は最終丸めの半刻みに近いだけでも生じ得るため、全演算の上界が同時に実現したとは言えない。** 最大比の面について、誤差・相対項・絶対項・中間値を分解すれば確認できる。

G1(iii) は**登録規則では棄却を維持**する。壁直近の蒸気方向は尺度以下であり、逆向きの機構が実証されたわけではない。帯内の内部面の射影和も閉じた収支ではないため、壁温低下や正味凝縮量の説明には使わない。

G2 は、G1 を閉じた後に次の条件で登録する。

- **IC:** `case/16.nozzle_wys/run_0482_passive_wys_s1_sfr2_c1/res_48000.h5` を共通起点とし、同一メッシュの `restart_field.py` でコピーする。これは旧作用素で得た場であり、**「OFF の場から始めない」という条件とは両立しない**。旧場から ON で開始してよい。ON 開始以降を別区間とし、過渡を判定から除く。
- **本数・長さ:** 各 ON レシピについて**独立プロセスで最低 3 本、初期予算 48000 step、400 step ごと保存**を推奨する。2 本では run 間分散の自由度が 1 しかない。48000 step は合格保証ではなく、未定常なら全反復を同じ長さだけ延長する。参照側も反復が必要で、`run_0521` 1 本を誤差ゼロの基準にしない。
- **判定区間:** ON 開始後の全系列に登録済み `--tail 0.5 --drift 0.0001 --osc 0.0001 --min-snaps 21` を適用し、48000 step 時は末尾約 24000 step を評価する。7 量すべて `STEADY`、全保存量の検査、全残差列の `RISING 0`、補正ゲートを確認する。親 #4j の受入を使う場合、残差停滞を許しても**「残差収束」とは書かない**。
- **σ と比較値:** 比較値は各 run の判定窓平均。σ は独立 run の窓平均の標本標準偏差とし、時系列の多数点を独立反復として数えない。生産・参照の差には両側のばらつきを伝播する。**`|ON−OFF| > 3σ` だけでは、その差の 10% を識別できるとは限らない。** 比較差の不確かさが 10% 許容幅をまたぐ場合は判定不能とする。ON−OFF 基準にも同じ集約方法と不確かさ評価が必要。
- **抽出照合:** 共通 IC の `g > 10⁻⁶` から M0 を一度だけ作る。全 run 間で節点順・座標、M0 の不一致数 **0**、中心線の両側節点と補間重み、上壁節点、出口列と双対長、体積重みを数値比較する。さらに**同一スナップショットを各 run の抽出設定で処理し、7 量が一致すること**を確認する。現行 [`twophase_ab_series.py:213`](/home/sano/work/forge-species/case/16.nozzle_wys/twophase_ab_series.py:213) のマスク長一致だけでは不十分。

結論: **G2 はまだ投入せず、G1 判定器の不整合を修正し、保存面入力の 0 step 演算監査で丸め上界の適用条件を閉じる。**

第 1 仮説: **旧 G1(ii) の FAIL は非正規化域の丸め尺度不足で説明でき、本番液流束の式を変更する必要はない。**　確度: **高**  
　根拠: `case/16.nozzle_wys/run_0560_g1_faces_off0520/` の保存判定では、指数スケール A/B が `SUPPORT`、改訂尺度は超過 0。コードの液流束は `twoPhaseDiffusion_d.cuh:114`。  
　反証条件: 実演算条件から導いた上界を本番評価が超える、または保存値を再現するはずの評価がビット一致しないこと。

第 2 仮説: **改訂尺度の説明が特定の FMA 縮約に依存し、縮約によらない上界としては不完全。**　確度: **中・未確認**。`ctg`・逆数の丸めと伝播係数の説明が不足している。

判別 A/B: **同じ保存面入力・同じ GPU で、変更点を試験プログラムの FMA 縮約許可だけにする。** A は本番と同じ条件、B は `--fmad=false`。時間更新 **0 step**、全評価面の液流束を固定した改訂尺度と比較する。A は保存値とのビット一致を前提とする。  
→ **両側が尺度内なら**縮約の有無で尺度が破綻する説を対象入力上で退け、導出・判定器を整えて G2 へ進む。  
→ **B だけ超過なら**縮約によらない上界という主張を撤回する。**A が超過・再現不一致なら**G1 を保留する。いずれも観測後に安全率を足さない。

やらない方がよいこと: **最大比 0.995 を導出の正しさの証明にすること、`PASS` の表示だけで exit 1 を無視すること、反復・マスク未照合の `run_0521` を確定基準にして長時間 G2 を先行すること。**

呼び出し側の前提への異議: **「改訂尺度での数値比較 PASS」と「G1 判定系の成立」は現状では別である。** また、OFF 流束の小ささだけから `∇Y_w ≈ 0` を逆算できず、係数と補正項の確認が要る。

不足情報: **`tp_faces.h5`、最大比を与える面の中間値、実ビルドの演算条件、参照 run の系列 CSV・判定区間付き VERDICT・マスク照合結果。** ローカルには判定テキストのみで、一次データの独立再集計はできていない。ファイル変更・forge 起動は行っていない。**plan 未反映**。呼び出し側で §5.1 #4sr・#4rr・#5、§6 G1・G2 に反映すること。
