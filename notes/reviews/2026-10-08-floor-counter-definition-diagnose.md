# codex 諮問 (diagnose): floor-counter-definition

- **brief**: [`notes/reviews/briefs/2026-10-08-floor-counter-definition.md`](../../notes/reviews/briefs/2026-10-08-floor-counter-definition.md)
- **plan**: [`plans/active/tooling-sern-te-wake-grid.md`](../../plans/active/tooling-sern-te-wake-grid.md)
- **date**: 2026-10-08
- **commit**: `30348fa4` (feature/sern-design)
- **codex**: effort `xhigh`, 7.4 min, rc=0
- **結論**: **A＋同じ前処理・判定関数による終了時監査を採用し、最後の更新だけに床未満状態を与える小型 A/B で計測の閉じ方を検証してください。**
- **extra**: `plans/active/convection-zero-thickness-edge-reconstruction.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0 / Major 4 / Minor 1）

| 重大度 | 対象・採否 | 根拠と対案 |
|---|---|---|
| **Major** | **C は却下。既存レビュー M1・M2 の趣旨は採用維持** | EOS は作業配列を変更しますが、定常 commit は `roN/roeN + dq` です（[dependentVariables_d.cu:207](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/dependentVariables_d.cu:207)、[update_d.cu:301](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/update_d.cu:301)）。**「commit に残った補正」ではなく「EOS の床を必要とした事象」**を数え、共通バイナリ・生産ゲートへの接続・欠測拒否を維持します。 |
| **Major** | **A＋終了時の読み取り判定を、時点の契約を修正して採用** | EOS 前には no-slip・軸・等温壁ピンがあります（[main.cpp:2012](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:2012)）。終了時だけ未射影の commit 状態を読むと、途中と判定対象が違います。**終了時はコピー上で同じ前処理を再現し、同じ判定関数を適用**してください。通常配列への EOS 再適用は禁止。入口状態・更新番号・終了時判定を区別し、同じ状態を二重計数しません。 |
| **Major** | **A の「再構成差が丸め上限を超えたら温度床」は却下** | 密度床は速度・内部エネルギーの算出より先、`roe` 再構成は温度床以外でも実行されます（[dependentVariables_d.cu:76](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/dependentVariables_d.cu:76)、同ファイル:209）。**床の種類別の述語で数え、補正量は付帯情報**にします。温度は EOS と同じ密度・組成での `e_in < e_mix(T_min)` と最終温度下限制約を確認。密度は `ρ_in < roMin`、圧力はクランプ直前の `P_raw < pMin`。TP 圧力床の量は **ΔP**、その操作による直接の **Δ(ρE) は 0** です。 |
| **Major** | **「dual-time では A＝C」は却下。適用拡大は要再検証** | in-place 加算は確認できますが（[update_d.cu:517](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/update_d.cu:517)）、TP 圧力床はそこでもエネルギー補正ではありません。また初期化では EOS 後に `updateVariablesOuter` が走るため、定常でも初期床補正は基準状態へ入ります（[main.cpp:1837](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:1837)、同ファイル:1866）。事象の意味は統一し、経路ごとに計測時点を登録してください。未検証経路をゼロ件で通してはいけません。 |
| **Minor** | **床近傍の別集計は採用** | 現行ゲートは保存された温度等の閾値判定です（[sern_gates.py:200](/home/sano/work/forge-sern-design/design/forge_design/metrics/sern_gates.py:200)）。`T ≤ T_min + 1 K` の件数・位置は補助情報として保持し、床事象件数と混ぜません。ちょうど床上で補正不要な状態も、こちらで可視化できます。 |

**丸めの扱いは、この生産受入れでは除外閾値 0 を推奨します。**  
これは「EOS 前後で `roe` が少しでも変われば失敗」という意味ではありません。上表の床述語が成立した事象だけを数え、補正量による足切りをしないという意味です。通常状態の再構成丸めは数えません。床作用が float32 の格納で消える場合も見逃しません。根拠のない固定許容値を「丸め上限」として導入する案は採用しません。これは元 plan の丸め除外条件を厳しくする変更です。

**判定区間 `(a, b]` は、区間入口と各更新結果を覆う契約にします。**  
入口 `Q_a` を別枠で記録し、`Q_(a+1)…Q_b` を次の EOS 入力、最後だけ終了時監査で確認します。すべて境界ピン後・EOS 床前という同じ状態の定義にそろえます。初期化・restart の EOS は別記録とし、リセットで痕跡を失わせません。定常の `nStepInner=5` は保存量を更新しない線形 sweep なので、**5 更新とは数えません**（[main.cpp:2177](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:2177)）。

結論: **A＋同じ前処理・判定関数による終了時監査を採用し、最後の更新だけに床未満状態を与える小型 A/B で計測の閉じ方を検証してください。**

第 1 仮説: 現行定義の見逃し原因は、床の使用を「後続 commit に補正が残ること」で判定している点にある。確度: **高**。  
  根拠: `dependentVariables_d.cu:203–209` と `update_d.cu:301–308` の処理順。保存記録でも `case/46.sern_design/run_1078_tewake_A0_m10/` は 500～20000 step の全40標本で床近傍1節点、節点517160は50 Kです。ただし毎更新の事象数は未取得です。  
  反証条件: 同一入力・同一増分を使う定常 commit 試験で、EOS の補正が `roN/roeN` またはアキュムレータへ伝播していることが示された場合。「全 run で非保持」という一般化は撤回します。

第 2 仮説: `roe_after − roe_before` だけの温度床判定には、密度床と通常の再構成誤差が混入する。確度: **中**。処理上の混入経路は確認済みですが、対象 run での誤分類数は未確認です。  
第 3 仮説: 無し。

判別 A/B: **本番流れ計算ではなく、2更新の小型カーネル試験**にします。変更するのは最後の commit の1内部節点の `ρE` だけ。密度・運動量・組成・種DB・前処理を固定し、A は床エネルギーより十分上、B は十分下へ置き、直後に終了します。差幅は `c_v(T_min) × 1 K` 相当以上かつ格納 `ρE` の64 ULP以上とします。終了時監査と、同じ状態のコピーへ既存 EOS を適用した結果を照合します。  
  → **結果A：A枝0件、B枝は最後の更新・指定節点に温度床1件、コピー上のEOSも床到達、監査前後の実配列はビット不変**なら、提案した計測契約を支持します。  
  → **結果B：B枝を見逃す、更新番号がずれる、EOSとの分類が異なる**なら、A＋終了時監査の実装を棄却し、時点・述語の不整合を修正します。A枝の誤検出も不合格です。

やらない方がよいこと: **この計測 plan の中で、EOS 床を `roN/roeN` に commit する修正を入れること。** 更新則を変え、既存 A/B と共通バイナリの前提を壊します。非保持そのものの是非は別 plan で扱ってください。本 plan では、床使用ゼロを実証できなければ受入れ不可とするので十分です。床使用ゼロでも、起動過渡の影響まで消えたとは主張できません。

呼び出し側の前提への異議: **「定常では構造上すべて0」「dual-timeではAとCが一致」「終了時のBは途中のAと同じ状態を見る」には同意しません。** 上表の初期化・圧力床・境界ピンの例外があります。また「費用ほぼ0」は未実測です。

不足情報: 対象 `run_1078/1079` の保存量 HDF5・残差CSVはローカルにありません。写しの収束記録は双方 **`NOT CONVERGED (stalled/plateau)`**。温度CSVは公式 `check_quasisteady.py --tail 0.5` で再判定し、**10000 < step ≤ 20000、双方 `OVERALL: ALL STEADY`**でした。これは毎更新の床使用ゼロの証拠にはなりません。

ファイル変更・`forge` 起動なし。**plan 未反映**。反映対象は `tooling-sern-te-wake-grid.md` §4・§5.1 #2・§6 #1 と、起点 plan §5 の3です。
