# codex 諮問 (diagnose): te-wake-blend-length

- **brief**: [`notes/reviews/briefs/2026-10-08-te-wake-blend-length.md`](../../notes/reviews/briefs/2026-10-08-te-wake-blend-length.md)
- **plan**: [`plans/active/convection-zero-thickness-edge-reconstruction.md`](../../plans/active/convection-zero-thickness-edge-reconstruction.md)
- **date**: 2026-10-08
- **commit**: `8c0870d1` (feature/sern-design)
- **codex**: effort `xhigh`, 6.8 min, rc=0
- **結論**: **L_b = 1.0H の既存 Hermite に固定し、上記の最終 node 格子検査を通してから、§6.0 を「局所変形だけで冷点を解消できるか」の A/B として行う。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0 / Major 4 / Minor 1）

| 重大度 | 対象・採否 | 根拠と対案 |
|---|---|---|
| **Major** | **3 次 Hermite・L_b = 1.0H を診断候補として採用。0.3H・0.5H は却下** | g3 の YAML から再生成すると、後縁／復帰部の最大折れは、0.3H が **7.03°／12.05°**、0.5H が **4.32°／8.03°**、1.0H が **2.20°／5.08°**。1.0H に固定し、下記の許容差で最終 node 格子を検査する。g0 の VERDICT は g3 の合格証拠にならない。[生成式](/home/sano/work/forge-sern-design/design/forge_design/meshing/mesh_sern3d.py:420) |
| **Major** | **「近壁間隔が不変だから解像度も同じ」は却下** | 層配置は帯の高さから再計算され、厳密な剛体移動ではない。1.0H の g3 では第一層の Δy はほぼ不変でも、格子線に垂直な間隔 Δy·cosθ は **A 比 0.876〜1.302**。対案は、Δy と垂直間隔を区別し、A/B を「局所格子変形の十分性」の試験とすること。[層生成](/home/sano/work/forge-sern-design/design/forge_design/meshing/mesh_sern3d.py:267)、[間隔の定義](/home/sano/work/forge-sern-design/case/46.sern_design/diag/te_wake_grid_check.py:262) |
| **Major** | **「旧中間線が噴流境界であり、B はそこから 0.11H 外れる」は要再検証** | θ_b は設計点 MOC の**終端角**。再計算では後縁直後 **−58.197°**、終端 **−40.235°**だった。さらに形状・格子の設計点と m10_on の作動点は別である。0.11387H は「旧中間線からの移動量」と記す。実際のせん断層との距離・解像度は、共通の物理断面で組成・速度分布から別途確認する。[角度の選択](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern.py:537)、[後縁の自由境界処理](/home/sano/work/forge-sern-design/design/forge_design/geometry/moc_sern.py:150) |
| **Major** | **`stage_key` 拡張の延期は条件付き採用。「別 run だから来歴確認は不要」は却下** | 現行 `stage_key` の格子署名追加は `w` 有効時だけで、今回の OFF 比較には効かない。今回に限り、各 run の実入力格子から再計算した署名、曲線版、L_b/H、生成コード版、設定、IC の転送記録を保存し、継続時も格子不変を確認する条件で延期できる。診断結果を設計 DB・学習へ取り込まない。生産化では属性の転記だけでなく、実格子の署名を識別に含める。[stage_key](/home/sano/work/forge-sern-design/solver_density_cuda/tools/stage_manifest.py:458)、[既存の署名関数](/home/sano/work/forge-sern-design/solver_density_cuda/tools/mark_zero_thickness_edges.py:105) |
| **Minor** | **「総回転 51°／station 数」を任意曲線の下限とする議論は却下** | 約 51° は、この Hermite の約 −47.8°への行き過ぎを含む値。位置・勾配を戻す制約だけから同じ総回転は決まらない。別曲線なら各 Δx を重みにした変位制約も必要。今回は新しい最適化曲線を追加せず、既存 Hermite を使う。[行き過ぎの式](/home/sano/work/forge-sern-design/design/forge_design/meshing/mesh_sern3d.py:434) |

事前登録する格子の許容差は次を推奨する。**診断への投入条件であり、温度やせん断層解像度の正しさを保証する基準ではない。**

| 確認量 | 許容差 |
|---|---|
| 最初の下流辺 | 下側の第一〜第10層・全 span で、下側接線 m0 の角度から **3°以内** |
| 後縁の折れ | 上下それぞれの上流辺を基準に、第一〜第10層・共有後縁点・全 span で **3°以下** |
| 復帰部の折れ | 接続先の最初の未変形辺まで含め、全 i 線・全 span で **6°以下** |
| 固定すべき量 | x・z、接続、節点数、上流座標、固体形状、外部境界形状は不変。対称面・側方遠方面上の節点の面内移動は区別して記録 |
| 層の Δy | 変換前の第一層は相対差 **1e−8以下**、第一〜第10層の各層厚は **5%以内**。変換後は座標丸めの上限を加える。第一層なら絶対差 ≤ 1e−8·Δy_A + 4ε32·max｜対象座標｜ |
| primal 品質 | 非正・非有限の頂点 Jacobian **0**、scaled Jacobian 最小 **0.65以上**、skew 最大 **0.9以下**かつ最大値の A からの増加 **0.10以下** |
| 品質・双対の判定 | 同じ閾値で `check_mesh_quality.py` が A と同等以上。変形領域に新しい AR 閾値超過を作らない。最終 node 入力で `check_dual_closure.py --tol 1e-5` が **PASS** |

今回、`HEAD 8c0870d1` の [g3 YAML](/home/sano/work/forge-sern-design/case/46.sern_design/problem_3d_prod_3op_wallres_lswx08.yaml:119) からメモリ上で再生成した 1.0H は、第一〜第10層の Δy 変化が下側 **0.179%**・上側 **1.938%**、scaled Jacobian 最小 **0.67245**、skew 最大 **0.44706 → 0.53048**、AR 最大 **4206.09 → 4206.09**、非正頂点 **0**だった。これは変換前の幾何検査で、**変換後の品質・双対 VERDICT は未確認**である。

結論: **L_b = 1.0H の既存 Hermite に固定し、上記の最終 node 格子検査を通してから、§6.0 を「局所変形だけで冷点を解消できるか」の A/B として行う。**

第 1 仮説: この限定した局所変形は、現行の再構成を変えずに m10_on の冷点を解消するのに十分である。　確度: **中**

  根拠: 保存診断集計では `run_1055_r7b_m10_A_c` の節点 517160 において、下流辺への Ux 外挿増分 17591.63 m/s のうち、壁法線側近傍からの寄与が 14967.6 m/s。格子変形はこの辺の向きを大きく変える。ただしこれは既に床へ到達した場の集計で、原因確定ではない。[診断集計:15](/home/sano/work/forge-sern-design/notes/investigations/2026-10-05-sern-r7a/TE_GRAD_LSQ_M10_CALL1.txt:15)

  反証条件: 対照 A が症状を再現し、品質・IC・定常性の条件を満たした B でも、判定領域の冷点が維持・移転する。その場合、**この 1.0H 変形の十分性**を棄却する。

第 2 仮説: 格子変形に伴うせん断層の解像度・数値拡散の変化が、温度改善または悪化を支配する。確度: **中、流れへの影響は未確認**。

第 3 仮説: 勾配係数・float32 適用・状態対応の不整合が残る。確度: **低、正式照合が未実施なので除外しない**。

判別 A/B: **A = `mesh3d.te_wake_blend_H: 0`、B = `1.0` の一対だけ**。共通バイナリ・設定、x 間隔固定、`w`・単面介入は無効。初期場は §4.2 の契約どおり、不動節点を index コピーし、移動節点のみ物理座標に基づいて転送する。

各 **20000 step、未定常なら +20000 を一度**。最初の100 callは早期異常検査、採否は末尾10000 step。固定領域 R と復帰区間で最低温度、低温体積・位置、床補正、4力量を追い、判定区間付きの `check_convergence.py`・`check_quasisteady.py` の VERDICT を残す。温度等の独自量は、既存の `--series-csv` 経路で判定できる。

→ **結果 A**：対照が症状を再現し、B が §6.0 の温度・床・定常性条件を満たすなら、局所変形の十分性を支持し、このケースでの `w` の必要性を棄却する。**折れ角単独の因果や物理解の正しさは確定しない。**

→ **結果 B**：対照が症状を再現し、B も定常化後に冷点を維持するなら、この変形の十分性を棄却する。格子原因全般の否定にはしない。

品質不合格、IC 不備、対照で非再現、延長後も未定常なら判定不能。

やらない方がよいこと: 結果を見ながら L_b を調整すること、数 station で変形を戻して折れを移すこと、位置差を遠方へ逃がして変形範囲まで同時に広げること。短い昇温だけで採用せず、既存の opt-in `w` は予備として保持する。

呼び出し側の前提への異議: **「変更キーが一つ」と「物理的な作用が一つ」は異なる。** 今回は折れ、辺長、層の向き、流れに対する節点配置が同時に変わる。また、A の他所の最大折れ 5°は、後縁でも安全という根拠にはならない。

不足情報: この checkout に `case/46.sern_design/run_1055_r7b_m10_A_c/` と `run_1061_r5h_m10_A_1step/` の原本がなく、流れ・収束・定常性の再検証はできていない。必要なのは、実際に投入する A/B の node HDF5、品質・双対 VERDICT、IC 転送記録、判定区間付きの流れの VERDICT である。

**plan 未反映。** ファイル変更・`forge` 起動は行っていない。呼び出し側で対象 plan の §4.2・§6.0・§5.1 #2c/#2e に反映すること。
