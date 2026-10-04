# codex 諮問 (diagnose): sern-next-priority

- **brief**: [`notes/reviews/briefs/2026-10-05-sern-next-priority.md`](../../notes/reviews/briefs/2026-10-05-sern-next-priority.md)
- **plan**: [`plans/active/tooling-nozzle-sern-chain.md`](../../plans/active/tooling-nozzle-sern-chain.md)
- **date**: 2026-10-05
- **commit**: `f97805ac` (feature/sern-design)
- **codex**: effort `high`, 9.1 min, rc=0
- **結論**: ?
- **extra**: `plans/active/tooling-nozzle-sern-3d.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：

| 重大度 | 指摘・採否 | 根拠と対案 |
|---|---|---|
| **Major** | **「MOO は 2D、3D は代表点確認」を却下** | [本体 plan:442](/home/sano/work/forge-sern-design/plans/active/tooling-nozzle-sern-chain.md:442) では、ユーザ判断で **3D SST を最適化ループに入れる**方針へ変更済み。一方、[現行 driver:33](/home/sano/work/forge-sern-design/design/forge_design/opt/driver_sern.py:33) は 2D runner を使用し、設計変数にも `L_sw` がない。対案は、まず既存の 3D runner による `L_sw` の直接比較。2D は形状生成・傾向把握に使う。 |
| **Major** | **R6(d) の先行を採用。ただし数値許容と設計制約を分ける** | [生産 YAML:132](/home/sano/work/forge-sern-design/case/46.sern_design/problem_moo_frozen_tp_cycle3op_transport.yaml:132) は加重平均の `cm_min: -7.0` のみ。作動点別 `cm_window` は未指定。**数値差の許容、検出したい設計差、実機のトリム許容窓は別物**。前二者を固定すれば判別試験には進めるが、最後が未定のまま最終 Pareto の実用成立性は認定できない。 |
| **Major** | **「R5h は局所帳簿以降未着手」を却下** | [3D plan:1933](/home/sano/work/forge-sern-design/plans/active/tooling-nozzle-sern-3d.md:1933) に速度のみの介入と、16000 step の力感度試験まで記録済み。恒久的な再構成変更を今回の試験の前提に戻さない。既存結果の適用範囲を明記し、新形状で低温域・床到達を監視する。 |
| **Major** | **`L_sw` をそのまま連続変数として探索する案は要再検証** | [メッシャ:148](/home/sano/work/forge-sern-design/design/forge_design/meshing/mesh_sern3d.py:148) は終端を最近傍の既存 station に丸める。現行 g3 問題の station 計算では、指定 `0.8` と `0.8001` がともに実位置約 `0.809103 H` になる。対案は、今回だけは**異なる既存 station を明示して比較**し、連続最適化への接続時に物理 station 化する。 |
| **Major** | **「farfield が片付いたので R4 全体も完了」を却下** | [farfield plan:165](/home/sano/work/forge-sern-design/plans/accepted/boundary-node-farfield-characteristic.md:165) の #4h は現形状・`m6_on`・g4・幅 2.50 H の LJ 比較。新 LJ の広幅応答は未測定。[3D plan:1942](/home/sano/work/forge-sern-design/plans/active/tooling-nozzle-sern-3d.md:1942) の上方・下方・出口距離、同時拡大、他作動点も別の残件。限定した判別試験と、最終順位の認定を分ける。 |
| **Minor** | **R4e-次・R4c・R5q の古い残作業記述を整理する案を採用** | [W2:1483](/home/sano/work/forge-sern-design/plans/active/convection-node-wall-reconstruction.md:1483) は実施済みで、両設定とも step 28 NaN。再投入を前提にしない。ただし後続レビューは「発散を防げなかった」までに限定しており、再構成一般の原因除外には使えない。旧トポロジ修正は [R4e:1914](/home/sano/work/forge-sern-design/plans/active/tooling-nozzle-sern-3d.md:1914) の実装・検証へ対応付ける。また新メッシャは [B1c:775](/home/sano/work/forge-sern-design/plans/active/tooling-sern-mesh-blocking.md:775) で 245 万節点案まで進んでおり、「1250 万で載らない」のみでは現況を表さない。 |

**結論:** R6(d) の判定条件を事前固定し、現行 3D・`m6_on`・g3 で MOC 輪郭を固定した `L_sw` 二水準×三反復の判別試験を、R7 の最初の一手にする。

**第 1 仮説:** 指定する二つの側壁長の間には、現行 3D 評価の再実行変動を超える、実務上意味のある `C_T_with_shear` 差がある。  
**確度:** 低。`L_sw` 感度そのものは未測定。

- **根拠:** `L_sw` は実際に側壁終端を変える変数で、採用済み方針でも最初の 3D 設計変数に挙げられている。ただし、既存の 2D–3D 差は `L_sw` 感度の証明ではない。
- 評価器側については、R5h の長時間比較が、局所冷点の変化と積分力への感度を分離している。対象は  
  `case/46.sern_design/run_0984_r5h_sens_alpha0/`  
  `case/46.sern_design/run_0985_r5h_sens_alpha1/`。
- 回収済み判定記録は両者とも **`GATES: PASS`、4 量 `STEADY`、`NOT CONVERGED (stalled/plateau)`**。収束解の一致とは扱わない。[証拠索引:100](/home/sano/work/forge-sern-design/notes/investigations/2026-10-03-farfield-evidence/INDEX.tsv:100)
- 最終 step 16000・末尾 20 標本の回収統計から、`|平均差|＋双方の振幅` は概算で `C_T_with_shear = 7.3×10⁻⁶`、`C_M = 6.9×10⁻⁴`。小さいが、**この介入・この形状の感度であり、新しい側壁長の誤差上限ではない**。[窓統計:1](/home/sano/work/forge-sern-design/notes/investigations/2026-10-03-farfield-evidence/V3_WINDOWS.txt:1)

**反証条件:** 下記の品質条件を満たした比較で、変動幅を加えても `|ΔC_T_with_shear| ≤ 0.002` なら、「この二水準に 0.002 を超える有用な差がある」という仮説を棄却する。

**第 2・第 3 仮説:**  
第 2 仮説：この二水準の差は 0.002 以下である。確度：低。未測定。これが支持されても、他の側壁長や 3D 設計変数まで無効とは結論しない。  
第 3 仮説：追加しない。

**判別 A/B:**

**一つの形状比較だけ**を行う。A は `L_sw` 指定 0.8、B は 1.0 に対応する既存 station。今回確認した g3 問題では実位置は約 **0.809103 H / 1.004328 H**。生成後の実位置を台帳に残す。

- **本数:** 2 形状×3 反復＝**6 本**。最初は `m6_on` のみ。多作動点 MOO の成立確認ではなく、固定条件での判別能力の試験とする。
- **固定するもの:** MOC 輪郭、幅、側壁以外の形状、z 座標列・`first_z_frac`・`sz` テーパ、物性、BC、リミッタ基準値、バイナリ、CUDA block size。診断介入は無効にする。
- **初期場:** 保存されている既存場から現行条件へ適応させる。同一メッシュは `restart_field.py`、異なるメッシュは `interp_field.py` を使い、同座標の壁内外ノードの取り違えと、新しい壁距離を確認する。削除済みの `run_1029` 全場を起点として指定しない。
- **長さ:** 適応後の同一設定区間を各 20000 step、500 step 間隔で保存。末尾 10000 step と直前 10000 step を比較し、平均差が `C_T`・`C_T_with_shear`・`C_L` で `5×10⁻⁵` 以下、`C_M` で `5×10⁻⁴` 以下。未達なら各 +20000、なお未達は判定不能。
- **必要条件:** メッシュ品質・有限値・床ゲートを通過し、4 量が `check_quasisteady.py` で `STEADY`。`check_convergence.py` の判定も区間付きで併記する。既存のユーザ判断に従い、残差プラトーだけで失格にはしない。

主判定量は `q = C_T_with_shear`。各群の三反復について、末尾平均の最大差を `r`、最大窓内振幅を `a`、最大前窓差を `d` とし、事前に

`N = max(3r, a, d, 10⁻⁶)`、`E = N_A + N_B`

と定義する。これは経験的な判別幅で、統計的な信頼区間ではない。**0.002 は今回提案する検出目標であり、既存の数値誤差許容と同一の意味ではない。**

→ **分岐 A:** `|平均_B − 平均_A| − E > 0.002`  
第 2 仮説を棄却。この二水準を識別できる。符号から優位側を記録し、作動点別制約と格子・領域感度を確認する段階へ進む。

→ **分岐 B:** `|平均_B − 平均_A| + E ≤ 0.002`  
第 1 仮説を棄却。この二水準では設定した規模の利得を示せない。大規模 DOE には進まない。

→ **中間／品質条件未達:** 判定不能。小数点の一致や平均値だけで採否を決めない。

**やらない方がよいこと:**

- R5h の全域一次化・局所帳簿、W2 の既実施比較を再び最優先にする。
- R5h の単一形状の感度を、設計箱全域の誤差上限として持ち越す。
- 現行 2D driver に `L_sw` を追記するだけで、3D 最適化になったと扱う。
- g3/g4 の二水準差を真値に対する誤差上限と呼ぶ。
- 新メッシャ完成や R5m 全件完了まで判別試験を止める。一方、試験の成功だけで側壁解像・他作動点・領域感度の確認を免除することもしない。

**呼び出し側の前提への異議:**

最大の問題は、**完了済み診断が未着手として残り、限定受理が全面完了として読まれていること**である。R7 を「限定した判別試験」と「MOO 本番」の二段に分ければ、優先順位を整理できる。

本体 plan §5.1 と 3D plan §5.1 は、次の内容へ書き直すことを推奨する。

| 項目 | 書き直す状態・依存関係 |
|---|---|
| **R6** | `(a)(b)(c)(e) 完了` を維持し、「残＝(a)(c)(e)…」という古い本文を整理。残件は数値許容・検出目標・作動点別 `C_M` 窓。前二者は R7a 前、実機制約値は MOO 本番前。 |
| **R7a：最優先** | 上記の 3D・単作動点・二水準×三反復。固定格子での判別能力のみを認定する。 |
| **R7b：R7a 後** | 3D runner と `L_sw` を最適化 driver に接続。対象作動点の評価成立、設計差の格子・領域感度、作動点別制約を確認してから MOO 再取得。 |
| **R4／R5** | 「現形状・m6_on の成立／側方境界の限定受理」と「未検証の方向・作動点・形状」を分割。2D 回帰を 3D 全作動点へ転用しない。 |
| **R5h** | 診断・現形状の力感度試験は完了。恒久対策は別件として保留。新形状で異常が増えれば再開。 |
| **R4e-次／R4c** | W2 比較済み、後流 station 対策済み等を履歴付きで閉じる。旧トポロジ修正の実装・試験への対応を明記し、有限厚・隅フィレットの R5q と分離する。 |
| **R5m／R5q** | R7a と並行・後続。ただし側壁解像不足のまま最終順位を認定しない。z 細分時は `sz` による形状変化を分離する。 |
| **完了行の整理** | R9/R9b・R10・R11 は完了表記と取り消し線を統一。R4b は完了を明示。R1–R3・R8 は履歴を保持。旧 #1 の再取得は R7b、#2 の剥離作動点追加は決着済み §8-6 へ対応付け、#2b の「最優先」は現行の起動 plan と整合させる。 |

**不足情報:**

- 作動点別の実機 `C_M` 許容窓と、その出典。
- 現行バイナリでの `L_sw` 感度、および新形状・他作動点への評価精度の持ち越し根拠。
- AWS 原本は今回再判定していない。数値は回収済み証拠と plan を照合したもので、`forge` は起動していない。

**ファイル変更なし。plan 未反映。** 上表は、依頼文の指定どおり呼び出し側が反映するための提案である。
