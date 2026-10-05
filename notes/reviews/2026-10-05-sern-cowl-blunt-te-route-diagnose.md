# codex 諮問 (diagnose): sern-cowl-blunt-te-route

- **brief**: [`notes/reviews/briefs/2026-10-05-sern-cowl-blunt-te-route.md`](../../notes/reviews/briefs/2026-10-05-sern-cowl-blunt-te-route.md)
- **plan**: [`plans/active/tooling-nozzle-sern-chain.md`](../../plans/active/tooling-nozzle-sern-chain.md)
- **date**: 2026-10-05
- **commit**: `4f895c95` (feature/sern-design)
- **codex**: effort `xhigh`, 8.4 min, rc=0
- **結論**: **T2 を採用し、最初の一手は既存接続模型の双対閉性 FAIL を、同一メッシュの幾何演算精度 A/B で切り分けることとする。**
- **extra**: `plans/active/tooling-sern-mesh-blocking.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：

| 重大度 | 対象・判定 | 根拠と対案 |
|---|---|---|
| **Major** | **T1：今回の経路として却下** | 実現不能ではないが、「ギャップ層追加」で閉じない。現行は後縁・自由側端の共有節点と側端テーパに依存する（[mesh_sern3d.py:271](/home/sano/work/forge-sern-design/design/forge_design/meshing/mesh_sern3d.py:271)、同ファイル:279）。接合面・側端壁・後流・タグを再設計するなら、既存の `SW`・`CW1/CW2` ブロックを再利用する **T2** に集約する。 |
| **Major** | **T3：記載された構成は却下** | `sz=0` の側端までギャップのヘキサ列を延ばすと辺が潰れる。今回、既存の [cornerJ:413](/home/sano/work/forge-sern-design/case/46.sern_design/cad/hex_junction_model.py:413) を用いたメモリ上の最小例で、側端を閉じたヘキサは **8 頂点中 4 頂点の Jacobian が 0**。プリズム化や専用遷移なら別設計だが、現行の共有処理を残すだけでは成立しない。 |
| **Major** | **「接続模型は全ゲート PASS」：却下** | 保存記録 [run_0446…/MESH_QUALITY.txt:7](/home/sano/work/forge-sern-design/case/46.sern_design/run_0446_3d_gs050_flag1/MESH_QUALITY.txt:7) は品質 `PASS` に続き、閉性 **max 7.41e−5、閾値超過 866 CV、`VERDICT: FAIL`**。245 万節点模型も閉性 max 1.7e−4・6,410 CV 超過と記録されている。**B1b を先に閉じる**。現ドライバの「閉性 FAIL でも続行」は受入経路では禁止する（[run_junction_model.py:72](/home/sano/work/forge-sern-design/case/46.sern_design/cad/run_junction_model.py:72)）。 |
| **Major** | **「B1d は旧診断のまま未解決」：要再検証・台帳訂正** | 後継の [accepted plan:834](/home/sano/work/forge-sern-design/plans/accepted/convection-slau-wall-normal-chi.md:834) には `case/46.sern_design/run_0439_3d_junction_outflow_steady/`、通算 66,000 step、要求系列 **`OVERALL: ALL STEADY`・床数 0** がある。ただし最後の 36,000 step の残差は **`NOT CONVERGED (stalled/plateau)`**。`slauWallNormalChi` も現在は auto＝実効 1（[recommended-settings.md:77](/home/sano/work/forge-sern-design/procedures/recommended-settings.md:77)）。既済みの排出診断を繰り返さず、この限定成立を B1d に戻し、新形状の成立確認へ進む。 |
| **Major** | **旧模型の寸法・細分則の流用：却下** | 端面間隔 0.25 mm は板厚 2 mm の模型用。後縁厚を生産板厚相当の **0.5 mm** にすると `t_te/5=0.1 mm` を超える。また [hex_junction_model.py:547](/home/sano/work/forge-sern-design/case/46.sern_design/cad/hex_junction_model.py:547) の `--scale` は、明示した **`H1_END` を細分しない**。端面解像を独立に指定・検査し、細分列でも確実に変える。 |
| **Major** | **帳簿・旧基準列の持ち越し：却下** | [forces3d:285](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern3d.py:285) は `cowl_base`・`cowl_side`・`sidewall_end` を集計しない。入口項も矩形面積の式（同:311）。**新端面をノズル力へ加え、入口実面積・流束を整合させ、新形状で基準列を取得する**。 |

結論: **T2 を採用し、最初の一手は既存接続模型の双対閉性 FAIL を、同一メッシュの幾何演算精度 A/B で切り分けることとする。**

T2 は即時の生産受入ではない。しかし T1 も接合部の再設計と受入試験を要し、短期の迂回路とは認められない。T2 なら既存の接続設計・生成器・検査を再利用でき、後の全ヘキサ移行との二重投資を避けられる。R7b 再開までの日数は、SERN 全体の節点数・変換メモリが未計測なので断定できない。

設計判断として、次を推奨する。

- **`t_te` は物理形状入力**とする。実機値が未定なら、初回候補は既存板厚と同じ `t_te/H=0.005`、H＝0.1 m なら **0.5 mm**。これは暫定モデル値であり、床が消えるまで厚くする調整値にはしない。厚さの定義を y 方向差か法線距離か明記し、保持する内面輪郭を固定して外側へ厚みを付ける。板厚分布・テーパ長・側壁厚・フィレットは格子数から決めない。初版でフィレットまで同時に変更しない。
- **後流格子は厚さ方向・流れ方向の両方でベースを分解する。** ベース直後の少なくとも `t_te` の範囲は Δx≤t_te/5、間隔比≤1.2 を満たす。それとは別に、端面法線の第一内部点距離と局所 y₁⁺ を検査する。`t_te/5` は壁解像の保証ではない。端面解像の緩和は、新形状での感度を根拠にする。
- **最小の予備受入は、新形状の 3 作動点×新 g3/g4。** メッシュは float32 後の全頂点 Jacobian 正、双対体積正、全 CV 閉性≤1e−5、境界タグ漏れ・重複 0、skew≤0.90、AR は既定範囲を必須とする。方式受入の 3 解像度試験は維持する。
- **CFD 判定は本段と同じ実効設定の 20,000 step**、前後各 10,000 step で行う。窓条件未達時のみ事前登録した追加 20,000 step を一度許す。4 係数の `STEADY`、窓平均差が `C_T/C_T_with_shear/C_L`≤5e−5、`C_M`≤5e−4、全残差の VERDICT と判定区間、有限性・化学種・各床ゲートを保存する。m10_on は末尾 10,000 step の全保存場で温度床到達 0 とし、**新ベース・側端・周囲の内部点への冷点移動も不合格**とする。
- **m6_on・m4_off の旧→新差は形状変更の影響として測る。** 同じ基準値・帳簿定義で、Dq＝|末尾平均差|＋両側の変動・ドリフト幅を報告する。「影響が小さい」の暫定閾値は `C_T/C_T_with_shear/C_L` 0.002、`C_M` 0.05。超過は新形状の即失格ではなく、従来結果を持ち越せないという判定である。**新形状自身の格子・領域誤差とは分ける。** 新 g3/g4 は取り直し、R7b の `L_sw` 二水準の設計差もその後に再評価する。

第 1 仮説: **接続模型の閉性 FAIL は、双対幾何の丸め誤差が主因であり、有限厚接合部の位相そのものの不成立ではない。** 確度: **中**

根拠: primal 品質 `PASS` と双対閉性 `FAIL` が併存する上記保存記録、および座標を `stof` で読む [gmshReader.hpp:504](/home/sano/work/forge-sern-design/solver_density_cuda/mesh/gmshReader.hpp:504)。ただし、これだけで演算精度を真因確定とはしない。

反証条件: 同じ座標・接続で双対幾何の中間演算だけを倍精度化しても、同じ CV 群に閾値超過が残ること。

第 2・第 3 仮説: **有限厚後縁により m10_on の近壁冷却が緩和する可能性はあるが、未確認。確度: 低。** [LEDGER100:5](/home/sano/work/forge-sern-design/notes/investigations/2026-10-05-sern-r7a/R5H_M10_LEDGER100.txt:5) が示すのは、介入後に冷点が移り、517199 が call 92 で床へ到達したことまでである。有限厚化の有効性を測った記録ではない。第 3 仮説は置かない。

判別 A/B: **既存の小さい接続模型を固定し、双対幾何の中間演算精度だけを変える。**

- A＝現行経路、B＝同じ float32 節点座標・同じ接続から、面ベクトル・重心・体積を float64 で計算し、最終保存形式は A と同じ。
- **CFD は 0 step、変換・検査を各 1 回。** 全 CV 閉性の最大値・閾値超過数・位置、双対体積、共有面の反対称性を見る。
- **A の FAIL を再現し B が全 CV≤1e−5** → 第 1 仮説を支持し、B1b を閉じる根拠にする。
- **B でも同じ超過が残る** → 演算精度だけで説明する仮説を棄却し、その CV の面集合・境界半割面を点検する。A が再現しなければ、保存記録との版・入力差を先に解消する。

やらない方がよいこと: `sz` を残した潰れヘキサ、閉性 FAIL のままの生産 CFD、新端面を帳簿外にした力比較、旧 g3/g4 の無条件流用、m10_on の 100 step 昇温だけでの受入。B1e の一般的な修正も今回の前提へ追加しない。

呼び出し側の前提への異議: **「有限厚にする」は形状方針として採用するが、「有限厚なら床が消える」は未検証。** またブリーフは、接続模型の閉性不合格を落とし、後継 plan の限定成立を B1d に戻していない。この二点は反対方向の誤差なので、両方訂正する必要がある。

不足情報: 実機の後縁厚・側壁厚・板厚分布、および新 SERN 全体の規模見積もり。対象 `run_0439`・`run_1055`・`run_1065/1066` の生 HDF5・残差系列は本 checkout に無く、保存診断記録とコードを照合した範囲の判断である。

**plan 未反映**（依頼どおりファイル無変更・`forge` 未実行）。呼び出し側の反映先は `tooling-sern-mesh-blocking.md` §4・§5.1・§6 と、`tooling-nozzle-sern-chain.md` §5.1 R7b 4-2。
