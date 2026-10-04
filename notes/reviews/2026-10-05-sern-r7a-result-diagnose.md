# codex 諮問 (diagnose): sern-r7a-result

- **brief**: [`notes/reviews/briefs/2026-10-05-sern-r7a-result.md`](../../notes/reviews/briefs/2026-10-05-sern-r7a-result.md)
- **plan**: [`plans/active/tooling-nozzle-sern-chain.md`](../../plans/active/tooling-nozzle-sern-chain.md)
- **date**: 2026-10-05
- **commit**: `bea7f70f` (feature/sern-design)
- **codex**: effort `high`, 5.0 min, rc=0
- **結論**: **R7a の登録判定は保持し、次は長側壁 B の g3 格子で、既存双子節点の内外を保存する初期写像への変更だけを試す。**
- **extra**: `notes/investigations/2026-10-05-sern-r7a/R7A_VERDICT.txt`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：

| 重大度 | 判断・採否 | 根拠と対案 |
|---|---|---|
| **Major** | **分岐 A は、登録条件内の判定として限定付き採用** | 集約値から再計算すると、Δ`C_T_with_shear`＝0.0020859、Δ−E＝0.00207901。算術は正しい。ただし、E は反復・窓内変動の経験的指標で、初期場依存性や格子・領域誤差を含まない。**「g3・m6_on・今回の初期化で登録基準を満たした」まで**とする。[判定原本:8](/home/sano/work/forge-sern-design/notes/investigations/2026-10-05-sern-r7a/R7A_VERDICT.txt:8) |
| **Major** | **「最近傍距離 0 なので、初期場の懸念は新設803節点だけ」を却下** | `interp_field.py` は座標だけで最近傍を選び、既存の同座標双子節点も内外を区別しない。B の反復は同じ入力の複製なので、共通する写像誤差は群内差 r に現れない。対案は**同じ B 格子で初期写像だけを変える対照**。[写像:212](/home/sano/work/forge-sern-design/solver_density_cuda/tools/interp_field.py:212)、[反復生成:67](/home/sano/work/forge-sern-design/case/46.sern_design/r7a_setup.py:67) |
| **Major** | **「ΔC_M＝−0.213 は許容0.05の4倍だから問題」を却下** | 0.05 は**同一形状の格子・領域感度**の許容で、異なる形状間の設計差や実機トリム限界ではない。大きな形状応答として記録し、B の数値感度と作動点別制約を別々に確認する。[3D plan:1970](/home/sano/work/forge-sern-design/plans/active/tooling-nozzle-sern-3d.md:1970) |
| **Major** | **現行 C_M をそのまま実機トリム判定に使う案を却下** | 実装の `C_M` はノズル面の**圧力モーメント**。摩擦モーメントと機体面は含まない。実機の許容窓とは、面集合・基準点・正規化・摩擦の扱いを合わせる必要がある。[力積分:259](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern3d.py:259)、[集計:309](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern3d.py:309) |
| **Major** | **R7a だけで MOO 本番へ進む案を却下** | 既存 G は短側壁形状の結果であり、長側壁形状への持ち越し根拠がない。設計差の格子依存性は **G_B−G_A** で決まる。初期写像の確認後も、設計差の格子・領域感度、対象作動点の評価成立、作動点別 C_M 窓が必要。[R7b:438](/home/sano/work/forge-sern-design/plans/active/tooling-nozzle-sern-chain.md:438) |

結論: **R7a の登録判定は保持し、次は長側壁 B の g3 格子で、既存双子節点の内外を保存する初期写像への変更だけを試す。**

第 1 仮説: **B の初期写像による過渡は十分減衰しており、写像を正しても `C_T_with_shear` は5×10⁻⁵以内で再現する。** 確度: **中・未確認**。  
根拠: 集約報告では延長後の B の窓内振幅は約3.6×10⁻⁶、前窓差は約6.5×10⁻⁷まで低下している。一方、これは同じ初期写像を共有する三反復である。[判定原本:8](/home/sano/work/forge-sern-design/notes/investigations/2026-10-05-sern-r7a/R7A_VERDICT.txt:8)  
反証条件: 下記対照が必要条件を満たし、変動幅を差し引いても写像間の差が5×10⁻⁵を超えること。

第 2 仮説: **同じ初期写像に由来する偏りが残り、三反復の小さな r では検出できていない。** 確度: **低・最終値への影響は未確認**。既存の `r4d_common_restart.py` も、この双子節点の誤写像を対処対象として明記している。[既存対策:2](/home/sano/work/forge-sern-design/case/46.sern_design/r4d_common_restart.py:2)

判別 A/B:

- **対照A**：既存の長側壁系列  
  `case/46.sern_design/run_1037_r7a_lsw10_1/` → `run_1040_r7a_lsw10_1_c40k/` と、その2・3反復。
- **試験B**：同じ長側壁格子・設定・バイナリで新規1本。起点も `run_1017_ff4f_g3_2p50/res_20000.h5` とし、**初期写像だけ変更**する。共通の既存節点は座標と境界所属・接続で内外を識別し、保存量を直接コピーする。新しく分裂した節点は従来と同じ初期値、新メッシュの `wall_dist` は保持する。
- **投入前確認**：共通領域の保存量一致と、変更節点数・位置・量別最大差を記録する。`r4d_common_restart.py` は幅追加用で一対一対応を要求するため、側壁延長へ無変更で流用しない。
- **長さ**：既存 B と同じ40000 step、500 step間隔。最後の20000 stepを判定区間とし、前後10000 stepの窓条件、4量の `STEADY`、有限値・床・品質ゲートを満たすこと。残差 VERDICT も併記する。未達なら判定不能。
- **見る量**：4係数の末尾平均・振幅・前窓差。既存B群の N と新規runの `max(a,d,10⁻⁶)` の和を比較幅 U とする。これは信頼区間ではない。
- **結果A**：全量で `|平均差|＋U` が、推力2量・C_Lは5×10⁻⁵、C_Mは5×10⁻⁴以下 → この写像対照では大きな初期場依存性を支持しない。短側壁群との比較でも分岐 A が維持されるか確認する。
- **結果B**：いずれかで `|平均差|−U` が同閾値を超える → 初期化に依存しないという仮説を棄却し、R7a の認定範囲を訂正する。境界をまたぐ場合は判定不能。

やらない方がよいこと: **同じ入力の4本目を追加して共通の初期場依存性を調べたことにすること**、検出余裕が小さいという理由だけで事前判定を撤回すること、逆に E を全数値誤差の上限と扱うこと。

呼び出し側の前提への異議:

- 検出目標を上回る余裕は **7.90×10⁻⁵**。問題は「4%だから不合格」ではなく、初期化・格子・領域への依存性をその余裕と比較していない点である。
- 報告平均からの算術では `|ΔC_M/ΔC_L| ≈ 23.75`。基準点は x_ref＝−20H なので、モーメント差の大きさは長い腕による増幅と整合する。基準点を x_ref＝0 に換算した差は約−0.03368になる。**これは異常の証明でも、圧力分布の正しさの証明でもない**。[基準点:42](/home/sano/work/forge-sern-design/case/46.sern_design/problem_3d_prod_m6on_wallres.yaml:42)
- 物理 station 化は連続最適化への接続前に必要だが、実位置を固定した今回の対照の前提にはしない。

不足情報: **R7a の生履歴・個別 VERDICT は本 checkout にない。** 確認できた原本は `R7A_VERDICT.txt` の集約値と「前提: OK」であり、残差・準定常性・有限値・品質を独立再判定していない。各runの `force_history.csv`、`residual_history.csv`、個別判定、実効config、初期写像照合記録が必要。run索引は [case README:373](/home/sano/work/forge-sern-design/case/46.sern_design/README.md:373)。

**ファイル変更・forge起動なし。plan 未反映。** 提案の反映先は `tooling-nozzle-sern-chain.md` §5.1 R7a／R7b。
