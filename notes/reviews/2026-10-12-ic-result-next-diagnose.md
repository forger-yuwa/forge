# codex 諮問 (diagnose): ic-result-next

- **brief**: [`notes/reviews/briefs/2026-10-12-ic-result-next.md`](../../notes/reviews/briefs/2026-10-12-ic-result-next.md)
- **plan**: [`plans/active/architecture-float-state-double-geometry.md`](../../plans/active/architecture-float-state-double-geometry.md)
- **date**: 2026-10-12
- **commit**: `92784219` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 3.5 min, rc=0
- **結論**: **M1 を採用し、G0c の継続場から、同一 FP64 バイナリで「保存量9本の末尾量子化 OFF／ON」を各2本・3,000 step 比較する。**
- **extra**: `case/45.isobutane_m6_d155/ric_an.py`, `case/45.isobutane_m6_d155/ric.sh`, `case/45.isobutane_m6_d155/_band_ab/cold_pair/ric_judge.json`, `solver_density_cuda/cuda_forge/quantizeStateDiag.cpp`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **M1 を採用し、G0c の継続場から、同一 FP64 バイナリで「保存量9本の末尾量子化 OFF／ON」を各2本・3,000 step 比較する。**

第 1 仮説: float では step 内の演算・格納の精度差が差を供給し、それが境界層の縁で増幅される。末尾の保存量9本の量子化だけでは、今回の大きな差も再現しない。確度: **中**。
  根拠: [plan §6.43](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:1024) では、G0・2,000 step の温度差 RMS が全量量子化で0.28〜0.31 K、float で3.7〜4.5 Kだった。ただし、**G0c・継続場への外挿は未検証**であり、残差演算や右辺への局在化まではできない。
  反証条件: 同じ起点・同じ FP64 バイナリで、量子化 ON が再実行差を十分超え、後述する float 応答の半分以上を両評価窓で再現すること。

第 2 仮説: G0c の継続場では末尾量子化だけでも十分な駆動になり、従来の G0・2,000 step の陰性結果はこの条件には当てはまらない。確度: 低。

第 3 仮説: 精度依存の初期化差が主要な種になっている。継続場でも `res_0` の ρE 最大相対差は4.90e-5で、初期化の寄与は除外されていない。確度: 低。

判別 A/B: **変更は `FORGE_DIAG_QUANTIZE_STATE` の OFF／ON だけ**とし、G0c・B0・上限50、共通入力 `_ic/cont_q32.h5`、各2本・3,000 step・250 step ごとに出力する。

- C＝FP64 OFF、Q＝同じ FP64 の全9量 ON。`FORGE_DIAG_QUANTIZE_STATE_ONLY` は未設定にする。実装上の介入は全 CV の保存量への丸めであり、float 演算全体の再現ではない（[quantizeStateDiag.cpp:171](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/quantizeStateDiag.cpp:171)、[同:275](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/quantizeStateDiag.cpp:275)）。
- **陽性対照 F を同じソース版の float で2本添える。** `fgeom11` の C を取り直すだけでは、旧 `fgeom7` の F との比較まで統制できない。合計6本だが、試験する介入は量子化の1点だけ。
- C/Q の `res_0` の9保存量・座標の一致、全出力の有限性、量子化対象と履歴の自己点検を前提にする。step 1・2,750・3,000 を自己点検し、化学種履歴は次 step が存在する点で確認する。
- 評価域を x/r_t ∈ [−9, −6]・壁から40〜90節点、評価窓を W1＝2,000〜2,500、W2＝2,500〜3,000 に固定する。各 run の温度と同時点の C 2本の平均との差について、空間・窓内の RMS を E とする。N は同じ定義による C の再実行差。
- 新試験用の判定を事前登録する。両窓で F の最小 E ≥ 10N を要求し、差の絶対値も K 単位で残す。
  - **A：** 両窓で Q の最小 E − N ≥ 0.5 × F の最大 E。第1仮説の「末尾量子化だけでは不十分」を退ける。この介入の十分性までで、float の実際の主因・同一モードの証明とはしない。
  - **B：** 両窓で Q の最大 E ＋ N ≤ 0.1 × F の最小 E。第2仮説を、この起点・3,000 step に限って退ける。第1・第3仮説は残る。
  - 中間、反復間の不一致、陽性対照不成立は判別不能。閾値や窓を事後変更しない。

やらない方がよいこと:

- **M2 を次に優先しない。** 上限1の効果は G0 ですでに観測されている。G0c の短期陰性を追加しても、増幅の抑制と進行の遅延を分けにくい。上限は平均流・SST の擬似 dt をともに変える（[solver-settings.md:76](/home/sano/work/forge-integ-1005/procedures/solver-settings.md:76)）。
- **M3 の「最初に差が見えた量」を原因としない。** 出力頻度を上げても、相互に結合した更新の因果順序は確定しない。
- **M4 を K1／K2 の判別としない。** 一度の乱数摂動は連続的な丸めと異なり、空間形状・変数間の相関にも依存する。陰性でも閾値や二状態の存在を否定できない。

呼び出し側の前提への異議:

- **採用：§6.53 の正式判定は「判定不能」のまま。** JSON から再計算すると、前提(a)の float 側は成立し、FP64 側だけが不成立。前提(b)の比は0.00679、(c)は1.063だった。「継続場の float 2本でも低い初期残差から再上昇した」という観測は記載できるが、正式な結果Aには変更しない。8本の記録された収束判定はすべて **`OVERALL: CHECK FAILURES ABOVE`**（[ric_judge.json:859](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/cold_pair/ric_judge.json:859)）。
- **Major／却下：K1 の「線形増幅なら約6,000 step 遅れて FP64 も離れるため、有限振幅閾値が必須」という推論。** 実効的な擾乱の比を機械精度比1e9と置く根拠がない。FP64 にも単精度経路があり（[plan:1009](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:1009)）、精度間で反復写像も同一とは限らない。対案は、**有限期間の増幅を作業仮説として残し、閾値を必須条件にしない**こと。
- **Major／要再検証：K4 の「同じモード」「飽和」と K2 の「二つの状態」。** 継続場の F1−D1 は3,000→6,000 stepで最大差194→278 K、RMS 3.05→8.12 K。後者は約2.66倍で、頭打ちの根拠にならない（[ric_judge.json:768](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/cold_pair/ric_judge.json:768)、[同:852](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/cold_pair/ric_judge.json:852)）。最大差の位置も動くため、似た成長率だけでは同一モードを同定できない。対案は「近い領域で差の増幅を観測」に限定し、飽和・二状態の確定を保留すること。対象量の準定常 VERDICT も提示されていない。
- **Major／要再検証：K3 の「系統的にずれた右辺」という局在化。** 末尾量子化の陰性結果は、右辺・左辺・原始量の格納・初期化を区別しない。対案は第1仮説のように範囲を広く保ち、今回も「十分／不十分」だけを判定すること。
- **Minor／修正：判定器の失敗時の扱い。** `res_0` が欠けると一致検査を飛ばす（[ric_an.py:56](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/ric_an.py:56)）。また、前提(a)不成立の文言が「旧起点で再上昇が再現しない」となり、実際に失敗した FP64 条件を隠す（[同:75](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/ric_an.py:75)）。次の判定器では必要ファイル欠落を判定不能にし、失敗した条件を個別に出す。

不足情報: 対象 `run_0725〜0736` の CSV/HDF5 はこの作業ツリーに存在せず、直接再集計できたのは保存済み JSON の値と判定式まで。同精度間の局所増幅はブリーフ・plan の記録に依存し、抽出の再検証は未実施。`fgeom11` の配置・実バイナリと、陽性対照に用いるソース版の確認が必要。

**plan 未反映。** 読み取り専用の診断依頼に従い、ファイル変更・forge 起動は行っていない。
