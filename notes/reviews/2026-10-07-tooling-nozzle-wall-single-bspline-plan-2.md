# codex レビュー: tooling-nozzle-wall-single-bspline (plan)

- **plan**: [`plans/active/tooling-nozzle-wall-single-bspline.md`](../../plans/active/tooling-nozzle-wall-single-bspline.md)
- **stage**: `plan`
- **date**: 2026-10-07
- **commit**: `04c130d7` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 3.6 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M1/m3
- **extra**: `notes/reviews/2026-10-07-wall-single-bspline-repr-diagnose.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
案 A の数学的構成と、W3 の「ソルバ入力がビット同一なら追加 CFD 不要」は支持します。  
実装前に、保存壁の復元仕様・検証用データ・ケース依存定数・STEP の誤差予算を明確にしてください。

コード基準は HEAD `04c130d7`。`plans/README.md`・関連 accepted plan・仕様・設定・検証手順を確認しました。過去の A14 は低自由度の近似壁を評価した計画であり、今回の表現統一とは重複しません。

`case/45.isobutane_m6_d155/run_0147_ns_mono_final/` の保存データを別ツリーから読み、MOC を再計算せず、保存された設計スプラインと δ_r から独立に試算しました。

| 確認項目 | 今回の実測 |
|---|---:|
| 再構成壁と保存 `wall_physical.csv` の最大差 | 3.33×10⁻¹⁶ m |
| 最小二乗行列のランク／係数数 | 1747／1747 |
| 行列の条件数 | 35.009 |
| 半径／r′／r″ の最大誤差 | 1.17×10⁻¹³／4.57×10⁻¹²／1.41×10⁻⁹ |
| 生成座標を直接 float32 化した不一致 | 0／582,000 成分 |
| 同じ座標を `.10g` で文字列化した不一致 | **3 成分** |

最後の結果は、[メッシュ出力時の直列化](/home/sano/work/forge-integ-1005/design/forge_design/meshing/mesh2d.py:237)まで含めて W3 を実施する必要性を裏付けます。文字列差をソルバ入力差とは判定していません。変換後 HDF5 の比較は未実施です。

1. **Major — 保存壁の復元仕様が不完全で、旧 run の報告と上流の差分図を壊す余地があります。**

   **根拠:** [plan:65](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-wall-single-bspline.md:65) は物理壁・設計壁の双方を保存スプラインから読む方針ですが、既存 `run_0147` には新しい物理壁の係数がありません。現行報告は [CSV を読む処理](/home/sano/work/forge-integ-1005/design/forge_design/report/nozzle_report.py:485)を、[報告生成から無条件に呼びます](/home/sano/work/forge-integ-1005/design/forge_design/report/nozzle_report.py:546)。

   また、`wall_fit.spline` は**下流の S のみ**です。[実際の設計壁](/home/sano/work/forge-integ-1005/design/forge_design/geometry/wall_axismach.py:85)は直管・上流 Hermite・S を切り替えます。保存 S を上流へ外挿すると、今回の実測で x＝−6 の半径は **10.0086 r_t**、正しい設計壁は **4.8675 r_t**でした。`BSpline` は既定で端区間を外挿するため、例外にもなりません。[SciPy 公式仕様](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.BSpline.html)

   **対案:** §4.2 に復元契約を追加してください。新形式では版・表現種別・有効域を保存し、設計壁には上流 Hermite と直管の復元情報も含める。新形式の欠損はエラー、旧形式の run は明示した旧経路で報告する。W4 に「旧 run」「新形式」「係数欠損」「上流を含む差分図」の検査を登録してください。

2. **Minor — W0・W4 の検証入力と実施方法が未確定です。**

   **根拠:** [W0](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-wall-single-bspline.md:113) は過去 run の座標との完全一致を要求しますが、参照できた `run_0147` には `nozzle.h5` がありません。また [W4](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-wall-single-bspline.md:121) の報告生成には、[入力 HDF5 と `res_*.h5` が必要](/home/sano/work/forge-integ-1005/design/forge_design/report/nozzle_report.py:60)です。`prepare_ns` までの CFD 0 step 検証だけでは、これらの結果ファイルは揃いません。

   **対案:** 実装前の固定入力・環境で `legacy` の基準成果物を確保し、W0 は変更前後を比較する。過去 run との再現検査は別項目にする。W4 は幾何報告部分の単体検査と、入力同一性を確認した既存結果を使う統合検査に分け、必要ファイル・取得元・ハッシュを §5.1 に登録してください。

3. **Minor — 対応範囲は一般の joint 壁ですが、構成規則は case/45 の定数で書かれています。**

   **根拠:** [plan:48](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-wall-single-bspline.md:48) の接合位置と [plan:71](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-wall-single-bspline.md:71) の `76.6539 mm` はケース固有です。一方、現行実装は [任意の `pw_ramp` と既定ランプ](/home/sano/work/forge-integ-1005/design/forge_design/geometry/wall_axismach.py:374)を扱い、[問題の `r_throat`](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner_axismach.py:1006)で実寸化します。

   **対案:** 接合位置は `−L_U`・実効ランプ両端・0、定義域は壁属性、mm 換算は `1000 × scale_m` から導くと明記してください。case/45 の数値は例として分離し、異なるスケールと既定ランプを使う試験を追加してください。

4. **Minor — W5 の「量ごとに誤差を足す」は、曲率について式が不足しています。**

   **根拠:** [W1](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-wall-single-bspline.md:114) の対象は r′・r″、[W5](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-wall-single-bspline.md:124) の対象は接線角・幾何曲率です。スロート以外では κ＝r″／(1＋r′²)³ᐟ² なので、r″ の許容差を単位換算しただけでは曲率の表現誤差予算になりません。

   **対案:** 同じ物理 x で、元壁と保存スプラインの接線角・曲率を直接比較し、その表現誤差に STEP 転送誤差を加える方式を事前登録してください。接線角の 10⁻⁹ rad 判定には、丸めに弱い `acos(dot)` ではなく `atan2` を使うことも明記すると確実です。

**推奨は、案 A を維持し、上記 1→2→3→4 の順で plan を補ってから実装することです。** 現在の精度と条件数から、ノット削減や別の近似方式への変更は不要です。メッシュ品質検査は既存 `prepare_ns` にあり、W3 が成立する範囲では、新たな IC・段階起動・node/cell・周期・block-DPLUR の流れ計算を追加する理由もありません。

なお、参照 run の保存判定は main 区間・最終 step 59999で **`NOT CONVERGED (stalled/plateau)`** です。本レビューは幾何表現の成立性を評価したもので、流れの収束や定常性を認定していません。ファイル変更・forge 起動はなく、提言は **plan 未反映**です。

指摘数: Critical 0 / Major 1 / Minor 3
