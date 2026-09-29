# codex 諮問 (diagnose): cary-origin-ab-result

- **brief**: [`notes/reviews/briefs/2026-09-26-cary-origin-ab-result.md`](../../notes/reviews/briefs/2026-09-26-cary-origin-ab-result.md)
- **date**: 2026-09-26
- **commit**: `f1e56cac` (feature/gap-heating-precision)
- **codex**: effort `high`, 3.4 min, rc=0
- **結論**: **遷移モデルや比較座標を変える前に、同一 checkpoint から「現行精度／全域 FP64」の継続 A/B を行い、残差プラトーの精度依存性を判別する。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **遷移モデルや比較座標を変える前に、同一 checkpoint から「現行精度／全域 FP64」の継続 A/B を行い、残差プラトーの精度依存性を判別する。**

第 1 仮説: **流束評価・乱流変数更新に残る float32 演算が、残差プラトーを制限している。** 確度: **中（実 run 未確認）**
  根拠: [gen_runs.py:48](/home/sano/work/forge/case/59.flatplate_cary_m6/gen_runs.py:48) は既に `qAccumulatorFP64: 1`。しかし FP64 正本の対象は流れの保存量 5 本だけで、`roK`・`roOmega` は対象外である（[main.cpp:2192](/home/sano/work/forge/solver_density_cuda/main.cpp:2192)）。通常の演算型も float のまま（[flowFormat.hpp:6](/home/sano/work/forge/solver_density_cuda/flowFormat.hpp:6)）。したがって「アキュムレータ有効だから丸めは除外済み」とはいえない。一方、提示された残差値だけでは丸め床と断定できない。
  反証条件: 同一条件の FP64 継続でも、停滞していた残差列が現行精度と同程度の床・振れ幅に留まる。ただし短期間で変化しないだけなら、遅い緩和との区別は未了。

第 2 仮説: **精度に依存しない反復停滞・緩和不足。** 確度: **中**。`ALL STEADY` は観測した St 系列の判定であり、全保存量の離散方程式が解けた証拠ではない。

第 3 仮説: **実験と完全乱流 CFD の発達履歴の違いが R の下流増加に寄与する。** 確度: **低〜中、未確認**。遷移後の overshoot 緩和という特定の説明までは、今回の A/B では検証できない。

判別 A/B: **変更点は演算精度だけ。**

- 共通の出発場は `case/59.flatplate_cary_m6/run_0002_tw02_re027_ext/res_40000.h5`。同一メッシュ・保存量の index コピーで、新しい 2 run に分岐する。
- A＝現行 float32＋`qAccumulatorFP64: 1`、B＝同じソース版の全域 FP64。CFL・緩和・SST・BC・比較点は固定し、メッシュも再生成しない。
- 各 **20,000 step**。全残差列と、既存の A/B 座標における St 34 列を取得する。壁出力は例えば 100 step 間隔とし、既存の 2,000 step 間隔より細かく振動を確認する。
- **判定を先に固定する**：末尾 4,000 step を比較し、B だけが停滞列で A より **1 桁以上低い水準**を持続すれば、精度依存を支持する。両腕とも停滞し、各列の末尾中央値が **±20 %以内**なら、この試験で精度変更が停滞を解消する仮説を退け、第 2 仮説を優先する。列ごとに割れる・まだ低下中なら判別保留。
- この閾値は**新しい診断用の提案**であり、収束合格基準ではない。両腕の `check_convergence` と St 系列の `check_quasisteady` の VERDICT は別に記録する。

やらない方がよいこと:

- **Major — `rounding_floor.py` だけで収束扱いにすること。** 同ツールは内部面の算術平均から作る目安で、再構成・リミッタを含まず、評価対象も質量とエネルギーだけ（[rounding_floor.py:5](/home/sano/work/forge/case/57.transition_flat_plate/tools/rounding_floor.py:5)、[同:25](/home/sano/work/forge/case/57.transition_flat_plate/tools/rounding_floor.py:25)）。対案は上記の精度 A/B。ツールは残差の位置と規模を見る補助に限定する。
- **Major — 今の結果に合わせて比較点を下流へ絞る、仮想原点を fitting すること。** 事前登録は閾値をまたぐ場合を判定不能としている（[acceptance.json:39](/home/sano/work/forge/case/59.flatplate_cary_m6/acceptance.json:39)）。対案は現行点を維持して予備結果として保存すること。`lm2009` の導入も今は保留する。現行資料の検証実績は T3A・C3X で、`dilatationCorrection` との組合せも未検証と明記されている（[recommended-settings.md:146](/home/sano/work/forge/procedures/recommended-settings.md:146)）。

呼び出し側の前提への異議:

- **Major — 「面積 0.4 %、前縁だけ」は未確認。** `check_wall_resolution.py` は **点数割合**を「面積割合」と表示している（[check_wall_resolution.py:368](/home/sano/work/forge/solver_density_cuda/tools/check_wall_resolution.py:368)）。また `index 0` は壁出力配列の添字であり、それだけでは位置の証明にならない。**対案:** 超過点の実座標・境界面積重み・未評価点数を出す。ゲートを満たす処置は、その分布に基づく壁法線細分化と再評価である。3.94 を固定した単純見積りでは y₁ は 1.5/3.94 ≈ **0.38 µm以下**だが、再計算後の壁応力で再判定が必要。満たせなければ **FAIL のまま予備扱い**とし、前縁除外を事後導入しない。

- **Major — 延長単独の低下桁数では、全体の収束ゲート未達を確定できない。** 判定器は入力系列のピークから低下桁数を測る（[check_convergence.py:182](/home/sano/work/forge/solver_density_cuda/tools/check_convergence.py:182)）。**対案:** 実効設定と継続を確認して、`run_0001` の ramp0→main と `run_0002` を接続した区間を再判定する。未収束ならその判定を維持する。未収束 run を `--from-floor` の参照にする方法は、実装も拒否する（[同:353](/home/sano/work/forge/solver_density_cuda/tools/check_convergence.py:353)）。

- **Major — D は遷移や実測偏差の原因を識別する指標ではない。** 定義から  
  **D = St_CFD(x_exp−22.22 cm)/St_CFD(x_exp) − 1**  
  となり、実測 St は約分される（[cary_compare.py:61](/home/sano/work/forge/case/59.flatplate_cary_m6/tools/cary_compare.py:61)）。**対案:** 記録は「**提示された未収束・壁解像 FAIL の場では、抽出座標変更により予測比が＋9.1〜25.8 %変化した。事前登録上は判別保留**」とする。この限定付きなら質問 1 の記述は採用できる。「比較座標が誤差を支配する」「B が物理的に正しい」はまだ書けない。

- **Minor — drift 0.0 % は厳密な不変を意味しない。** 表示は小数 1 桁である（[check_quasisteady.py:291](/home/sano/work/forge/solver_density_cuda/tools/check_quasisteady.py:291)）。**対案:** 丸め前の drift・fluctuation・判定窓を併記する。

不足情報: ローカルに対象 run はなく、提示された `NOT CONVERGED`／`ALL STEADY`／壁解像 `FAIL` は独立再検証できていない。必要なのは、実行時 YAML・`stage_manifest.json`・バイナリのソース版と精度設定、接続区間の全残差判定、St の判定窓と fluctuation、壁超過点の座標・面積重みである。**ファイル変更・forge 起動はしていない。plan 未反映。**
