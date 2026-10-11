# codex 諮問 (diagnose): after-g9-next

- **brief**: [`notes/reviews/briefs/2026-10-11-after-g9-next.md`](../../notes/reviews/briefs/2026-10-11-after-g9-next.md)
- **plan**: [`plans/active/architecture-float-state-double-geometry.md`](../../plans/active/architecture-float-state-double-geometry.md)
- **date**: 2026-10-11
- **commit**: `183e973a` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 5.3 min, rc=0
- **結論**: **案 A の `pRef` 試験を先に行い、FP64 の同等性を確認したうえで、登録した温度差成分を抑えられるか判別する。**
- **extra**: `case/45.isobutane_m6_d155/rg_an.py`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **案 A の `pRef` 試験を先に行い、FP64 の同等性を確認したうえで、登録した温度差成分を抑えられるか判別する。**

第 1 仮説: `space.pRef: 5.499e6` による圧力項の組み方の変更だけで、2,000 step の登録成分を十分に抑えられる。確度: **低**
  
根拠:

- `pRef` は、計算済みの `p_tilde` から差し引かれ、運動量流束に入る。[内部面](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:618)、[境界面](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/convection/convectiveFlux_boundary_d.inc.cuh:197)、[hoop ソース](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/axisymmetricSource_d.cu:52)が同時に変わる。圧力項の積・残差への加算の精度を変える介入として合理的だが、効果は未測定。
- 保存済みの [rg_judge.json](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/cold_pair/rg_judge.json) の `a` から再計算すると、G9 は e = −0.012453、区間 [−0.023641, −0.001264]、F は e = 0.912544、区間 [0.574291, 1.250798]。登録規則の「結果 B」は再現できる。ただし、これは**末尾の保存量量子化だけでは不十分**という証拠であり、圧力項が真因という証拠ではない。

反証条件: FP64 の同等性・陽性対照が成立しているのに、下記の結果 B となること。この特定の `pRef` 変更の十分性を退ける。

第 2 仮説: 圧力などの格納・EOS・再構成・面値計算で生じる丸めが重要で、後段の `pRef` 減算では十分に抑えられない。未確認。
  
第 3 仮説: その他の演算・幾何の丸め、または精度別の初期化差との相互作用が必要。未確認。G9 の陰性でも除外されない。

判別 A/B: **変更は `space.pRef: 0 → 5499000` の一つだけ。**

- 新規腕は Cp＝FP64 ×3、Fp＝float ×3。各 2,000 step 固定。予約済みの `case/45.isobutane_m6_d155/run_0700〜0705` を使用する。
- 対照 C は同 case の `run_0595_rq_c1`〜`run_0597_rq_c3`、F は `run_0621_rq_f1`〜`run_0623_rq_f3`。各精度内でバイナリ・入力・環境・初期化後の9保存量が一致し、設定差が `space.pRef` だけと確認できる場合に限り再利用する。
- 評価域・参照 `v_s`・時点 1,000/1,500/2,000・指標 `a` は変更しない。まず Cp を評価し、次の前提が通ってから Fp に進む。

**FP64 の同等性は二段で確認する。**

1. **作用素の確認。** 同じ状態でのゲージ変更による運動量残差の差は、丸めを除けば  
   `ΔR_i = pRef × (Σ_f W_if − A_i)`  
   となる。`W_if` は実際に使う外向きの r 重み付き面ベクトル、`A_i` は実効ソースの閉性補正。境界半割面と拘束行を含めて照合する。`hoopAreaFromClosure: 1` では x 成分の補正も有効になるため、省かない。[実装根拠](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/axisymmetricSource_d.cu:109)

   FP64 で `‖ΣW − A‖ / (Σ‖W‖ + ‖A‖) ≤ 1e−12` を暫定ゲートとする。分母ゼロは絶対値で扱う。これは幾何・ゲージ整合の確認であり、軌跡比較で代用しない。

2. **有限期間の同等性。** 次を結果を見る前に固定する。
   - `mean(a_Cp) − mean(a_C)` の95% Welch区間全体が **[−0.002, +0.002]** に入ること。既存 C の `a` の最大差は 0.0001742、その10倍は 0.001742。
   - 保存時点 500/1,000/1,500/2,000 の Cp–C の全組合せで、温度差が **RMS ≤ 0.05 K、最大絶対差 ≤ 1 K**。登録上流域と全域を別々に確認する。射影と直交する差や、下流の悪化を見落とさないためである。

これらは診断用の提案許容幅で、実証済みの物理精度ではない。外れた場合は「前提不成立」とし、幅を広げて救済しない。

**主判定の式。**

`e_Fp = mean(a_Fp) − mean(a_Cp)`  
`d = [mean(a_F) − mean(a_C)] − e_Fp`

新しく判定する `e_Fp` と `d` の2比較について、各区間の t 分位を `1 − 0.05/4` とする。`d` は4腕の平均の対比として計算し、`v_g = s_g²/n_g` とすると、

`SE_d² = Σ_g v_g`  
`ν_d = (Σ_g v_g)² / Σ_g[v_g²/(n_g−1)]`

独立標本は run ごとの `a`。節点や保存時点を反復に数えない。

- **結果 A:** `CI(e_Fp) ⊂ [−0.5, +0.5]` かつ `L_d ≥ 0.5`。登録参照の尺度で0.5以上削減し、残りも±0.5以内になった。第1仮説を支持する。
- **結果 B:** `L_Fp > 0.5`、`U_Fp < −0.5`、または `U_d < 0.5`。第1仮説の「この変更だけで上記の抑制目標を満たす」を退ける。逆符号の大差は抑制と呼ばない。
- それ以外は判別不能。延長・窓変更はしない。

既存 F の陽性判定は過去の確認として明記する。新規2比較の補正と、探索全体の95%保証を混同しない。4比較へ機械的に補正を増やすと、既存 F の下限は **0.430** となり、0.5の陽性条件を満たさなくなる。

やらない方がよいこと: 今回同時に Gp・Gp9 を実装・投入すること。`pRef` の陽性を「残差の加算が真因」、陰性を「圧力の丸めは無関係」と読むこと。既定値の変更や本番採用まで進めること。

呼び出し側の前提への異議:

- **Major — 「Cp と C の `a` が近い＝作用素が同じ」は採らない。** `a` は一方向への射影である。[rg_an.py:86](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/rg_an.py:86)  
  対案は上記の閉性確認と、射影以外の場の比較。結論も登録した有限期間・成分に限る。

- **Major — 案 B の挿入位置では、列挙した原始量・物性の格納精度を一括して試せない。** `gasProperties` の後に境界条件が入り、さらに `turbulent_viscosity` が実行されるため、そこで丸めた境界値・`vis_turb` は後から書き換わる。[main.cpp:1913](/home/sano/work/forge-integ-1005/solver_density_cuda/main.cpp:1913)、[main.cpp:1961](/home/sano/work/forge-integ-1005/solver_density_cuda/main.cpp:1961)  
  対案は、将来この診断を設計するときに、各量の生成・上書き・消費位置を定義し、消費直前の実効性を検査すること。

- **Major — 「G9 対 Gp で原始量の格納精度が要るか分かる」は却下。** G9 は保存量の量子化あり、Gp はそれを外して原始量の量子化を加えるため、二因子が同時に変わる。現行の介入対象は保存量の配列である。[quantizeStateDiag.cpp:173](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/quantizeStateDiag.cpp:173)  
  対案は一因子だけ異なる比較に直すこと。陰性でも初期化・境界・途中の上書きとの相互作用は残る。

不足情報: 対象 run の HDF5・残差 CSV・実効設定はローカルに無く、生データからの再検証はできなかった。確認できたのはコードと保存済み判定 JSON。C・E・F の収束記録は `OVERALL: CHECK FAILURES ABOVE`、C・E・F・G9 の準定常記録は `OVERALL: NOT ALL STEADY`。G9 の `rg_judge.json` は `records` が空で、同ファイルでは収束ツールの記録を確認できない。次の判定では実効設定・入力・座標対応と各 VERDICT を保存すること。

ファイル変更・forge の起動は行っていない。**plan 未反映**。上記は呼び出し側が事前登録するための提案である。
