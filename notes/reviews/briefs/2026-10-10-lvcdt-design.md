# ブリーフ: 値 3・マスク 7 (全部入りの薄層 Jacobian) の破綻の切り分け — 次の判別 A/B (方向別 dt の 1/0) の設計

諮問の目的: ユーザは次の一手として「マスク 7 のまま方向別 dt の 1/0 だけを変える A/B」(codex 諮問の推奨) で進めることを選んだ。
そのうえで、ユーザの指示により上位 (diagnostician) に設計を諮る。AGENTS.md のエスカレーション条件 1 (§4・§6 の新規の設計) と 4 (承認済みの手順に無い run) に当たる。
ユーザの方針は「速度でなく筋のいい手法」。熱伝導も温度をライン内で結合する全部入り (マスク 7) を成り立たせたい。マスク 5 は回避策で、全部入りの代わりにはならない。

## 1. 観測事実 (これまでの連鎖。plan `plans/active/time_integration-line-viscous-jacobian-faceh.md`)

共通の条件:
- case/45 (M6 ノズル、冷却等温壁 300 K、TP 多成分、SST、node、FP64 のビルド)。
- 出発の場は `run_0183` の res_100000 (sha256 207d39f0…、親 plan で NOT CONVERGED と記録)。
- 値 3 (`lineViscCoupling: 3` = 薄層の D/K + 従来のスカラーの対角を全行に重ねる)、`implicitThermalJacobian` 5、`lineDtDirectional` 1、上限なし、cfl 4、緩和 0.7、sweep 5、ISP 0 (LHS は float)。最大 2000 step。

| 節 | 介入 (各腕 2 本) | 結果 | 判定 |
| --- | --- | --- | --- |
| §6.2 | 面エンタルピーを double で評価 (`FORGE_DIAG_FACE_H_DOUBLE`) | 4 本とも 122〜123 step で非有限 | 棄却 |
| §6.3〜§6.6 | 共通関数 `accumulate_thinlayer_visc_jacobian` の U-J の列ごとの照合 (host) | 多倍長の参照で PASS | 丸めの仮説を支持 |
| §6.7〜§6.8 | 製品の経路の照合 (監査用のビルド) | 係数の再現は PASS。壁際のライン面の β・κ が float の座標の差で double の係数と最大 10.8 % ずれる | FAIL (段 C)・T 保留 |
| §6.9〜§6.10 | LHS の座標の差を double の差 e に (段 ②/③ のバイナリ) | 4 本とも 199〜357 step で非有限。序盤の増え方は段 ③ が −1.8 % | 棄却 |
| §6.11〜§6.13 | マスク 7/5 (熱伝導の近傍 K の有無だけ、段 ③) | マスク 7: 348・280 step で非有限。マスク 5: 2 本とも 2000 step 有限・未収束 | 支持 (この条件・期間に限る) |

§6.13 の詳細 (run は AWS `~/forge-faceh-audit/case/45.isobutane_m6_d155/`。転送が終わると手元の `~/forge-evidence/2026-10-10-faceh/` にも置く):

| run | マスク | `rms_ro`: 0 / 20 / 50 / 100 / 200 / 300 / 500 / 1000 / 1999 |
| --- | --- | --- |
| `run_0570_lvc75_m7_a1` | 7 | 7.97e-6 / 5.72e-5 / 3.83e-3 / 1.38e-2 / 1.13e-2 / 1.20e-2 / — (348 で非有限) |
| `run_0571_lvc75_m5_b1` | 5 | 7.97e-6 / 5.74e-6 / 4.81e-6 / 5.93e-6 / 2.17e-4 / 3.22e-3 / 4.47e-3 / 2.07e-3 / 1.47e-3 |
| `run_0572_lvc75_m7_a2` | 7 | 同上で 280 で非有限 |
| `run_0573_lvc75_m5_b2` | 5 | 7.97e-6 / 5.74e-6 / 4.81e-6 / 5.90e-6 / 2.18e-4 / 2.74e-3 / 4.09e-3 / 2.02e-3 / 6.6e-6 |

- マスク 7 は step 16 で出発の 10 倍を超えた。マスク 5 は step 100 まで下がった後、step 170 前後から上がり、step 316 付近で最大になって下がった。
  マスク 5 の 2 本とも `check_convergence --segment` は NOT CONVERGED で、末尾 500 step の `rms_roK` は上昇していた。
- 場の変化が大きい位置は、どちらのマスクでも縮流部 (x ≈ −0.90〜−0.73) の冷却壁際 (wall_dist 1e-6〜3e-5) だった。
  マスク 7 の非有限は 2 本とも、隣り合う 11 本の壁法線のライン (列 41〜51、x −0.817〜−0.772) の全 1331 節点。sweep 5 回で両側 5 列に広がった形と数は合うが、確認していない。
- 事前のゲート (§6.12) で、マスク 7/5 の差は次のとおりだった。D と K の行 0〜3 はビット一致した。
  K の行 4 の変化は、熱伝導の近傍 K の式 `κ (γ_j/c_p,j)/ρ_j · [−(e_j − ½|u_j|²), −u_j, −v_j, −w_j, 1]` と 1200 面すべてで一致した。

## 2. 期待値と出典

- 共通関数の式と符号系: `solver_density_cuda/cuda_forge/block_dplur_jacobian_d.cuh` の `accumulate_thinlayer_visc_jacobian` (94〜175 行付近)。
  `D ΔQ_i − K ΔQ_j = rhs`。熱伝導の D は常に入り、マスクのビット 2 は熱伝導の近傍 K だけを切り替える。前処理行列としての近似で、実残差の厳密な微分ではない (コメントに明記)。
- 時間項: `timeIntegration_d.cu` 868 行付近で、D に `V/Δτ` (`v / dt_local`) を全行の対角として入れる。
- 方向別 dt: `procedures/solver-settings.md` の `lineDtDirectional` (75 行付近)。内部のライン面の音響の制約を局所 dt から外し、壁法線のラインでは Δτ が大きくなる。`lineDtDirectionalCap` で上限を付けられる (今回は上限なし)。
- codex 諮問 (`notes/reviews/2026-10-10-lvc75-result-diagnose.md`) の推奨: マスク 7 のまま `lineDtDirectional` 1/0 だけを変える A/B。
  - 介入のゲート: 残差・物性・拘束・K が不変で、`dt_local` が変わったこと。D の差は、拘束の前では時間項 `V/Δτ` の変更から丸めの範囲で再現し、拘束の行では上書きされること。
  - 第 1 仮説 (中): 方向別 dt の大きな刻みが、全部入りの更新で過大な補正を許している。
  - 反証条件: 方向別 dt なしも各 2 本とも非有限。

## 3. 再現条件

- 段 ③ の FP64 (`~/forge-fgeom3-fp64`、sha256 129de3f4…、元のセッションのバイナリを読むだけで使う)。
- 台本・ゲート・判定の雛形: `case/45.isobutane_m6_d155/lvc75.sh`・`lvc75_pregate.py`・`lvc75_judge.py` (commit 511cc614 以降)。
- run の作成: `cold_cfl.py prep <src> <run> --line dir|... --itj 5 --lvc 3 ...` (`--line` の値が方向別 dt の切替に当たるかは要確認)。

## 4. 実施済みの操作と結果

- 上の表のとおり。各節の VERDICT・ゲートの記録は `case/45.isobutane_m6_d155/_band_ab/cold_pair/` の `lvcfh_judge.json`・`lvcgeom_pregate.json`・`lvcgeom_judge.json`・`lvc75_pregate_r1.json`・`lvc75_pregate.json`・`lvc75_judge.json`。
- 諮問の記録: `notes/reviews/2026-10-10-lvc-faceh-result-diagnose.md`・`2026-10-10-uj-mp-result-diagnose.md`・`2026-10-10-lvcaudit-result-diagnose.md`・`2026-10-10-lvcgeom-result-diagnose.md`・`2026-10-10-lvc75-pregate-rhs-diagnose.md`・`2026-10-10-lvc75-result-diagnose.md`。
- 親 plan (`plans/accepted/time_integration-line-viscous-jacobian.md`) の観測 (別のバイナリ):
  - 値 2・マスク 7 は 29 step で非有限。
  - 値 3・マスク 15 (熱伝導の K の密度の列だけ外す) は 566 step で非有限。
  - 親 §6.16: 面エンタルピーを float で評価すると、エネルギーの行の実残差の方向微分が double の約 28 倍に見えた (局所、S0・p7・ライン 2183)。

## 5. 仮説 (呼び出し側のもの。確かめていない)

- **H1 (減衰の不足)**: 方向別 dt では、壁法線のラインの `V/Δτ` が小さく、熱伝導のライン内の結合 (D と近傍 K) がほぼ Newton の補正になる。
  前処理行列 (薄層・凍結物性) と実残差の熱流束の微分が少しでも違うと、滑らかなモードの補正が過大になって増える。
  マスク 5 は近傍 K が無いので、熱伝導の補正が `rhs/(V/Δτ + Σκ)` に縮み、減衰が残る。dt の A/B は、この「時間項の減衰が全部入りを救うか」を見る。
- **H2 (LHS と実残差の不整合)**: 冷却壁際の熱流束について、実残差の T の Q への依存 (TP の c_p(T)・組成・壁の弱形式・面の k_f の補間) が、薄層の K と食い違っている。H1 と併存しうる。
- **H3 (ライン外との釣り合い)**: ライン方向だけ熱を陰的に強く結合し、横方向 (ライン外のキー 5 のスカラー + 温度の項) は遅れるので、2 次元の熱のモードで釣り合いが崩れる。

## 6. 聞きたいこと

1. ユーザは案 1 (方向別 dt の 1/0) を選んだ。この A/B は、全部入りを成り立たせる方向で次に回す判別として妥当か。
   もっと情報の多い単因子の介入があるなら挙げてほしい (例: `lineDtDirectionalCap` で Δτ に上限を付ける、cfl を下げる)。ただしユーザの選択を覆すほどの理由がなければ案 1 で設計する。
2. 案 1 の設計:
   - `lineDtDirectional 0` で変わるもの (`dt_local` だけか、ライン外の節点・SST・化学種の dt も変わるか)。
   - 介入のゲートで確かめるべきこと。
   - 分岐と言える範囲。
   - 見る量 (壁際の `dt_local` の比、補正量など)。
   方向別 dt を外すと同じ step 数でも擬似時間の進み方が変わるが、2000 step の窓のままでよいか。
3. 結果の分岐ごとに、全部入りに向けた次の手 (H1〜H3 のどれに進むか) をどう決めるか。
4. この診断は faceh plan の中で続けるべきか、後継の active plan を作って移すべきか。faceh は result 段のレビューを受けて閉じる予定。

## 読んでよいファイル (リポジトリは `/home/sano/work/forge-faceh`、ブランチ `feature/faceh-audit-viscjac-close`)

- `plans/active/time_integration-line-viscous-jacobian-faceh.md` (全体、特に §1・§5.1・§6.9〜§6.13)
- `plans/accepted/time_integration-line-viscous-jacobian.md` (§4.1・§6.9・§6.16・§6.17)
- `solver_density_cuda/cuda_forge/block_dplur_jacobian_d.cuh`・`solver_density_cuda/cuda_forge/timeIntegration_d.cu`・`solver_density_cuda/cuda_forge/setDT_d.cu`
- `procedures/solver-settings.md` (`lineDtDirectional`・`lineViscCoupling`・`implicitThermalJacobian` の節)
- `case/45.isobutane_m6_d155/cold_cfl.py`・`lvc75.sh`・`lvc75_pregate.py`・`lvc75_judge.py`・`_band_ab/cold_pair/*.json`
- `notes/reviews/2026-10-10-lvc75-result-diagnose.md` ほか上の諮問の記録
- 手元の証拠 `~/forge-evidence/2026-10-10-faceh/` (run の CSV・HDF5。転送中のものは不完全なことがある)
