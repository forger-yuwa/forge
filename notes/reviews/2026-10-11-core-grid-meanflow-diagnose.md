# codex 諮問 (diagnose): core-grid-meanflow

- **brief**: [`notes/reviews/briefs/2026-10-11-core-grid-meanflow.md`](../../notes/reviews/briefs/2026-10-11-core-grid-meanflow.md)
- **plan**: [`plans/active/tooling-nozzle-core-grid.md`](../../plans/active/tooling-nozzle-core-grid.md)
- **date**: 2026-10-11
- **commit**: `c72ad4cf` (feature/nozzle-core-grid)
- **codex**: effort `high`, 3.5 min, rc=0
- **結論**: **M-a を採用し、B0 の `time.deltaT.implicitThermalJacobian` だけを `5 → 0` にした 20,000 step の試験を、既存の `run_0619` の同区間と比較する。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **M-a を採用し、B0 の `time.deltaT.implicitThermalJacobian` だけを `5 → 0` にした 20,000 step の試験を、既存の `run_0619` の同区間と比較する。**

第 1 仮説: **追加した温度の対角項と、壁側で拡大した擬似 dt の組合せが、境界層のエントロピー成分を含む連成反復の減衰を悪くしている。** 確度: 中。機構は未確定。

根拠:

- B0 の熱伝導項は、`diag_block[4][0..4]` に加わる**自節点の微分**である。この分岐には近傍の温度増分との結合がない。`5 → 0` では従来の全 5 行のスカラー対角が残り、追加した温度項を外せる。等温壁の拘束を変えるビット 2 は、5・0 の双方で無効である。[timeIntegration_d.cu:1026](/home/sano/work/forge-coregrid/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1026)、[同:1045](/home/sano/work/forge-coregrid/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1045)、[同:1163](/home/sano/work/forge-coregrid/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1163)
- `case/45.isobutane_m6_d155/run_0626_cg_gref_freeze/` では、SST 保存量の更新を止めても縮流部の残差比は ρ が 1.00・1.12、ρE が 0.94・1.00。記録された `check_convergence` は **NOT CONVERGED**。平均流側の反復を次に調べる根拠になる。[plan §4.39:776](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:776)
- 上限 1 では、残差集中域の dt 比は中央値 0.977 とほぼ不変なのに残差が大幅に減った。したがって、**集中域自身の dt だけで説明する仮説は不十分**で、壁側からの結合を含める必要がある。[plan §4.36:746](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:746)

反証条件: ITJ の無効化が実際に効き、対照が停滞している条件で、後述の両窓・両成分が変更側／対照 ≥ 2/3、残差の集中場所も維持されること。「追加温度項を外せば、この停滞が大きく緩和する」という第 1 仮説を退ける。熱拡散一般の関与までは否定しない。

第 2・第 3 仮説:

- 第 2: ITJ への依存が小さく、ライン外の結合を 5 回の内部反復で解き残すことが支配的。確度: 中。ライン外では前 sweep の増分を読み、各 outer step で増分をゼロから解く実装は確認できるが、**5 回での線形解き残しは未測定**。[timeIntegration_d.cu:899](/home/sano/work/forge-coregrid/solver_density_cuda/cuda_forge/timeIntegration_d.cu:899)、[main.cpp:2052](/home/sano/work/forge-coregrid/solver_density_cuda/main.cpp:2052)
- 第 3: ρk・ρω の更新を止めても残る、密度・ひずみを介した μt と熱伝導係数の再計算が関与する。確度: 低。凍結試験でも μt は 1.8〜4.2 % 動いており、未除外。ただし変化率だけでは原因とは言えない。

判別 A/B: **M-a の 1 組だけを行う。**

- 起点: `case/45.isobutane_m6_d155/run_0608_cg_gref_cont_k1/res_60000.h5`。新しい run に `restart_field.py` で移し、9 保存量のビット一致を確認する。
- 対照: `case/45.isobutane_m6_d155/run_0619_cg_gref_b0_ctl/` の 0〜20,000 step。変更側は ITJ だけ 0。ライン・方向別 dt・上限 50・実効 CFL・緩和・内部反復 5・SST 更新・バイナリを固定する。`FORGE_FREEZE_TURB` は無効。
- 指標: 縮流部 x/r_t ∈ [−8, −2] の `q = √(ΣV(R/V)² / ΣV)`。窓は **10k〜15k、15k〜20k**、各窓の出力平均を使う。ρ・ρE それぞれの変更側／対照を取る。対照の停滞条件も、両成分で後窓／前窓 ≥ 1/2 と事前登録する。
- **結果 A:** 両窓・両成分で ≤ 1/3 → ITJ 追加項への強い依存を支持。「ITJ にほぼ依存しない解き残しだけが支配する」という第 2 仮説の限定形を退ける。ただし、熱対角単独と内部反復との相互作用は分離できない。
- **結果 B:** 両窓・両成分で ≥ 2/3、集中場所も維持 → 第 1 仮説の上記予測を退け、第 2 仮説を優先する。第 2 仮説の確定ではない。
- 中間・成分間の不一致・対照の停滞不成立は判定不能。非有限は発散として別扱いにし、自動延長しない。全保存量の `check_convergence` も残し、SST が悪化した場合は「平均流のみ改善、連成系は未解決」とする。θ_r・Q_w の参照値への到達は、この試験の合否に使わない。

やらない方がよいこと: M-b・M-c・低マッハ前処理を同時に足すこと。特に M-b は、20 回でも線形反復が収束した保証がなければ、陰性結果から内部反復不足を除外できない。さらに、反復回数を増やしても、B0 の行列に存在しない近傍熱伝導項は追加されない。

呼び出し側の前提への異議:

- **[Major] 「流れ方向に細かいからライン外結合が強い」は、まだ機構の説明ではない。** 方向別 dt 自体もライン外面の CFL から決まるため、絶対的な係数の増大だけでは判断できない。[setDT_d.cu:135](/home/sano/work/forge-coregrid/solver_density_cuda/cuda_forge/setDT_d.cu:135)  
  対案: 今回は ITJ への依存だけを判別し、機構確定には時間項・対角・近傍項の相対強度と線形解き残しを求める。「追加対角が大きすぎるだけ」という説明も、dt を縮めると改善する理由まで示せていない。
- **[Major] 他セッションの「値 3・マスク 7」を ITJ の比較と混同しない。** コードには `lineViscCoupling == 3` と `FORGE_LVC_TERMS` の別経路があり、近傍熱伝導 K の追加は B0 の ITJ 対角項とは異なる操作である。[timeIntegration_d.cu:978](/home/sano/work/forge-coregrid/solver_density_cuda/cuda_forge/timeIntegration_d.cu:978)、[同:1613](/home/sano/work/forge-coregrid/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1613)  
  対案: 当該 run の実効キー名を明示し、今回は補助情報に留める。
- **[Major] F1 で除外した範囲を「乱流との連成全体」に広げない。** 除外を支持するのは「SST 保存量の更新を止めれば停滞が崩れる」という予測であり、μt の再計算は残っている。[plan §4.39:780](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:780)  
  対案: F1 の結論は維持し、再計算される輸送係数の経路だけを未除外として記録する。

不足情報: 内部 sweep ごとの線形残差、熱対角と時間項の相対強度、実行バイナリと今回読んだソースの対応。これらは機構の確定に必要だが、M-a の選択を妨げない。run 数値・VERDICT は指定された plan の記録を引用し、禁止された生データは読んでいない。ファイル変更・forge 実行は行っておらず、**plan 未反映**。呼び出し側で §4.39 に続く試験登録へ反映すること。
