# codex レビュー: discretization-node-face-weight-midpoint (result)

- **plan**: [`plans/active/discretization-node-face-weight-midpoint.md`](../../plans/active/discretization-node-face-weight-midpoint.md)
- **stage**: `result` (diff base `18916ed9`)
- **date**: 2026-09-22
- **commit**: `8ebe72f8` (feature/sern-design)
- **codex**: effort `high`, 3.5 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M4/m1
- **focus**: 実装 diff は solver_density_cuda/cuda_forge/calcStructualVariables_d.cu の nodeMode 1 行 (fx=0.5 を node 内部双対面の固定スキームに) と、診断の追加 (ifaceRconv/ifaceRpre/ifaceRro, wi_eheat/wi_ework, iface_q_eff から e_w*R_rho を引く)。見てほしい点: (1) §5.1 #2/#2b/#3 の回帰と回転試験の結果が、#4 恒久実装と accepted への移動を裏付けるか。#7 (格子系列での次数測定) を未実施のまま accepted にしてよいか、active に残すべきか、推奨を 1 つに絞ってほしい。(2) §4.3 の表と §5.1 #5 の数字の主張に過大なものが無いか。読んでよいのは本 plan、git diff 18916ed9 -- solver_density_cuda (grep で該当箇所のみ)、各 run の CONVERGENCE_VERDICT.txt。禁止: *.log, residual_history.csv, res_*.h5, *.vtu, plans/README.md, case/*/README.md の全文を読まない (grep は可)。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**NO-GO**  
`fx=0.5` の恒久実装は設計どおりですが、定常回帰は plan 自身の合格条件を満たしていません。  
#7 の次数測定を実施し、回転試験・翼の数値主張の裏付けを整えるまで、`active` に残すことを推奨します。

指定の差分、保存済み `CONVERGENCE_VERDICT.txt`、許可された README の該当行を確認しました。ファイル変更・禁止された生データの読み込みは行っていません。

実装については、[calcStructualVariables_d.cu:59](/home/sano/work/forge/solver_density_cuda/cuda_forge/calcStructualVariables_d.cu:59) の `nodeMode == 1 && ip < nNormalPlanes` と同ファイルの wrapper により、node 内部面だけが `0.5` になります。cell・境界半割面への適用拡大や、新しい絶対閾値ガードはありません。[conjugateWall.cpp:228](/home/sano/work/forge/solver_density_cuda/conjugateWall.cpp:228) の補正も `roe/ro × Rro` を引いており、係数・符号・単位は plan の半離散式と整合しています。ただし、これだけで周期 seam・軸・全拡散経路の数値的妥当性まで確認できたとはいえません。

1. **Major — 定常回帰の必須 `PASS` がなく、#2 の完了扱いが不適切。**

   根拠：[plan §6:119](/home/sano/work/forge/plans/active/discretization-node-face-weight-midpoint.md:119) は平板・スラブに `check_convergence.py` の `PASS` を要求しています。実際の保存判定は以下です。

   | run（リポジトリ相対パス） | 保存済み判定 | 最終 `rms_roe` |
   |---|---|---:|
   | `case/48.flat_plate_cooled_m4/run_0020_fx_ctrl/` | `NOT CONVERGED (stalled/plateau)` | 0.214 |
   | `case/48.flat_plate_cooled_m4/run_0021_fx_half/` | 同上 | 0.375 |
   | `case/52.conjugate_slab/run_0006_ctrl_newbin/` | 同上 | 2.23e−5 |
   | `case/52.conjugate_slab/run_0007_fxhalf/` | 同上 | 2.17e−5 |

   各数値の出典は当該 run の `CONVERGENCE_VERDICT.txt:7`。平板では `rms_roK` も `still converging` です。[case/48 README:41](/home/sano/work/forge/case/48.flat_plate_cooled_m4/README.md:41) の「収束場からの継続で低下桁が出ないだけ」は、要求された判定の代わりになりません。これは新実装が悪化した証明ではなく、**回帰合格が未証明**という指摘です。

   **対案：** #2 を未完了に戻す。適格な参照 run を使った `--from-floor` で判定を保存するか、収束可能な定常検証で所定の `PASS` を取得する。比較差が小さいことだけで完了にしない。

2. **Major — #2 の追加回帰は、#7 の次数測定を代替していない。**

   根拠：[case/09 README:198](/home/sano/work/forge/case/09.Taylor-Green/README.md:198) は、`run_0169`–`run_0172` が一様直交格子で「幾何 `fx` も 0.5」と明記しています。[case/44 README:759](/home/sano/work/forge/case/44.vitiated_air_wt/README.md:759) も、`run_0505`–`run_0508` は Euler で変更対象の粘性・拡散経路をほぼ通らないとしています。最大比 1.39／1.88 の `check_field_regress` 合格は無影響確認として有用ですが、**重みが実際に変わる格子での精度検証ではありません**。

   [plan:108](/home/sano/work/forge/plans/active/discretization-node-face-weight-midpoint.md:108) の #7 は採用済みレビュー対応として未実施です。辺中点での値補間の性質から、双対面上の流束積分や解の収束次数は保証できません。

   **対案：** 本 plan を `active` に保持し、#7 の伸長・曲面高 AR・回転格子で三格子検証を実施する。温度・速度・壁熱流束の次数を記録し、`p≥1.8` を満たした範囲だけに二次精度の主張を限定する。なお、非定常の `case/09` に定常残差の `PASS` がないこと自体は問題にしていません。

3. **Major — 回転試験の「残りは float32 丸め」という原因断定は裏付け不足。**

   根拠：[plan:103](/home/sano/work/forge/plans/active/discretization-node-face-weight-midpoint.md:103) の数値から、差の縮小率は `0.085/0.024 ≈ 3.54`、`0.254/0.062 ≈ 4.10` で、算術は正しいです。しかし、`case/48.flat_plate_cooled_m4/run_0022_rot30_fx_ctrl/` と `run_0023_rot30_fx_half/` はともに `NOT CONVERGED`。各 `CONVERGENCE_VERDICT.txt:9` の `rms_roOmega` は、それぞれ **1.12e4／1.99e4、低下 0.7／0.6 桁、stalled** です。

   座標の丸め幅と第一層厚の比を示すだけでは、残る差を丸めに帰属できません。反復誤差や他の離散化誤差が分離されていません。

   **対案：** 現状は「4000 step 時点で回転前後の差が約 3.5／4.1 倍縮小」と限定する。収束・対象量の定常性を確認した比較と、座標精度または原点移動の対照実験で、丸めへの帰属を検証する。

4. **Major — §4.3 の定量表を支える準定常判定が、別の run の判定になっている。**

   根拠：§4.3 の対象は `case/53.c3x_vane_cht/run_0122_fxhalf/` と `run_0123_fxgeom_ctrl/`。両者の保存判定は `NOT CONVERGED` です。一方、[plan:104](/home/sano/work/forge/plans/active/discretization-node-face-weight-midpoint.md:104) が `ALL STEADY` の根拠として挙げるのは、生産設定の `run_0125/0126` と Mark II の `run_0024/0025` です。壁温設定も熱流束定義も異なるため、この判定を §4.3 に転用できません。[case/53 README:743](/home/sano/work/forge/case/53.c3x_vane_cht/README.md:743) にも、旧 A/B の対象量について準定常判定は示されていません。

   表の算術自体は整合しています。`1.22→0.46` は **約62.3%減**、最大の領域別偏差変化は **2.7 pt** で、記載された閾値内です。しかし、許可された証拠だけでは、その値が定常化した量であることを独立確認できません。

   **対案：** `run_0122/0123` の対象量・判定窓・`check_quasisteady` の実出力を対応づけて提示する。それまでは表を「未収束試行の暫定値」とし、「25 節点周期のうねりは消えた」という確定表現を弱める。

5. **Minor — 計画一覧が実装済みの状態に追随していない。**

   根拠：[plans/README.md:31](/home/sano/work/forge/plans/README.md:31) は今も `draft`・環境変数 A/B の説明です。一方、[plan:105](/home/sano/work/forge/plans/active/discretization-node-face-weight-midpoint.md:105) は恒久実装済みです。

   **対案：** 一覧を「恒久実装済み、回帰合格・次数測定・result レビュー対応待ち」に同期する。今回の閲覧制限では `methods/`・`procedures/`・`methods/index.md` の整合性までは確認しておらず、文書整合全体を合格とはしていません。

§5.1 #5 の更新 run は、挙げられた全対象で最終 step **59999** を確認できました。ただし、全対象の残差判定は `NOT CONVERGED` です。README は C3X 層流の遷移後を `DRIFTING`、Mark II 層流の一部を `DRIFTING`／`OSCILLATING` と区別しており、この限定は維持すべきです。また、[Mark II README:75](/home/sano/work/forge/case/54.markii_vane_cht/README.md:75) の層流域偏差は **75.9→79.7%、+3.8 pt**。§6(a) は C3X 対象なので直ちに基準違反とはしませんが、「翼全般で変化3 pt以内」と一般化できません。

**推奨は、恒久実装を維持したまま本 plan を `active` に残すことです。** 優先順は、定常回帰の合格証拠取得、#7 の次数測定、回転・翼の主張と判定記録の対応付け、文書同期です。現状の #2/#2b/#3 は、`accepted` 移動を裏付けるには不足しています。本回答はレビュー提案であり、ファイル変更禁止のため **plan 未反映**です。

指摘数: Critical 0 / Major 4 / Minor 1
