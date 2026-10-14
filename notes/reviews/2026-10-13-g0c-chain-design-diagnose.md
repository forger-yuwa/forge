# codex 諮問 (diagnose): g0c-chain-design

- **brief**: [`notes/reviews/briefs/2026-10-13-g0c-chain-design.md`](../../notes/reviews/briefs/2026-10-13-g0c-chain-design.md)
- **plan**: [`plans/active/tooling-nozzle-isothermal-wall-chain.md`](../../plans/active/tooling-nozzle-isothermal-wall-chain.md)
- **date**: 2026-10-13
- **commit**: `a8f66843` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 3.7 min, rc=0
- **結論**: **メッシャ内の `x_respace` を採用し、上記の再現性・段分離を先に確かめたうえで、固定した `run_0179/res_20000` から段階起動した G0c・FP64・B0 の検証を1本だけ事前登録して行う。**
- **extra**: `design/forge_design/meshing/mesh2d.py`, `case/45.isobutane_m6_d155/make_g0c.py`, `case/45.isobutane_m6_d155/cold_cfl.py`, `case/45.isobutane_m6_d155/cold_pair.py`, `case/45.isobutane_m6_d155/problem_d155_ns_prod_coldmesh_tw300.yaml`, `plans/active/architecture-float-state-double-geometry.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（**現案は修正後に採用**。Critical の指摘はありません）。

| 重大度 | 対象・採否 | 根拠と必要な修正 |
|---|---|---|
| **Major** | メッシャ内の `x_respace`：**採用。ただし端点の定義を修正** | 旧実装 `96197821:case/45.isobutane_m6_d155/core_grid_mesh.py:458–465` は、区間端を最寄り station に固定し、ブレンド内端は **−7.5、−2.5 のまま**にする。読み取り専用の再計算では区間端が −7.9904141995、−2.0011935537、旧方式は **4467 点**。一方、固定した区間端から幅 0.5 を取る解釈では **4468 点**になる。`blend` の内端を `x0 + blend`、`x1 − blend` と明記し、局所積分の 20,001 点、丸め、端点の上書きも旧方式に合わせる。 |
| **Major** | 本段だけ B0、S1・S2 は従来どおり：**採用。ただし現行の設定生成では実現しない** | [runner_axismach.py:1585](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner_axismach.py:1585) は本段設定から CFL・次数・内部反復数だけを変更する。関数を読み取り専用で実行すると、S1 にも `lineImplicit: 1`、`lineDtDirectional: 1`、上限 50、ITJ 5 が残った。**従来の前段設定を別に保持し、本段へだけ B0 を適用する**設計にする。S1・S2・本段それぞれの実効設定を検査する。 |
| **Major** | `evaluate.main_line` の5項目で B0 を再現：**要再検証** | B0 の系譜にはリミッタ基準値の固定もある。[hp7.sh:21](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/hp7.sh:21)、[cold_cfl.py:129](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/cold_cfl.py:129)。チェーンの設定にはこの配管がない。`main_line` という明示的な opt-in はよいが、**5項目が一致したことを全設定一致の代わりにしない**。`run_0717` の設定実物を取得し、リミッタ基準値・出力量・BC・物性ファイル・実バイナリ・環境変数まで照合する。既存の `cfl_main`／`implicit_relax` と競合した場合の扱いも決める。 |
| **Major** | 見張りを本段だけに適用：**採用。ただし起動境界を実装仕様にする** | 前段は同じディレクトリに `res_3000.h5` を出し、次段前に削除する。[runner_axismach.py:1566](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner_axismach.py:1566)。一方、[m9_watch.py:136](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9_watch.py:136) は最初の出力から 2,500 step 間隔を要求する。起動時から監視すると前段出力の混入で `DATA_ERROR` になり得る。**本段設定への切替・前段出力の処理後に見張りを開始する明示的な通知点**を設ける。見張りによる正常な停止と、予期しないプロセス終了も区別する。 |
| **Major** | `run_0179/res_20000`＋段階起動：**採用。ただし入力を固定** | [runner_axismach.py:1427](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner_axismach.py:1427) は指定ディレクトリの**最新**出力を選び、補間後には ω フロアも適用する。したがって `ic_from` の run 名だけでは開始場が固定されない。既存の [cold_pair.py:154](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/cold_pair.py:154) と同様に `res_20000.h5` を検査し、その SHA、補間・ω フロア適用後の入力、S2 終了後の本段入力を記録する。 |
| **Major** | M1〜M7 と停止規則による判定：**採用。ただし運用上の比較に限定** | [float plan:1493](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:1493) は収束合格を前提から外し、停止時の出力比較と定義している。この判断は維持する。対案は、停止規則・M1〜M7・収束・準定常を独立して記録し、合格名を「**登録した開始手順で停止規則に到達し、T2 両参照との比較が許容内**」とすること。 |
| **Minor** | 一般化・バイト一致・AR 緩和：**採用範囲を限定** | 密度を [mesh2d.py:113](/home/sano/work/forge-integ-1005/design/forge_design/meshing/mesh2d.py:113) と同じ式にするのは妥当。ただし非自明なガウス細分との組合せは旧 G0c の検証範囲外。バイト一致は固定入力の回帰条件として採用し、キーなしの出力不変も確認する。複数区間を許すなら重なり・順序・`ni_ref` の意味を定義する。AR ≤10,000 は今回の条件つき例外に限り、各生成メッシュで品質と高 AR セルの位置・直交性を再確認する。 |

結論: **メッシャ内の `x_respace` を採用し、上記の再現性・段分離を先に確かめたうえで、固定した `run_0179/res_20000` から段階起動した G0c・FP64・B0 の検証を1本だけ事前登録して行う。**

第 1 仮説: **G0c・FP64・B0 の停止時出力の互換性は、登録 Q32 起点から、`run_0179/res_20000 → 補間 → S1 → S2` という開始手順へ拡張できる。** 確度: **中**

根拠: [float plan §6.58・§6.61](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:1490) と保存済みの `sp_judge.json` は、登録 Q32 起点での比較要件充足を記録している。ただし、新しい開始手順の結果はない。**生産系譜が共通であることは、今回の開始手順が合格する証拠ではない。**

反証条件: 入力・設定・監視の事前条件を満たした試験で、本段 200,000 step までに停止規則へ到達しない、発散する、または T2a・T2b のいずれかに対して M1〜M7 が許容を超えること。

第 2・第 3 仮説: 今回は追加しない。失敗前から物理的な真因を決める根拠がない。

判別 A/B: **変える因子は、本段へ渡す開始場の生成手順だけ。**

- **対照**: `case/45.isobutane_m6_d155/run_0804_sp_g0c/` を再利用する。再投入しない。
- **試験側**: SHA を固定した `run_0179_ns_n012_N2_ext/res_20000.h5` と `delta_r_initial.csv` を使い、指定 YAML・`offset="radial"` で準備する。S1・S2 は従来設定、本段は対照と同じ G0c・FP64・B0 とする。`.msh` だけでなく、変換後の幾何・接続・壁距離も照合し、開始場以外の差を残さない。
- **予算**: **本段だけで最大 200,000 step**。S1 3,000＋S2 3,000 は別枠で、最大 206,000 step。現行 `full` は前段分を本段から差し引かない。[runner_axismach.py:1592](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner_axismach.py:1592)。準備・段階起動を含む総時間と本段時間を分ける。
- **停止**: 本段だけに `--phase line --budget 200000 --consec 2 --interval 2500 --min-points 9 --keep-tail 9`。停止条件は変更しない。準定常判定を理由に事後延長しない。
- **比較**: 停止出力を T2a・T2b **それぞれ**と、§6.58 の M1〜M7 全項目・同じ式・同じ許容で比較する。`run_0804` との近さだけでは合格にしない。
- **段のゲート**: S1・S2 は正常終了、必要保存量の存在・有限性・ρ>0、P・T の妥当性、全残差履歴の非有限を確認する。現行 `stage_gate` は最終場の必要保存量・ρを確認する機構で、収束認定ではない。[runner_axismach.py:1520](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner_axismach.py:1520)。前段に収束合格まで要求してレシピを変えることは、今回の試験へ追加しない。
- **判定区間**: `stage_manifest.py <run> --segments` で最後の区間が本段だけであることを確認し、`check_convergence.py <run> --segment` を実行する。S1・S2 の残差を本段の低下桁数へ混ぜない。準定常判定は本段末尾 20,000 step の9出力を使い、§6.11 の閾値を維持する。

**A＝到達し、T2 両参照との全比較が許容内**なら、この開始手順への適用拡大を支持する。**B＝未到達・発散・許容超過**なら拡大を退ける。入力不一致、出力欠損、監視の誤作動は判定不能とし、数値手法の失敗に数えない。

やらない方がよいこと:

- 今回に `ic_from` なしの準一次元開始を追加すること。今回は既存 NS ドナーを使う経路を検証する。新規設計をドナーなしで開始できるという主張には、別の検証が必要。
- 生成後の `.msh`／HDF5 へ別工程で格子変更を掛けること。station 配置をメッシャに統合し、壁法線層・接続・変換を同じ生成経路で扱う。
- 合格後に、壁温変更、δ* ループの全 pass、別形状、float ビルドへ一括適用すること。δ* pass は壁も変わるため、今回の固定壁の試験では保証できない。[deltastar_loop.py:369](/home/sano/work/forge-integ-1005/design/forge_design/feedback/deltastar_loop.py:369)。
- 起点の違う対照との時間比を高速化率とすること。

呼び出し側の前提への異議: **「停止時の比較要件を満たす」と「収束解が一致する」は分ける。** 保存記録で `run_0804_sp_g0c` の収束判定は `OVERALL: CHECK FAILURES ABOVE`、本段末尾 20,000 step の θ_r 3量・Q_w はすべて `DRIFTING`。[sp_judge.json:146](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/cold_pair/sp_judge.json:146)。M7 は腕間の流量差であり収支閉合ではなく、熱流束・せん断には代理物性による後処理値が含まれる。今回 A でも、定常解・物理精度・壁温感度の保証にはならない。

不足情報: ローカルには `run_0179`、`run_0717`、T2 両 run、`run_0804` の実ディレクトリがなく、指定の `delta_r_initial.csv` も確認できなかった。したがって入力 SHA、完全な設定一致、HDF5・残差履歴からの判定再実行は未確認。既存結果は保存された判定 JSON と plan を根拠にした。事前登録前に AWS 上の実物を照合する必要がある。

**plan 未反映。** 呼び出し側の反映先は `tooling-nozzle-isothermal-wall-chain.md` §4.8・§5.1 #29〜#31・§6 V-g0c-chain。依頼どおり、ファイル変更・forge 起動は行っていない。
