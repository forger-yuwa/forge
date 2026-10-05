# codex 諮問 (diagnose): m6-euler-recalibration

- **brief**: [`notes/reviews/briefs/2026-10-06-m6-euler-recalibration.md`](../../notes/reviews/briefs/2026-10-06-m6-euler-recalibration.md)
- **plan**: [`plans/active/tooling-nozzle-cfd-pinned-initial-line.md`](../../plans/active/tooling-nozzle-cfd-pinned-initial-line.md)
- **date**: 2026-10-06
- **commit**: `ccacd780` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 4.8 min, rc=0
- **結論**: **Euler の格子引数の受け渡しを修正・確認したうえで、凍結線と設計壁を固定し、新バイナリ・同一実効設定による G0/G1 の Euler A/B を先に行う。**
- **extra**: `case/45.isobutane_m6_d155/README.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（**草案のままの実行は不可。Euler 再較正の方針は採用し、最初の比較を修正する**）

コードは直接確認した。AWS の対象 run・残差・VERDICT 原本は今回確認できていないため、以下の実測値は plan の記録に基づく。

| 重大度 | 対象・採否 | 根拠と対案 |
|---|---|---|
| **Major 1** | G1 の生成方法：**要修正** | Euler の [`runner_axismach.py:539`](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner_axismach.py:539) が渡す格子引数は `ni/nj/wall_first_frac/throat_refine/scale` だけ。NS 経路の同ファイル:878 が渡す `throat_width`、`wall_first_frac_throat`、前後のブレンド等は無視される。**NS の YAML をコピーしても指定どおりの G1 にならない。** 両経路で格子パラメータの受け渡しを揃え、生成後の座標から実効配置を検査する。 |
| **Major 2** | 「Euler 差 ≈ NS 差 −0.00093 ±50%」：**判定根拠として却下** | `run_0094→0109` は格子だけでなく、`k_f 1.03736→1.055734`、`r_t 76.7531→76.6715 mm` と物理壁も変わっている。[plan:173](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-cfd-pinned-initial-line.md:173) がその解き直しを記録している。さらに出口平均は自格子の節点平均であり、評価作用素も変わる（[`nozzle_report.py:214`](/home/sano/work/forge-integ-1005/design/forge_design/report/nozzle_report.py:214)）。**同じ Euler 設計壁・共通標本で格子差を測り、NS の差との等値を要求しない。** |
| **Major 3** | 旧 G0 の数値をそのまま対照にする：**却下** | 旧較正の実行経路は [`prep_wallfit_euler.py:14`](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/prep_wallfit_euler.py:14)・:46 により **CFL 6、緩和 0.7**。YAML の `cfl_main: 2` が実効値ではない。新旧バイナリの確認も NS の δ_E・1000 step に限られる（[plan:170](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-cfd-pinned-initial-line.md:170)）。**G0 は新バイナリで再計算し、G1 と実効設定を揃える。** |
| **Major 4** | 閾値と NS 成功予測：**要再検証** | [草案 #11f](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-cfd-pinned-initial-line.md:98) の予測帯は補正量で `+0.000465〜+0.001395`。NS に係数 1 で移ると仮定しても、小さい側では `5.998329+0.000465=5.998794` と下限 `5.9988` を満たさない。G2 差が縮小するだけでも `1e-4` 精度の裏付けにはならない。**E2 の ±1e-4 は較正残差の基準として維持し、時間変動を含めて判定する。NS 成功、格子独立性、`r_t/k_f ±0.3%` の保証には使わない。** |
| **Major 5** | 再較正後の Euler 参照：**処置を追加** | [`run_finemesh_final.sh:10`](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/run_finemesh_final.sh:10) と [`exitM_sampling_ab.py:12`](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/exitM_sampling_ab.py:12) は旧 `run_0086` を固定参照し、後者はその場で δ_E も抽出する。**新設計へ進む際は合格した E2 を新しい固定 Euler 参照にし、抽出・較正・報告の参照先を同期する。** 旧参照の継続を無自覚に混ぜない。 |
| — | 凍結線 `run_0062` の維持：**採用** | 今回は固定された設計の評価格子依存を測るため、線まで変えると原因を分離できない。ただし、凍結線自身の格子依存が消去されたわけではない。[plan:111](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-cfd-pinned-initial-line.md:111) でも未確認事項である。 |
| — | 1 係数の線形補正：**初回更新として採用** | [V3′ の記録](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-cfd-pinned-initial-line.md:149) は粗格子の近傍で応答係数 ≈1 を支持する。細分格子・逆向き・より大きな補正での保証ではない。**初回は係数 1 で更新し、E2 を独立の検証にする。外れたら補正を重ねず再諮問する。** |

結論: **Euler の格子引数の受け渡しを修正・確認したうえで、凍結線と設計壁を固定し、新バイナリ・同一実効設定による G0/G1 の Euler A/B を先に行う。**

第 1 仮説: 粗い Euler 格子で決めた出口較正が、生産格子相当の配置では負側へずれる。 **確度: 中。ただし Euler の格子差は未測定。**

根拠: 粗格子 Euler の較正値は `case/45.isobutane_m6_d155/run_0083〜0088` で約 `6.000000`。一方、plan に記録された `case/45.isobutane_m6_d155/run_0109_ns_finemesh_final_ext` の出口平均は `5.998326`、共通の粗格子標本に補間しても `5.998520`。記録上の `check_quasisteady` は両者 **STEADY**。標本変更だけでは不足を消せない。ただし、これ自体は Euler 起因の証明ではない（[plan:176](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-cfd-pinned-initial-line.md:176)）。

反証条件: 同一条件の Euler A/B で、共通標本による差の評価幅全体が `|ΔM| < 2e-4` に入るなら、**今回の不足を主に説明する規模の Euler 格子持越し説**を棄却する。明確な正符号なら、負側の不足を説明する仮説として棄却する。

第 2 仮説: 上流を含む δ 閉包・物理壁変更が、NS の出口コア M に残差を作っている。 **確度: 中。** 出口一点の `δ_E/δ_C≈1` は、上流の有効壁形状や流れの一致を保証しない。Euler A/B が小差でも、この仮説を確定したとは扱わない。

第 3 仮説: バイナリ・実効 CFL・未収束床の差が較正値に寄与する。 **確度: 低、出口 M について未確認。** 新 G0/G1 間ではこれらを固定して交絡を防ぐ。

判別 A/B: **変更因子は評価格子だけ。**

- **腕 A＝新 G0**：旧較正の格子配置。
- **腕 B＝G1**：NS 生産格子の全分布パラメータを反映した配置。**Q2 はこの配置を推奨する。** `2000×97・wall_first_frac=5e-3` は別の格子なので、今回の生産配置への持越しを直接検証できない。薄セルの slip 壁誤差も、この配置の感度に含めて測る。ただし設計壁と NS 物理壁は異なり、同じパラメータでも物理座標・離散作用素まで同一とは呼べない。
- 両腕とも同じ設計壁・`r_t`・凍結線・ガス・BC・バイナリ。旧較正済み Euler 場を共通の初期場の出所とし、同一格子は `restart_field.py`、変更格子は `interp_field.py` を使う。本段は両腕とも **CFL/CFL_pseudo 2、緩和 0.7、2 次 SLAU、limiter 2** を明示し、その他の実効設定も揃える。
- 長さは **soft 3000 → 本段 12000 step、1000 step ごと保存**を初回予算とする。判定未達なら本段を **6000 step、1 回だけ延長**し、なお未達なら保留。薄セルで必要な長さは未確認なので、完走だけで判定しない。
- 生産指標の「自格子の節点平均」に加え、**G1 の出口帯内 η 節点列を固定標本として、両腕を同じ方法で補間・平均する**。この共通標本での差を D とする。自格子平均との差も残し、標本依存を分ける。
- `check_convergence.py` は同一設定の本段区間、`check_quasisteady.py --series-csv` は出口平均そのものに適用する。**STEADY に加え、末尾 5 枚の幅と直前 5 枚との平均差を各 5e-5 以下**とする提案。補間感度・単調時の漸近値も含め、閾値をまたぐ場合は判定しない。`NOT CONVERGED` は未収束のまま記録する。

**結果 A：D の評価幅全体が −2e-4 未満**  
→ 負側の Euler 格子依存を支持し、「格子依存は無視できる」を退ける。ただし NS 不足の原因確定ではない。

**結果 B：D の評価幅全体が ±2e-4 内**  
→ 第 1 仮説の主要因説を棄却する。境界層側を確定したとは言わない。

境界をまたぐ・準定常未達なら保留。**G2 はこの判別後に回し、E2・NS 連鎖へ自動続行しない。**

やらない方がよいこと: 旧 G0 の `6.000000` を新環境の無誤差な基準にすること、各格子固有の平均だけで離散化差を論じること、ピン壁から初期線を取り直すこと、NS の不足をそのまま `Md_moc_offset` に足すこと。旧 YAML は `spec.M_design=5.999584` なので、生産経路へ移す際は **`spec.M_design=6` と `geometry.Md_moc_offset` に統一し、二重補正を避ける**（[`runner_axismach.py:247`](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner_axismach.py:247)）。

呼び出し側の前提への異議: **「NS の細−粗＝純粋な格子差」「他ゲートは全て合格」は受け入れない。** 後者は、plan の最新記録で `overshoot01` が **DRIFTING** のため、値の閾値内と準定常合格を分ける必要がある。`run_0107/0109` の残差記録も **NOT CONVERGED (plateau)**。出口 M の FAIL と④保留を維持する。

不足情報: 対象 run の実効 config、生成座標、バイナリ・ガスのハッシュ、出口平均の時系列、判定区間付き VERDICT 原本。今回 `forge` は起動せず、ファイル変更もしていない。**plan 未反映**。呼び出し側で `plans/active/tooling-nozzle-cfd-pinned-initial-line.md` の **§5.1 #11f・§6** に反映する。
