# codex 諮問 (diagnose): m6-finemesh-pred-out

- **brief**: [`notes/reviews/briefs/2026-10-05-m6-finemesh-pred-out.md`](../../notes/reviews/briefs/2026-10-05-m6-finemesh-pred-out.md)
- **plan**: [`plans/active/tooling-nozzle-cfd-pinned-initial-line.md`](../../plans/active/tooling-nozzle-cfd-pinned-initial-line.md)
- **date**: 2026-10-05
- **commit**: `2d3a4e07` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 3.2 min, rc=0
- **結論**: **③の前に、未緩和 δ の入力だけを替える CFD 0 step の A/B で②の計算を確認し、再現できれば `k_f=1.055734`・`r_t=76.6715 mm` の③単独へ進む。**
- **extra**: `case/45.isobutane_m6_d155/README.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0 / Major 2 / Minor 1）

| 指摘 | 採否・根拠・対案 |
|---|---|
| **Major 1：Q1「全量を格子差と形状差で説明できる」** | **要再検証**。旧値 0.722471 は、計画の訂正履歴では未緩和 `delta_E` ではなく緩和後の値に相当する。旧 `k_f`・`r_t` が緩和値から計算されていたことは [plan:157](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-cfd-pinned-initial-line.md:157) に明記され、現行コードも [deltastar_loop.py:155](/home/sano/work/forge-integ-1005/design/forge_design/feedback/deltastar_loop.py:155) で訂正されている。**格子差・列の訂正・形状変更を分けて説明すること**。 |
| **Major 2：Q2「③④をそのまま続行」** | **③は下記確認後に採用、④への自動続行は却下**。現行 [run_finemesh_pin.sh:44](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/run_finemesh_pin.sh:44) は③の終了・報告から④へ進み、③の生産ゲートによる停止がない。報告失敗も同ファイル18行で握りつぶす。**③の全ゲート判定を④の前提にする**。 |
| **Minor 1：割線法の終了条件と報告残差** | **要再検証**。提示残差 −1.4e−6 は [c2pin_solve.py:43](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/c2pin_solve.py:43) の終了条件 `abs(f2)<1e-7` を満たさない。反復上限到達でも結果を出力する実装なので、`k_f_hist` と実行コードを照合する。相対残差は約 0.00019 % であり、これ自体が今回の大きな予測外れを説明するとは考えない。 |

結論: **③の前に、未緩和 δ の入力だけを替える CFD 0 step の A/B で②の計算を確認し、再現できれば `k_f=1.055734`・`r_t=76.6715 mm` の③単独へ進む。**

第 1 仮説: **予測外れの主因は較正対象 δ の増加であり、その増加には格子差と旧来の緩和値使用の訂正が含まれる。** 確度: **中**。

根拠:
- 提示値の算術を再計算すると、細格子／旧入力は **+1.47065 %**、細格子／同形状の粗格子未緩和値は **+1.07821 %**、粗格子未緩和値／旧入力は **+0.38825 %**。最後の項を「形状・`k_f` の差だけ」とする根拠はない。
- [solve_rt の式](/home/sano/work/forge-integ-1005/design/forge_design/feedback/deltastar_loop.py:166) に提示値を代入すると、Re 補正後 δ は **0.73325190**、出口半径は **0.77499991 m**。半径変化の説明は成立する。
- 一方、`p=0.8` の近似からは `k_f≈1.05646`。方向と大きさの説明にはなるが、**実装による全量再現の証明ではない**。実装は `Cf` を倍率変更して運動量厚さを積分するため、固定べき乗則ではない（[deltastar_integral.py:168](/home/sano/work/forge-integ-1005/design/forge_design/feedback/deltastar_integral.py:168)、[同:202](/home/sano/work/forge-integ-1005/design/forge_design/feedback/deltastar_integral.py:202)）。

反証条件: 設計壁・ガス・基準半径・平滑化・コードを固定し、提示された未緩和 δ を与えても今回の `k_f`・`r_t` を再現できない場合、「δ の変更だけで説明できる」は棄却する。

第 2 仮説: AWS 側のコード・入力・平滑化条件と、レビュー対象との不一致。確度: **低・未確認**。割線法の終了条件との不整合が照合理由になる。

判別 A/B: **変更するのは `solve_rt` と `dC(k_f)` に渡す未緩和 δ(x_F) だけ**。

- 腕A：`case/45.isobutane_m6_d155/run_0105_ns_coarse_cfl1` の **0.725276**。
- 腕B：`case/45.isobutane_m6_d155/run_0103_ns_finemesh_pass_cfl1` の **0.733096**。
- 共通：`r_t,prev=76.7531 mm`、同じ MOC 壁・熱力学条件・CONTUR・平滑化。**forge は起動せず**、半径反復と積分法の割線求解だけを行う。
- 半径の予測は、提示された丸め済み `r_F=9.3748` なら **A≈76.73179 mm、B≈76.67155 mm**。`k_f` は腕Bが大きくなるはず。
- **結果A：** 腕Bが今回値を `|Δk_f|≤1e−4`、`|Δr_t|≤0.001 mm` で再現し、両腕の較正残差が相対 `1e−5` 以下 → 入力 δ の変更による説明を支持し、③へ進む。
- **結果B：** 再現しない → 「予測が古いだけ」を棄却し、入力列・尺度・実行コードの不一致を調べる。③は保留。

やらない方がよいこと: 予測範囲を書き換えること、`k_f` を手で戻すこと、NS の出口 M に合わせて `M_design` を再較正すること、③の判定前に④へ進むこと。

呼び出し側の前提への異議:
- **Q1:** 検算は変化量の整合確認として有効だが、旧 δ の定義が混ざっているため原因説明としては不足する。
- **Q2:** `k_f=1.055734` という数値だけを不自然として棄却する根拠はない。ただし、これは CONTUR の局所閉包に掛ける倍率であり、実際の壁面摩擦が一律 5.6 % 増えることや、物理的妥当性を保証する値ではない。
- **Q3:** 既存の数値閾値は維持する。追加すべきは、①生成した実際の物理壁で出口半径と δ_C を確認、②③の未緩和 δ_E・出口コア M・波・オーバーシュート自身に準定常判定、③段階起動を混ぜず同一設定区間で全残差を判定、④その結果で④への進行を止められる仕組み。`R_exit_pred` は目標 δ を使った代数的予測なので、実際の壁の検査を代用できない（[c2pin_solve.py:55](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/c2pin_solve.py:55)）。

不足情報: 対象 run と `c2pin_solve_{pass2,fine}.json` はローカルに存在しない。AWS の原 CSV、求解履歴、実行コード、各 VERDICT は直接確認できていない。計画に記録された判定は **`NOT CONVERGED (plateau)`／δ_E `STEADY`** であり、収束済みとは扱わない。旧 0.722471 の実際の列名確認も必要。

**ファイル変更・forge 実行なし。plan 未反映。呼び出し側が本判断を対象 plan §5.1 #11・§9 に反映すること。**
