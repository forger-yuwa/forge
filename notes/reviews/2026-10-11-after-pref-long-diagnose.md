# codex 諮問 (diagnose): after-pref-long

- **brief**: [`notes/reviews/briefs/2026-10-11-after-pref-long.md`](../../notes/reviews/briefs/2026-10-11-after-pref-long.md)
- **plan**: [`plans/active/architecture-float-state-double-geometry.md`](../../plans/active/architecture-float-state-double-geometry.md)
- **date**: 2026-10-11
- **commit**: `0e99e66b` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 4.7 min, rc=0
- **結論**: **案 P を、全格子・point 固定で精度だけを変える再現試験として先に行い、「B0 の構成がこの大きな精度差に必要か」を調べる。**
- **extra**: `case/45.isobutane_m6_d155/_band_ab/cold_pair/lp_diffmap.txt`, `case/45.isobutane_m6_d155/_band_ab/cold_pair/lp_judge.json`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：

| 重大度 | 対象 | 採否・根拠・対案 |
|---|---|---|
| Major | 「pRef は下流の差を減らさなかった」 | **却下。** 保存済み `lp_judge.json` では、W2 の θ_r の差は F0−Ap が +1.450〜2.215%、Bp−Ap が +1.223〜1.939%。対応する断面では約8〜17%小さい。ただし F0 は1本で、正式判定も前提不成立。**「登録した抑制目標の達成は確認できず、今回の2本には大きな差が残った」**と書く。[plan §6.47](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:1087) |
| Major | 初期の型とは「別の通り道」がある | **要再検証。** 初期の射影の縮小と後期の差の残存だけでは、別成分の成長と、残った同じ成分の遅れた増幅を分けられない。`lp_diffmap.py` は領域別の絶対差の統計で、モードや因果経路を同定していない。**場所・時点・正規化を付けた観測として記録する。** [抽出処理](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/lp_diffmap.py:25) |
| Major | forge-25 が H-slow を裏付ける | **仮説として採用、原因認定は要再検証。** 記録は Gref の観察で、対象 G0 とは壁法線解像度が異なる。また point 化は `lineImplicit`・`lineDtDirectional`・ITJ の組を変える。**「B0 の参照が未収束の途中状態かもしれない」まで採用する。** [plan §5.1 #22](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:184) |
| Major | 案 P | **修正して採用。** point を両精度に共通で固定する「精度の A/B」とする。B0 との比較から個別設定の原因を特定しない。方向別 dt は内部 line 面のスペクトル半径を除外するため、同じ30,000 step は同じ進行量を意味しない。**差が小さいだけで抑制判定しない。** [設定の意味](/home/sano/work/forge-integ-1005/procedures/solver-settings.md:75) |
| Minor | Q_w の差の表記 | **訂正を採用。** ブリーフの +1.8〜2.5% は主に Bp 自身の起点からの移動。Ap の移動を引いた差は全窓で **+1.648〜2.304ポイント**。θ_r と同じく、ΔBp と ΔBp−ΔAp を区別する。[指標の定義](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/lp_an.py:80) |

結論: **案 P を、全格子・point 固定で精度だけを変える再現試験として先に行い、「B0 の構成がこの大きな精度差に必要か」を調べる。**

第 1 仮説: **point に替えても、十分に進行させれば大きな精度差が残る。** B0 の停滞だけでは説明できず、初期化や step 内の演算の精度差が関与する。確度: **低**。

根拠: pRef の下でも、`case/45.isobutane_m6_d155/run_0710_lp_bp1/`・`run_0711_lp_bp2/` の起点からの変化は、`run_0712_lp_ap/` に対して θ_r で +0.905〜1.939ポイント離れている。また、保存量9つの末尾量子化だけでは初期の型を再現しない。ただし、これらは **point での再現の証拠ではない**。[量子化試験の結論](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:1028)

反証条件: 下記の進行・定常性の条件を満たした point の比較で、3断面とも精度差が変動幅込みで 0.1|D| 以下になること。棄却するのは、この条件での「大差が残る」という予測であり、残差演算の丸めへの関与全般ではない。

第 2 仮説: **B0 の反復特性が精度差を大きく増幅しており、point では大差を維持しない。** 確度: 低。forge-25 の観察が動機だが、G0 で未確認。

第 3 仮説: **精度別の初期化差の増幅が一部を担う。** 未確認。共通の入力 Q32 と、初期化後の保存量の一致は別である。[plan §6.31](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:800)

判別 A/B:

- **変える因子はビルド精度だけ。** 両腕とも `time.deltaT` の `lineImplicit: 0`・`lineDtDirectional: 0`・`implicitThermalJacobian: 0` を固定する。共通の `_tr/q32_init.h5`、G0 全格子、同一 BC・CFL・緩和・空間離散化を使う。pRef は §6.26 と同じ0。FP64・float を各2本とし、FP64 側の再実行幅も測る。
- **期間は最大120,000 stepを事前登録する。** まず30,000 stepまで実施し、20,000〜25,000／25,000〜30,000の両窓で明瞭な大差を再現した場合だけ早期終了できる。差が小さい場合は抑制と読まず、登録済みの100,000〜110,000／110,000〜120,000の窓まで進める。12万 step 自体を収束保証にはしない。
- 主量は θ_r(40/70/94)。各 run 自身の起点からの変化を取り、**point の FP64 の変化を差し引く**。D = (1.892, 1.647, 1.431)% は尺度としてのみ使用する。Q_w、上流・スロートの温度差、全残差、局所擬似時間刻みは別に記録する。
- **A：大差が残る。** 両窓・3断面・全再実行の組合せで、変動幅を差し引いても D の向きの差が 0.5|D| 以上。→ 第1仮説を支持し、**「B0 の構成が大差に必須」を退ける**。面エンタルピーが原因とは言わない。
- **B：大差が残らない。** 最終両窓・全組合せで、変動幅込みの絶対差が 0.1|D| 以下。さらに全残差の収束判定と、対象量の厳しい準定常判定が通る。→ 第1仮説の予測を退け、第2仮説を支持する。B0 のどの設定が効いたかは未確定。
- 進行不足、未収束・ドリフト、反復間のばらつき、中間的な差は**判別不能**。窓内の保存点を独立反復に数えない。各精度内の res_0 の一致と、精度間の初期化差も記録する。

やらない方がよいこと: **案 S・E・低マッハ補正を同時に始めない。** S の「出口が超音速」は中心流だけでは足りず、壁際の亜音速域を含めた切断面の確認が必要。E は float ビルドで実際にどの演算・格納まで double になるか未確認である。今回の試験を理由に opt-in を削除したり既定値を変えたりしない。

呼び出し側の前提への異議: **H-slow と H-res は排他的ではない。** 残差の偏りが減衰の遅い方向を継続的に励起することは両方に該当する。「FP64 の揺れは10⁹倍小さい」も実測ではなく、演算精度の比を軌跡の揺れの比に置き換えている。case/48 についても、書けるのは登録した有限期間の再現条件を満たさなかったことまでで、「収束解が一致」「line 系は除外」とはしない。

既存結果の根拠となる保存済み判定は、`run_0709_lp_f0/`〜`run_0712_lp_ap/` の4本について次のとおり。

```text
check_convergence: OVERALL: CHECK FAILURES ABOVE
check_quasisteady: DRIFTING
判定区間: 20,000〜30,000 step、drift/osc = 1e-4
§6.46: 判別不能
```

不足情報: このワークツリーには上記 run の元の CSV・HDF5 が無く、今回は [保存済み判定記録](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/cold_pair/lp_judge.json) とコードを照合した。有限性・実効設定・収束ツールの再実行は独立確認できていない。forge-25 の根拠 run、全残差の VERDICT、準定常判定の閾値・区間も必要。

**plan 未反映。** 指示どおりファイルは変更していない。呼び出し側で §6.47 の解釈、§5.1 #22、次の試験の事前登録へ反映すること。
