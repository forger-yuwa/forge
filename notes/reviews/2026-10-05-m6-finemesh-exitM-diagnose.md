# codex 諮問 (diagnose): m6-finemesh-exitM

- **brief**: [`notes/reviews/briefs/2026-10-05-m6-finemesh-exitM.md`](../../notes/reviews/briefs/2026-10-05-m6-finemesh-exitM.md)
- **plan**: [`plans/active/tooling-nozzle-cfd-pinned-initial-line.md`](../../plans/active/tooling-nozzle-cfd-pinned-initial-line.md)
- **date**: 2026-10-05
- **commit**: `35ca6472` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 3.1 min, rc=0
- **結論**: **④は保留し、既存場を使う CFD 0 step の「出口コア平均の標本位置だけを替える A/B」を先に行う。**
- **extra**: `case/45.isobutane_m6_d155/README.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **④は保留し、既存場を使う CFD 0 step の「出口コア平均の標本位置だけを替える A/B」を先に行う。**

採否表:

| 対象 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| A：細分 Euler 1本で H1/H2 を判別 | **次手として却下・Major** | `exit_core` は帯内節点の単純平均で、格子間で標本位置・重みが変わる。[nozzle_report.py:214](/home/sano/work/forge-integ-1005/design/forge_design/report/nozzle_report.py:214)。さらに「Euler 用第1セル」の格子は NS の半径方向配置を再現しない。同じ `ni/nj` だけでは対照にならない。まず集計の寄与を分離する。 |
| 「δ_E/δ_C≈1 なので境界層側では説明できない」 | **却下・Major** | 出口1点の閉包成立から上流全域の成立は導けない。過去の1D換算も訂正後なお実測との差が残る。[plan:159](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-cfd-pinned-initial-line.md:159)。境界層・物理壁側を候補に戻す。 |
| B：登録 FAIL を維持し④保留 | **採用** | 現行条件は出口 M ±0.02%、④は③のゲート判定後。[plan:92](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-cfd-pinned-initial-line.md:92)、[plan:96](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-cfd-pinned-initial-line.md:96)。 |
| C：今回の値を理由にゲートを見直す | **却下・Major** | 現行登録に、格子変更を理由とする免除規定はない。診断用の集計変更も、既存 FAIL の合格への読み替えには使わない。仕様変更が必要なら、旧 FAIL を保持した別の設計判断とする。 |

第 1 仮説: **格子による節点平均の重みの変更が、出口 M の差に無視できない寄与を混ぜている。** 確度: **中**。FAIL 全体を説明する大きさかは未確認。

根拠: [nozzle_report.py:214](/home/sano/work/forge-integ-1005/design/forge_design/report/nozzle_report.py:214) は、η∈[0.05, 0.7] の `Mx.mean()` を使う。共通ηへの補間や面積重みはない。第1セル変更で軸側間隔も 0.0752→0.0888 r_w に変わることは [plan:117](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-cfd-pinned-initial-line.md:117) に記録済み。提示された run_0107 の帯内 M 幅は 0.00274 であり、集計位置の影響を無視できる根拠はない。

反証条件: 下記 A/B で、末尾5枚すべての集計差が |ΔM|<1e−4 なら、「標本位置変更が今回の差の主要因」という仮説を退ける。

第 2 仮説: **共通の指標でも実際のコア場が格子に依存する。** 確度: 低〜中、未確認。ブリーフの H1 と H2 は、どちらも非粘性部分の格子依存を述べており、独立した対立仮説ではない。

第 3 仮説: **上流を含む δ 閉包・物理壁形状の差が出口に残っている。** 確度: 低〜中、未確認。比較元と比較先で `k_f`・`r_t`・壁形状も変わっているため、粗細差を格子だけには帰属できない。

判別 A/B: **変えるのは後処理の標本位置だけ。追加計算は 0 step。**

- 対象は `case/45.isobutane_m6_d155/run_0109_ns_finemesh_final_ext/` の末尾5枚。同一スナップショットの出口 M(η) を固定する。
- 腕A：現在のコードどおり、細分格子の帯内節点を単純平均する。
- 腕B：同じ M(η) を、`case/45.isobutane_m6_d155/run_0094_ns_c2pin_pass2_ext6k/` の帯内η節点位置へ線形補間して単純平均する。補間方式・参照節点列を結果を見る前に固定する。
- 見る量は各時刻の ΔM_sampling = M_B−M_A、その末尾平均・幅、および粗格子との差がどれだけ変わるか。各腕の時系列には `check_quasisteady.py --series-csv` の VERDICT を付ける。

→ **B−A が正で、末尾5枚を通じて +0.000471 以上なら**、標本位置だけで今回の下限超過量を説明できる可能性を支持する。  
→ **全枚で |B−A|<1e−4 なら**、標本位置変更の主要因説を棄却する。ただし、その結果だけで Euler 格子依存を確定しない。  
→ 中間、逆符号、または準定常未達なら、寄与だけ記録して原因判定は保留する。**腕Bは診断値であり、生産ゲートの代替値ではない。**

やらない方がよいこと: NS の出口 M に合わせた `Md_moc_offset` の再較正、無条件の再延長、±0.02% の事後緩和。`Md_moc_offset` は表示値の補正ではなく MOC・軸 law・壁へ渡る設計入力である。[runner_axismach.py:247](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner_axismach.py:247)。

呼び出し側の前提への異議:

- **Major：末尾幅が小さいことを準定常 VERDICT の代わりにしている。** `nozzle_report.py` は末尾5枚の最大−最小を出すだけで、準定常判定を実行していない。[nozzle_report.py:237](/home/sano/work/forge-integ-1005/design/forge_design/report/nozzle_report.py:237)。今回の出口 M・波・オーバーシュート・δ 比の判定記録が必要。
- バイナリ変更の同等性確認は、1000 step の δ_E 比較である。[plan:167](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-cfd-pinned-initial-line.md:167)。これを出口 M の同等性証明へ拡張できない。旧 Euler と新バイナリの細分 Euler だけを比較すると、この交絡も残る。
- **④へ進む条件**は、③の登録性能・壁解像ゲートの成立、対象量の `STEADY`、同一設定区間の全残差 VERDICT の記録が揃うこと。`NOT CONVERGED` はそのまま明記する。診断 A/B の成功だけでは④へ進まない。
- **ユーザへの位置づけ**は「資料記載値では出口コア M=5.998329、目標比−0.02785%。許容下限を0.000471、目標比0.00785ポイント下回り、登録 FAIL。η=0.1 の波・オーバーシュートは別指標として閾値内」。一般的な「試験部全体の一様性合格」や「数値誤差なので実質合格」には広げない。

不足情報: 対象 run はローカルに存在せず、AWS の場・実効設定・VERDICT 原本は未確認。資料記載の収束判定は **run_0107／0109 とも `NOT CONVERGED (plateau)`**。今回必要な準定常 VERDICT、出口η節点列と M 分布、Euler 較正時との評価器の同一性が不足している。

ファイル変更・forge 起動は行っていない。**plan 未反映**。依頼どおり、反映は呼び出し側が行う。
