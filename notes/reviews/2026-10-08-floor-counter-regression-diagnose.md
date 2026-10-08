# codex 諮問 (diagnose): floor-counter-regression

- **brief**: [`notes/reviews/briefs/2026-10-08-floor-counter-regression.md`](../../notes/reviews/briefs/2026-10-08-floor-counter-regression.md)
- **plan**: [`plans/active/tooling-sern-te-wake-grid.md`](../../plans/active/tooling-sern-te-wake-grid.md)
- **date**: 2026-10-08
- **commit**: `9c930e0a` (feature/sern-design)
- **codex**: effort `high`, 4.8 min, rc=0
- **結論**: **現状だけで TP の数値非干渉を受け入れず、同一の凍結入力に対する「旧版 EOS 対カウンタ有効の新版 EOS」の全出力ビット比較を追加する。**
- **extra**: `procedures/verification/README.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **現状だけで TP の数値非干渉を受け入れず、同一の凍結入力に対する「旧版 EOS 対カウンタ有効の新版 EOS」の全出力ビット比較を追加する。**

採否表:

| 重大度 | 判断対象 | 採否・根拠・対案 |
|---|---|---|
| Major | 既存のビット比較＋CPG 回帰＋小型 A/B で TP も受け入れる | **要再検証**。TP の `res_1` 比較では `roe` などが除外され、`gamma`・`cp`・`Rmix` は保存されていない。これらは EOS の更新対象である（[dependentVariables_d.cu:230](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/dependentVariables_d.cu:230)、同:251）。**EOS 直後・残差組立て前で全出力を比較する。** |
| Major | TP の時間発展の超過を変更の悪影響と認定する | **却下**。保存データから追加旧版の超過も再現した。下記の数値は、現行の物差しが旧版自身を通せない証拠である。既存結果は「判定不能」のまま保存し、別の比較を事前登録する（[verification/README.md:89](/home/sano/work/forge-sern-design/procedures/verification/README.md:89)）。 |
| Minor | 比較スクリプトが非有限を必ず拒否するという前提 | **却下**。`regress_floor_events.py:156` で校正幅が `inf` になり得る一方、:186 は `worst <= W` だけなので `inf <= inf` が通る。幅・差の有限性を明示的に必須化する。今回再確認した評価標本には非有限はなく、今回の判定不能の原因ではない。 |

第 1 仮説: **TP の超過の主因は、旧版にもある実行間変動とその時間発展での増幅であり、EOS の数値写像の変更ではない。** 確度: **中**  
  根拠: 以下を HDF5 と残差 CSV から独立に再計算した。パスの基点は [scratchpad](/tmp/claude-1000/-home-sano-work-forge/60a080db-fad5-4e19-bdaa-3590a1035176/scratchpad)。

| 評価標本 | 量 | 最大絶対差 | 凍結幅 |
|---|---|---:|---:|
| `regress_tp/sN_head_r4`、200 step | `T` | 41.5969 K | 32.5696 K |
| `regress_tp/sN_new_on_r1`、200 step | `Ux` | 36.5339 | 34.6726 |
| `regress_tp_dev/sN_head_r6`、200 step | `P` | 1132.7109 Pa | 1005 Pa |
| `regress_tp_dev20/sN_head_r6`、20 step | `rms_roUx` | 0.01306749 | 0.01237048 |

  既存の比較対象についてはビット一致を再確認した。小型試験の記録も TP は `59 OK / VERDICT: PASS`、CPG は `50 OK / 1 SKIP / VERDICT: PASS`。ただし、これで未比較の EOS 出力まで保証したことにはならない。  
  反証条件: **完全に同じ EOS 入力から、旧版と新版で流れの配列に再現性のあるビット差が出ること。** その場合、「時間発展の非決定性だけ」という説明を棄却する。

第 2 仮説: カウンタ追加に伴うコンパイル結果の変化が、未比較の EOS 出力を変えている。確度: **低・未確認**。`feOn` の分岐追加は確認できるが、機械語由来の数値差は測定されていない（[dependentVariables_d.cu:80](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/dependentVariables_d.cu:80)）。

判別 A/B: **変えるのは EOS の実装だけ：旧版／新版カウンタ有効。時間積分は行わず、各凍結入力から EOS を１回実行する。**

- 入力集合を先に固定する。TP 初期場・発達場に加え、対象の g3 継続場と既存小型試験の床上下の状態を使う。保存量だけでなく、温度反転の初期推定に使う入力 `T`、組成、物性 DB、設定も同一にする。
- 境界ピン後・EOS 直前の状態を複製し、**壁・ghost を除外せず**比較する。`ro`・`roe`・速度・`T`・`P`・`Ht`・音速・`k`・`omega`・`gamma`・`cp`・`Rmix` を含む全書込み先と、変更されないはずの入力を対象にする。比較対象を反復結果から選別しない。
- 合格条件は、有限性・配列構造の一致と**差分バイト数 0**。各腕２回を同じ入力から実行して比較器の再現性も確認する。

→ **結果 A：全比較がビット一致**なら、試験した状態での第 2 仮説を退ける。既存の監査不変・計数試験と合わせ、対象経路に限定して共通バイナリによる g3/g4 検証へ進めてよい。TP 時間発展は「判定不能」のまま、代替検証と適用範囲を plan に明記する。  
→ **結果 B：再現する差がある**なら、第 1 仮説の「非決定性だけ」を退け、最初に異なる配列・節点・演算を調べる。非有限・入力不一致・腕内の非再現は試験不成立として進めない。

やらない方がよいこと: 幅・本数・時間窓を後付けで変更して合格を得ること。現段階で残差の `atomicAdd` を決定化すること。**SASS 全体一致を受入れ条件にすること**も勧めない。分岐・計数命令の差は不具合を意味せず、無効経路の一致だけでは有効経路を保証できない。SASS は上記 B になった場合の差分定位に使う。

呼び出し側の前提への異議: 「旧版５回で一致した出力＝EOS 全体の決定性」は受け入れない。また、旧版の追加３本が幅内だった CPG についても、確認できたのはその標本での受入れであり、検出力が定量的に保証されたわけではない。共通バイナリは g3/g4 間の版の交絡を防ぐが、カウンタの数値非干渉の証明にはならない。

不足情報: EOS 直前・直後の全配列比較、および検証バイナリのハッシュとソース・ビルド条件の対応。**plan 未反映**。呼び出し側で `plans/active/tooling-sern-te-wake-grid.md` §5.1 #2・§6 に上記の代替検証と受入れ範囲を反映する。
