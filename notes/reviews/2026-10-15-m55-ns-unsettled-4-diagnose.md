# codex 諮問 (diagnose): m55-ns-unsettled-4

- **brief**: [`notes/reviews/briefs/2026-10-15-m55-ns-unsettled-4.md`](../../notes/reviews/briefs/2026-10-15-m55-ns-unsettled-4.md)
- **date**: 2026-10-15
- **commit**: `9edb227f` (feature/card-local-mesh)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 別 (依頼側 `claude-opus-5-5`、実効 `gpt-6-astra`)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 3.4 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-15-m55-ns-unsettled-4-diagnose.codex.log` (codex: ok)
- **結論**: 長時間延長とδ*補正を保留し、既存の後期保存場から各1更新だけ再生して、`output.interfaceDiag` による保存形壁熱量と現在の `qwall` 積分を比較する。
- **extra**: `case/71.m55_combustion_wt/_runs/run_0003_ns_cont/judge/check_quasisteady_monitor.txt`, `case/71.m55_combustion_wt/_runs/run_0003_ns_cont/judge/check_quasisteady_series.txt`, `plans/active/tooling-nozzle-card-local-mesh.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表:

| 対象 | 判断 | 重大度・根拠・対案 |
|---|---|---|
| 直ちに約6万 step延長する | **要再検証** | **Major**：停止規則を満たしても、現在の熱収支の定義問題は解消しない。先に下記の短い診断を行う。なお、新系列の初出力が2,000なら窓が揃う最速は22,000、保持込みで**最短62,000 step**（`stop_rule.py:570`, `:626`）。 |
| θ_r・出口Mの `STEADY` を理由にδ*補正へ進む | **却下** | **Major**：汎用の `STEADY` は停止規則の代替にならない。さらに、手元の記録を既存の `two_window_check` に渡すと、160k〜200kのδ記録不足で **pass=False、判定不能**。δ_E・δ_targetの2窓判定を完了させること（`deltastar_loop.py:698`）。 |
| 収支差1.9%をソルバの非保存と解釈する | **却下** | **Major**：現在の `Q_wall` は温度拘束の反力を含まない。収支用には既存の `iface_Qf_eff` を使う（下記）。非保存の可能性自体は未否定。 |
| double延長で `thermoFloat` 起因を除外する | **要再検証** | **Major**：複数経路を同時変更し、旧経路の保存場から継続している（plan `tooling-nozzle-card-local-mesh.md:73`）。「現在のdouble経路でも超過が残る」は採用するが、旧経路の寄与まで除外しない。追加の精度A/Bは今回の優先事項にはしない。 |

結論: **長時間延長とδ*補正を保留し、既存の後期保存場から各1更新だけ再生して、`output.interfaceDiag` による保存形壁熱量と現在の `qwall` 積分を比較する。**

第 1 仮説: **報告されたエネルギー収支差の主成分は、壁熱量の定義の不一致である。** 確度: 中  
  根拠:
  
- `run_0003_ns_cont/judge/monitor_series.csv` の200kを再集計すると、差は **102,647.97 W**、入口エネルギー流束比 **7.37056e−4**。
- [`monitor_0003.py:41`](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/_runs/run_0003_ns_cont/monitor_0003.py:41) は `qwall` を壁節点の折れ線で台形積分している。
- [`viscousFlux_d.cu:683`](/home/sano/work/forge-cardmesh/solver_density_cuda/cuda_forge/viscousFlux_d.cu:683) の `qwall` は壁面伝導寄与を面積で割った値。一方、[`nodeWallDirichlet_d.cu:94`](/home/sano/work/forge-cardmesh/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:94) は等温壁ノードのエネルギー残差をゼロ化する。この拘束を収支から落とせない。
- 既存実装の保存形熱量は  
  **`iface_Qf_eff = ifaceRraw − ifaceFw − (roe/ro)·ifaceRro`**、固体向き正。最後の項は壁CVの半離散的な蓄積補正である（[`conjugateWall.cpp:277`](/home/sano/work/forge-cardmesh/solver_density_cuda/conjugateWall.cpp:277)）。

  反証条件: 同じ場で保存形熱量へ置き換えても、後述の収支差の大きさが半分以下にならなければ、「壁熱量の定義差が主成分」を棄却する。**定義が違うことは確定しているが、102.6 kWを説明する量的寄与は未測定。**

第 2 仮説: 開境界の節点積分と実際の境界流束の差、および残存過渡が収支差に寄与している。確度: 中、寄与量は未確認。`monitor_0003.py:91` は節点のρh₀Uを積分するが、node境界対流は境界状態を使う（`convectiveFlux_d.cu:500`）。

第 3 仮説: 全圧超過には冷却によるエントロピー変化が寄与している。確度: 低、未確認。固定集合Sの平均h₀不足と平均P₀超過だけでは、物理的な説明にも数値欠陥の証明にもならない。

判別 A/B: **変更する設定は `output.interfaceDiag` の1点、A=0／B=1。**

- `case/71.m55_combustion_wt/run_0003_ns_cont/` の **140k・160k・180k・200k** をそれぞれ共通入力とし、同一メッシュのindex restartで各腕 **1更新だけ**。本段の数値設定・精度・BCは維持する。これは収束比較ではなく、熱量の測定経路の比較である。
- A/Bの保存量が実行ばらつき内で同じことを確認する。診断は既存実装上、残差を変更しない（`nodeWallDirichlet_d.cu:268`）。
- 各入力について、入口・出口積分を共通に固定し、次を比較する。  
  **A:** `D_A = E_in − E_out + Q_wall_trap`  
  **B:** `D_B = E_in − E_out − 2πΣ iface_Qf_eff`
- **4入力の ‖D_B‖₂ ≤ 0.5‖D_A‖₂ なら**、第1仮説を支持し、「元の差を非保存の証拠とする解釈」を棄却する。**それ以外なら**、第1仮説の「主成分」を棄却する。どちらでも、ソルバ全体の保存性が証明されたとは扱わない。

壁積分の照合では、Bの同じ出力から **`2πΣ ifaceFw`** も集計する。これは静止壁で実際に加算した壁面寄与であり、台形積分との差が幾何・積分の差、`iface_Qf_eff` との差が拘束・蓄積補正になる。`qwall`から復元する場合は `BCONDS/3/iPlanes` に対応する実効 `ss` を使う。区間r重みが有効なら **`ss = ‖PLANES/rSurfVect‖`** である（[`variables.cpp:795`](/home/sano/work/forge-cardmesh/solver_density_cuda/variables.cpp:795)）。**積分済みの `iface_Qf_eff` に面積やrを再度掛けない。**

やらない方がよいこと: 停止閾値を緩めること、未収束runを参照床にすること、`qwall`との差を見て直ちに壁BCを変更すること、全圧超過の最大値だけでδ*補正への影響を小さいと決めること。

呼び出し側の前提への異議:

- **Major：「θ_rがSTEADYだから受け渡せる」は成立しない。** 再実行した `check_quasisteady` はθ_rを `STEADY` としたが、既定のドリフト許容は末尾窓で5%（`check_quasisteady.py:543`）。停止規則は別定義で、再生結果も **D=6.3163e−4・5.3347e−4、閾値5e−4超過、`first_stop=None`**。対案は両判定を別々に記載すること。
- **Major：「δ*補正はθ_rと出口Mで判断する」という枠組みが不足している。** 実装はδ_E・δ_targetの全使用域について、隣接2窓の変動を **0.01 mm以下**で確認する。今回の写しではこの判定が不能。Euler参照もplan上で未適格の記録が残る（plan:71）。対案は正式な受け渡し条件を揃えること。
- **Major：P₀超過とh₀不足の同居は、それだけでは矛盾ではない。** `total_quantities.py:147` のP₀はエントロピー差から計算される凍結組成全圧であり、冷却流に入口Pt以下という一般的な上限を置けない。逆に、同居だけで「冷却だから正常」とも言えない。壁熱量A/BはP₀の原因を直接判別しないため、超過は未解明として残す。
- **Minor：終了状態はまだ確定記録ではない。** 手元の `stop_rule.json` は `RUNNING`。`CENSORED` の保存はδ抽出待ちの後である（`stop_rule.py:1003`）。対案は「予算上限終了、停止条件未到達、終了記録確定待ち」と書くこと。

結果は次のように並べて報告するのが妥当:

> `case/71.m55_combustion_wt/run_0003_ns_cont/`：本段0〜199,999の `CONVERGENCE_VERDICT` は **NOT CONVERGED**。停止規則は200kで未到達。θ_rの3断面は汎用準定常判定で **STEADY・単調増加**、推定漸近値は0.0536868／0.0798884／0.0902438。ただし、r・S平均h₀・S平均P₀・現行収支差は **DRIFTING**。δ*補正への受け渡しは未承認。

不足情報: HDF5と全残差履歴は手元になく、場の数値は回収済みJSON・CSV、収束は保存された判定原文に依拠した。停止規則・準定常・δの2窓判定は既存ツールで再評価した。上記の短い診断を選ぶ根拠は揃っているが、ソルバの保存性や全圧超過の真因は未確定。

証拠: 足りる
