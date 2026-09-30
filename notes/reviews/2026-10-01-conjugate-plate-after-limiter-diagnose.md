# codex 諮問 (diagnose): conjugate-plate-after-limiter

- **brief**: [`notes/reviews/briefs/2026-10-01-conjugate-plate-after-limiter.md`](../../notes/reviews/briefs/2026-10-01-conjugate-plate-after-limiter.md)
- **plan**: [`plans/active/boundary-cht-conjugate-flat-plate.md`](../../plans/active/boundary-cht-conjugate-flat-plate.md)
- **date**: 2026-10-01
- **commit**: `ee3a8b47` (feature/cht-conjugate-benchmarks)
- **codex**: effort `high`, 3.2 min, rc=0
- **結論**: **(a) を選び、既存 A を固定対照として B の同一設定継続を追加 36000 step・累計 48000 step に限定して行い、延長区間の判定を事前登録する。**
- **extra**: `plans/accepted/boundary-cht-conjugate-benchmarks.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0 / Major 3 / Minor 1）

| 重大度 | 論点 | 採否・根拠・対案 |
|---|---|---|
| Major | H1「係数の切り替わりが停滞の原因」 | **要再検証**。`limiter: 0` は係数をすべて 1 にするため、時間変動と平均的な再構成の両方を変える（[limiter_d.cu:450](/home/sano/work/forge-cht/solver_density_cuda/cuda_forge/limiter_d.cu:450)）。現状は「再構成への依存を示唆する」と記載し、切り替わり単独の因果とはしない。 |
| Major | (a) PASS まで延長 | **固定期限・固定区間を付けて採用**。B は保存された判定で全能動列が下降中。ただし現行判定器は先頭 12000 step 固定で、延長しても後半を評価しない（[ab_judge_lim.py:72](/home/sano/work/forge-cht/case/65.conjugate_flat_plate/ab_judge_lim.py:72)）。延長判定を別登録し、旧「中間」を保持する。 |
| Major | (c) 直ちに本番 6 本を `limiter: 0` で再走 | **却下**。壁温固定の未収束 run では、連成を戻した際の G-if・G-cons と精度を保証できない。[子 plan §6](/home/sano/work/forge-cht/plans/active/boundary-cht-conjugate-flat-plate.md:149) の全条件を維持し、まず今回の診断を完結させる。 |
| Minor | 「未達は `rms_roUx` のみ」 | **訂正を採用**。第 2 区間の `rms_roUy` も B/A = **0.109 > 0.1**（[判定出力:9](/home/sano/work/forge-cht/case/65.conjugate_flat_plate/AB_JUDGE_LIM_run0017_0018.txt:9)）。「最後の 2 区間」の条件ではこちらも未達。 |

結論: **(a) を選び、既存 A を固定対照として B の同一設定継続を追加 36000 step・累計 48000 step に限定して行い、延長区間の判定を事前登録する。**

第 1 仮説: 当該条件の停滞にはリミッタを含む再構成が関与し、B の残る `rms_roUx` は再構成変更後の緩和過渡である。確度: **中**。係数の時間変動そのものが真因という部分は未確認。

  根拠:
  
  - A = `case/65.conjugate_flat_plate/run_0017_lim_c1_n64_lim2/`、B = `case/65.conjugate_flat_plate/run_0018_lim_c1_n64_lim0/`。手元の YAML 差分は `space.limiter` のみ。関連するリミッタ実装・設定読み込みは、申告された `86115cb1` と現行で差分なし。
  - 保存判定の区間は再開後 0〜11999 step。A は **`NOT CONVERGED (stalled/plateau)`**、B は **`NOT CONVERGED (still converging)`**。B の `rms_roUx` 区間平均は **5.511e-7 → 9.067e-8 → 5.765e-8**。一方、A は約 **4.1e-8** で横ばい。
  - A の残差二乗和の後縁帯割合は **94.9〜98.7 %**。これは局在と係数変動の併存を示すが、因果の方向は示さない。
  
  反証条件: 延長後も `rms_roUx` が A の 1/10 を超えて横ばいになるなら、「残る未達は延長で抜ける過渡」という部分を棄却する。再構成の関与まで否定する結果ではない。

第 2 仮説: B は別の離散定常解へ移動しており、後縁の変動は減っても別領域の残差床が残る。確度: **低・未確認**。

判別 A/B: **比較因子は引き続き `space.limiter: 2 / 0` の 1 点**。A の既存最終区間を固定し、B だけを新しい `run_NNNN_*` へ FP64 のまま継続する。CFL・壁温・BC・物性・参照尺度は変更しない。

- **延長前に固定する区間**: B の累計 step `[36000,40000)`、`[40000,44000)`、`[44000,48000)`。比較対象は A の `[8000,12000)`。旧 12000 step の判定は上書きしない。
- **見る量**: 能動残差全列の区間平均、下流 slip の規格化 P 振幅。旧登録の「最後の 2 区間が A の 1/10 以下、3 区間で減少」を使う。`rms_roUx` の閾値は約 **4.095e-9**。
- **結果 A：全指標が登録条件を満たす** → 「同程度の停滞が無制限再構成でも残る」を棄却し、再構成依存を支持する。係数の切り替わり単独の証明にはしない。
- **結果 B：未達列が残り、3 区間の max/min ≤ 1.1** → 「累計 48000 step までの延長で過渡が抜ける」を棄却する。なお下降中なら引き続き判定不能で、成功扱いも自動延長もしない。
- **収束・準定常は別判定**: B の同一設定区間を継いで `check_convergence.py` を実行する。`limiter: 2` の履歴は混ぜず、未収束の A/B を `--from-floor` の参照にしない。登録尺度で P・界面 q/Tw の `check_quasisteady.py` の VERDICT も併記する。**PASS 単独で A/B 成功にはしない**。

やらない方がよいこと: **今は (b) の `venkatK` 調整と (c) の本番 6 本へ進まない。** `venkatK` 増加も制限の強さを変えるので、「時間的な切り替わりだけを抑える操作」ではない（[methods/limiter.md:83](/home/sano/work/forge-cht/methods/limiter.md:83)）。

`limiter: 0` は、将来 **C 専用の事後改訂した検証設定**として採用する余地はある。ただし現時点では未採用とする。発注元 §4.5 は物理条件の登録であり、これを維持しても数値レシピの変更は消えない。採用時は子 plan に変更理由・適用範囲を明記し、連成を有効にした 6 本すべてで従来の前提ゲート・主判定・格子差・固体効果を評価する。旧結果は「判定不能」のまま残す。

呼び出し側の前提への異議: **「B の残りは前縁帯 91〜93 %」を現在の残差位置とは読めない。** この割合は全 12000 step の二乗和で、上位節点表示は最後の 4000 step に限る（[ab_judge_lim.py:89](/home/sano/work/forge-cht/case/65.conjugate_flat_plate/ab_judge_lim.py:89)、[ab_extract.py:112](/home/sano/work/forge-cht/case/65.conjugate_flat_plate/ab_extract.py:112)）。実際、末尾の `res_roUx` 上位点は x/L = **0.064〜0.071** で登録前縁帯の外にある。対案は、領域別の**絶対二乗和と割合を同じ 4000 step 区間ごと**に出すこと。H2 の「衝撃なしだからリミッタの役割が薄い」も、本番採用の根拠には不足する。

不足情報: 手元には両 run の `residual_history.csv`・`ab_series.npz`・`res_*.h5`・準定常 VERDICT がない。`check_convergence.py` の今回の実行結果は両方 **`NO residual_history.csv`**。上記数値は保存された判定出力との照合であり、一次データからの再計算、NaN/Inf、初期場・固定壁温の同一性は独立確認できていない。呼び出し側でこれらを照合してから延長すること。run 索引は [case README](/home/sano/work/forge-cht/case/65.conjugate_flat_plate/README.md:9)。

**plan 未反映**。本回答は読み取り専用の諮問であり、反映先は呼び出し側の `plans/active/boundary-cht-conjugate-flat-plate.md` §4.4・§5.1。
