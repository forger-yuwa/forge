# codex 諮問 (diagnose): settle-float-long-span

- **brief**: [`notes/reviews/briefs/2026-10-16-settle-float-long-span.md`](../../notes/reviews/briefs/2026-10-16-settle-float-long-span.md)
- **plan**: [`plans/active/tooling-settle-convergence-and-monitor.md`](../../plans/active/tooling-settle-convergence-and-monitor.md)
- **date**: 2026-10-16
- **commit**: `e311025c` (feature/card-local-mesh)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 不明 (`--requester-model` 未指定)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 2.9 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-16-settle-float-long-span-diagnose.codex.log` (codex: ok)
- **結論**: 事前登録に上記の実行条件を補い、`case/71.m55_combustion_wt/run_0017_ns_float_free130/` を途中 restart なし・停止規則なし・本段130,000 stepで1本回す。

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0 / Major 3 / Minor 1）

| 重大度 | 論点・採否 | 根拠と対案 |
|---|---|---|
| **Major** | #14 を省略して閉じる案は**却下** | `case/71.m55_combustion_wt/_runs/settle_evidence/v3b_replay/` の原本では、FP64 の `run_0012_ns_free` は **A**、32,000〜112,000 の41点が SETTLED。一方、float の `run_0015_ns_float_a`・`run_0016_ns_float_b` は **未成立**、32,000〜60,000 の15点だけ。60,000での精度間比較は、その後の時間変化を拘束しない。**新しい float の自由継続1本で埋める**。 |
| **Major** | 130,000 step・t_s ≤ 50,000・閾値固定は**採用。ただし実行条件を補う** | `settle_virtual_stop.py:99` の t_s 上限の既定は120,000。`:954` では `--criteria` 使用時に別の `--out` が必須。**`--span 80000 --ts-max 50000 --criteria … --out …` を事前登録**する。さらに `p5_all.sh:12` は run_0017 を受け付けず、`:18` のビルド既定は FP64。既存 `p5_ns.py` の `--precision float --stop-rule off --nsteps 130000` を使い、float バイナリの sha256・実効設定・初期場を固定して照合する。 |
| **Major** | 現行ツールの A だけで検証完了とする扱いは**却下** | `settle.py:1199`〜`:1218` の `residual_finite` は、途中の行欠けや必要な列の欠けを保証しない。plan #10 の既知の未修正項目。**新しい run の合格確定前に #10 の修正・回帰試験を通し、完全な履歴で再生する**。残差の振幅条件を戻す必要はない。今回の既存 float 2本は、実データで `outer_end` の欠測・重複0、全 `rms_*` の非有限0を確認した。 |
| **Minor** | restart で延ばさない判断は**採用**。「必ず跳ねる」は**要再検証** | `procedures/recommended-settings.md:173` は SST restart の非忠実性を記載するが、この2本の再開直後の跳ねを示す実測は提示されていない。理由は**「同じ連続軌道を検証するため」**とする。既存 plan の `tooling-settle-convergence-and-monitor.md:172` も、restart 系列の継ぎ足しを除外している。 |

結論: **事前登録に上記の実行条件を補い、`case/71.m55_combustion_wt/run_0017_ns_float_free130/` を途中 restart なし・停止規則なし・本段130,000 stepで1本回す。**

第 1 仮説: float でも監視対象量は t_s 後80,000 stepで固定許容内に留まり、現在の未成立は観測期間の不足による。確度: **中**。  
　根拠: v3b 再生の `run_0015` / `run_0016` は、観測済み15点すべて SETTLED。t_s からの差／許容の最大は、出口コア M **0.06853 / 0.04981**、線 **0.14077 / 0.11959**、流れ角 **0.00827 / 0.00709**。ただし原本の **`VERDICT: 未成立`** は維持する。  
　反証条件: 新しい連続軌道の対象区間で、累積差／許容が1を超える、または全条件の持続が破れる。

第 2 仮説: float の遅いドリフトまたは間欠的な変動が60,000 step以降に現れる。確度: **低・未確認**。現状の最大比が FP64 より大きいことだけでは、増大し続ける証拠にならない。

判別 A/B: **追加する検証は自由継続1本だけ**。数値条件を固定し、最初の停止候補 t_s を再生で求め、[t_s, t_s＋80,000] の2,000 stepごとの41点を評価する。

- **A**: 累積差／許容 ≤ 1、全点 SETTLED、欠測なし。この条件・観測間隔・検証長について第2仮説を退け、#14を完了できる。
- **B**: 累積超過または持続違反。「観測期間不足だけ」を退ける。**累積超過と持続違反は別に報告**し、後者だけで「目的量を早く止めすぎた」と断定しない。
- **未成立**: t_s > 50,000、候補なし、または反例はないがデータ不足。合格するまで上限を延ばさない。t_s は既存 run の32,000に固定しない。

やらない方がよいこと: 結果を見て閾値を緩める、restart 後の系列を連続軌道として継ぐ、A を局所壁熱流・δ_E・別条件のノズルまで一般化する。

呼び出し側の前提への異議: 「60,000まで精度間で近い」と「停止後80,000 stepで変化しない」は別の検証である。また、`v3_replay` は壁熱流を含み、`v3b_replay` は除外後なので、**プリセット名だけでなく基準ファイルとハッシュを固定**する。

不足情報: 手元の float の写しには `nozzle.h5` と全場がなく、抽出・再生全体の独立再実行はしていない。今回は保存済みの判定・基準ハッシュ・格子点別記録と、残差実データを照合した。これは次の検証を選ぶ判断を妨げない。

証拠: 足りる
