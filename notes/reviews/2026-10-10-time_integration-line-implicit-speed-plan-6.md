# codex レビュー: time_integration-line-implicit-speed (plan)

- **plan**: [`plans/active/time_integration-line-implicit-speed.md`](../../plans/active/time_integration-line-implicit-speed.md)
- **stage**: `plan`
- **date**: 2026-10-10
- **commit**: `d2a6644d` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `medium`, 1.7 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M4/m2
- **focus**: §6.16 (1 step の時間のふるい: FORGE_LINE_MAXLEN 100/64/48・nStepInner 3/2、部分被覆での LAYOUT2 のビット一致) と case/45.isobutane_m6_d155/scr.sh・scr_judge.py だけを、投入前の不具合 (起動できない・誤合格・設定が効いていないのに通る・GPU 専有でない・ディスク) と、ふるいの規則 (5 % 以上速い案だけ総時間へ) が誤った案を通す・落とす穴に絞って短く点検してほしい。部分被覆で LAYOUT2 の経路が point の節点を正しく扱うか (mesh.cpp:1090 の打ち切り、timeIntegration_d.cu の LP カーネル) も見てほしい
- **extra**: `case/45.isobutane_m6_d155/scr.sh`, `case/45.isobutane_m6_d155/scr_judge.py`, `case/45.isobutane_m6_d155/time_pairs_judge.py`, `case/45.isobutane_m6_d155/line_cmp_judge.py`, `solver_density_cuda/mesh/mesh.cpp`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

部分被覆時の point 節点を LP カーネルが誤って上書きする問題は見つかりませんでした。ただし、比較不合格の見逃しと計測条件の穴があるため、そのままの投入は勧めません。以下は静的確認と模擬入力による確認で、GPU 実測はしていません。

1. **Major — 部分被覆の照合に失敗しても M 系が進出できる。**  
   根拠: [scr.sh:23](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/scr.sh:23) は `c1` をログに書くだけで、[scr_judge.py:16](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/scr_judge.py:16) は比較結果を読みません。§6.16 の「FAIL なら M 系を進めない」が実装されていません。  
   **対案:** `cmp_judge_LAYOUT2.json` の `PASS` を M100/M64/M48 共通の必須条件にし、欠落・判定不能・判定器の異常終了も進出不可にする。

2. **Major — GPU 照会の失敗を「0 本」と記録できる。**  
   根拠: [scr.sh:35](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/scr.sh:35)。`nvidia-smi` が標準出力なしで失敗すると `wc -l` は `0` を保存します。`pipefail` はありますが、終了コードを検査せず続行します。また、投入前だけの確認では実行中の競合を検出できません。  
   **対案:** 照会成功とプロセス 0 本を別々に検査してから起動し、計測中も対象 GPU の他プロセスを監視する。競合・照会失敗のある測定は無効にする。

3. **Major — 性能判定が残差の NaN/Inf を受理する。**  
   根拠: [time_pairs_judge.py:19](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/time_pairs_judge.py:19) が確認する有限性は時間だけです。実際に `check()` へ「101〜999 の全行で `rms_ro=nan`、時間は 30 ms」という模擬ログを渡すと、正常な測定値として受理しました。[scr.sh:42](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/scr.sh:42) は最終場も検査せず削除します。  
   **対案:** 削除前に全残差列と最終場の有限性を確認し、異常な run は性能候補から除外する。短期ふるいに収束合格まで要求する必要はありません。

4. **Major — 「交互計測」になっておらず、順序効果で 5% を越え得る。**  
   根拠: [scr.sh:32](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/scr.sh:32) は毎巡 `B0→M100→M64→M48→S3→S2` の固定順です。[time_pairs_judge.py:27](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/time_pairs_judge.py:27) の B0 の幅では、毎巡繰り返す温度・クロック等の順序効果を除けません。仮に順序効果だけで B0 が毎回 32 ms、候補が毎回 30 ms なら、幅 0・差 −6.25% で通ります。  
   **対案:** 各候補を直近の B0 と組にし、AB/BA の順序を均衡させる。判別不能は「捨てる」でなく再測定対象にする。

5. **Minor — ディスク検査が開始時のホーム側だけ。**  
   根拠: [scr.sh:17](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/scr.sh:17)。出力先と別ファイルシステムなら無効で、`-BG` の切り上げ表示では実空きが 3 GiB 未満でも通ります。  
   **対案:** 各 `prep` 前に出力先の空きをバイト単位で検査し、コピー・入力 HDF5・同時に残る出力のピーク容量を確保する。

6. **Minor — 5% 未満を捨てても「総時間で効果なし」とは判定できない。**  
   根拠: [plan:452](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-implicit-speed.md:452)。M/S 系は必要 step 数も変えます。例えば単価 −4%・必要 step 数 −20% なら総時間は約 −23% ですが、このふるいでは落ちます。  
   **対案:** 5% は計算予算上の優先基準と明記し、未通過案は「総時間未評価」と残す。進出も速度改善の採用判定とは区別する。

部分被覆の実装については、`mesh.cpp:1086,1129,1137` で打ち切り先の `prev/next=-1` が保たれ、`timeIntegration_d.cu:831,1191` で point 解へ分岐します。LP は `line_cells` 内だけを処理し、比較側も `timeIntegration_d.cu:3278` で point の補正を先にコピーしています。この経路の扱いは整合しています。

**推奨:** 上記 1〜4 を優先して判定・計測処理を修正してから、§6.16 の短期ふるいを実施すること。`bash -n` は通りました。ファイル変更なし、plan 未反映です。

指摘数: Critical 0 / Major 4 / Minor 2
