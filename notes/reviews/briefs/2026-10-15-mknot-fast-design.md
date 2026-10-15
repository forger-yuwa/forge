# 諮問ブリーフ: M_knot の自動選択を「選ぶ網」と「作る網」に分けてよいか (2026-10-15)

- 依頼者: 主セッション `M6ノズルの設計` (エスカレーション条件 1: §4・§6 を新規に書く)
- plan: `plans/active/tooling-nozzle-mknot-fast.md` (§3 実測、§4.1〜§4.3、§6 V1〜V4)。前身 `plans/accepted/tooling-nozzle-mknot-auto.md` (規則・ゲート・試験 (n)(o))
- コード: `design/forge_design/evaluate/mknot_select.py` (`select_mknot`・`resolution_gate`・`evaluate_mknot`)、`design/forge_design/intake/wt_card.py` (`mk_select_first`・`mk_gate`・`mk_finalize`、`mknot_workers`・`_mknot_mem_check`)、`design/forge_design/evaluate/runner_axismach.py` (`from_length`、`L_total_tol`)。HEAD 9441b630 (main 4b1708ac + plan)

## 観測事実

- 所要: 要求カード 1 枚の M_knot の自動選択が 45〜97 分 (AWS、並列 1/1)。大半が 2N = 4800 の評価。
- 1 評価 (AWS A、2026-10-15): 600 で約 3 s・0.13 GB、1200 で 7〜12 s、2400 で 9〜23 s・0.8〜1.3 GB、4800 で 45〜81 s・3.9〜5.6 GB (表は plan §3)。
- プロファイル (2400): 6〜7 割が逆 MOC の内部点 `moc_kernel.interior_vec`。
- 600 では `L_total_tol` をカードの値 (0.005〜0.020 r_t) のままだと全点 `from_length 未収束`。case/70 で 0.1 r_t にすると境目 1.97 (不成立)/1.98 (成立)、最小 2.04。2400 は境目 1.971875・最小 2.07、4800 は境目 1.970625・最小 2.05。θ_max の値は 600 で約 0.023° 高い。2400 と 4800 の θ_max の差 0.001〜0.002°。谷 2.00〜2.10 の中の差は 0.006°。
- (追加の測定を実施中: 600・1200 を `L_total_tol` 0.05 で、case/71 の 2400 の掃引。諮問の時点で結果が無ければ、それを前提にしない判断をお願いしたい)
- 生産の 2400 は壁 (出口の特性線の交点・壁の末尾) の精度のために決まった値 (1200 → 2400 で x_F が 0.05 r_t 動く)。M6 (case/42) の M_K の掃引は 600 で行っていた。

## 期待値と出典

- ユーザ (2026-10-15): 「劇的に早くしてほしい」「そもそもが、もっと荒くていいのでは」「規則というかメッシュのつもり」。
- 前身の plan の規則 (制約・ΔM_safe 0.05・同点の規則・N/2N のゲート・窓の外の文脈・settled) は変えない。試験 (n)(o) (2N にだけある境目) は通す。

## 案 (plan §4)

1. 選ぶ網 N_s = N_p/4 (600) と 2N_s (1200) で今の `select_mknot`・`resolution_gate` をそのまま回す。長さの許容差は L_s = max(カード, 0.05 r_t) (0.05 で未収束なら 0.1)。
2. 作る網 N_p (2400、カードの許容差) で `chosen_N ± 0.05` を 0.01 刻み 11 点評価し、全点成立を確かめる。1 点でも不成立なら、その組だけ今の方式 (N_p/2N_p) に落ちる。2N_p は既定では評価しない。
3. r_t の固定点の最終の組でも 1・2 をやり直す (最初の r_t の結果を使い回さない)。
4. 並列はメモリ・CPU から自動 (上限 4、環境変数で下げられる)。
5. MOC の JIT 化は後回し (生産の壁の数値が変わるので別の受入)。

## 問い

1. 案 1・2 の「選ぶ網で選び、作る網では近傍だけ確かめる」で、前身の plan のゲートが守っていたもの (2N にだけある境目・選択の網への依存) を失わないか。失うなら最小で何を足すか (例: 作る網の確認を ±ΔM_safe より広く、2N_p を選んだ点 1 点だけ評価、など)。
2. 選ぶ網で長さの許容差を緩めることの影響: θ_max の順位と境目の位置への効きをどう扱うか。L_s の決め方 (固定値か、N_s の x_F の離散化の大きさから決めるか)。
3. §6 V3 の合格条件 (θ_max の差 ≤ 0.005°、旧の評価表の全部の境目から 0.05 以上、フォールバックなし、寸法ゲート) は妥当か。0.005° の根拠 (2400 と 4800 の差 0.001〜0.002° の約 3 倍) で足りるか。
4. 見落としている危険 (例: 600 では成立するが 2400 で不成立になる M_K の型、from_length の許容差を緩めた評価が別の分枝に入る、など)。

## 読んでよいファイル

- `plans/active/tooling-nozzle-mknot-fast.md`、`plans/accepted/tooling-nozzle-mknot-auto.md`
- `design/forge_design/evaluate/mknot_select.py`、`design/forge_design/intake/wt_card.py`、`design/forge_design/evaluate/runner_axismach.py`
- `design/tests/run_mknot_select_tests.py`
- `notes/reviews/2026-10-14-mknot-auto-result-diagnose.md`、`notes/reviews/2026-10-14-tooling-nozzle-mknot-auto-result.md`
