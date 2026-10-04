# 諮問ブリーフ: 凝縮二相拡散カーネル (#4) — 未確定 5 点の決定 (2026-10-02)

作業ツリー `/home/sano/work/forge-species` (`feature/species-transport`)。plan: `plans/active/condensation-two-phase-transport.md` §4.2・§6・§5.1 #1/#4。設計メモ全文: `notes/investigations/2026-10-02-twophase-diffusion-kernel-design.md`。AGENTS.md エスカレーション 6 (カーネル編集の前)。


**問い**: plan `condensation-two-phase-transport` §4.2 の二相拡散カーネル (§5.1 #4) に着手してよいか。§4 の未確定 1–5 をどう決めるか。
特に (2) 相変化ソース・更新クランプと蒸気/液の増分の関係、(3) 総水分の FCT。

**観測事実** (ホスト参照実装のみ。forge の run・カーネルは無い):
- `tests/unit/test_twophase_diffusion_harness.py` (float64 参照, 2026-10-02, HEAD `33520773` + 未 commit のハーネス) が §6 の単体条件を全件 PASS (§3 の表)。
- 3 セル判別: A (点対角 1 回) 総液量 +3.3333333 % / +2.8571429 %、B (解き切った BE) +2.2e-9 (1000 更新)。上位検算 (+3.3333346 % / +2.8571441 %, float32) と一致。
- 2 セル非負: 改訂前の更新で ρv −0.016667 (上位検算と一致)、新 (蒸気/液変数) で 0。
- 補正の面 z: D 比 100・段差 1 の 2 セルで、算術平均は ρv −6.2e-3、風上は +3.2e-4 (点対角 1 回)。通常の組成のランダム試験では両者とも負にならない。
- float32: (i) の許容は dx 1 mm・Δg 1e-4/面で 98 倍超過。保存は停止の丸め床 16ε で 7.2e-6、4ε で 5.0e-7 (1000 更新)。
- エネルギー: 接線投影 6e-9、潜熱流束なし・エネルギーだけ潜熱はそれぞれ 27・33 倍ずれる。面の L を T^(n+1) で評価すると節の近くで 4.4 % (J_l·∇L)。

**期待値と出典**: plan §6 の事前固定の数値 (単体・保存・非負・エネルギー)。判別規則は diagnose 2026-09-27 (`notes/reviews/2026-09-27-condensation-transport-redesign-diagnose.md`)。

**再現条件**: `python3 solver_density_cuda/tests/unit/test_twophase_diffusion_harness.py` (約 9 分, `--only S3` などで個別)。CFD・GPU 不要。

**実施済みの操作と結果**: methods §7c に §4.2 の仕様を記述 (未実装と明記)。ハーネスで決めた解釈は本メモ §3 末尾 (二元 D の混合平均、解き切りの定義と丸め床、1D の前処理、熱伝導なし)。
現行コードの該当箇所は本メモ §1。

**仮説** (検証していない):
- H1: 契約 (保存形 BE の固定点 + 点対角の前処理 + 蒸気/液変数) は、輸送だけなら §6 を満たす (ハーネスが支持)。
- H2: 実装上の主なリスクは輸送ではなく、(a) 相変化ソースと θ・床が蒸気の増分にどう入るか、(b) 水の commit を液の後に遅らせる順序変更 (`main.cpp:2127-2152`, `:2355-2385`) と
  coupling 2 との非互換、(c) dual-time で総水分が FCT を通らないこと。
- H3: 補正の面 z は風上にすると非負が構造で保たれるが、算術平均 (現行) からの変更は g=0 のビット一致を崩さない範囲 (TP carrier 凝縮の新カーネル内だけ) に限られる。

**呼び出し側の前提で疑わしいもの**: 「解き切った」の定義 (相対 1e-7 + 丸め床) と、float32 の受け入れ許容を float64 の条件から分けるべきか。

## 禁止事項 (厳守)

- ファイルを変更しない。`*.log`, `residual_history.csv`, `res_*.h5`, `*.vtu`, `plans/README.md` を読まない。推奨は論点ごとに 1 つ。根拠は `ファイル:行` か本ブリーフ。

## 読んでよいもの

- `notes/investigations/2026-10-02-twophase-diffusion-kernel-design.md` 全体、`plans/active/condensation-two-phase-transport.md` 全体、`methods/condensation.md` の §7c (grep '7c')
- `solver_density_cuda/tests/unit/test_twophase_diffusion_harness.py`
- メモが `ファイル:行` で挙げた範囲 (speciesTransport_d.cu・condensationUpdateLimiter_d.cuh・main.cpp・passiveFct_d.cuh)
