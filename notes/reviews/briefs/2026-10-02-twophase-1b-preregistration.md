# 諮問ブリーフ: 凝縮二相拡散 #1b (case/16 A/B) の事前登録 — 判定条件の確定 (2026-10-02)

作業ツリー `/home/sano/work/forge-species` (HEAD は git log -1)。plan `plans/active/condensation-two-phase-transport.md` §5.1 #1b・#4e・#4f (受入済み、#1b の条件に両 run の `[twophase-audit] VERDICT: PASS`)。準備物: `case/16.nozzle_wys/prepare_twophase_ab.sh`・`twophase_ab_series.py`・`twophase_ab_compare.py` (docstring に定義)。設計メモ `notes/investigations/2026-10-02-twophase-diffusion-kernel-design.md` §14.4・§15。AGENTS.md エスカレーション 1・3 (§6 の合否を新規に書く; 結果を見てから作らない)。

## 問い

下の事前登録文案の [ ] を埋め、穴があれば直して確定してほしい。特に:
1. 収束判定: B は `rms_roYv` 列が増え `--from-floor run_0482` は列不一致で必ず FAIL。A も run_0482 以後の既定変更 (limiterScaled/venkatK 0.05、ljSource、#3/#3b の D_mix、内蔵 H2O) で run_0482 の不動点ではなく再開の過渡から始まる。ピークから 3 桁落ちは過渡が小さいと届かない可能性。判定区間・方式を 1 つ。
2. 監査 `[twophase-audit]` の基準 max|r| ≤ max(1e-7·r0, 6ε·max A) は r0 が小さいと実質 6ε·max A。check_convergence PASS でも監査 NOT CONVERGED が続きうる — どう扱うか。
3. 「補正が収束時 0」の定義: step 1 で Qcut 2.3e-44 (非正規数)、モーメント射影 788 節点、再正規化の相対量 6.0e-5 (#2 で予想された常時 WARN)。各しきい値。
4. 定常判定の許容: 既定 (drift 5 %/osc 10 %) は A/B の差に粗い。`Tw_mean_x10_K` は run_0482 で 20k→48k step に 283.55→283.19 K と動き続けた。許容と、延長の規則 (再開でつなぐと limiter 基準値が取り直され作用素が変わる → 同じ初期場から nStepOuter を増やして回し直すか、`limiterRoRef 0.5279554756`/`limiterPRef 43069.4137`/`limiterARef 329.7385022` を A・B に固定)。
5. 交絡: B は拡散作用素と同時に蒸気・液・Q の前処理も DPLUR→点対角。両者収束なら場の比較は作用素の差として読める — この読み方でよいか。

## 事前登録文案 (implementer 作成)

- 共通初期場: run_0482 `res_48000.h5` を `prepare_twophase_ab.sh` で 2 run に index コピー (流れ・乱流・化学種・液・Q2/Q1/Q0; 液とモーメントは空データセットを作ってから写す — 14 量 SRC とビット一致を確認済み)、species_db から H2O を外す、`--force-species`、起動は `FORGE_ALLOW_UNVERIFIED_SPECIES=1` (A・B 同じ)。
- A = run_0482 の config (nStepOuter・outStepInterval のみ変更)、B = A + `condensation.condTwoPhaseDiffusion: 1` (relax 既定 1)。同じバイナリ。初回 4000 step・200 step ごと出力 (1 run 約 0.4 GB)。
- 報告量: `onset_c_g1e4_mm, g_exit_mw, pw_mean_x10, pw42, pw52, Tw_mean_x10_K, T_cond0_vmean_K` (壁は断熱で熱流束 0 → `Tw_mean_x10_K` で代替)。
- 使う前の条件: (a) 両 run `check_convergence` PASS [区間・方式]、(b) 報告量すべて `check_quasisteady --series-csv --drift [ ] --osc [ ]` で STEADY、(c) B の終了時 `[twophase-audit] VERDICT: PASS`、(d) B の末尾 10 % で θ<1 のセル 0、Qcut [ ]、vround [ ]、`[cond-corr]` の射影・負値床・増分制限・floor [ ]、(e) NaN なし。
- 延長: (a)/(b) 未達なら同じ共通初期場から nStepOuter を 2 倍で A・B とも回し直し (または基準値固定で継続)、上限 [例 32000]。上限で未達、または check_convergence PASS でも (c) NOT CONVERGED が 2 回続けば許容を緩めず上位に諮る。
- 記録: 共通凝縮域 (両者 g>1e-6) の体積重み p95 |ΔT| (p50/p99/max/符号付き平均)、報告量の B−A と相対差、片方だけ凝縮の体積。δṪ は省略 (EOS 接線をソルバと同じに組む C++ ツールが要り重い)。合否にしない (1 K は大小の分類)。

## 禁止事項 (厳守)

- ファイルを変更しない。`*.log`, `residual_history.csv`, `res_*.h5`, `*.vtu`, `plans/README.md` を読まない。推奨は問いごとに 1 つ。根拠は `ファイル:行` か本ブリーフ。

## 読んでよいもの

- plan 全体、設計メモ §14–§15、`case/16.nozzle_wys/prepare_twophase_ab.sh`・`twophase_ab_series.py`・`twophase_ab_compare.py`、`solver_density_cuda/tools/check_convergence.py` (`--segment`・`--from-floor` 付近)、`check_quasisteady.py` (`--series-csv` 付近)、`procedures/recommended-settings.md` の凝縮の節
