# 諮問ブリーフ: 推力ノズル (③ベル・④デュアルベル・⑤SERN) の現状方針の横断レビュー

- 日付: 2026-10-11
- 依頼元: 主セッション (Opus)。ユーザ指示「現状方針、diagnostician にも聞いて」(切り替えスイッチは codex だが、ユーザがその場で diagnostician を指定)
- 作業ツリー: `/home/sano/work/forge-thruster` (ブランチ `feature/thruster-nozzle-design`)。**ここのファイルを読む** (`/home/sano/work/forge` は別ブランチで古い)
- 段階: 計画は書けた。実装はベルの共通評価器の TP 化 (§5.1 #2a) を implementer が実施中。他は未着手

## 1. 経緯 (この日の対話で決まったこと)

ユーザの問い: 「推力装置向けにベル・SERN・デュアルベルを開発する。特性曲線法でできた形状をベースとするか、形状そのものを動かすか、悩み中」。
対話で決まった方針は親 plan `plans/active/tooling-nozzle-design-tool.md` §3 の確定事項 5〜15 (46〜57 行)。要約:

- 5〜7: スラスタは海面作動を含む / ③④共通で外気圧を与えて剥離を評価 / ノズル剥離の文献データを検証に
- 8: 全機種 thermally perfect (MOC と CFD の両方)
- 9: ③は「特性曲線で辿るモード (TP の Rao)」と「形状を直接振るモード」を作って比べる
- 10: ③の基準ケース LOX/LH2・10 MPa・O/F 6・ε 50・凍結 TP・断熱
- 11: ③の壁は両モードとも C2 (スロートの曲率 1 つ)
- 12〜13: 剥離評価 (外部領域・DLR TC1A 検証) は**別セッションの担当**で、設計側は使える前提で進める
- 14: ④は高度で指定 (h1 海面モード、h2 = 50 km、使えなくてよい高度帯)、遷移高度は定義せずモードの保証を制約に
- 15: 剥離・遷移の判定は公知の試験データ (経験式) を正、CFD は検証のみ (CFD で係数を較正しない)
- MOC は①風洞で検証済みの配列版 `interior_vec` (analytic・converge・対の分類) に乗せる (同日ユーザ指示)

子 plan は 3 本。いずれも codex の diagnose と plan 段のレビューを受け、指摘は全件採用で改訂済み:

| plan | 内容 | codex |
| --- | --- | --- |
| `plans/active/tooling-nozzle-bell-two-mode.md` | ③の 2 モード比較 (TP Rao と C2 B-spline)、共通評価器、比較の不確かさ | diagnose NO-GO → 改訂、plan GO-with-changes (M8) → 改訂 |
| `plans/active/tooling-nozzle-dual-bell.md` | ④ 高度指定、全長 CP 延長、方向別の遷移閾値、保証の状態、P6.2 の作り直し | diagnose GO-with-changes (M5)、plan GO-with-changes (M8) → 改訂 |
| `plans/active/tooling-nozzle-ambient-separation.md` | 外部領域と剥離評価 (別セッション向け) | diagnose NO-GO → 改訂。外気 BC の方式 §8-1 は未決 |

⑤ SERN は別セッションが `feature/sern-design` で進めている (`plans/active/tooling-nozzle-sern-chain.md` は本ツリーのものが古い可能性あり。本諮問では親 plan の扱いだけ見ればよい)。

## 2. 問い (横断的な判断。個々の plan の細部は codex が見ているので重複不要)

1. **ユーザの元の問い (MOC ベースか形状を直接動かすか) に、この計画群は答えられるか**。③の 2 モード比較の設計 (`tooling-nozzle-bell-two-mode.md` §4.5・§6 V3) で、比較結果がどちらに出ても④・⑤の手法選択に使える形になっているか。比べ方に構造的な偏り (M は 1 変数で S1 ゲートを通す、D は 4 変数で MOO) が残っていないか。
2. **順序と費用**: 計画が保証・ゲート・検証で重くなっている。最初に価値の出る成果 (ユーザが判断に使えるもの) までの経路が長すぎないか。削る・後回しにするなら何か。1 つに絞って推奨してほしい。
3. **plan 間の整合**: 推力の定義 (CV の 2 会計)、TP の共通契約、C2 と変曲点の例外、MOC の経路、作動点の書き方 (高度のリスト) が 3 本で食い違っていないか。
4. **見落としている大きな危険** があれば 1〜2 個 (例: 剥離評価が別セッションで遅れたときの影響、TP の Rao の導出が成立しない場合の代替)。

推奨は 1 つに絞り、根拠は `ファイル:行` で示してほしい。決まったこと (§1) を覆す提案は、覆すべき強い根拠があるときだけ、その根拠と一緒に。

## 3. 読んでよいもの

- `plans/active/tooling-nozzle-design-tool.md` (全体 344 行。特に 37〜60 行 §3、270〜300 行 §5.1)
- `plans/active/tooling-nozzle-bell-two-mode.md` (161 行)
- `plans/active/tooling-nozzle-dual-bell.md` (158 行)
- `plans/active/tooling-nozzle-ambient-separation.md` (160 行)
- codex の記録: `notes/reviews/2026-10-11-thruster-nozzle-child-plans-diagnose.md`、`notes/reviews/2026-10-11-tooling-nozzle-bell-two-mode-plan.md`、`notes/reviews/2026-10-11-dual-bell-plan-diagnose.md`、`notes/reviews/2026-10-11-tooling-nozzle-dual-bell-plan.md`
- Phase 2 の実績: `plans/accepted/tooling-nozzle-moo-loop.md` (変更ログ)
- 必要ならコードを grep・部分読み (`design/forge_design/`)

## 4. 禁止事項

- ファイルを変更しない。
- `*.log`, `residual_history.csv`, `res_*.h5`, `*.vtu`, `*.pdf`, `plans/README.md` を読まない。
- 両論併記で逃げず、推奨は 1 つに絞る。
