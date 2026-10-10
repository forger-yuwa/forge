# 諮問: g3 継続場の EOS ビット比較 — 全配列 0 バイト差だが v2 で INVALID (遠方境界の ghost の roe が NaN) (2026-10-11)

関連 plan: `plans/active/tooling-sern-te-wake-grid.md` §5.1 #2「規則の改訂と追加の試験」(4) と #8 (ghost の初期化の欠落)。判断の記録 `notes/reviews/2026-10-08-eos-dump-rule-diagnose.md` (v2 = 読む入力の pre・全書き込み先の post の有限性、壁・ghost を除外しない、差分 0 バイト)。エスカレーション条件 3。
記録の写し `notes/investigations/2026-10-08-eos-dump/g3/` (SUMMARY.json・cmp_*.txt/json・prepare.log・INPUT.json)。

## 実施 (AWS-B、2026-10-10 UTC)
- 入力: `case/46.sern_design/run_1079_tewake_B10_m10` の res_20000 を `restart_field` (9 量ビット一致) で複製、nStepOuter 1・`output.floorEvents: 1`。harness `eos_ab_dir.py run` (旧 2 回・新 2 回、5 組)。
- バイナリ: 新 = HEAD a80600d0 のクリーンビルド (sha256 61cf8f4f…)、旧 = 9f35e3e7 + 同じフック/再生パッチ (ea2d101d…、eosDump.cpp/hpp の sha256 が新と一致 dd37b919…/f0caddc4…)。RelWithDebInfo・sm_86。

## 結果
- 5 組すべて: **EOS の全配列 post 差分 0 / 175,204,008 バイト、pre 差分 0、その他の配列も pre/post 差分 0**。新版の floor_events: init 行・step 1 の eos 行とも温度・密度・圧力の床事象 0、床近傍 0。
- **VERDICT_v2 INVALID (全組)**: 「EOS が読む入力 roe・T の pre に非有限」「書き込み先 T・P・Ht の post に非有限」(各 51,143、A・B とも同数)。
- 位置の特定 (新版のダンプを 1 回残して調べた): **非有限 51,143 点はすべて ghost** (index ≥ nCells = 1,920,103、nCells_all = 2,085,762、ghost 165,659)。実節点の非有限は 0。該当 ghost は ro = 1e-4 (= roMin)、pre の roe = NaN。gamma・cp の pre の NaN は ghost 全 165,659 点 (2D でも見た書き込み専用の配列)。
- 51,143 は **BCONDS/10 (bcondKind `farfield`、`side_far`) の境界面の数とちょうど一致**。他の境界 (inlet・outflow・wall_isothermal・slip) の境界面数とは一致しない。2D の比較 (TP) には farfield が無く、非有限は書き込み専用の gamma・cp だけだった。
- 推定 (未検証): node モードの farfield は境界節点の値で流束を作り (`farfieldFlux_d.inc.cuh` に再構成なし — 実装担当の確認)、ghost の値を埋めないので、初期化の EOS で作られた NaN が残る。実節点は 20000 step の run を通じて有限 (detectNaN 無し、GATES PASS)。

## 問い
1. g3 の比較を、EOS の数値非干渉の判定として「合格」と読んでよいか。v2 の有限性の要求 (読む入力の pre) に ghost の遠方境界の NaN が当たったのは、結果を見た後の規則の変更になる。認めるなら条件 (例: 実節点の全配列が有限・ghost を含む全配列が NaN のビット列まで一致・非有限が特定の境界の ghost に限られることを index で確認) を書く。認めないなら、何を追加で示せばよいか (例: farfield の ghost を有限値で埋めた入力で再比較する — 再生機能で pre の ghost を書き換えれば出力専用の範囲で作れる)。
2. 遠方境界の ghost の NaN を #8 (ghost の初期化の欠落) に含めるか、別件にするか。g3/g4 の時間積分 (#4) を止める理由になるか (実節点は 20000 step 有限、run_1079 は GATES PASS)。
