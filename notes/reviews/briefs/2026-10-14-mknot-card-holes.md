# 諮問ブリーフ: M_knot の自動選択をカードに組み込んだ後の穴 (plan tooling-nozzle-mknot-auto §4.2〜§4.4、§9 の 2026-10-14 の最後の項目)

AGENTS.md 条件 1 (§4 の方針の穴) と 7 (result の解釈の前)。判断役のスイッチが codex。実装 commit 1e178365。

## 観測事実
- 実装担当の返却 (plan §9 の最後の項目に要約): V2(k)・V1a・V6・V3・V4 final 合格。V4 explore は未了 (手元の PC がメモリ不足で落ちた; 以後の重い計算は AWS)。
- 穴 (1)〜(7) は plan §9 の最後の項目のとおり。(1) は合成の試験 V2(k) で、ゲートが入力 2.05 から −0.1 の窓の端 1.95 を選び直して PASS を出した。
- `resolution_gate` は `design/forge_design/evaluate/mknot_select.py:602` 付近、`wt_card` の流れは `design/forge_design/intake/wt_card.py:1528` 付近。

## 私の案
(1) ゲートの判定に「選び直した N・2N の値が詰めの窓の端に無い」を足し、端なら「未確定」(例外)。カードの流れでは差 > 0.01 の組み直しで窓が動くので、その後のゲートで判定される。
(2)(3)(4) は実装担当の選択を採る。(5) 所要は AWS で回す前提で受け入れ、短縮は別の項目。(6) 台本は数値の M_knot を明示して直す。
(7) YAML には要約 (規則の名前・選んだ値・距離・境目・ゲートの判定・card_log の sha256) だけを置き、掃引の表は card_log の JSON だけに置く。Euler の写しには要約だけ。

## 問い
**私の案のうち、間違っている・足りないものを 1 つ挙げ、直し方を 1 つ推奨してほしい** (無ければ「無し」)。特に (1) の端の判定が、平らな谷 (θ_max の差 1e-3° 以内が広い) で正しい設計を「未確定」にしてしまわないか。

## 禁止事項 (厳守)
- ファイルを変更しない (read-only サンドボックスで動いている)。
- **`*.log`, `residual_history*.csv`, `res_*.h5`, `*.vtu`, `*.msh`, `plans/README.md`, `case/70.m5_air_wt_demo/_mknot_*/*.json` を読まない** (巨大)。
- 下に列挙した範囲以外のファイル読みをしない。grep は可。
- 推奨は 1 つに絞る。根拠は `ファイル:行` か提示した数値で示す。

## 読んでよいもの
- `sed -n '54,100p' plans/active/tooling-nozzle-mknot-auto.md` と plan の末尾 30 行 (`tail -30`)
- `sed -n '590,740p' design/forge_design/evaluate/mknot_select.py`
