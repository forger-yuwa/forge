# 諮問ブリーフ: 熱物性の経路の止まる先の試験 (§4.97、判定不能) の解釈と、ユーザの問いへの答え

- **plan**: `plans/active/tooling-nozzle-core-grid.md` (§4.94〜§4.97)
- **date**: 2026-10-13
- **諮問の条件**: AGENTS.md「エスカレーション条件」3 (判定不能)・7 (解釈の確定の前)
- **ユーザの問い**: 「既定の熱物性の評価の経路によって、収束先が異なる事象が生じていた可能性はある?」

## 禁止事項 (厳守)

- ファイルを変更しない (read-only サンドボックスで動いている)。
- **`*.log`, `residual_history.csv`, `res_*.h5`, `*.vtu`, `*.f64`, `plans/README.md` を読まない** (巨大)。
- 下に列挙した `sed -n 'A,Bp' <file>` 以外のファイル読みをしない。grep は可。
- 両論併記で逃げず、**推奨は 1 つに絞る**。根拠は `ファイル:行` か下に示した数値で示す。

## 読んでよいもの

```
grep -n '^### 4\.9[4-7]' plans/active/tooling-nozzle-core-grid.md    (§4.94〜§4.97 を sed で 25 行ずつ)
python3 -c "import json;j=json.load(open('case/45.isobutane_m6_d155/_band_ab/core_grid/tfs_judge.json'));print(j['qualification'],j['result']);[print(q,{w:{k:(round(v,4) if isinstance(v,float) else v) for k,v in d.items() if k in ('D','F')} for w,d in r['windows'].items()},r.get('H'),r.get('verdict')) for q,r in j['q'].items()]"
sed -n '1,30p' case/45.isobutane_m6_d155/cg_tfs_judge.py
```

## 1. 結果 (§4.97)

- Gref の B0 の停滞の場 (`run_0608` res_60000) から、旧バイナリ 7e95e5bb で B0 のまま `thermoFloat` 1 (T1a・T1b) / 0 (T0a・T0b)、各 6 万 step を 1 本ずつ。入力は 4 本で同じ sha。
- θ_r[65,75] の窓平均の point の止まった値 (`thermoFloat: 1` で求めた値、μ = 0.113931) からの差: W1 (20k〜40k] T1a +0.569・T1b +0.608・T0a +0.317・T0b +0.213 %、W2 (40k〜60k] T1a +0.593・T1b +0.719・T0a +0.285・T0b +0.226 %。D = −0.323・−0.401 %、F = 0.104・0.126 %、U = 0.067 %、H = 0.192 %。θ_r[84,94] D −0.286・−0.365 % (H 0.197)、Q_w D −0.219・−0.305 % (H 0.215)。3 量・2 窓のすべてで T0 の 2 本 < T1 の 2 本 (重ならない)。各 run の各窓は `classify_series` で STEADY。
- 登録の資格「主指標で H ≤ 0.025 %」を満たさず判定不能。量ごとの判定式 (|D| − H ≥ 0.1 % かつ符号が同じ) は θ_r[65,75] で満たす (0.13・0.21 %)。
- `check_convergence`: 4 本とも NOT CONVERGED (プラトー)。

## 2. 判断してほしいこと

推奨を 1 つに絞り:
1. §4.97 の判定不能のまま、ユーザの問いに今どこまで答えてよいか (言える・言えない)。2〜4 文で、ユーザに伝える言い方で。
2. 次の一手が要るなら 1 つ (事前の読み方つき)。候補 (縛られなくてよい): (a) 複製を増やす・期間を延ばして H を下げる、(b) point も `thermoFloat: 0` で回して「両方の止まる先が経路でどう動くか」を見る、(c) ここで区切って記録だけ残す。
3. 実行間差 (同じ入力・同じ設定で θ_r の窓平均が 0.04〜0.13 % 違う) をどう記録すべきか。
