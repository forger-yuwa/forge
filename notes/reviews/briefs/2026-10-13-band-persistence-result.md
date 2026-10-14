# 諮問ブリーフ: 持続 K = 4 の検証 (S2b) が FAIL・確認不能 — 帯選択の頑健化の方針を決め直す (plan tooling-nozzle-cooled-wall-ns-settling §4.2b・§6 S2b)

AGENTS.md 条件 3 (事前に書いた判定が FAIL・判定不能)。判断役のスイッチが codex。前回の諮問 `notes/reviews/2026-10-13-band-persistence-diagnose.md` は K = 4 を採用 (条件つき) としていた。

## 観測事実 (`case/70.m5_air_wt_demo/band_persist_verify.py`・`_band_persist_verify/verify.json`、forge 0 step)

- 実装: `edge_band_positions(persist)`・`deltastar_from_core_matched_euler(edge_persist)`。persist 1 は旧の抽出とバイト同一 (3 ファイル)、persist 4 は測定用の差し替えとバイト同一。
  **作業ツリーでは既定を 1 に戻してある (未 commit)**。
- (a) 窓をまたぐ再現性 (run_0010b、W1 = 通算 70k〜90k の 5 枚、W2 = 100k〜120k の 5 枚、Euler run_0012b、使用域 [0, 48.03]):

| persist | 帯の再選択による差 = (独立の帯の W1−W2) − (共通帯の W1−W2) | 独立の帯の W1−W2 | 共通帯の W1−W2 (= 場の変化) | 判定 (≤ 1e-4) |
| --- | --- | --- | --- | --- |
| 4 | 1.42e-4 (x 4.91; x 3〜10 に 1.2〜1.4e-4) | 2.27e-4 | 9.5e-5 | FAIL |
| 1 | **6.1e-5** (x 11.14) | 1.60e-4 | 9.9e-5 | PASS (記録として回した) |

- (b) 非劣化 (case/45 の 4 対、persist 4 − 1): [0, x_e − 3] の最大 3.45e-4・3.84e-4・6.90e-4・6.02e-4 r_t (全部 FAIL)、出口の端 +0.063 %・+0.048 %・−0.24 %・−0.11 %。
  使用域の y* の NaN は K1・K4 ともほぼ 0 で、差は NaN の埋め方でなく、**K4 で y* が外側の「4 段続く平らな部分」へ移る**ため (例: run_0038 x 80 で y_b 1.75 vs 1.48、+19 %)。
- (d) NaN の区間の埋め方の A/B: case/70 の 10 枚は**確認不能** — K4 では x < 3.9 (直管・縮流部・スロート) の y* が全部 NaN になり、hold の基準位置 x = 3.0 に左側の支持点が無い。
  両側に支持点がある区間だけなら max|B − A| 7.1e-7〜5.5e-6。case/45 は PASS (差 0)。
- 1 枚ごとの切り替わり (run_0010 の 50k→55k で 5.4e-4 r_t) は、K1 でも**窓の 5 枚の平均では窓をまたいで 6.1e-5 に収まる** (上の (a) の K1 の行)。

## 当方の案 (確定前)

1. K = 4 は**不採用**。既定は K = 1 (今と同じ) のまま。`edge_persist` の引数は残す (既定 1 で旧の抽出とバイト同一) か、外す。
2. 帯選択の頑健化は「**整定した場の 5 枚の窓で、各枚は帯を独立に選び、δ_E を点ごとに平均したものを設計の入力にする**」で行う (K1 の (a) が 6.1e-5 で合格)。
   帯を窓で固定する案 (§4.2) は、窓をまたぐ切り替わりを解かない (plan 段のレビュー M3) ので採らない。

## 問い

**上の案 (K = 4 不採用・既定 K = 1・窓の 5 枚の平均を入力にする) に反対の点があれば 1 つ挙げ、採るなら追加で要る確認を 1 つ推奨してほしい。**
候補 (縛られなくてよい): (a) K1 の (a) の合格は 1 組の窓 (run_0010b の W1・W2) だけで、別の run・別の時刻で再現するか未確認、(b) 窓の平均は「整定した場」が前提で、遅いモード (保留中) が解けるまで実データで使えない、
(c) 1 枚ごとの 5.4e-4 の切り替わりが窓の中に 1 回入ると平均を 1e-4 動かしうる (5 枚なら 1/5)、(d) `edge_persist` の引数を残すことの是非 (使われない選択肢)。

## 禁止事項 (厳守)

- ファイルを変更しない (read-only サンドボックスで動いている)。
- **`*.log`, `residual_history*.csv`, `res_*.h5`, `*.vtu`, `*.msh`, `plans/README.md`, `delta_r_equiv_diag.json`, `_ladder_probe/*.json`, `_band_persist_verify/` の CSV を読まない** (巨大)。
- 下に列挙した `sed -n 'A,Bp' <file>` 以外のファイル読みをしない。grep は可。
- 両論併記で逃げず、**推奨は 1 つに絞る**。根拠は `ファイル:行` か提示した数値で示す。

## 読んでよいもの

- `sed -n '1,200p' plans/active/tooling-nozzle-cooled-wall-ns-settling.md`
- `python3 -c "import json;d=json.load(open('case/70.m5_air_wt_demo/_band_persist_verify/verify.json'));print(json.dumps({k:v for k,v in d.items() if k!='raw'},indent=1)[:6000])"`
- `sed -n '330,380p' design/forge_design/metrics/deltastar.py`
- `sed -n '1,40p' notes/reviews/2026-10-13-band-persistence-diagnose.md`
