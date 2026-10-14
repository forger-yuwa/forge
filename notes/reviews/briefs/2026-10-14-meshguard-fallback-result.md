# 諮問ブリーフ: result 段の解釈 (plan tooling-nozzle-mesh-size-guard と tooling-review-fallback-backend)

AGENTS.md 条件 7 (result 段の解釈を確定する前)。判断役のスイッチが codex。実装は commit `7e7c529e` (基準 `88675245`)。

## 観測事実 (各 plan の §6「結果 (2026-10-14)」と §5.1 に詳細)
### mesh-size-guard
- V3: `node_count` は 5 経路で実際の座標数と一致、HEAD 版 mesh2d と 9 通りでビット一致。A/B (上限 25 万に差し替え) A 194,000 通過・B 402,000 で生成・設計チェーン・run dir の前に停止。
- V2: c2b カードは変更前の生成器と、節点数のコメント行・時間の値を除いて一致。NS ni 6000 のカードは `_mesh_ar_full` の前に停止。明示の `mesh_euler.ni` は保持。
- V4: 文書。case/70 の Euler を `--prepare-only` で n_nodes 174600・source default。
- 振る舞いの変化: 既存の YAML 7 本 (case/44 ×4・case/45 ×3、43.5〜51.4 万点) が prepare_ns で止まるようになる。
- 残作業: 他の runner 4 本に判定が無い、生成器と prepare_ns の NS nj の既定 (81 / 97) の食い違い。
### review-fallback-backend
- V1 (同調の対照試験、opus、旧版の抜粋): 6/6 正答、B 3/3 で誤った前提を明示的に退ける。A3 は 2 回 600 s の timeout (出力なし)、3 回目 8 s で正答。
- V2: 偽の CLI の試験 71/71、既存の記録 321 件が ok に分類。V3: haiku で写し 9006 ファイルのハッシュ一致、init の道具 [Glob, Grep, Read]、渡していない道具の拒否は permission_denials に残らない。
- V4: haiku の diagnose・plan の記録が ok、表示された行で check_plans PASS。V5: 文書。
- 残作業: 判定の語の正規表現が本文で最初の語を取る (既存の欠陥)、子が先行の記録を読める。

## 問い
**2 つの plan について、「事前に登録した §6 に照らして合格」と解釈して accepted にしてよいか、それとも解釈を変える・追加の確認が要る点があるか。最も重要な 1 点を挙げ、推奨を 1 つ。**
候補 (縛られなくてよい): (a) V1 の A3 の timeout 2 回を「判定役が使えない」とみなし再試行で埋めたのは、事前の 6/6 の条件の読み替えか、(b) 既存の YAML 7 本が止まるのは他のセッションの作業を壊すか、
(c) V3 の「または permission_denials に残る」が当たらなかったのに合格としてよいか、(d) mesh-size-guard の V1b を判定不能のまま閉じたこと。

## 禁止事項 (厳守)
- ファイルを変更しない (read-only サンドボックスで動いている)。
- **`*.log`, `residual_history*.csv`, `res_*.h5`, `*.vtu`, `*.msh`, `plans/README.md` を読まない** (巨大)。
- 下に列挙した `sed -n 'A,Bp' <file>` 以外のファイル読みをしない。grep は可。
- 両論併記で逃げず、**推奨は 1 つに絞る**。根拠は `ファイル:行` か提示した数値で示す。

## 読んでよいもの
- `sed -n '1,170p' plans/active/tooling-nozzle-mesh-size-guard.md`
- `sed -n '1,170p' plans/active/tooling-review-fallback-backend.md`
- `git diff --stat 88675245 7e7c529e`
