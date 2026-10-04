# codex 諮問 (diagnose): gradient-scalar-lsq-4a-closeout

- **brief**: [`notes/reviews/briefs/2026-09-26-gradient-scalar-lsq-4a-closeout.md`](../../notes/reviews/briefs/2026-09-26-gradient-scalar-lsq-4a-closeout.md)
- **plan**: [`plans/active/gradient-scalar-lsq-unification.md`](../../plans/active/gradient-scalar-lsq-unification.md)
- **date**: 2026-09-26
- **commit**: `eb8fb6a4` (feature/sern-design)
- **codex**: effort `high`, 3.0 min, rc=0
- **結論**: **既存のダンプ ON/OFF・1 step 出力を、決定的配列の厳密比較と NS 勾配の順列集合検査で再判定することを、次の一手とする。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0）

| 対象 | 採否・重大度 | 理由と対案 |
|---|---|---|
| **1. #4a の commit** | **採用：commit 可。Minor あり** | 差分は添付 patch と SHA256 が一致した。`preGatherDump` は env 無しで GPU 操作前に戻り、有効時も D2H 読み取りである（[periodicNode_d.cu:304](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/periodicNode_d.cu:304)）。`S1_hook4a.txt:32,62,92,136` の S1(1) は全4ケース PASS。数値更新を壊す欠陥は確認できず、診断 hook として commit してよい。ただし「ダンプ非干渉検証済み」とは記録しない。Minor：null 成分をゼロに置換する処理（同ファイル:312）は欠測と真のゼロを混同するため、診断失敗として明示する方がよい。また毎回の vector 構築（`speciesTransport_d.cu:732–738` 等）は呼び出し側の env 判定で避けられる。「無視できる性能差」は未測定。 |
| **2. H2 の判定規則** | **原案は却下・修正版を採用〔Major〕** | `#4b` の6本は**すべてダンプ有効**であり、ダンプ無しとの pre-gather 一致は測っていない（[pregather_check.py:49](/home/sano/work/forge-sern-design/case/09.Taylor-Green/_g0_lsq_seam/pregather_check.py:49)、同:79–90）。「両側ともビット一致なら厳密比較」という条件だけでは、hook が新たに生んだ揺らぎまで許容できてしまう。対案は、**決定的な配列はコード上の依存関係で指定して厳密比較し、周期 NS 勾配は部分和の順列集合で判定、その他の既知の非決定的経路の出力だけ両側ノイズ規則を使う**こと。事後規則であることと初回 FAIL を保存する。 |
| **3. S0/S1 を閉じて S2 へ** | **現時点では要再検証〔Major〕** | `#4b=A` を条件とする S1(2) 訂正の採用は支持する。しかし現在の集合検査はダンプ有効 run だけで、元のダンプ無し S1 の出力を検査していない（[pregather_check.py:99](/home/sano/work/forge-sern-design/case/09.Taylor-Green/_g0_lsq_seam/pregather_check.py:99)）。まず既存のダンプ無し出力も同じ部分和集合に照合し、新規則の VERDICT を残す。加えて `S0_axi.txt`・`S0_axi_m1.txt` には a/b の PASS しかなく、plan §6 が要求する S0-c の結果または適用除外理由が必要。S2 投入前には #2e の起点固定と #2f の判定スクリプトを完了させる。 |

判定規則には、さらに次の訂正が必要です。plan §6 S1(2) の「**勾配・リミタを gather 前で比較**」は実装順序と合いません。`limiter_d_wrapper` は gather **後**です（[main.cpp:1479](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:1479)、同:1485）。pre-gather の対象を **NS 18 勾配配列**に限定し、リミタは評価後の配列で比較してください。今回のリミタは `S1_hook4a.txt:139–144` の全15対一致側に含まれており、リミタ差が観測されたわけではありません。これは **Minor の判定仕様訂正**です。

結論: **既存のダンプ ON/OFF・1 step 出力を、決定的配列の厳密比較と NS 勾配の順列集合検査で再判定することを、次の一手とする。**

第 1 仮説: tgv の NS 勾配差は、周期 gather の `atomicAdd` の加算順序差で説明できる。 確度: **高（ダンプ有効の測定範囲）**  
　根拠: [periodicNode_d.cu:65](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/periodicNode_d.cu:65) は member を root に `atomicAdd` する。`PREGATHER_4b.txt:15–52` では pre-gather NS 18 配列が6本でビット同一、post-gather の相違は `res_0` 124節点・`res_1` 232節点、すべて4 member、順列集合外は0。記録の判定は **`VERDICT #4b: A`**。  
　反証条件: 同じ入力・バイナリのダンプ無し出力に、共通部分和から作る順列集合の外の NS 勾配が1成分でもある、または非周期節点・2 member 以下に厳密不一致がある。

第 2 仮説: case48 の `roe` 差も、既存の非決定的演算の影響が伝わったもの。 確度: **中**。`DUMPCHECK_case48.txt:27` の差は最大 `0.0156`、無し–有りの不一致数は2、有り同士は3。旧バイナリ同士にも同じ最大差がある（`S1_hook4a.txt:29`）。ただし、この値だけで発生経路までは確定できない。

判別 A/B: **変更点は `FORGE_DUMP_PREGATHER` の OFF/ON だけ**。tgv、同一 `gg` バイナリ・入力、1 step × 各3本の既存データを使い、追加計算なしで再解析する。OFF 側は `/home/ubuntu/sglsq/s1h/s1_tgv_gg_{a,b,c}`、ON 側は `/home/ubuntu/sglsq/pg/pg_tgv_gg_{a,b,c}`。`res_0/1` の NS 勾配を ON 側の共通部分和に照合し、非周期節点は部分和そのもの、2 member は唯一の和、3 member 以上は順列集合への所属を要求する。  
→ **A：全点が適合**なら、観測された NS 勾配差について順序差の説明を採用する。  
→ **B：集合外または決定的部分の不一致あり**なら、順序差だけという仮説を棄却し、S2 移行を保留する。

やらない方がよいこと: FAIL の3配列だけを例外登録すること、両側ノイズの2倍以内という理由だけで全配列を「数値不変」とすること、今回の判定を通すために gather 自体を決定的実装へ変更すること。

呼び出し側の前提への異議: 「gather 以外の非決定源を除外済み」は範囲が広すぎる。確認済みなのは、**ダンプ有効 run の NS 勾配が順列和で説明できること**であり、ダンプ無しの局所入力、保存量、GG スカラー勾配まで除外した証拠ではない。また、ダンプ比較表は両側内部対の**最大差**を出力していない（`pregather_check.py:156–157`）。したがって「新規則なら3件とも PASS」は、不一致数だけでなく最大差も元データから再計算して記録する必要がある。

不足情報: AWS の元配列はこの環境に無く、上記再判定は未実施。軸対称2変種の S0-c 記録も不足している。今回は1 step の作用素・非干渉レビューであり、収束・準定常性の合格を意味しない。**plan 未反映**。呼び出し側で `plans/active/gradient-scalar-lsq-unification.md` §5.1・§6・§6.1 に反映すること。
