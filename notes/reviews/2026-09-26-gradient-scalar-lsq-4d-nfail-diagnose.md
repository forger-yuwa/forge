# codex 諮問 (diagnose): gradient-scalar-lsq-4d-nfail

- **brief**: [`notes/reviews/briefs/2026-09-26-gradient-scalar-lsq-4d-nfail.md`](../../notes/reviews/briefs/2026-09-26-gradient-scalar-lsq-4d-nfail.md)
- **plan**: [`plans/active/gradient-scalar-lsq-unification.md`](../../plans/active/gradient-scalar-lsq-unification.md)
- **date**: 2026-09-26
- **commit**: `eb8fb6a4` (feature/sern-design)
- **codex**: effort `high`, 1.8 min, rc=0
- **結論**: **(b) を一度限りの追加試験として採用し、`lsq` のダンプ無し・有りを各 +5 本、計8本で、全 N 配列の2倍規則と E/P を再判定する。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0 / Major 2）

| 対象 | 採否 | 根拠と対案 |
|---|---|---|
| `#4d = A` だからダンプ非干渉を記録する | **却下・Major** | [dump_recheck.py:176](/home/sano/work/forge-sern-design/case/09.Taylor-Green/_g0_lsq_seam/dump_recheck.py:176) は `N` の成否を総合判定に使っていない。plan の [#4d:99](/home/sano/work/forge-sern-design/plans/active/gradient-scalar-lsq-unification.md:99) と [§6:120](/home/sano/work/forge-sern-design/plans/active/gradient-scalar-lsq-unification.md:120) の結論範囲が不整合。**「E/P は A、N は FAIL、ダンプ全体の非干渉は未確認」**と分けて記録する。 |
| (a) `gg/lsq` のノイズ対をプールする | **却下** | [実測:167](/home/sano/work/forge-sern-design/case/09.Taylor-Green/_g0_lsq_seam/DUMP_RECHECK_tgv.txt:167) の `gg` 側の最大不一致数 45 を採れば通るが、同一ノイズ分布の証拠にはならない。演算式が同じでも、前段カーネルが違えば実行順序の分布まで同じとは限らない。`lsq` 内で比較する。 |
| (b) 各 +5 本、計8本で再判定する | **修正採用** | 一度限りの標本数感度確認として妥当。ただし **FAIL 継続から「ダンプが分布を変えた」と確定する解釈は却下・Major**。[判定式:154](/home/sano/work/forge-sern-design/case/09.Taylor-Green/_g0_lsq_seam/dump_recheck.py:154) は最大値の比較であり、分布差の検定ではない。PASS でも「固定した許容規則を満たした」まで。 |
| `N FAIL` が残ってもダンプ無しで S2 へ進む | **条件付き採用** | [着手条件:101](/home/sano/work/forge-sern-design/plans/active/gradient-scalar-lsq-unification.md:101) の `E/P=A` と `#2e・#2f` を維持する。追加試験でも E/P が通り、未解決事項が N の比較に留まるなら、ダンプ無しの S2 を保留する根拠にはしない。`#4a` は「診断専用・非干渉未確認」のまま残す。 |

結論: **(b) を一度限りの追加試験として採用し、`lsq` のダンプ無し・有りを各 +5 本、計8本で、全 N 配列の2倍規則と E/P を再判定する。**

第 1 仮説: atomic 加算の順序変動に対し、各3本から作った群内最大不一致数が小さく出たための閾値超過。確度: **中**
  
  根拠: [DUMP_RECHECK_tgv.txt:168](/home/sano/work/forge-sern-design/case/09.Taylor-Green/_g0_lsq_seam/DUMP_RECHECK_tgv.txt:168) の `res_1/res_ro` は、群内最大が11・14、交差最大が40で、`40 > 2×14`。最大絶対差は三組とも `2.98e-8`。さらに [periodicNode_d.cu:32](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/periodicNode_d.cu:32) に **`res_ro` 自体の周期集約の `atomicAdd`** がある。面ループだけをノイズ源とする説明は不十分。
  
  反証条件: 実行条件を揃えて順序を混ぜた追加試験でも、ダンプ有無に対応する分離が再現されること。その場合、単なる小標本の偶然を主説明として採らない。ただし2倍規則の再 FAIL だけでは、統計的な反証確定にはならない。

第 2 仮説: ダンプが実行タイミングを変え、許容される加算順序の出現頻度を変えている。**未確認**。E/P の集合検査は、値の到達可能性を検査しており、出現頻度の同一性は保証しない。

判別 A/B: 変更点は **`FORGE_DUMP_PREGATHER` の無効／有効だけ**。同一バイナリ・入力・GPU・BLOCKSIZEで各1 step。追加5組の実行順を事前固定して有無を混ぜ、計8本ずつで判定する。

- 群内28対＋28対、交差64対について、配列ごとに **最大絶対差・最大不一致数が両方とも群内最大の2倍以下**。`res_ro` だけでなく `res_0/res_1` の全 N 配列を対象とし、E/P も再確認する。
- **A：E/P 通過・全 N 通過** → 「固定8本の規則では超過しない」。今回の FAIL を持続的干渉の証拠とする解釈を退ける。初回 FAIL は保存する。
- **B：E/P 通過・N FAIL 継続** → 小標本の偶然として閉じる案を退け、ダンプ非干渉は未確認とする。追加本数をさらに増やして合格を探さない。S2 は上表の条件で進める。
- E/P が失敗した場合は既存の停止条件を適用し、S2 を保留する。

やらない方がよいこと: `gg` の45を借りて合格にすること、PASS するまで反復を追加すること、N FAIL を既存の「P 集合外」の **B′** と同一視すること。

呼び出し側の前提への異議: 「3本では検出力不足」は未立証。観測された11・14・45はそれぞれ**対比較の最大値**であり、ノイズ分布の範囲ではない。また、28対は8本を共有するため、28個の独立標本ではない。「質量残差へスカラー勾配が入らない」も、今回の設定と呼び出し順を確認できておらず、確定根拠には採用しない。

不足情報: 元の AWS run はこの環境に存在せず、HDF5 の再集計はできなかった。確認できたのは結果テキストと判定コード。全 run の入力・実効設定・バイナリ同一性も未検証で、スクリプトが表示する provenance は [基準 run のみ:82](/home/sano/work/forge-sern-design/case/09.Taylor-Green/_g0_lsq_seam/dump_recheck.py:82)。追加試験では各 run の同一性を記録すること。今回は1 step の非干渉診断であり、収束・定常性の判断はしていない。

**plan 未反映**（read-only 指示による）。呼び出し側で §5.1 #4a・#4d・#5、§6 に判定範囲と固定8本の停止規則を反映すること。
