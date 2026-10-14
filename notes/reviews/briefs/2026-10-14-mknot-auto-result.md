# 諮問ブリーフ: M_knot の自動選択の result 段の解釈 (plan tooling-nozzle-mknot-auto)

AGENTS.md 条件 7 (result 段の解釈を確定する前、予想どおりに見えるときも)。判断役のスイッチが codex。コードの最終 b2528f7f、検証の記録 5c885f30。

## 観測事実 (plan §5.1 #2・#2b・#2c・#3・#3b・#3c、§9 の 2026-10-14 の項目)
- 事前登録の §6: V1a (固定の入力の再現)・V1b (M6 の目的の回帰)・V2 (a)〜(k) (合成)・V3 (数値の M_knot の後方互換)・V4 (統合、雛形 final と explore)・V6 (N/2N の A/B)。
- 結果: V1b 2.6/2.7 (規則を切った設定でも既定でも)。V1a 2.069270614・距離 0.0975・照合 37/37・ゲート PASS。V6 2.07/2.07 (素の最小は 2.07 → 2.05 に動き同点の規則で 2.07)。
  V2 89 ok (手元)。V3 PASS (AWS、m5 2.25・c2b 2.069270614 の数値のカードで YAML・card_log が一致、探索を呼ばない)。
  V4 final・explore (AWS) とも 2.069270614、ゲート PASS (窓を広げず)、寸法 PASS (拡大部長 +0.37 mm / 公差 1 mm、物理スロート径 −2〜−4e-8 m / 2e-5)、YAML の要約 2.6 kB・card_log と照合 ok、
  所要 74 分・67 分 (§6 V4 の目安「手元で 15 分以内」を超える; §5.1 #6 に記録)。
- 諮問・レビューの採否はすべて plan §6.1 に記録 (plan 段 GO-with-changes 全件採用、諮問 5 回)。

## 問い
**事前登録の §6 に照らして「合格」と解釈して accepted にしてよいか。最も重要な 1 点を挙げ、推奨を 1 つ。**
候補 (縛られなくてよい): (a) V4 の所要 (目安の 4〜5 倍) を合否に入れるか、(b) V4 で窓を広げる経路が実データで一度も通っていない (合成だけ)、(c) V6 の合格が同点の規則に支えられている、
(d) M5 の 1 形状 (c2b) だけで既定を自動にしてよいか (V1b は M6 の記録の再生)。

## 禁止事項 (厳守)
- ファイルを変更しない (read-only サンドボックスで動いている)。
- **`*.log`, `residual_history*.csv`, `res_*.h5`, `*.vtu`, `*.msh`, `plans/README.md`, `case/70.m5_air_wt_demo/_mknot_*/*.json` を読まない** (巨大)。
- 下の範囲以外のファイル読みをしない。grep は可。
- 推奨は 1 つに絞る。根拠は `ファイル:行` か提示した数値で示す。

## 読んでよいもの
- `sed -n '1,160p' plans/active/tooling-nozzle-mknot-auto.md`
- `tail -40 plans/active/tooling-nozzle-mknot-auto.md`
- `cat case/70.m5_air_wt_demo/_mknot_card_v4_3b/check_final.txt case/70.m5_air_wt_demo/_mknot_card_v4_3b/check_explore.txt`
