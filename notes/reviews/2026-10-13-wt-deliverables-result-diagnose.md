# codex 諮問 (diagnose): wt-deliverables-result

- **brief**: [`notes/reviews/briefs/2026-10-13-wt-deliverables-result.md`](../../notes/reviews/briefs/2026-10-13-wt-deliverables-result.md)
- **plan**: [`plans/active/tooling-nozzle-wt-deliverables.md`](../../plans/active/tooling-nozzle-wt-deliverables.md)
- **date**: 2026-10-13
- **commit**: `0785d720` (feature/nozzle-deliverables)
- **codex**: effort `high`, 1.6 min, rc=0
- **結論**: `accepted` は保留し、**総熱量を一定に保ったまま局所熱流束だけが変動する負例で、`thermal_use` が誤って PASS しないか**を一度検証する。

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: `accepted` は保留し、**総熱量を一定に保ったまま局所熱流束だけが変動する負例で、`thermal_use` が誤って PASS しないか**を一度検証する。

第 1 仮説: **［Major］`thermal_use` は総熱量の定常性しか検査せず、用途である局所の壁熱伝導に必要な熱流束分布の定常性を保証していない。** 確度: 高（検査対象の欠落はコードで確認、実 run での発生は未確認）。
  
  根拠: `design/forge_design/export/wall_tables.py:649–650` は時系列を作るが、準定常判定に渡すのは `--series-cols Q_total_W` だけ。壁の一部で熱量が増え、別の部分で同量減れば、この検査には現れない。一方、用途の式は局所の `q` を使う `T_back = T_w − q t/k` である（`design/tests/run_wall_tables_tests.py:420–421`）。現在の `run_0006` が FAIL になることは、この偽陽性を排除する証拠ではない。
  
  対案: 総熱量に加え、壁ごとの局所 `q` 分布の変化を検査し、その合格を `thermal_use` の条件にする。
  
  反証条件: 総熱量・収束・壁解像・有限性をすべて合格に保っても、局所熱流束の変動を理由として現行 `thermal_use` が FAIL を返すこと。

第 2・第 3 仮説: 追加しない。

判別 A/B: **変更点は熱量の空間配分の時系列だけ。forge は 0 step、既存の合成 PASS 試験を使う。** A は各節点の熱量を一定にする。B は離れた二領域へ、各時刻で総和が厳密にゼロになる増減を加え、局所 `q` が基準値の ±20 % で交互に変動する系列にする。残差・幾何・壁応力・総熱量は同一とし、判定窓を満たす十分な snapshot 数を用意する。  
→ **A が PASS、B も PASS なら仮説を支持**し、局所定常性の条件を追加する。**A が PASS、B が局所変動を理由に FAIL なら仮説を棄却**する。入力不備による FAIL は反証にならない。

やらない方がよいこと: V5 の単位検査のためだけに M5 run を延長すること。今回の穴は長時間計算では閉じない。また、`ΣQ = ΣqA` は `q = Q/A`（`wall_tables.py:358`）から成立する自己整合であり、独立な熱収支検証として数えない。

呼び出し側の前提への異議: **「未収束だから V5 の単位検査が無効」は却下。** plan §6 V5 は桁・符号の検査で、熱設計への使用可否は別条件である。合成試験には既知の `q` と `Q` の照合もある（`run_wall_tables_tests.py:438–450`）。また、`−qwall` と `iface_q_compact` の近さだけで非独立とは断定できない。コード上は再構成勾配（`viscousFlux_d.cu:615–617,654`）と第一内部点の温度差（`conjugateWall.cpp:264`）で式が異なる。単位検査の支持と、熱設計用ゲートを含む書き出し器全体の受入れを分けるべきである。

不足情報: 上記負例の結果。生の run データと VERDICT 原本は指定に従い読んでおらず、実測値は提示記録として扱った。ファイル変更なし、**plan 未反映**。
