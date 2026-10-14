# 引き継ぎ: ノズルの外部領域と過膨張剥離の評価 (別セッション向け)

- 日付: 2026-10-11
- ブランチ: `docs/thruster-nozzle-design` (作業ツリー `/home/sano/work/forge-thruster`。別セッションは自分の作業ツリーを切ること — 1 作業ツリー 1 セッション)
- **正本は plan**: [`plans/active/tooling-nozzle-ambient-separation.md`](../../plans/active/tooling-nozzle-ambient-separation.md)。残作業・判定基準・採否はそこを見る。本メモは写しとポインタだけ。

## 何をするか

③ベル・④デュアルベルの設計評価に、ノズル出口より外の外部領域と外気圧を入れ、過膨張剥離を計算できるようにする。
forge の SST の剥離予測を DLR TIC test case 1A で検証し、同じ格子の SU2 SST と比べて forge 固有の誤差を切り分ける。
設計側 (③④⑤、親 plan [`tooling-nozzle-design-tool.md`](../../plans/active/tooling-nozzle-design-tool.md)) はこの機能が使える前提で進めている (親 §3-13)。

## 最初にやること

1. plan §8-1 の外気 BC の方式 (a 軸対称の双方向開放境界を作る / b 現行 BC + 境界の流入出診断と領域感度で受け入れ) をユーザに確認する。現行 BC の契約 (逆流は内部の等エントロピー、`inlet_Pressure` は流出で速度 0) は plan §3。
2. 決まったら codex `--stage plan` (`python3 solver_density_cuda/tools/codex_review.py plans/active/tooling-nozzle-ambient-separation.md --stage plan`)。
3. plan §5.1 #2 (TC1A の表の転記)。原典 PDF は `/home/sano/work/forge/papers/nozzle_separation/` (git 管理外)。

## 参照

- 検証データの調査: [`notes/investigations/nozzle-flow-separation-validation-survey.md`](../investigations/nozzle-flow-separation-validation-survey.md)
- codex diagnose (2026-10-11、NO-GO → 改訂済み): [`notes/reviews/2026-10-11-thruster-nozzle-child-plans-diagnose.md`](../reviews/2026-10-11-thruster-nozzle-child-plans-diagnose.md)
- 推力の会計の定義は兄弟 plan [`tooling-nozzle-bell-two-mode.md`](../../plans/active/tooling-nozzle-bell-two-mode.md) §4.1 と共通 (変えるときは両方)
