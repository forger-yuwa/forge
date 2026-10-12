"""入力層: 実寸の要求カード → 問題 YAML (plans/accepted/tooling-nozzle-wt-deliverables.md §4.1)。

- `wt_card`: 風洞ノズルの要求カード `forge.wt_card/v1` の読み込み・検査と、問題 YAML (`wind_tunnel_axisym_axismach`) の生成
  (長さの r_t 単位への写し、物理スロート径からの r_t の固定点、壁解像の第 1 層と ni の初期見積もり、空気の液化の事前警告)。
  `python -m forge_design.intake.wt_card CARD --out problem.yaml`。
- `recipes/*.yaml`: 生成器が写す構成キーの雛形 (探索段階 `axismach_explore_v1`・最終の NS 用 `axismach_final_v1`)。
- `examples/*.yaml`: カードの例。

`wt_card` はここで import しない (`python -m forge_design.intake.wt_card` の実行で同じモジュールを 2 回読む警告を避ける)。
"""
