# codex 諮問 (diagnose): wall-temperature-profile

- **brief**: [`notes/reviews/briefs/2026-10-13-wall-temperature-profile.md`](../../notes/reviews/briefs/2026-10-13-wall-temperature-profile.md)
- **plan**: [`plans/active/tooling-nozzle-wall-temperature-profile.md`](../../plans/active/tooling-nozzle-wall-temperature-profile.md)
- **date**: 2026-10-13
- **commit**: `138e83b2` (feature/nozzle-deliverables)
- **codex**: effort `high`, 1.2 min, rc=0
- **結論**: plan §4.2・§6 に「壁温表の無次元化には寸法逆算中の実効 `scale` を使う」という契約と、その A/B 検証を追加することを推奨する。

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: plan §4.2・§6 に「壁温表の無次元化には寸法逆算中の実効 `scale` を使う」という契約と、その A/B 検証を追加することを推奨する。

第 1 仮説: **[Major・要再検証] 寸法逆算中に、積分法の長さ換算には更新後の半径、`Tw_table` の換算には元の `spec.r_throat` が使われ、積分法だけ壁温分布がずれる。** 確度: 中（経路の存在は確認済み、分布対応は未実装）
  
  根拠:
  
  - `design/forge_design/feedback/deltastar_loop.py:134` は、反復中の `rt` を `integral_delta_r(..., scale=rt)` に渡す。`:185` で反復ごとに呼ばれ、`:197` で半径が更新される。
  - `design/forge_design/evaluate/runner_axismach.py:1161` は明示された `scale` を採用する一方、`:1168` の熱条件取得は現在、半径を受け取らない property。`:1172` では積分法へ実効 `scale` を渡す。この接続点が改修の要所である。
  - `design/forge_design/feedback/deltastar_integral.py:200` は無次元座標 `x` で温度を評価するため、表の換算半径が違えば、そのまま異なる壁温を閉包へ渡す。
  - plan `tooling-nozzle-wall-temperature-profile.md:104` の V1 は表の換算を検査するが、**`scale != spec.r_throat` の呼び出し経路**を指定していない。`:106` の V2 は積分法について記録確認に留まり、NS 側の温度一致ではこの誤配管を排除できない。

  対案: `integral_delta_r` 内で解決した同じ `scale` を `p.wall_thermal_bc_integral(scale)` と `integral_bl(..., scale, ...)` の両方へ渡すことを契約化する。分布表を元の半径で換算してキャッシュしない。plan §4.2 の方向性は採用するが、実効半径の指定と検証が不足している。

  反証条件: 元の問題定義を固定し、実効 `scale` だけ変えても、積分法が実際に使用した全評価点の `Tw` が、原表を実寸座標 `scale × x` で評価した温度と一致すること。

第 2・第 3 仮説: 追加しない。確認した `run_staged_ns` の段切替処理は同じ run 内で動き、削除対象も `res_*`（`runner_axismach.py:1555`–`:1557`）。CSV 脱落を第一候補にする証拠はない。

判別 A/B: **変更点は `integral_delta_r` の `scale` 引数だけ**とする。`p.spec.r_throat = S`、非一様の実寸温度表、設計壁、積分設定を固定し、A は `scale=S`、B は `scale=1.2S`。両方の壁範囲を覆う表を使い、換算を間違えた場合に内部の評価点で 10 K 以上の差が出る勾配を選ぶ。実行量は積分法を各 1 回、forge は不要。戻り値の `res_init["Tw"]`（実使用温度、`deltastar_integral.py:223`）について、独立に原表から計算した `T_ref(scale × res_init["x"])` との最大差を測る。

  → **A は最大差 ≤ 1e−8 K、B だけ外れ、旧半径による分布に一致するなら仮説を支持。A/B とも ≤ 1e−8 K なら、この換算漏れ仮説を棄却。** それ以外の失敗を、この仮説の裏付けにはしない。δ_r の一致は判定対象にしない。

やらない方がよいこと: V2 の最終壁ノード温度が最大差 ≤ 0.1 K だったことだけで、積分法・NS・成果物の配管全体を合格にすること。NS へ正しい表が渡り、積分法だけ誤る経路を見逃す。

呼び出し側の前提への異議: 「単一ソース」は同じ実効温度分布を保証しない。表が同じでも換算半径が違えば別条件になる。case/48 の最大差 0.068 K は提示された実測値であり、今回直接検証していない。また、短い run の温度拘束検査と、流れの収束・派生量の定常性は別である。

不足情報: 分布対応後の実装と、上記 A/B の実使用温度。現時点で実装済みの不具合とは断定しない。ファイル変更は禁止されているため **plan 未反映**。呼び出し側で §4.2・§6 に反映すること。
