# codex 諮問 (diagnose): mesh-size-guard

- **brief**: [`notes/reviews/briefs/2026-10-14-mesh-size-guard.md`](../../notes/reviews/briefs/2026-10-14-mesh-size-guard.md)
- **plan**: [`plans/active/tooling-nozzle-mesh-size-guard.md`](../../plans/active/tooling-nozzle-mesh-size-guard.md)
- **date**: 2026-10-14
- **commit**: `88675245` (feature/nozzle-cooled-wall-settling)
- **codex**: effort `high`, 1.8 min, rc=0
- **結論**: §4.2 の事前判定を、指定された `ni × nj` ではなく、生成器と共通の配点処理で求める実効節点数に変更することを推奨する。

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: §4.2 の事前判定を、指定された `ni × nj` ではなく、生成器と共通の配点処理で求める実効節点数に変更することを推奨する。

第 1 仮説: **[Major] `axis_gap_frac` 指定時に節点数を過小評価し、上限超過のメッシュを生成してから停止する。** 確度: 高  
  根拠: `plans/active/tooling-nozzle-mesh-size-guard.md:57` は生成前の `ni × nj` 判定と生成後の照合を規定している。しかし `design/forge_design/meshing/mesh2d.py:170–173` は `nj = len(s)` に置き換え、同ファイル:189 で実効 `ni × nj` の座標配列を確保する。`axis_gap_frac` は `design/forge_design/evaluate/runner_axismach.py:718` で許可されている。  
  例えば `ni=2000, nj=97, axis_gap_frac=0.005` なら指定値では **194,000 点**だが、半径方向の間隔が最大 0.005 なので実際には **nj ≥ 201、402,000 点以上**となる。これはコードからの下限計算で、run の実測値ではない。生成後の照合では「大きなメッシュを作らせない」という目的に間に合わない。  
  対案: 径方向の実効点数を求める処理を生成器と事前判定で共有し、座標・接続配列の確保と run directory 作成より前に上限を適用する。  
  反証条件: `axis_gap_frac` を含む実効点数が事前判定に入り、上記入力が大規模配列の確保・run directory 作成より前に拒否されること。

第 2・第 3 仮説: 挙げない。

判別 A/B: Euler の `ni=2000, nj=97, wall_first_frac=0.005`、上限 250,000 点を固定し、**`axis_gap_frac` だけを A: 未指定／B: 0.005** に変える。実装後の事前判定試験を各 1 回、forge は **0 step**。見る量は判定点数・停止位置・メッシュ生成処理の呼出し有無。  
→ A が 194,000 点で許可され、B が 402,000 点以上として生成前に拒否されれば仮説は棄却。B が事前判定を通り、生成後の照合で初めて止まれば仮説を支持。

やらない方がよいこと: YAML の `ni × nj` だけを検査して、生成後の点数照合を事前ガードの代わりにすること。

呼び出し側の前提への異議: 「構造格子だから指定 `ni × nj` が実際の節点数になる」は成立しない。既存コード自身が `runner_axismach.py:862–863` で、`axis_gap_frac` 使用時の `nj` を導出値として記録している。

不足情報: ガードは未実装なので、これは設計上の失敗予測である。禁止された run データは読んでおらず、提示された V1a の数値や収束・準定常性は独立検証していない。
