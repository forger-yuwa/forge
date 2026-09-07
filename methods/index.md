# forge ドキュメント目次

forge の理論的背景と実装解説の索引。運用方針は [README.md](README.md) を参照。

新規ファイル追加・改名・削除時は本目次を必ず更新すること。

## 機能ディレクトリ

小〜中規模の領域は 1 ファイル (`<area>.md`) に理論+実装をまとめ、大規模な領域のみ
`<area>/theory.md` + `<area>/implementation.md` に分割する (方針は [README.md](README.md) 参照)。

| 領域 | ドキュメント | 状態 |
| --- | --- | --- |
| アーキテクチャ全体 | [architecture/overview.md](architecture/overview.md) | 整備済み |
| 離散化レイアウト (cell/node) | [discretization.md](discretization.md) | 整備中 (median-dual 両対応) |
| 勾配再構成 | [gradient.md](gradient.md) | 整備済み |
| リミッタ | [limiter.md](limiter.md) | 整備済み |
| 対流フラックス | [convection/theory.md](convection/theory.md) ・ [implementation.md](convection/implementation.md) | 整備済み (分割) |
| 粘性・熱伝導 | [diffusion.md](diffusion.md) | 整備済み |
| 時間積分 | [time_integration/theory.md](time_integration/theory.md) ・ [implementation.md](time_integration/implementation.md) | 整備済み (分割) |
| 境界条件 | [boundary.md](boundary.md) | 整備済み |
| Poisson 求解 | [poisson.md](poisson.md) | 整備済み |
| 軸対称 | [axisymmetric/theory.md](axisymmetric/theory.md) ・ [implementation.md](axisymmetric/implementation.md) | 整備済み (分割) |
| 乱流モデル | [turbulence/theory.md](turbulence/theory.md) ・ [implementation.md](turbulence/implementation.md) | 整備済み (分割) |
| 多成分熱物性 (TP gas) | [thermophysics.md](thermophysics.md) | M1-M4 + TP境界条件 完了 |
| 有限速度化学 | [chemistry.md](chemistry.md) | 理論・設計確定、Phase 0 (DB ツール・機構ファイル) 完了、ソルバ実装は未着手 |
| 非平衡凝縮 (4 モーメント) | [condensation.md](condensation.md) | 整備中 (Phase 1 受動スカラー骨格) |
| ノズル設計ツール | [design/overview.md](design/overview.md) | 整備中 (Phase 0/2 完了・①風洞 = axis-Mach チェーン A0–A5 完了) |

新規領域を追加した場合は本表に行を追加すること。統合済みの領域は 1 リンク、分割の領域は
theory ・ implementation を併記する。
