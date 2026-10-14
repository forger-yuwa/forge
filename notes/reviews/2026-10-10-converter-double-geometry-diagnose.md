# codex 諮問 (diagnose): converter-double-geometry

- **brief**: [`notes/reviews/briefs/2026-10-10-converter-double-geometry.md`](../../notes/reviews/briefs/2026-10-10-converter-double-geometry.md)
- **plan**: [`plans/active/architecture-float-state-double-geometry.md`](../../plans/active/architecture-float-state-double-geometry.md)
- **date**: 2026-10-10
- **commit**: `f704c1de` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 2.7 min, rc=0
- **結論**: **案 (b) を採用し、変換器専用の double 幾何を入力から HDF5 出力まで保持する設計に §4.6・§6.14 を補強してから、float/FP64 ビルド間の変換比較を行う。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **案 (b) を採用し、変換器専用の double 幾何を入力から HDF5 出力まで保持する設計に §4.6・§6.14 を補強してから、float/FP64 ビルド間の変換比較を行う。**

確認対象は HEAD `f704c1de`。以下はコードに基づく設計レビューであり、既存 run の結果は再検証していない。ファイル変更・変換器・`forge` の実行は行っていない。**plan 未反映**。

採否表:

| 重大度 | 対象・採否 | 根拠と対案 |
|---|---|---|
| **Major** | 案 (b)：**採用。ただし局所変数の型置換だけでは不十分** | `gmshReader` は `mesh` を継承し、読み取った座標を `node` に格納する。双対構築後も共用の `cell`・`plane` に戻し、出力時はそれらを読む。根拠: `solver_density_cuda/mesh/gmshReader.hpp:34,516,2224,2240,2343,2423,2465`。**座標・primal 幾何・dual 幾何・境界半割面の正本を変換器専用の double 配列／構造体に置く。** 初期値設定に必要な共用 `mesh` への変換は一方向とし、その丸め済みデータから幾何を再計算・出力しない。 |
| **Major** | 壁距離も型置換で済む：**却下** | 壁点・評価点を `msh.nodes/cells/planes` から取り、最後に `var.c["wall_dist"]` へコピーしている。書き出しも `var.c` 経由。根拠: `solver_density_cuda/input/calcWallDistance_kdtree.cpp:137,159,175`、`solver_density_cuda/mesh/gmshReader.hpp:2525`。**double の点群を受け、double 配列を返す計算部分を切り出し、`/VALUE/wall_dist` はその配列から直接書く。** ソルバ用のコピーと HDF5 用の正本を分ける。 |
| **Major** | 案 (a)：**今回の推奨として却下** | 変換器は `variables.cpp` 等と `cuda_forge` をリンクし、実際に `allocVariables`・H2D コピーを呼ぶ。根拠: `solver_density_cuda/CMakeLists.txt:114`、`solver_density_cuda/mesh/convertGmshToForge.cpp:57`、`solver_density_cuda/input/calcWallDistance_kdtree.cpp:187`。変換器側だけ typedef を変えると、共有型・関数の整合を保証できない。**共用型と `cuda_forge` のビルドは維持し、幾何だけ独立した型で扱う。** CUDA 依存の完全除去は別変更とする。 |
| **Major** | 座標 dtype だけで警告：**要再検証・条件追加** | ソルバは座標とは別に面ベクトル・面積・面重心・セル重心を読む。座標だけ double の混在ファイルを見逃す。根拠: `solver_density_cuda/mesh/mesh.cpp:238,261,265,313`。**使用する幾何データセット全体を検査し、不足する精度を一つの警告にまとめる。** §6.14 には混在 dtype の検査を追加する。double への後付け変換は元の精度を回復しないため、dtype 確認を生成精度の証明とは扱わない。 |
| **Major** | §6.14 のビット一致：**条件付き採用。許容差は緩めない** | 幾何は CPU 上で生成され、壁距離は `HAVE_KDTREE` と環境変数で実行経路が変わる。根拠: `solver_density_cuda/mesh/gmshReader.hpp:1864`、`solver_density_cuda/input/calcWallDistance_kdtree.cpp:165`。**同じソース・コンパイラ・浮動小数点フラグ・依存ライブラリ・探索経路・入力で固定する。** 「`forge-wallfit-bin-fp64` か段④」という参照の選び方では不足。参照バイナリと入力をハッシュで一意に固定する。 |
| **Minor** | 「float32 格子では段①〜④の利点が無い」：**却下、表現修正** | float32 座標を double で読み直しても失われた座標情報は戻らない。しかし、それだけで double による後続演算・累積の改善まで否定できない。根拠: `solver_density_cuda/mesh/mesh.cpp:237`。警告は **「入力時に失われた幾何精度は回復できない」** とする。また §2 の「変換器は対象外」は §4.6 と矛盾するため同期する。 |

§6.14 の合格条件は、現状では**不足**している。次を明文化すべきである。

- **新 float 対新 FP64**：幾何と壁距離は binary64、ビット不一致 **0 件**。接続・番号・境界所属・可視化接続も完全一致。float の初期状態や BC の実数値まで、ビルド間で同一にする要求は加えない。
- **新 FP64 対旧 FP64**：同じ計算順序を維持する今回の変更では、既存データセットの型・形・値・属性の一致を維持する。不一致が出たら、先に変更箇所を特定する。根拠なく `1e-6` 等へ緩和しない。
- **比較器**：必要データの欠損、形状不一致、非有限は FAIL。符号付きゼロまで含むビット比較を行い、両側に同じ NaN があることを合格にしない。`/BCONDS/*`・`/VIZMESH/*`・存在時の `/MESH/RENUMBER_PERM` も対象を明記する。
- **読み込み確認**：旧 float32 格子での20 stepは互換性試験として残す。ただし、それだけでは新変換器の出力を検証しない。**新しい double 格子をソルバが読み、double の写しに保持できること**も別途確認する。20 stepを収束の証拠にはしない。

ビット一致は、完全に独立した同一の double 幾何経路なら**成立を見込める**。現時点で量別の非ゼロ許容差を設ける根拠はない。異なるコンパイラ等の比較が必要になった場合は、今回のビルド型独立性試験と分けて登録する。

第 1 仮説: **案 (b) の主な実装漏れは、double で計算した幾何が共用 `mesh` または `variables` を通って float に戻ることになる。** 確度: 高  
　根拠: `solver_density_cuda/mesh/gmshReader.hpp:2229,2248,2423`、`solver_density_cuda/input/calcWallDistance_kdtree.cpp:177` に、その受け渡しが実在する。これは将来の実装に対する仮説であり、修正後の欠陥を確認したものではない。  
　反証条件: 修正後の経路監査で中間の float 格納がなく、固定入力の float/FP64 変換比較でも、壁距離を含む全幾何のビット不一致が 0 件となること。

第 2 仮説: 残る不一致が壁距離だけなら、探索バックエンド・コンパイル条件の不一致が候補。確度: 中。実際のビルド条件は未確認。  
第 3 仮説: なし。

判別 A/B: **案 (b) 実装後、case/48 の `fp_y1_3um.msh` と両 YAML を固定し、変更点をビルドの状態型 float/FP64 だけにして各1回変換する。ソルバは0 step。** 上記の幾何 dtype・ビット差・整数表を比較する。  
→ **結果 A：全件一致**なら、この入力での「ビルド型が幾何に漏れる」仮説を棄却し、§6.14 の残りへ進む。  
→ **結果 B：1件でも不一致**なら、ビルド型独立性を棄却する。最初に異なる生成段階を追い、許容差の緩和で通さない。

やらない方がよいこと: `geom_float` の一括置換だけで完了扱いすること、float 配列を double に広げて精度を保証したとすること、CUDA の終了エラーを「HDF5 が存在する」だけで成功扱いすること。

呼び出し側の前提への異議: **「案 (b) は2ファイルの型変更で閉じる」「座標 dtype が double なら幾何全体を保証できる」は受け入れない。** また、`cuda_forge` の二重ビルドで時間が「倍」になるという見積もりは未計測であり、案 (a) を退ける主根拠にはしない。

不足情報: 比較に使う3ケースの正確な run パス・入力ハッシュ、旧 FP64 変換器のソース revision／ビルド条件、共用 `node/plane/cell`・`Point`・`setInitial` の変更案。今回指定された閲覧範囲では、それらの内部まで含む実装の安全性は確定できない。
