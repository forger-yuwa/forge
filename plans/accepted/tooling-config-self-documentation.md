# solverConfig.yaml を読めるようにする (注釈生成ツール + 既定値ログ)

## メタ

- **area**: `architecture`
- **status**: `done`
- **related_docs**:
  - `procedures/solver-settings.md` (数値設定リファレンス)
  - `procedures/recommended-settings.md` (推奨レシピの正本)
- **related_plans**: `plans/active/architecture-bndfirstorder-removal.md` (廃止キーの扱い)
- **created**: `2026-09-13`
- **owner**: `Claude`

## 1. 目的

`solverConfig.yaml` にはコメントが一切なく、キーの意味・既定値・廃止済みかどうかが config を見ても分からない。
一方で config に手でコメントを書くと、run ディレクトリの複製 (case/16 だけで 100 以上) で古い注釈が運ばれ、
既定値や意味が変わっても更新されない。**正本はコードに置いたまま、必要なときに注釈を生成する**形で読めるようにする。

## 2. スコープ

- **やる (A)**: `solver_density_cuda/tools/config_doc.py` — コードからキー・型・既定値・説明を抽出し、
  `check` (未知/廃止/節違い/非既定/必須欠落の点検)、`annotate` (注釈つきの写しを別ファイルに出力)、
  `template` (全キー入り雛形)、`list` (Markdown 表) を提供する。
- **やる (A)**: 説明の無いメンバに `solverConfig.hpp` の行末コメントを追加し、抽出率を 100 % にする。
- **やる (B)**: 省略されたキーを既定値のまま使ったことを起動ログに `[default]` 行として残す。
- **やらない**: 実行される `solverConfig.yaml` そのものへのコメント書き込み (複製で腐るため)。
- **やらない**: ソルバ側での未知キー検出 (キー一覧を C++ に二重持ちすることになるため、ツール側で見る)。
- **やらない**: 条件付き必須 (`physProp.chemistry.mechanismFile` は `enabled: 1` のときだけ要る 等) の検査。
  ソルバ自身が起動時に理由付きで弾く (`solverConfig.cpp` の `requires` / `is required` の throw) ので、
  同じ条件をツールに写すと二重の正本になって腐る。

## 3. 関連 docs と前提

- 推奨値の正本は `procedures/recommended-settings.md`。ツールが出すのは**コード上の既定値**であり、推奨値ではない。
- 使用禁止キー (`mesh.bndFirstOrder`) は `AGENTS.md` のルール由来。コードは受け付けるのでツール側の内蔵リストで持つ。

## 4. 設計方針

**抽出元をコードに限定する (腐らせないための中心判断)**。

- キー・型・既定値: `input/solverConfig.cpp` の読み出し式から正規表現で取る。4 通りの読み方に対応する。
  1. `getOptionalValidatedValue<T>(node, "key", default, "section")` → 省略可能キー + 既定値
  2. `getValidatedValue<T>(node, "key", "section")` → 必須キー
  3. `node["key"].as<T>()` (if-else で既定を書く古い形) → 既定は `(コード側で分岐)`
  4. `for (... : node["key"])` → 配列キー
- 節名は第 3 引数 (`"time.deltaT"`) を最優先し、無ければ `auto deltaT = config["time"]["deltaT"];` の束縛から復元する。
  config 側もドット付きフルパスに平坦化して突き合わせるので、`time.deltaT.cfl` と `time.cfl` を取り違えない。
- 説明: `input/solverConfig.hpp` のメンバ宣言に付いたコメント。**行末コメント + その桁に揃った継続行**を第一候補、
  無ければ**宣言と同じ字下げの直前ブロック**を使う。桁で見分けるのは、右端に流した継続行と「次のメンバの前置き」が
  見た目で区別できないため (これを誤ると隣のキーの説明が混ざる)。
- 廃止キー: コード自身の拒否メッセージ (`"Key 'x' in 'y' is no longer supported"`) から取る。
  `recommended-settings.md` §9 の表は自然文なので機械抽出しない (「`thermoHrefTemp` 未指定が発散要因」のような行を
  「キーが廃止」と誤読した)。ルールで禁止されたキーだけツール内の内蔵リストに置く。

**抽出漏れを「無い」ことにしない**: 解析できなかった読み出しは行番号付きで `coverage` に出し、`list` / `check` も
件数を表示する。**0 件でないかぎり一覧は完全ではない**という扱いにする (codex 指摘 M1/M6)。

**照合はフルパス完全一致**: `space.keepDissCoeff` のように節を間違えたキーは、たとえ同名キーがトップレベルに
あっても `[節違い]` として弾く。`detectNaN` のように複数箇所で読まれるキーは、実在する読み出し位置を
それぞれ登録することで両方を通す (親をたどる曖昧照合はしない。codex 指摘 M2)。

**説明の正しさはコードが保証しない**: 説明は人が書いたコメントなので、追加するときは**値の消費箇所まで当たる**。
実際に `gpu` を「GPU 番号」と誤記していた (正しくは 0/1 フラグ。`boundaryCond.cpp` が 1 以外をエラーにする)。
同様に `initial` は変換器専用、`time.last.control` は現在どこからも参照されていない (codex 指摘 M3)。

**生成物の安全性**: `annotate` は**元のテキストに一切触れず、キーの行の前にコメント行を挿入するだけ**にする。
値を読み書きし直す方式は捨てた (`species: [N2, NO]` の `NO` が YAML 1.1 の真偽値として `false` になり、
空の節 `output: {}` が null に潰れた)。挿入位置は正規表現ではなく **YAML 構文木のキー位置** (`yaml.compose` の
`start_mark`) から決める。行を見るだけだと、ブロックスカラー (`|-`) の中の `drive: mesh.h5` をキーと誤認して
文字列の内側に `#` を挿し込み、複数行のフロー形式の入れ子を取り違える。
検査は 2 本立て: (1) 挿入行を外すと元テキストに 1 文字違わず戻る、(2) YAML として読んだ値が元と一致する。
(1) だけではブロックスカラー内への挿入を見逃す。出力先が入力と同一ファイルなら拒否する。
`template` は廃止・使用禁止キーを出さず、既定値が無いキーは `REQUIRED` / `SEE_CODE` の明示プレースホルダにする
(codex 指摘 M4/M5)。

**(B) 既定値ログ**: `getOptionalValidatedValue` が既定に落ちたとき `[default] 'key' in 'section': value` を出す。
明示されたキーの行は従来書式のまま (ログを読む既存手順を壊さない)。`FORGE_CONFIG_LOG_DEFAULTS=0` で抑止できる。
**これはヘルパーを通った省略の記録であって最終的な実効値ではない**: `discretization` のようにヘルパーを通らない
既定、トップレベル `detectNaN` による上書き、node での `gradLSQ` 強制はここに出ない (codex 指摘 m7)。

## 5. 実装ステップ

1. `solver_density_cuda/tools/config_doc.py` を追加 (抽出・4 サブコマンド)。
2. `solver_density_cuda/input/solverConfig.hpp` に行末コメントを追加 (説明の無い 28 メンバ)。
3. 配列読み・ローカル変数受けで宣言にコメントを置けない 4 キーはツール内の `KEY_DESC` に書く。
4. `solver_density_cuda/input/solverConfig.cpp` の `getOptionalValidatedValue` に `[default]` ログを追加。
5. `procedures/solver-settings.md` にツールの使い方を追記。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 |
| --- | --- | --- |
| 1 | 死にキーの掃除 | 実測で見つかった `turbulence.kInf` / `turbulence.omegaInf` / `time.last.time` / `time.implicit.nLoop` / `mesh.nodeWallViscGradFlux` はソルバが読まない。生産 config から落とす (別コミット) |
| 2 | 旧 config の必須欠落 | `case/04.laval_nozzle/run_slau` 等の古い run は現ソルバの必須キーを欠く。参照用に残すなら README に「現ソルバでは起動しない」と明記 |
| 3 | CI 的な常用 | 新規 run 投入前に `config_doc.py check` を回す運用を `calculation-workflow.md` に入れるか検討 |
| 4 | 実効値のログ | `[default]` はヘルパー経由の省略のみ。**`detectNaN` のトップレベル上書きは最終値がログに一切残らない**。解決後の最終値を 1 箇所にまとめて出す (`[config effective]`) のは別途 |
| 5 | `time.last.control` | 必須なのに未参照 (`endTimeControl` がどこからも読まれない)。廃止するか終了条件として実装するか決める |

## 6. 検証

- **網羅性 (M1/M6)**: `config_doc.py coverage` が「解析できなかった読み出し **0 件**」であること。
  件数が 0 でなければ一覧は不完全として扱う。
- **期待値つき試験 (M6)**: 仕込みの config 1 本で次を同時に確認する
  — 節違い (`space.keepDissCoeff`) / 未知キー (`turbulence.kInf`) / 配列キー (`mesh.wallDistExtraPhysIDs`, 誤検知しない) /
  廃止キー (`mesh.bndFirstOrder`) / 引用符・バックスラッシュを含む文字列の往復。
- **生成物 (M5)**: `annotate` の出力から挿入コメント行を外すと**元テキストに 1 文字違わず戻る**こと、かつ
  **YAML として読んだ値が元と一致する**こと。出力先に入力自身を指定したら拒否 (exit 2) すること。
  `template` の出力が YAML として読め、廃止キーを含まないこと。
- **注釈の中身 (result-2 M2)**: フロー形式を含む実 config で、正しいキーに正しいフルパスの注釈が付き、
  正しいキーに `[節違い]` を付けないこと。
- **説明の正しさ (M3)**: 追加した説明は値の消費箇所 (`variables.cpp` / `main.cpp` / `setInitial.hpp` 等) まで当たって確認する。
- **C++ 側 (B)**: `g++ -std=c++17 -fsyntax-only` を通し、`solverConfig::read()` だけを呼ぶ最小ハーネスを
  実 config に対して実行して「明示キーの行が従来書式のまま」「省略キーが `[default]` として出る」
  「`FORGE_CONFIG_LOG_DEFAULTS=0` で `[default]` が 0 行になる」を確認する。数値計算経路は触らない。
- **既存 config の一掃**: 既存 run の config に `check` をかけ、出た指摘がすべて実在の問題であること。
- **記録**: 実行結果・対象一覧・指摘の内訳は
  [`notes/investigations/2026-09-13-config-doc-verification.md`](../../notes/investigations/2026-09-13-config-doc-verification.md) に残す。
  試験入力と期待値はツール内の `SELFTEST` / `SELFTEST_TEXT` にあり、`config_doc.py selftest` でいつでも再実行できる。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| result (2) | `2026-09-13` | [`notes/reviews/2026-09-13-tooling-config-self-documentation-result-2.md`](../../notes/reviews/2026-09-13-tooling-config-self-documentation-result-2.md) | GO-with-changes, C0/M2/m2 | **全件採用** (下記) |
| result | `2026-09-13` | [`notes/reviews/2026-09-13-tooling-config-self-documentation-result.md`](../../notes/reviews/2026-09-13-tooling-config-self-documentation-result.md) | GO-with-changes, C0/M2/m3 | **全件採用** (下記) |
| plan | `2026-09-13` | [`notes/reviews/2026-09-13-tooling-config-self-documentation-plan.md`](../../notes/reviews/2026-09-13-tooling-config-self-documentation-plan.md) | GO-with-changes, C0/M6/m3 | **全件採用** (下記) |

**採否** (すべて採用・対応済み):

| 指摘 | 内容 | 対応 |
| --- | --- | --- |
| M1 | 抽出が現行コードを網羅していない (`mesh.wallDistExtraPhysIDs` / `output.extraFields` / `physProp.chemistry.mechanismFile`、入れ子型、`LESmodel`/`DESmode` の動的廃止列挙) | 束縛追跡 (`const YAML::Node ch = physProp["chemistry"]` を含む) と入れ子型に対応。`for (const char* old_key : {...})` から廃止キーを抽出。`coverage` で未解析 0 件を保証 |
| M2 | 節の照合が親をたどるため `space.keepDissCoeff` を見逃す。条件付き必須も見逃す | 照合をフルパス完全一致に変更。条件付き必須は**ツールでは検査しない**方針を §2 に明記 (ソルバが起動時に弾くため二重の正本にしない) |
| M3 | `gpu` の説明が誤り (GPU 番号ではなく 0/1 フラグ) | 訂正。あわせて `initial` (変換器専用)、`convMethod` (0/1/2 の対応)、`cfl` (定常では効かない)、`nStepOuter`、`dt_min`/`dt_max`、`time.last.control` (未参照) を消費箇所まで当たって訂正。§4 に「説明の正しさはコードが保証しない」を明記 |
| M4 | `template` が使用禁止キー (`mesh.bndFirstOrder`) を再投入する | 廃止・禁止キーを除外。既定値の無いキーは `REQUIRED` / `SEE_CODE` のプレースホルダに |
| M5 | `annotate` が文字列をエスケープせず値を壊す。`-o` に入力自身を指定すると元 config を上書き | YAML シリアライザ + 往復検証 (一致しなければ出力しない)。入出力が同一ファイルなら拒否 |
| M6 | §6 の検証が不具合を素通しする | §6 を網羅性 0 件・期待値つき試験・往復検証・C++ ハーネス実行に差し替え |
| m7 | `[default]` は実効設定の完全な記録にならない | §4 とツール/文書の説明を「ヘルパー経由の省略記録」に限定。最終値の記録は §5.1 #4 へ |
| m8 | `recommended-settings.md` 自体が死にキー `kInf`/`omegaInf` を推奨している | 推奨レシピを `kInit`/`omegaInit` に訂正し、適用条件 (IC に `roK`/`roOmega` が無いときだけ効く) を追記。`methods/turbulence/implementation.md` の旧キー列にも注記 |
| m9 | plan/result レビュー統合は免除条件に当たらない。`methods/` 更新も無い | 本レビューを `plan` 段として記録し `result` 段は別途実施。`methods/architecture/overview.md` §7.1 を更新 |

**result 段 (2 回目) の採否** (すべて採用・対応済み):

| 指摘 | 内容 | 対応 |
| --- | --- | --- |
| M1 | 行を正規表現で見る注釈挿入が、ブロックスカラー (`\|-`) の中身をキーと誤認して**文字列の内側に `#` を挿す**。原文復元検査は通ってしまう | 挿入位置を **YAML 構文木のキー位置** (`yaml.compose` の `start_mark`) から決める方式に変更。検査に「YAML として読んだ値が一致するか」を追加 (原文復元だけでは見えないため) |
| M2 | 複数行のフロー形式の入れ子を取り違え、正しい `time.deltaT.dt_min` に誤った `[節違い]` を付ける。フロー形式内のキーに注釈が付かない | 同じ構文木からフルパスを取得。1 行に複数キーが並ぶ行には、その行の前にフルパス付きの注釈をまとめて置く |
| m3 | `detectNaN` の最終値をログで確認できるという説明が誤り (上書き時は何も出さない) | `procedures/solver-settings.md` と `methods/architecture/overview.md` を「上書き後の最終値はログに残らない」に訂正。実効値ログは §5.1 #4 |
| m4 | plan §4/§6 が旧方式 (YAML シリアライザ + 読み戻し) のまま | §4 と §6 を実際の方式 (コメント挿入 + 2 本立て検査) に更新 |

**result 段 (1 回目) の採否** (すべて採用・対応済み):

| 指摘 | 内容 | 対応 |
| --- | --- | --- |
| M1 | `annotate` が化学種名 `NO` を YAML 真偽値として読み `false` に書き換える。往復検証も変換後どうしを比べていて通ってしまう | **値を読み書きし直す方式をやめた**。元テキストの行の上にコメントを挿入するだけにし、検証も「挿入した行を外すと元テキストに 1 文字違わず戻るか」に変更。`NO` / 空節 / 引用符・バックスラッシュを `selftest` に追加 |
| M2 | `check` が節をドット連結するため、トップレベルの `"time.deltaT"` という 1 キーと入れ子 `time: {deltaT: ...}` を同一視。`mesh.bndFirstOrder: {}` のような空辞書も見逃す | パスをタプルで保持する方式に変更。辞書そのものも点検対象にし、`[型違い]` (スカラーのキーに辞書) と未知の節を検出。廃止・型違いの下は掘らず、未知の節の下は「正しい節が分かる子」だけ出す |
| m3 | 空の節 (`output: {}`) を含む config を注釈生成できない (exit 2) | M1 の方式変更で解消 (テキストを保つので空節も原文のまま) |
| m4 | `keepDissCoeffMax` の説明が適用領域を逆に書いている | `convectiveFlux_keep_d.inc.cuh` の `max(keepDissCoeff, ransFrac*keepDissCoeffMax)` を確認し「RANS 帯で使う上限」に訂正 |
| m5 | 検証の実測・対象一覧・期待値の記録が無い | `config_doc.py selftest` (期待値つき 12 項目) を追加し、実行結果・対象 config 一覧・指摘の内訳を `notes/investigations/2026-09-13-config-doc-verification.md` に保存。§6 から参照 |

## 7. 影響範囲

- 追加: `solver_density_cuda/tools/config_doc.py`
- 変更: `solver_density_cuda/input/solverConfig.hpp` (コメントのみ)、`solver_density_cuda/input/solverConfig.cpp` (ログ追加)
- 変更: `procedures/solver-settings.md` (使い方)
- 計算結果への影響なし (既定値ログは標準出力のみ)。

## 8. 完了条件

- [x] 関連 `procedures/` を更新済み (`solver-settings.md` に「config を読む・点検する」節、`recommended-settings.md` の死にキー訂正)
- [x] `methods/architecture/overview.md` §7.1 と `methods/turbulence/implementation.md` の旧キー列を更新
- [x] 実装・検証完了 (§6。記録は `notes/investigations/2026-09-13-config-doc-verification.md`)
- [x] codex レビュー 3 回 (plan / result / result-2) を §6.1 に記録し、Critical 0 / Major 10 / Minor 8 をすべて採用
- [x] `status` を `done` に変更し §9 に変更ログ
- [x] `plans/accepted/` へ移動、`plans/README.md` を同期

## 9. 変更ログ

- `2026-09-13` — 初稿 + 実装 (A: 抽出ツールと hpp コメント、B: `[default]` 起動ログ)。
- `2026-09-13` — codex plan レビュー (C0/M6/m3) を全件採用。抽出を束縛追跡・入れ子型・動的廃止列挙に対応させ
  `coverage` を追加 (未解析 0 件)、照合をフルパス完全一致に、`annotate` を往復検証 + 上書き拒否に、
  `template` から禁止キーを除外。誤説明 (`gpu` ほか 6 件) を消費箇所まで当たって訂正。
  `recommended-settings.md` の死にキー推奨 (`kInf`/`omegaInf` → `kInit`/`omegaInit`) と
  `methods/architecture/overview.md` §7.1 を更新。
- `2026-09-13` — codex result レビュー 2 回目 (C0/M2/m2) を全件採用。注釈の挿入位置を YAML 構文木から決める方式に
  変更し (ブロックスカラー内への `#` 挿入・複数行フロー形式の取り違えを根絶)、検査を「原文復元 + YAML 値一致」の
  2 本立てに。`detectNaN` 上書きの最終値がログに残らない件を文書に明記。`selftest` を 16 項目へ。
- `2026-09-13` — codex result レビュー (C0/M2/m3) を全件採用。`annotate` を「元テキストにコメント行を挿入する」
  方式へ作り直し (`species: [N2, NO]` の `NO` が false になる事故を根絶)、config の照合をパスのタプルに変更して
  入れ子の取り違えと空辞書の見逃しを解消。`selftest` (期待値つき) と検証記録を追加。`keepDissCoeffMax` の説明を訂正。
