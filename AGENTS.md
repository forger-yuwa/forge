# Forge Agent Rules

このファイルは forge リポジトリの開発・運用ルールの正本であり、Claude (CLAUDE.md 経由) と GitHub Copilot (copilot-instructions.md 経由) の双方から参照される。

## 文書構成

各文書の役割は次の通り。詳細手順は本ファイルに重複記載せず、該当文書を参照する。

| 文書 | 役割 |
| --- | --- |
| `AGENTS.md` (本ファイル) | 全ルール・方針の正本。常時参照される基準 |
| [`procedures/`](procedures/README.md) | forge を使う/開発するための**運用・開発ハンドブック**。以下が主要文書 |
| [`procedures/calculation-workflow.md`](procedures/calculation-workflow.md) | 計算ケース準備・メッシュ生成/変換・`forge` 実行の標準手順 |
| [`procedures/divergence-and-startup.md`](procedures/divergence-and-startup.md) | 発散の主因と安定起動の手順 (新規計算 / 発散時に必ず参照) |
| [`procedures/solver-settings.md`](procedures/solver-settings.md) | `convMethod` / `limiter` などの数値設定リファレンス |
| [`procedures/su2-cross-check.md`](procedures/su2-cross-check.md) | 同一メッシュ・同一 BC で SU2 と比較し forge 固有の問題を切り分ける手順 |
| [`procedures/development-environment.md`](procedures/development-environment.md) | 開発環境とビルド (Docker / WSL native) の方針 |
| [`procedures/coding-conventions.md`](procedures/coding-conventions.md) | ソース構成・C++/CUDA 命名規約・ビルド/テスト実行手順 |
| [`procedures/verification/`](procedures/verification/README.md) | 検証ケース選定 (`README.md`) と各標準検証ケースの個別手順 |
| [`plans/`](plans/README.md) | 変更単位の設計判断文書。`active/` (検討中・進行中) / `accepted/` (現役の設計判断) / `archived/` (superseded・終了)。着手前に参照する基準文書 |
| [`notes/`](notes/README.md) | 調査メモ・作業ログ。`investigations/` (技術調査・サーベイ) / `sessions/` (使い捨て作業プロンプト) |
| [`methods/`](methods/index.md) | 現在の仕様と解説 (機能単位 `<area>/`)。「なぜそうしたか」は `plans/` 側 |
| `papers/` | 参照文献 (PDF)。理論・スキームの一次資料。本文では引用元として参照 |

## 言語方針

- 文書本文・コード内コメントは日本語で記述する。
- 識別子・関数名・スキーム名・ファイルパスなどコード由来の語は原語 (英語) のまま `` `code` `` 表記とする。
- git commit メッセージは英語の命令形で記述する (識別子・コード由来語は原語のまま)。

## 計算・実行ルール

計算の実行方法、ケース準備、メッシュ生成、メッシュ変換、`forge` の起動、Docker 経由の Gmsh/ParaView 利用について回答するときは、まず `procedures/calculation-workflow.md` を参照し、その手順に合わせて案内すること。計算手順の本文はこのファイルに重複記載しない。

`solverConfig.yaml` の `convMethod`・`limiter` などの数値設定を変更・確認するときは、必ず `procedures/solver-settings.md` を参照すること。設定値の意味を記憶や推測で判断しないこと。

**`mesh.bndFirstOrder` は使用禁止**。新規 config に書かないこと。既存 run から複製するときも必ず落とすこと。このフラグは (1) 対流再構成だけでなく `viscousFlux` が読む速度勾配も 0 にするため**粘性応力を破壊**し、(2) `bnode_flag` が「いずれかの bcond の CV」なので疑似 2D (1 セル厚) メッシュでは spanwise 境界が全ノードを覆い**全域に効いてしまう**。したがって安定化しても交絡した副作用であり、切り分けにも対策にも使えない。詳細と削除計画は [`procedures/solver-settings.md`](procedures/solver-settings.md) の該当節と [`plans/active/architecture-bndfirstorder-removal.md`](plans/active/architecture-bndfirstorder-removal.md) を参照。

**新規計算の投入時、および計算が発散 (NaN / 残差爆発) したときは、まず [`procedures/divergence-and-startup.md`](procedures/divergence-and-startup.md) を参照し、そのチェックリストと段階起動手順に従うこと。** 発散の大半は物理やメッシュでなく投入設定 (初期値が入口流れと不整合 / 超音速向け BC を亜音速に使用・出口逆流 Pt/Tt 未設定 / 初手から 2 次移流・乱流・no-slip / 非直交での free-stream 桁落ち) が原因であり、易しい条件で収束させてから引き継ぎ計算で段階的に上げるのが基本。

エージェント自身が計算検証を実行する場合も、既存の `run_*` ディレクトリをそのまま使い回さず、必ず複製した新しい `run_*` ディレクトリで実行すること。
また、計算を実行した場合は `residual_history.csv` から `residual_history.png` も生成して残すこと。

**計算結果ディレクトリの明示 (必須)**: 後からどの検証がどこにあるか追えるよう、計算の「投入時」と「結果まとめ時」の両方で、結果が出力される `run_*` ディレクトリを**リポジトリルートからの相対パスで明示**すること。

- **投入時**: 計算を起動する応答内で、各 run の出力先パスを列挙する (例: `case/05.sod_shock_tube/run_0004_slau_tracer_n2n2/`)。複数 run を投入する場合は run ごとに何を変えたか (スキーム・設定差分) を 1 行で添える。
- **結果まとめ時**: 結論や比較表とあわせて、その結論の根拠となった `run_*` パスと主要成果物 (`res_*.h5` / `residual_history.png` / 後処理図 `*.png` 等) を明示する。複数 run を比較したときは「どの図/数値がどの run のものか」を取り違えないよう run パスを併記する。
- 一時 run や破棄予定の run を作った場合も、その旨とパスを明記する (放置して所在不明にしない)。

**case の run 索引 (必須・恒久)**: 応答での明示はスクロールで流れて後から辿れないため、**各 case ディレクトリの `README.md` に「## 計算 run 一覧」表を維持する**。これが「どの検証がどの run か」の恒久的な一次情報であり、結果をまとめる応答ではこの表へのポインタも示す。

- **対象 case に `README.md` が無ければ新規作成**し、ケース概要に続けて run 一覧表を置く。
- **run_* を作成・破棄したら同表を必ず同期**する。列は最低限: `run_*` ディレクトリ名 / 目的・主要設定差分 / 主要結果・成果物 / 状態 (`active` / `破棄予定` / `ref`=入力リファレンスとして保持)。
- **命名規約**: `run_NNNN_<slug>` とし、`<slug>` は目的が一目で分かる語にする (例 `run_0015_prod_double`)。**既存 run と番号が衝突しない連番**を使い、衝突や系列違いが起きたら表の備考でグルーピングを明記する。
- 表の本文は簡潔に (1 run = 1 行)。詳細な考察は `plans/` の該当計画に書き、表からはその計画へリンクする。

**NaN / 発散チェック (必須)**: 計算が「完走 (exit 0)」しても妥当とは限らない。**初期ステップから発散している**ことが多いため、結果を使う前に必ず NaN/Inf の有無を確認すること。

- **確認タイミング**: ① 投入直後 (序盤数十ステップで NaN/Inf が出ていないか早期確認)、② 結果まとめ前 (最終 `res_*.h5` と全 `residual_history.csv` を確認)。長時間 run では途中でも定期的に見る。
- **確認内容**:
  - `residual_history.csv` の `rms_*` 列に NaN/Inf が無いか、**最初に NaN になったステップ** (`step 0`/`step 1` からなら初手発散) を特定する。
  - 最終および中間 `res_*.h5` の `VALUE/*` (特に `ro`,`P`,`T`,`roe`,勾配 `dUxdx` 等) に NaN/Inf が無いか、値が物理的か (例: 静圧 ≤ 全圧、`ro>0`、`T>0`) を確認する。
- **発散していた場合**: 「収束した」と報告しない。最初に NaN が出たステップ・場所 (どの境界・どのセル) を切り分け、原因 (BC のゼロ割・IC とBC の不整合・CFL 過大・メッシュ不良 等) を特定してから対処する。`max cfl` のログだけ見て安定と判断しないこと (CFL が有限でも残差が NaN のことがある)。
- 早期確認の例: `python3` で `residual_history.csv` の先頭数十行の `rms_ro` と `res_0/res_1` 相当出力の `isnan` を見る。専用ツールがあれば優先する。

**収束確認 (必須・NaN チェックとは別)**: 「NaN が無い」=「収束した」ではない。結果を使う/報告する前に、各保存量の残差が**実際に下がりきっているか**を必ず確認すること。

- **判定は手作業の目視ではなく `solver_density_cuda/tools/check_convergence.py <run_dir> ...` を実行して行う** (本ルールの実体化ツール)。全残差列の低下桁数・トレンド・NaN を判定し `PASS / NOT CONVERGED / DIVERGED` を返す。**「収束した」「一致した」と報告する応答には、このツールの VERDICT を根拠として必ず貼ること** (これを怠り `rms_ro` だけで「収束」と誤報告した事例があるため、ツール経由を必須とする)。未収束なら「未収束」と明記し、未収束のトランジェント同士の場の比較を「一致」と表現しない。
- **`rms_ro` だけで判断しない**。`residual_history.csv` の **全列** (`rms_ro`,`rms_roUx`,`rms_roUy`,`rms_roUz`,`rms_roe`、RANS 時は `rms_roK`,`rms_roOmega`) のトレンドを見る。`rms_ro` が低くても、運動量 (特に軸対称の `rms_roUy`) や乱流 (`rms_roK`/`rms_roOmega`) が**下げ止まり・横ばい・上昇**していれば未収束。実例: 軸対称 SST で `rms_ro`≈3e-5 でも `rms_roUy`≈1e-2 停滞・`rms_roK` 増大=近軸が未収束だった ([architecture-axisym-axis-singularity.md](plans/accepted/architecture-axisym-axis-singularity.md))。
- **残差プラトーは「収束」ではない**。下げ止まる場合は、積分量 (massflux/推力/出口諸量) が**定常化**しているか、場が**発達しきっている**かを併せて確認する (リミットサイクルの可能性)。
- **場の発達も確認する** ([develop-flow-before-reporting] と同趣旨): 残差が下がっていても、境界層・乱流・衝撃などが発達途中なら結果は使えない。中間 `res_*.h5` を時系列で見て、注目量が定常化したことを確認する。
- 外部ソルバ (SU2 等) でクロスチェックする場合の収束確認も同様。手順は [`procedures/su2-cross-check.md`](procedures/su2-cross-check.md) を参照。

**準定常確認 (必須・収束確認とは別)**: 残差が下がっていても (またはプラトーでも)、**報告する派生量そのもの** (衝撃位置・上下非対称・CL/CD・massflux・推力・peak μt/μ・出口諸量 等) が**定常化 (頭打ち) しているか**は別問題である。**残差プラトー ≠ 量の定常化**であり、**過渡ピーク ≠ 飽和値**である。これを怠り、過渡ピークの量を定常値として報告した事例があるため、量の報告にはツール経由の確認を必須とする。

- **判定は `solver_density_cuda/tools/check_quasisteady.py <run_dir> [--quantity shock,asym,...]` を実行して行う** (本ルールの実体化ツール)。全 `res_*.h5` スナップショット時系列から対象量を計算し、末尾の頭打ちを `STEADY / DRIFTING / OSCILLATING / TRANSIENT-UNSETTLED` で判定する。**衝撃位置・非対称・CL/CD・massflux 等の派生量を「○○だ」と報告する応答には、このツールの VERDICT を必ず貼ること**。
- **単一スナップショット・短窓で量を報告しない**。`DRIFTING` / `TRANSIENT-UNSETTLED` を「定常」「飽和」と表現しない。`OSCILLATING` (リミットサイクル) は瞬時値でなく**平均±振幅**で報告する。
- **過渡が減衰しきる長さまで回す**。短い run では非対称・衝撃位置などの**過渡ピークを定常的な値と誤認する** (例: 擬似衝撃波の上下非対称は visc=0 Euler で過渡 0.25→減衰 0.05 だが、~12k step では 0.25 を定常偏りと誤判定した。十分長く=量が頭打ちするまで回す)。量が `DRIFTING` なら「未だ動いている」と明記し run を伸ばす。
- 強い偏り流など中心線量が破綻する場合は対象量を選ぶ (`--quantity asym` 等)。詳細閾値は `--tail/--drift/--osc`。
- 外部ソルバ比較や「一致」の主張でも同様: 比較する両者がともに `STEADY` であることを確認してから「一致」と述べる。

**メッシュ変更後の restart (必須)**: メッシュを変えた (quad↔tri↔構造化, 解像度変更) ときは **uniform 初期値から計算を始めない** (超音速/衝撃波/SST は uniform IC から step 数回で発散する)。`solver_density_cuda/tools/interp_field.py SRC.h5 新メッシュ入力.h5` で過去の収束済み場を**最近傍interpolateして cross-mesh restart** する (保存量+roK/roOmega+スカラー輸送を移植、wall_dist は移植せず新メッシュの値を使う)。同一メッシュの restart は `restart_field.py`。

**メッシュ品質チェック (計算前・必須)**: メッシュを HDF5 化したら計算投入前に必ず `solver_density_cuda/tools/check_mesh_quality.py <mesh.h5>` で品質を確認する。**アスペクト比 ≤ 1000、スキューネス ≤ 0.9 を目標**とし、`VERDICT: FAIL` のメッシュは投入しない。近壁細分化で AR が増えやすいので接線長と第一セル厚のバランスを取る (高 Re では y+~1 と AR≤1000 が両立しないことがあり、その場合 y+~30-80 + `wallTreatmentSST=1` を選ぶ)。詳細は [`procedures/calculation-workflow.md`](procedures/calculation-workflow.md) の「メッシュ品質チェック」。「メッシュできた/収束した」と報告する応答には品質 VERDICT も併記する。

forge の結果が「軸対称・乱流・近軸で forge だけ妙な値になる」ようなときは、推測で結論づけず [`procedures/su2-cross-check.md`](procedures/su2-cross-check.md) の手順で **同一メッシュ・同一 BC の SU2 と比較**して切り分けること。

開発環境に関する既定ルールは `procedures/development-environment.md` を参照すること。通常の開発は Docker コンテナを基本とし、NVIDIA 提供の GPU 速度計測・プロファイリングツールを使う場合は、Docker 内でうまくいかないケースを想定して WSL native または Linux native の手順を優先する。さらに、**計算速度の評価・比較・プロファイル作業は原則 native で行う** (Docker は計測値が不安定でツールも揃わないため、速度評価の基準環境としない)。詳細は `procedures/development-environment.md` の「速度評価・ベンチマークは native を既定とする」を参照。

コード変更時の検証ケースの選び方と個別の確認手順は `procedures/verification/README.md` を参照すること。

`plans/` 配下の `.md` は、変更単位の設計判断文書 (実装計画) として扱うこと。対象タスクに対応する計画が存在する場合、実装や調査を進める前にまずその `.md` を参照すること。計画一覧は [`plans/README.md`](plans/README.md) を使う (`active/` が進行中、`accepted/` が現役の設計判断、`archived/` が superseded・終了)。コード変更を伴わない技術調査・サーベイは `notes/investigations/`、使い捨ての作業プロンプトは `notes/sessions/` に置く。

実装方針を変更する場合は、先に対応する `plans/active/*.md` を更新してから実装に移ること。計画未更新のまま実装だけを先行させないこと。**実装を伴わない方針決定も同じ**で、決まった時点で plan に反映する (「[開発フロー](#開発フロー-新規機能設計変更)」の「方針の反映トリガは『実装』ではなく『決定』」節を参照)。

ユーザーが solver のコード構造、アーキテクチャ、モジュール責務について質問した場合は、`methods/architecture/overview.md` を既定の参照先とすること。

## ドキュメント整備方針

forge の理論的背景と実装解説は `methods/` 配下に機能単位 (物理プロセス / 数値処理) のサブディレクトリで蓄積する。運用方針の本文は `methods/README.md` を参照し、本ファイルには重複記載しない。

- 新規に理論背景や実装解説を書くときは `methods/<area>/` (例: `gradient/`, `limiter/`, `convection/`, `diffusion/`, `time_integration/`, `boundary/`, `poisson/`, `architecture/`) の配下に追加する。スキーム名 (Roe / SLAU / KEEP など) は対応する機能ディレクトリ内のセクションまたはファイルとして扱い、新たなトピック単位ディレクトリは作らない。
- 既存ドキュメントがある話題は新規ファイルを作らず、該当ファイルを更新する。
- 理論・実装に関する質問に答えるときは、該当する `methods/<area>/*.md` を既定の参照先とする。
- ファイルの追加・改名・削除を行った場合は `methods/index.md` の目次を必ず同期させる。
- 数式は KaTeX 記法 (インライン `$...$`、ブロック `$$...$$`)、図は mermaid フェンスまたは画像で統一する。
- `methods/` は現在の仕様・解説 (理論+実装を機能単位 `<area>/` に集約。`theory.md`/`implementation.md` の強制分割は廃止し、規模が大きい領域のみ分ける)、`plans/` は変更単位の設計判断、`notes/` は調査・作業ログ、`procedures/` は運用・開発の手順/ルール (計算手順・発散対処・設定・SU2比較・開発環境・コーディング規約・検証ケース)、という棲み分けを守る。`.github/` は GitHub プラットフォーム設定 (`copilot-instructions.md` 等) のみ。

## 開発フロー (新規機能・設計変更)

バグ修正と軽微な変更を除く「新規機能追加」「数値スキームや設計方針の変更」を行うときは、必ず次の 4 ステップを順に踏むこと。

1. 該当する [`methods/<area>/`](methods/) の現在仕様ドキュメントを更新する (理論・実装の両面。領域が小さければ 1 ファイル内の節で分ける)。
2. [`plans/_template.md`](plans/_template.md) を雛型に `plans/active/<area>-<short-slug>.md` を作成し、`related_docs` で 1. のファイルをリンクする。[`plans/README.md`](plans/README.md) の `active` 一覧にも追記する。
3. 以上 2 つが揃ってから実装に着手する。

完了時には次を行うこと。

- 計画の `status` を `done` に更新し、`## 変更ログ` に実装・検証結果を追記する。**ファイルを `plans/active/` から `plans/accepted/` へ移動する** (superseded になった場合は `plans/archived/`、後継を本文に明記)。ファイル名は変えない (リンク維持)。
- [`plans/README.md`](plans/README.md) の一覧表も同期させる (移動元・移動先の両セクション)。
- 1. で触った docs と [`methods/index.md`](methods/index.md) の整合性を確認する。

例外 (本フローを要しないもの):

- typo やコメント修正のみ
- 1 ファイル内で閉じるバグ修正 (振る舞い変更なし)
- リネーム/象徴のない単純リファクタ (振る舞い同一)
- docs のみの記述補正

判断がつかないときはフローを踏む側を選ぶこと。

### 方針の反映トリガは「実装」ではなく「決定」

上のフローは実装に着手するときの手順だが、**実装を伴わない設計判断が対話の中だけで完結して失われる**事故が
実際に起きた (⑤ SERN の作動点定義: 方針を対話で確定させたのに plan には一行も入らず、残作業リストも
`notes/sessions/` にしか無かった)。原因は「実装に移る前に更新する」という条件が **実装しない turn では発火しない**
ことにある。守られているルール (run パス明示・VERDICT 貼付・case README の run 一覧) は、いずれも
**応答に痕跡が残る形**をしている。plan 反映も同じ形に揃える。

- **決まった時点で書く**: 方針・設計判断・優先順が決まったら、実装の有無を問わず**その応答の中で**
  対応する `plans/active/*.md` に反映する。**対話の中だけで決めて plan に書いていないものは「決めていない」**
  とみなす (次のセッションには残らない)。
- **応答に反映先を明示する** (run パス明示・VERDICT 貼付と同じ扱い): 方針・優先順・設計判断を述べた応答には、
  反映した plan のファイルと節番号を書く。反映していないなら「plan 未反映」と明記する。
- **残作業の正本は plan の残作業表** ([`plans/_template.md`](plans/_template.md) §5.1): 各 `plans/active/*.md` は優先順つきの残作業表を持つ。
  「未解決・次にやること」の実体はここに置き、`notes/sessions/` の引き継ぎ文書には**写しとポインタだけ**を置く
  (`notes/sessions/` は使い捨ての置き場なので、恒久的な残作業を置くと腐る)。
- **決着した未確定事項は消さない**: `## 未確定事項` の項目が決着したら削除せず、取り消し線 +
  「決着 (日付, 参照節)」を付ける。判断の履歴を残す。

構造の存在確認は `python3 solver_density_cuda/tools/check_plans.py [PATH ...]` で行う (残作業表・変更ログ・
`plans/README.md` 記載の有無のみを見る機械的な lint)。`plans/active/*.md` を編集すると PostToolUse フック
(`hook_plan_todo_gate.py`) が編集したファイルだけを検査して不足を知らせるので、**触った plan から順に直す**。
フックの登録は `.claude/settings.json` (git 追跡対象。`.claude/` の他は追跡外)。これは Claude Code のみが読む
仕組みなので、Copilot 側には本節の記述そのものが規範として効く。

## コミット・push 運用

実装または修正が一段落し、検証 (ビルド成功 + 該当検証ケースの結果が妥当) まで確認できたら、その機能・計画単位で git commit し、現在の feature ブランチへ push すること。本節は検証済みマイルストーンでの commit/push を事前承認するものとして扱い、都度の確認は不要とする。

- **commit タイミング**: 機能・計画単位でまとめて commit する。段階検証ごと・小さな編集ごとには commit しない。plan を伴う作業では plan の `status: done` 化と整合させる。
- **commit 対象**: `solver_density_cuda/` などのソース、`methods/`、`plans/`、`notes/` といった意味あるソース・文書変更を対象とする。`case/` 配下の検証 run 成果物 (`res_*.h5`, `*.xmf`, `residual_history.csv` / `.png`, `forge_run.log` など) は commit に含めない。検証 run の入力 config (`solverConfig.yaml`, `bcondConfig.yaml` 等) は、リファレンスとして意図的に残したい場合のみ明示的に `git add` する。
- **push**: commit 後は現在の feature ブランチへ push する。`main` へ直接 commit / push しない。`main` 上にいる場合は先に feature ブランチを切ってから作業する。
- **commit メッセージ**: 英語の命令形で記述する (例: `Add Menter SST turbulence model with dilatation correction`)。
- 無関係な既存の作業ツリー変更を巻き込まないよう、commit 前に `git status` / `git diff` で対象を確認すること。
