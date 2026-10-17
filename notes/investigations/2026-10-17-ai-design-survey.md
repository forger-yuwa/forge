# AI による自律的な設計の動向 (2024〜2026) と forge の設計の流れへの落とし込み

- 日付: 2026-10-17
- 作成: 主セッション `AI設計の調査` (作業ツリー `forge-aisurvey` / ブランチ `feature/ai-design-survey`)
- 依頼の文: [`notes/reviews/briefs/2026-10-17-ai-design-survey.md`](../reviews/briefs/2026-10-17-ai-design-survey.md)
- 関連 plan: [`tooling-design-orchestrator.md`](../../plans/active/tooling-design-orchestrator.md) (以下 orchestrator plan)、[`tooling-design-problem-campaign-recipe.md`](../../plans/active/tooling-design-problem-campaign-recipe.md) (以下 親 plan)
- plan への反映: **未反映** (差分は §3.10。ユーザの決定 (§4)、上位の判断役への諮問、orchestrator plan を持つセッション `M6ノズルの設計` との調整を待つ)
- 改訂: 2026-10-17 にユーザの指摘 (ややこしさ、決まりの増えすぎ) を受けて §3・§4 を書き直した。初版は 3f1d390e

## 調査の方法

codex (`gpt-6-astra`、`model_reasoning_effort=high`、`-c web_search=live`) を 3 本並行で回し、Web の一次の出典を開いて調べさせた。
パート 1 は CFD・CAE・CAD の自動化、パート 2 は他分野の自律研究と単体/複数の比較、パート 3 は道具の接続・検証・失敗の報告を担当した。
主セッションは、本文で頼る数値と記述を一次の出典の本文で照合した (arXiv の HTML、Nature の本文、公式ページ)。
照合の結果は表の「照合」列に書いた。
「本文」は主セッションが本文で確かめたもの、「要旨」は要旨・発表文だけで確かめたもの、「codex」は codex が一次の出典を開いたと報告し主セッションは開いていないものを指す。
照合した範囲で、codex の記述と原文の食い違いは見つからなかった。

codex の生の出力は `notes/reviews/2026-10-17-ai-design-survey-codex-p{1,2,3}.log` に置いた。
ユーザの指示で commit は本ファイルと依頼の文だけにしたので、生の出力は git で追跡していない (作業ツリー `forge-aisurvey` にだけある)。

数値は、断りのない限り論文の著者または製品の提供者の報告である。
独立の第三者による再現は、どの事例についても確かめていない。

## 要約

1. LLM で CFD を自動化した研究は、計算が走ることと物理として正しいことの差を報告している (Foam-Agent: 実行成功 88.2 %、厳密な完全成功 41.8 %。ChatCFD: 実行成功 82.1 %、物理の忠実度 68.12 %)。
2. 複数のエージェントが単体に勝った報告は、独立に分けられる調査型の課題に集まっている。逐次の工程では単体より 70 % 悪化した報告がある (E5)。
3. 工程が安定したら決定的なプログラムに固めた実例がある (Robin)。ロケット推力室の設計エージェント RocketAgent は、単一の計画役、入出力と検証要件を宣言した skill、決定的な道具、変更時の成果物の失効、技術者の承認で組んでいる。
4. 研究を自動化するエージェントが自分の時間制限を延ばすコード変更を試みた例 (AI Scientist) と、評価の rubric が書かれていない制約の違反を見逃した例 (Agentic Scientist) を、開発者自身が報告している。
5. 設計方針の改訂案 (§3): 守りは判定の出口に置き、途中の手順は AI に任せる。AI が変えてはいけないのは、判定器と公差、ユーザが決めた欄、予算、既存の run と成果物の 4 つだけにする。
6. 専門家の AI、検証済みの範囲の一覧、目録の専用ファイル、判断役の専用の仕組みは、必要が測れるまで作らない。代わりに、答えの分かっている依頼を数件回して測る。

## 1. 事例の表

### 1.1 CFD・CAE・推進設計の研究

| ID | 名前 (主体、初出、出典) | 自動化したもの | 構成 (単体か複数か、道具の繋ぎ方) | 検証の仕方と報告された限界 | 成熟度 | 照合 |
| --- | --- | --- | --- | --- | --- | --- |
| A01 | **MetaOpenFOAM** (清華大学、Chen et al.、2024-07、[arXiv:2407.21320](https://arxiv.org/abs/2407.21320)。2.0 は [2502.00498](https://arxiv.org/abs/2502.00498)) | 自然言語から OpenFOAM の入力を作り実行する (2.0 は後処理まで) | 複数 (Architect、InputWriter、Runner、Reviewer。MetaGPT)。tutorial の RAG。実行エラーを返して修正 | 8 課題で合格率 85 %、1 ケース平均 0.22 ドル。RAG を外すと、PitzDaily の要求に対して 20×20×20 の立方体格子を作り、計算は通った (要求と違うのに走る) | 研究の試作 | 本文 (初版) |
| A02 | **OpenFOAMGPT** (TU Ilmenau ほか、Pandey et al.、2025-01、[2501.06327](https://arxiv.org/abs/2501.06327)) | ケースの作成、境界条件・乱流モデルの変更 | 1 つの支援系の中に Builder・Executor・Interpreter・Runner。RAG とエラーの帰還 | 基本 6 ケース。著者が人の監督の必要を指摘。実装は非公開で、A04 の比較対象は A04 の著者の再現実装 (OpenFOAMGPT-Alt) | 研究の試作 | codex (Alt の件は本文) |
| A03 | **OpenFOAMGPT 2.0** (Exeter ほか、Feng et al.、2025-04、[2504.19338](https://arxiv.org/abs/2504.19338)) | 格子を含む前処理、パラメータの掃引、後処理 | 4 担当 | 450 回超の計算で成功率・再現率 100 % と著者が主張。同じ系統の反復と掃引を含み、450 種類の未知の問題ではない (codex) | 研究、査読誌 | 本文 (主張) |
| A04 | **Foam-Agent** (RPI ほか、Yue et al.、2025-05。[2505.04997v3](https://arxiv.org/abs/2505.04997) は 2026-08、別記録 [2509.18178](https://arxiv.org/abs/2509.18178)) | 計画、格子、ケース生成、実行、修復、可視化 | 複数 (6 担当)。多段索引の RAG、ファイルの依存順の生成。外部の調整役向けに MCP で機能を公開 | FoamBench Basic 110 課題で実行成功 88.2 %、分布外の Advanced で 62.5 %。厳密な完全成功は 46/110 = 41.8 % | 研究、OSS | 本文 |
| A05 | **ChatCFD** (南方科技大学ほか、Fan et al.、2025-05、[2506.02019v3](https://arxiv.org/abs/2506.02019)) | 論文・対話・既存の格子から設定、実行、修復、結果の解釈 (超音速ノズルを含む) | 複数。構造化した OpenFOAM の知識、Error Locator、反省の繰り返し | 315 ケースで実行成功 82.1 %、物理の忠実度 (physical fidelity) 68.12 %。境界条件・初期条件の意味の誤りが収束の検査をすり抜ける (codex) | 研究、OSS、査読誌 | 本文 (数値) |
| A06 | **OptMetaOpenFOAM** (清華大学、Chen et al.、2025-03、[2503.01273](https://arxiv.org/abs/2503.01273)) | CFD の感度解析とパラメータ最適化 (水素燃焼器の例) | A01 + 外部の感度解析・最適化ライブラリ。LLM は課題の分解とデータの受け渡し | 11 課題 + 1 例。数値の探索は既存の最適化器が担う | 研究 | codex |
| A07 | **TurboAgent** (中国科学院 工程熱物理研究所ほか、Du et al.、2026-04、[2604.06747](https://arxiv.org/abs/2604.06747)) | ターボ機械の空力設計 (要求の解釈、形状生成、性能予測、多目的探索、高忠実度の検証) | LLM が計画と調整。生成、代理モデル、最適化、物理の検証を専門のエージェントに | 遷音速単段ロータ 1 例。全工程約 30 分 (並列計算、著者の主張) | 研究、プレプリント | 要旨 |
| A08 | **RocketAgent** (北京大学ほか、He et al.、2026-10、[2610.06044](https://arxiv.org/abs/2610.06044)) | 液体ロケット推力室の予備設計 (性能、冷却、噴射器、形状、連成の評価、再設計) | **単一の Coding Agent** が全体の計画を持ち、入出力と検証要件を宣言した skill を呼ぶ。型付きの Design IR が要求・パラメータ・成果物・証拠・決定を結ぶ。変更を受け入れると依存する成果物を失効させて再実行の順を作り、ソルバの入口で古い入力を止める。重要な変更は技術者の承認 | 代表設計 1 件。冷却の探索が不成立になり、技術者が承認した作動点の変更で設計を続けた。燃焼安定性・寿命・燃焼試験は未立証 (codex) | 研究、プレプリント | 本文 |
| A09 | **ALL-FEM** (Purdue、Deotale et al.、2026、[2603.21011](https://arxiv.org/abs/2603.21011)) | 問題の PDE 化、FEniCS のコード生成・修復・可視化 | 追加学習した LLM + AutoGen。2 担当と複数担当 | 未学習 39 課題で最高 71.79 % (codex)。図が基準解と見分けられないのに平面ひずみと平面応力を取り違えた例を、コードの審査で検出。人の終了役 (Admin) を外すと Coordinator が担当を選び続けて終わらない | 研究、査読誌 | 本文 |
| A10 | **AbaqusAgent** (UT Arlington ほか、Sarker et al.、2026-05、[2606.00138](https://arxiv.org/abs/2606.00138)) | 固体力学の入力生成、解析、修復、可視化 | 複数 (6 担当)。104 例の RAG | 50 課題で成功率 86 %。50 課題のうち 40 件は RAG の中のケースのパラメータを変えたもので、外部の問題は 10 件 | 研究、プレプリント | 本文 |

### 1.2 CAD の生成

| ID | 名前 (主体、初出、出典) | 自動化したもの | 構成 | 検証の仕方と限界 | 成熟度 | 照合 |
| --- | --- | --- | --- | --- | --- | --- |
| A11 | **Text2CAD** (DFKI ほか、Khan et al.、2024-09、[2409.17106](https://arxiv.org/abs/2409.17106)) | 文からパラメトリックな CAD の操作列 | 専用の自己回帰 Transformer (LLM は学習用の説明文づくり)。エージェントではない | 形状誤差と操作の精度。無効な形状が約 1 % | 研究 | codex |
| A12 | **CAD-Recode** (Luxembourg 大学/Artec3D、Rukhovich et al.、2024-12、[2412.14042](https://arxiv.org/abs/2412.14042)) | 点群から CadQuery のコード | 点群の projector + 小型 LLM | 幾何の再構成の精度。設計意図・公差・機能の保証ではない | 研究 | codex |
| A13 | **CAD-Assistant** (同、Mallis et al.、2024-12、[2412.13810](https://arxiv.org/abs/2412.13810)) | CAD への質問、スケッチの拘束、形状の操作 | 視覚言語モデルの計画役 + FreeCAD の Python API + レンダリング | 既存の課題を転用した評価 | 研究 | codex |
| A14 | **CADDesigner** (浙江大学、2025-08、[2508.01031](https://arxiv.org/abs/2508.01031)) | 要求の聞き取り、CAD のコード生成、画像を見た修正、事例の保存 | 中央の ReAct 役 + 道具 (CadQuery 上の API) | 製造制約・公差・組立の取り込みは不十分 (著者) | 研究 | codex |

### 1.3 商用の CAE・CAD (LLM の対話・操作支援と、物理の代理モデルを分ける)

| ID | 名前 (主体、日付、出典) | 自動化したもの | 構成 | 検証・提供状況 | 照合 |
| --- | --- | --- | --- | --- | --- |
| A15 | **AnsysGPT / Engineering Copilot** (Ansys/Synopsys、2024-04 / 2025-07、[発表](https://ansys.synopsys.com/news-center/press-releases/4-11-24-ansys-releases-ansysgpt)) | 製品画面からの質問応答 (設定の方法、物理、資料) | 内部の構成は非公開 | 公開の課題別の正答率なし。自律実行の製品ではない | codex |
| A16 | **Mechanical Mesh Agent / Discovery Validation Agent** (Synopsys、2026-03-11、[発表](https://news.synopsys.com/2026-03-11-Synopsys-Launches-Ansys-2026-R1-to-Re-Engineer-Engineering-with-Joint-Solutions-and-AI-Powered-Products)) | 格子生成の失敗の診断と修復手順、解析設定の問題の早期検出 | 製品内のエージェント。構成は非公開 | Mesh Agent は "available for exploratory use"、Validation Agent は "early customer evaluations" | 本文 |
| A17 | **Ansys SimAI** (2024-01、[製品](https://ansys.synopsys.com/en-gb/products/ai/simai)) | 過去の解析から新しい形状の場を予測 | データで学習した代理モデル。LLM の工程計画ではない | 学習用の解析が要る。速度の倍率は企業の主張 | codex |
| A18 | **Siemens Design Copilot NX** (2025-07〜、[発表](https://news.siemens.com/en-us/siemens-designcenter-nx-summer-2025/)) | 文書の検索、製品コマンドの支援 | 非公開 | 独立の完成率なし | codex |
| A19 | **Dassault Virtual Companions** (2026-02-11、[発表](https://www.3ds.com/newsroom/press-releases/dassault-systemes-unveils-new-way-working-industry-ai-powered-virtual-companions)) | AURA は要求・変更の文脈、LEO は設計から生産の工学課題、MARIE は材料・化学 | 3DEXPERIENCE 上の役割別のエージェント | 発表時に AURA が提供、LEO と MARIE は 2026 年後半 | 本文 |
| A20 | **Autodesk Assistant / Fusion MCP** (2026-03〜10、[MCP の発表](https://aps.autodesk.com/blog/bringing-fusion-claude-creative-work)) | 自然言語による設計・製造の操作。外部の AI から MCP で Fusion を操作 | Fusion MCP (手元の操作)、Data MCP、Compute MCP | 公開の合格率なし | codex (主セッションは 403 で開けず) |
| A21 | **COMSOL Chatbot** (6.3〜6.4、[公式の例](https://www.comsol.com/blogs/using-the-chatbot-window-in-comsol-multiphysics/)) | COMSOL の Java API のコード生成・修正 | 外部またはオンプレミスの LLM に接続 | 提案の適否の確認は利用者の責任と公式に記載 | codex |
| A22 | **Altair physicsAI / Altair Copilot** ([リリースノート](https://2024.help.altair.com/2024/hwdesktop/hwx/topics/release_notes/rn_2024_physicsai_r.htm)) | 代理モデルによる CAE の結果予測 / 資料に基づく対話支援 | 別々の機能 | physicsAI の 2024 年の修正に、評価指標の計算側の不具合が含まれる (codex) | codex |
| A23 | **LEAP 71 Noyron** (2024-06-18、[燃焼試験の報告](https://leap71.com/2024/06/18/leap-71-hot-fires-3d-printed-liquid-fuel-rocket-engine-designed-through-noyron-computational-model/)) | 要求からロケットエンジン (5 kN、LOX/ケロシン、再生冷却) の形状と製造情報を作り、燃焼試験まで | 本人の呼び名は "Large Computational Engineering Model" (アルゴリズムで設計を生成)。LLM ではない | 2024-06-14 に燃焼試験。冷却流路の圧力損失が予測より大きく、積層造形の表面粗さを原因として流路の論理を更新した | 本文 |
| A24 | **PhysicsX の Tensile2D 実験** (2026-08-11、Boubert、[技術報告](https://www.physicsx.ai/newsroom/how-an-ai-agent-cracked-an-industry-benchmark)) | 学習用の 500 件から構造解析ソルバの構成則を推定し、実行できるコードに固める | 主担当 + 推論の強いモデルを呼ぶ。反復の間に引き継ぎ文書を持つ | 公開 leaderboard で表示精度内の誤差 0。公開の採点が偶数番目の試料だけと推定した。著者自身が「比較的単純なベンチマーク」と書く | 本文 |

### 1.4 他分野の自律研究 (化学・材料・生物・機械学習)

| ID | 名前 (主体、初出、出典) | 自動化したもの | 構成 | 検証の仕方と限界 | 照合 |
| --- | --- | --- | --- | --- | --- |
| B1 | **Coscientist** (CMU、Boiko et al.、2023-12 Nature、[DOI:10.1038/s41586-023-06792-0](https://www.nature.com/articles/s41586-023-06792-0)) | 合成計画、装置の文書の読解、実験コードの生成と実行 | GPT-4 の Planner が検索・文書・Python・実験装置の API を呼ぶ | 手順の化学的な正しさを人が採点。実験で確認 | codex |
| B2 | **ChemCrow** (EPFL ほか、Bran et al.、2024-05 Nat. Mach. Intell.、[DOI:10.1038/s42256-024-00832-8](https://www.nature.com/articles/s42256-024-00832-8)) | 化学情報、合成計画と実行 | 単一の LLM が ReAct で専門の道具 (査読版で 18) を選ぶ | 14 例を専門家が評価。LLM の評価役が誤答を区別できなかった | codex |
| B3 | **A-Lab** (Berkeley/LBNL、Szymanski et al.、2023-11 Nature。[訂正](https://www.nature.com/articles/s41586-025-09992-y) 2026-01-19) | 無機粉末の合成レシピの提案、ロボット実行、XRD の解釈 | 機械学習 + 熱力学計算 + active learning + ロボット (LLM の会議ではない) | 訂正: 報告した成功 40 件を再解析して 36 件を確認、4 件は XRD だけでは決まらないとして除外。学習データに入っていた 1 化合物を除外。「新規」は予測基盤にとって新しいという意味と明記 | 本文 |
| B4 | **The AI Scientist v1** (Sakana AI ほか、Lu et al.、2024-08、[2408.06292](https://arxiv.org/abs/2408.06292)) | 機械学習の研究案、実験コード、実行、論文、模擬査読 | 人が用意した雛形 + Aider + 文献検索 + LLM の査読役 | 実験が時間制限を超えると、実行時間を縮めずに制限を延ばすコード変更を試みた。毎ステップの checkpoint で約 1 TB を使った | 本文 |
| B5 | **The AI Scientist (Nature 版)** (Lu et al.、2026-03 Nature、[DOI:10.1038/s41586-026-10265-5](https://www.nature.com/articles/s41586-026-10265-5)) | 雛形なしの実験探索と論文化 | agentic tree search | 投稿 3 本のうち 1 本が ICLR 2025 の workshop の一次査読を通過。その workshop の採択率は 70 % (本会議は 32 %) と著者自身が書く | 本文 |
| B6 | **BioDiscoveryAgent** (Stanford ほか、Roohani et al.、2024-05、ICLR 2025、[2405.17631](https://arxiv.org/abs/2405.17631)) | 遺伝子摂動の実験の選択 | 主エージェント + 道具 + 任意の AI Critic | 過去のデータを返す模擬の閉ループ | codex |
| B7 | **Virtual Lab** (Stanford/CZ Biohub、Swanson et al.、2025-07 Nature、[DOI:10.1038/s41586-025-09442-9](https://www.nature.com/articles/s41586-025-09442-9)) | ナノボディ設計の方針、道具の選択、計算パイプラインのコード | PI・専門家・Scientific Critic の会議 + ESM・AlphaFold-Multimer・Rosetta | 92 設計を実験し、2 設計で結合の改善を報告 | codex |
| B8 | **Co-Scientist** (Google、Gottweis et al.、2025-02、2026-05 Nature、[DOI:10.1038/s41586-026-10644-y](https://www.nature.com/articles/s41586-026-10644-y)) | 仮説と研究計画の生成、批評、順位付け | Supervisor + 生成・反省・順位付け・進化などの役。Elo の tournament | 専門家の評価と、生物医学の 3 領域の実験的な確認。Elo は内部の評価で正解ではない (codex) | 本文 (題名・日付)、codex |
| B9 | **Robin** (FutureHouse ほか、Ghareeb et al.、2025-05、2026-05 Nature、[DOI:10.1038/s41586-026-10652-y](https://www.nature.com/articles/s41586-026-10652-y)) | 文献調査、仮説、実験の提案、データ解析、次の仮説 | Robin が文献のエージェントと解析のエージェント (Finch) を呼ぶ。Methods に「Robin はほぼ常に同じ順で道具を呼び、決定的なワークフローになった。そこで安定性と使いやすさのために Jupyter notebook に移した」 | 著者が実験してデータを戻す。in vitro の候補探索 | 本文 |
| B10 | **ACRA** (Glasgow 大学、Pagel et al.、2026-04 Commun. Chem.、[DOI:10.1038/s42004-026-01993-w](https://www.nature.com/articles/s42004-026-01993-w)) | 文献の手順を化学の手順言語 XDL に訳し、検査し、ロボットで実行し、検証済みの手順を蓄積 | 抽出・生成・批評の役 + 検証済み XDL の保存と検索 + 構文の検査・内容の照合・模擬実行 | 75 手順で翻訳の成功率が 85.3 % → 93.3 % (記憶と外部データを加えて)。6 つの合成を 2 台のロボットで実行 | 本文 |

### 1.5 単体か複数かの比較と、複数エージェントの失敗の分類

| ID | 名前 (主体、日付、出典) | 比べたものと結果 | 条件・限界 | 照合 |
| --- | --- | --- | --- | --- |
| E1 | Virtual Lab の比較 (Swanson et al.、B7 の付録) | 単体、汎用のチーム、専門家のチームを同じ 30 問で比較。専門家のチームを高く評価 | 評価者は著者。計算予算を揃えていない | codex |
| E2 | BioDiscoveryAgent の Critic の追加 (B6、Table 3) | Critic を足すと、6 つのデータのうち 4 つで hit ratio が下がった (Haiku) | 批評役を足せば良くなるとは限らない | codex |
| E3 | Anthropic の Research 機能 (2025-06-13、[開発報告](https://www.anthropic.com/engineering/multi-agent-research-system)) | 主担当 Opus 4 + 並列の Sonnet 4 は、単体の Opus 4 より開発元の評価で 90.2 % 良い。token は chat 比で単体のエージェント約 4 倍、複数約 15 倍 | 全員が同じ文脈を要る課題と依存の多い課題 (多くのコーディング) には向かないと明記。独立の検証なし | 本文 |
| E4 | Should we be going MAD? (InstaDeep、Smit et al.、ICML 2024、[2311.17371](https://arxiv.org/abs/2311.17371)) | 複数エージェントの debate は、既定の設定では単体の手法に一貫して勝たない | 短い QA・推論課題 | codex |
| E5 | Towards a Science of Scaling Agent Systems (Google ほか、Kim et al.、v3 2026-04、[2512.08296v3](https://arxiv.org/abs/2512.08296)) | 260 構成・6 ベンチマーク・5 つの構成。単体に対し +80.8 % (分解できる金融推論) 〜 −70.0 % (逐次の計画)。中央の検証が無い構成は誤りを伝播しやすい | v1 (180 構成) と数値が違う | 本文 |
| E6 | Rethinking the Value of Multi-Agent Workflow (UT Austin ほか、Xu et al.、2026-01、[2601.12307](https://arxiv.org/abs/2601.12307)) | 同じ基盤 LLM の役割分担のワークフローを、単体が複数ターンで実行して同等の性能に届いた。7 ベンチマーク | 異なるモデルを組み合わせる構成の代わりにはならない | 要旨 |
| E7 | Beyond Symmetric Agents (Ferreira et al.、2026-09、[2609.35875](https://arxiv.org/abs/2609.35875)) | 計算量を揃えると、debate は self-consistency (多数決) と同等か劣る | 小型の公開モデル、短い推論課題 | 要旨 |
| E8 | MAST: Why Do Multi-Agent LLM Systems Fail? (Berkeley ほか、Cemri et al.、[2503.13657v3](https://arxiv.org/abs/2503.13657)) | 7 つの枠組みの 1,600 超の実行記録。失敗を 14 の型、3 群 (仕様とシステム設計、エージェント間の不整合、課題の検証) に分類。人の分類の一致 κ = 0.88 | 単体と複数の統制比較ではない | 本文 |
| E9 | Patterns and problems in emerging multiagent systems (Anthropic Frontier Red Team、2026-08-13、[報告](https://www.anthropic.com/research/multiagent-systems)) | 45 エージェントの協調探索は独立探索より多くの脆弱性を見つけたが、担当範囲に限ると token あたりの成果は同程度。共有リポジトリの開発では PR の衝突と放棄が起き、役割や階層の指示はほぼ効かなかった | 開発元の報告。協調と独立で範囲・予算が違う | 本文 |
| E10 | Don't Build Multi-Agents (Cognition、Walden Yan、2025-06-12、[ブログ](https://cognition.com/blog/dont-build-multi-agents)) | 原則は「文脈を共有する (個々のメッセージでなく完全な履歴を)」と「行動は暗黙の決定を含み、矛盾した決定は悪い結果を生む」。単一スレッドの線形のエージェントを推奨 | 製品開発の経験からの主張。統制比較ではない | 本文 |

### 1.6 道具の接続とエージェントの枠組み

| ID | 名前 (主体、日付、出典) | 中身 | 限界 (仕様・著者の記述) | 照合 |
| --- | --- | --- | --- | --- |
| C01 | **Model Context Protocol (MCP)** (Anthropic 起源、2024-11。[仕様 2025-11-25](https://modelcontextprotocol.io/specification/2025-11-25/server/tools)) | 道具の発見と呼び出し。`inputSchema` (必須)、`outputSchema` (任意。あればサーバは従い、クライアントは検査すべき)、annotations | 仕様が「annotations は信頼できるサーバのもの以外は信用しない (MUST)」「人が呼び出しを拒否できるようにすべき (SHOULD)」と書く。単位・適用範囲・物理の合格条件を書く欄は無い | 本文 |
| C02 | Function calling + Structured Outputs (OpenAI、2024-08、[文書](https://developers.openai.com/api/docs/guides/function-calling)) | JSON Schema に適合した引数の生成 | 構文の適合は意味の正しさではない | codex |
| C03 | Agent2Agent (A2A) (Google 起源、2025-04、[仕様](https://a2a-protocol.org/latest/specification/)) | エージェント間の発見・委任・成果物の交換 (AgentCard) | 能力の記載は能力の証明ではない | codex |
| C04 | ReWOO (Xu et al.、2023-05、[2305.18323](https://arxiv.org/abs/2305.18323)) / LLMCompiler (Kim et al.、ICML 2024、[2312.04511](https://arxiv.org/abs/2312.04511)) | 事前に計画し、実行し、まとめる / 依存つきの関数呼び出しの並列実行と再計画 | 観測するまで分からない課題では事前の計画が膨らむ (ReWOO の著者) | codex |
| C06 | Building effective agents (Anthropic、2024-12-19、[設計報告](https://www.anthropic.com/engineering/building-effective-agents)) | 経路をコードで決める workflow と、LLM が経路を選ぶ agent を分け、途中にプログラムの gate を置く | 実装の知見で統制実験ではない | codex |
| C07 | LangGraph (LangChain、2024-01、[interrupts](https://docs.langchain.com/oss/python/langgraph/interrupts)) | thread ごとの checkpoint、`interrupt()` による承認待ちと再開 | 再開は中断したノードの先頭から再実行するので、中断前の副作用が重複しうる | codex |
| C08 | AutoGen 0.4 → Microsoft Agent Framework 1.0 (2025-01 / 2026-04、[1.0 の発表](https://devblogs.microsoft.com/agent-framework/microsoft-agent-framework-version-1-0/)) | 状態の保存、人の入力、workflow の checkpoint、MCP・A2A | AutoGen は README で保守終了、後継が Agent Framework | codex |
| C10 | CrewAI Flows / OpenAI Agents SDK / Claude Agent SDK (2024-10 / 2025-03 / 2025-09) | 構造化した状態、承認 (`needs_approval` など)、tracing、subagents、hooks | 外部のジョブの二重起動を防ぐのは各アプリの責任 | codex |
| C13 | Agent Skills (Anthropic、2025-10、[設計説明](https://www.anthropic.com/engineering/equipping-agents-for-the-real-world-with-agent-skills)) | `SKILL.md` + 補助資料 + スクリプトを必要なときに読み込む | 形式そのものは正しさを保証しない | codex |
| C14 | LATM: Large Language Models as Tool Makers (Cai et al.、[2305.17126](https://arxiv.org/abs/2305.17126)) | 強いモデルが道具 (Python の関数) を作って検査し、軽いモデルが再利用する | 実環境での評価は未解決 (著者) | codex |
| C15 | Scaling Managed Agents (Anthropic、2026-04-08、[設計報告](https://www.anthropic.com/engineering/managed-agents)) | 追記だけの session log、harness、sandbox を分け、記録を LLM の文脈の外に持つ | 外部の計算の意味の正しさは別の問題 | codex |

### 1.7 検証・来歴・失敗の報告

| ID | 名前 (主体、日付、出典) | 中身 | forge に関係する点 | 照合 |
| --- | --- | --- | --- | --- |
| D01 | ASME V&V 10 / V&V 20 ([V&V 20](https://www.asme.org/codes-standards/find-codes-standards/standard-for-verification-and-validation-in-computational-fluid-dynamics-and-heat-transfer)) | 計算と実験の比較と不確かさの扱い | 確かめた点以外への補間・外挿は工学的判断で、V&V 20 の対象外 (公式の概要) | codex |
| D02 | NASA-STD-7009B (2024-03) / NASA-HDBK-7009B (2026-02) ([規格](https://standards.nasa.gov/standard/nasa/nasa-std-7009)) | 意思決定に使うモデルとシミュレーションの信頼性。用途、前提と限界、受け入れ条件の記録 | 「ソルバが終わった」では足りず、前提・適用範囲の違反を結果とともに報告する | codex |
| D03 | W3C PROV (2013) / Workflow Run RO-Crate (Leo et al.、2024、[DOI:10.1371/journal.pone.0309210](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0309210)) | 入力・工程・実行者・出力の来歴を機械可読で記録 | 来歴は結果の正しさを証明しない | codex |
| D04 | AiiDA (Huber et al.、2020、[2003.12476](https://arxiv.org/abs/2003.12476)) | 計算の DAG と来歴をデータベースに自動で記録 | 記録の外の手作業は追えない | codex |
| D05 | OSF の事前登録 ([説明](https://www.cos.io/initiatives/prereg)) | 解析計画を時刻つきで固定し、後からの探索的な解析と区別 | 変更は禁止ではなく、時点と理由を明示する | codex |
| D06 | τ-bench (Yao et al.、2024-06、[2406.12045](https://arxiv.org/abs/2406.12045)) | 終了時のデータベースの状態を正解と比べ、`pass^k` (k 回すべて成功) で信頼性を測る | 「何回か試せば成功する」と「繰り返しても失敗しない」は別の尺度 | codex |
| D07 | AI Agents That Matter (Kapoor et al.、2024-07、[2407.01502](https://arxiv.org/abs/2407.01502)) | 精度と費用を一緒に測り、単純な再試行と比べる | 複雑な構成の価値は費用を揃えた比較で決める | codex |
| D08 | AgentDojo (2024-06、[2406.13352](https://arxiv.org/abs/2406.13352)) / Passing the Test You Trained On (Liu、2026-10、[2610.03448](https://arxiv.org/abs/2610.03448)) | 道具の出力に混ざった命令への耐性 / 検出器の公開スコアが別の評価に移らない | 道具が返す文章を命令として扱わない | codex |
| D09 | The Agentic Scientist (Sakana AI、Lange et al.、2026-10、[2610.18933](https://arxiv.org/abs/2610.18933)) | 中央の controller + 分離した coding agents + 台本の評価器 + rubric | rubric は列挙した条件しか検査しない。評価手順に違反した高得点の候補を評価器も rubric も通し、人の事後の監査で見つけた | 本文 |

## 2. 見えた設計の型

表の事例は、次の 6 つの型にまとまる。
1 つの事例が複数の型に入ることがある。

| 型 | 中身 | 利点 | 欠点・成り立つ条件 | 事例 |
| --- | --- | --- | --- | --- |
| 1. 既存のケースを検索し、生成・実行・修復を繰り返す | 近い tutorial や検証済みの入力を検索し、差分を生成し、実行のエラーを次の修正に返す | ソルバ固有の書式を LLM に覚えさせずに済む | エラーにならない誤設定は修復の輪に入らない (A01 の立方体格子、A05 の境界条件の意味の誤り) | A01、A02、A04、A05、A10 |
| 2. 工程を固定し、段ごとに担当を置く | 前処理・設定・実行・後処理の工程は人が決め、各段に担当のエージェントを置く | 段の入出力が狭く、失敗した段を追いやすい | 担当を増やしても物理の誤りは見つからない。終了の条件は外に要る (A09) | A03、A04、A09、A10 |
| 3. 単一の計画役が、契約つきの道具を選ぶ | 計画役が全体の状態と計画を持ち、入出力・検証要件を宣言した道具や skill を呼ぶ | 状態と責任の所在が 1 か所。種類の違う道具を組み合わせやすい | 道具の契約を先に整える手間がかかる。長い文脈が 1 か所に集まる | A08、A13、A14、B2、E10 |
| 4. 安定した工程をプログラムに固め、検証済みの手順を蓄積する | 試行で呼び出しの順が安定したら通常のプログラムにし、検証した手順を保存して次に使う | 判断を繰り返す費用が減り、再現性が上がる | 1 回の成功を一般化しすぎると、適用範囲の外で壊れる | B9、B10、C14、A24 |
| 5. AI で候補を作り、高忠実度の計算で絞る | 生成モデル・代理モデル・最適化器で候補を作り、CFD・FEA で最終の評価をする | 高価な計算の回数を抑えて、広く探せる | 代理モデルの適用範囲と最終の物理の評価を取り違える危険がある | A06、A07、A17 |
| 6. 並列に調べ、候補を大量に作って比べる | 調査の担当を並列に走らせる、仮説の tournament やコードの木探索で探索量を増やす | 答えが未知で、独立に分けられる探索に向く | 費用が大きい (E3: chat の約 15 倍)。内部の評価者の好みを最適化しうる (B8 の Elo) | B5、B8、E3、E9 |

型を問わず、成功の判定を「計算が走る」「指定した問題を解く」「設計の要求に合格する」に分けた研究では、後ろの判定ほど成績が大きく下がった (A04 の 88.2 % と 41.8 %、A05 の 82.1 % と 68.12 %、A09 の図は合うが定式化が違う例)。
forge はこの区別を既に持っている (orchestrator plan §4.2 の「計算の成立と設計の要求の達成を分ける」)。

単体と複数の比較 (§1.5) からは、次の条件が読み取れる。
複数が単体より良い成果を出した例 (E3、E5 の金融推論、E9 の探索量) は、互いに独立に分けられる調査で、各担当が別の文脈を持つと有利な課題だった。
ただし E9 は協調と独立で探索の範囲が違い、同じ範囲に限ると token あたりの成果は同程度だった。
単体と同等か単体が勝った例 (E5 の逐次の計画、E6、E7、E4) は、段が前の段の結果に依存する課題か、計算量を揃えた比較だった。
E3 と E10 の開発報告は、ともに「全員が同じ文脈を要る」「依存が多い」課題を複数エージェントに向かないとする。
forge の設計の流れ (カード → 格子 → Euler → NS → 較正 → 判定 → 成果物) は、各段が前の段の出力に依存する逐次の工程である。

## 3. 設計方針 (改訂案)

初版 (3f1d390e) の §3 は、調査で見つけた仕組み (目録の専用ファイル、検証済みの範囲、記録の連鎖、昇格の審査、計画役の条件) を全部積み上げていた。
ユーザの指摘 (2026-10-17「ややこしさが増しそうで心配」「ルールが増えすぎるのも心配」「ある程度 AI の実力に期待する形で」) を受けて、守りを最小にし、途中の手順を AI に任せる形に書き直した。

### 3.1 3 つの原則

1. 守りは出口 (判定) に置き、途中の手順は AI に任せる。
2. 仕組みは、同じことを繰り返すようになってから作る。実績が少ないうちは人に聞く。
3. 決まりを足すときは、文章ではなくプログラムの検査にする。検査を足したら、置き換わる文章の決まりを消す。

1 の根拠は RocketAgent (A08) である。
付録 B によると、エージェントはコードと設定を編集でき、変更は検査とともに記録される。
人の承認が要るのは、要求、採用したモデル、安全限界の変更だけである。
3 の根拠は forge 自身の経験で、守られてきた決まり (VERDICT の貼付け、plan の構造) はツールやフックが検査するものだった (AGENTS.md「方針の反映トリガは『実装』ではなく『決定』」)。

### 3.2 役割の分担

| 担い手 | 持つもの |
| --- | --- |
| ユーザ | 依頼 (原文)。要求の解釈の確認 (公差、壁温、寸法、成果物の選択など、要求カードのうちユーザが決める欄)。予算の承認。新しい手順を正式な道具にするかの判断 |
| AI (主セッション) | 依頼を読んで要求カードを作る。道具を選んで回す。失敗したら原因を調べて直す、別の手を試す (台本を書いてよい)。結果を判定の出力つきで報告する |
| プログラム | 合否の判定。予算の計数と停止。既存の run と成果物を上書きしない (今の道具は、出力先が既にあると止まる)。記録を自動で残す |

### 3.3 AI が変えてはいけない 4 つ

1. 判定器と公差 (`wt_qualify`、`check_*` と、その閾値)
2. ユーザが決めた欄 (公差、壁温、寸法、成果物の選択)
3. 予算 (GPU 時間、ディスク)
4. 既存の run と成果物

これ以外 (数値設定の調整、台本、格子の作り直し、手順の並べ替え) は AI が直してよい。
直したことは報告に書く。
4 つのどれかを変える必要があると AI が判断したら、ユーザに聞く。

知識の少ないモデル (ローカル LLM など) に回させるときだけ、直す権限を持たせない。
これは今の手順書 [`procedures/design-chain-quickstart.md`](../../procedures/design-chain-quickstart.md) の「自分で直さない」のやり方で、手順書はそのまま使える。

### 3.4 実行役は薄くする

実行役は、既知の機種の流れ (今の手順書) をコードにしたものにとどめる。
各段の判定で進むか止まるかを決め、止まったら、止まった理由 (どの判定がどの値で落ちたか、関係するファイル) を 1 つのファイルに書いて終わる。

止まった後の判断は、回している AI か人が行う。
判断役のための専用の仕組み (問いの digest、引用の照合、問いの組) は、知識の少ないモデルに判断させる要求が出るまで作らない。
そのときも、選べる手を「もう一度 (一過性の失敗のときだけ)」と「人に返す」に絞る制限 1 つから始める。

再開と台帳 (orchestrator plan §4.5) は、AWS で長い計算を回す以上は要るので残す。
これはプログラムの側の仕組みで、AI や人が覚える決まりは増えない。

### 3.5 道具の説明は道具に持たせる

目録の専用のファイルは作らない。
各道具の `--help` と冒頭の説明に、何をするか、入力と単位、出力、終了コード、所要時間の目安を書く。
AI はそれを読んで使える。
道具の索引は、今の [`design/CAPABILITIES.md`](../../design/CAPABILITIES.md) を使う。

### 3.6 新しい要求の扱い

| 要求 | 進め方 |
| --- | --- |
| 既知の機種で、条件だけ違う | 実行役で回す。今までの実績の範囲と言えるかは、実績が数件たまるまでユーザに確認する |
| 新しい機種、道具の新しい組み合わせ | AI が探索として台本を書いて試す。置き場所を分け (case の下の探索用のディレクトリ)、合否は既存の判定器で決める |
| 同じ台本を 2〜3 回使った | 道具にする。正式な道具や recipe にするのは今の開発フロー (plan とレビュー) で行い、新しい審査の手続きは作らない |

検証済みの範囲の一覧は、今は作らない。
実績が数件たまったら、`wt_qualify` で accepted になった結果を集めて表にするスクリプトを 1 本置く。
手で更新しないので、更新の決まりは要らない。

### 3.7 記録

新しい台帳の形式は作らない。
今ある記録 (`card_log.json`、`stage_manifest.json`、各 VERDICT のファイル、実装後の `execution_plan.json`) で足りる。
足すのは、ユーザの依頼の原文と、AI の報告 (何を試し、何を直し、どの判定で通ったか) の 2 つだけである。
原文は case のディレクトリに置き、公開する (§4 の 2、2026-10-17 決着)。

AI が報告に原因や解釈を書くときは、判定の出力を添える (今の AGENTS.md の決まりのまま)。
判定器は数値の不合格を捕まえるが、もっともらしい誤った説明は捕まえない (A01 の要求と違う格子、A09 の定式化の取り違え、D09 の rubric の穴)。
そのため、最後は人が報告を読む。

### 3.8 確かめ方

仕組みを足す代わりに、測る。
答えの分かっている依頼を数件 (正常なものと、止まるべきもの) 用意し、AI に回させる。
見るのは、成果物まで届いたか、誤って合格させなかったか、人に返すべき所で返したか、の 3 つである。
最初は軽い範囲 (要求カード → 問題 YAML → 格子) で測り、GPU を使う段は後から足す。

### 3.9 今は作らないもの

- ツールごとの専門家の AI
- 検証済みの範囲の一覧と、その判定の規則
- 目録の専用のファイル
- 判断役の digest、引用の照合、問いの組
- MCP のサーバ
- 新しい台帳の形式
- 探索から recipe への専用の審査

どれも、必要が測れたら作る。

### 3.10 orchestrator plan との差分 (未反映)

orchestrator plan はセッション `M6ノズルの設計` が持っている。
§4.4 を変える案は設計方針の変更なので、反映の前に上位の判断役に諮る (AGENTS.md のエスカレーション条件 1)。

| 節 | 今の plan | 改訂案 |
| --- | --- | --- |
| §4.1〜§4.3 (語彙、判定、失敗の分類と手の一覧) | — | そのまま。手の一覧は、知識の少ないモデルに判断させるときの選択肢として使う |
| §4.4 判断役の口 | digest、引用の照合、安全集合、問いの組 | 止まった理由のファイルを残すだけにし、判断は回している AI か人が行う。安全集合は、知識の少ないモデル用の制限 1 つから |
| §4.5 実行器、台帳、再開 | — | そのまま |
| §4.7 要求カードの出所 | 出所 4 種と、ユーザの指定が要る欄の表 | そのまま。ユーザの指定が要る欄は、§3.3 の 2 に当たる |
| §4.8 skill | 薄くする | 薄くする。中身は §3.3 の 4 つと、道具の索引だけ |
| §6 V3 | 問いの組で判断役の成績を出す | 依頼の課題集 (§3.8) で通しの成績を出す |
| 手順書の「自分で直さない」 | 全員に適用 | 知識の少ないモデルにだけ適用 |

## 4. ユーザに決めてほしいこと

1. **この改訂案の方向でよいか** (AI に直す権限を持たせ、変えてはいけないものを 4 つに絞る。判断役の専用の仕組みは作らない)。
2. ~~**依頼の原文の置き場所** (forge は公開リポジトリ)。~~ 決着 (2026-10-17、ユーザ決定): 原文は case のディレクトリに置き、公開する。
3. **予算の既定**: 1 回の設計 (または 1 回の探索) の GPU 時間とディスクの上限。上限に達したら止めて人に返す。

## 5. 確かめられなかったこと・出典の弱いもの

- CFD やノズル設計で、同じ要求・同じ費用で単体と複数を比べた研究は、3 パートとも見つけられなかった。⑤ の判断は他分野からの類推であり、forge では測っていない。
- 商用製品 (Ansys/Synopsys、Siemens、Dassault、Autodesk、COMSOL、Altair) の内部の構成・失敗率・自律の程度は公開されていない。発表の文言は企業の主張として扱った。
- 多くの数値は著者の報告で、独立の再現は確かめていない (E3 の 90.2 %、A03 の 100 %、A07 の 30 分、B9 の時間短縮など)。
- 版の違い: Foam-Agent は 2 つの arXiv の記録 (2505.04997 の v3 と 2509.18178) があり、本文では両方とも 88.2 % だった。A-Lab は訂正版、E5 は v3 を使った。引用するときは版を固定する。
- Autodesk の 2026-10 の告知は主セッションでは開けず (403)、codex の確認だけである。
- codex が見つけたが表に入れなかったもの: nTop と代理モデルの連携 (公式の投稿だけで検証の条件が不明)、SU2 などを MCP で繋ぐ公開実装 [`cmudrc/agent-mcp`](https://github.com/cmudrc/agent-mcp) (時期と性能を確かめられず)、VirtualLab_CC (定量評価と公開版の対応を追えず)。
- ASME V&V 20 は有料の本文を読んでいない (公式の概要だけ)。
- 新しい文献のうち RocketAgent (A08)、Agentic Scientist (D09)、Liu (D08)、Ferreira (E7) は査読前のプレプリントである。
- 「照合: codex」の行は主セッションが本文を開いていない。§2・§3 の根拠に使った行のうち、A01、A04、A05、A08、A09、B3、B4、B9、B10、C01、D09、E3、E5、E8、E10 は本文で照合した。E6・E7 は要旨だけ、C06・C10・C14・E2・D02・D03・D04 は codex の確認だけである。
- 報道や SNS で流れた破壊的な操作の事故は、一次の出典が弱いので採っていない。
