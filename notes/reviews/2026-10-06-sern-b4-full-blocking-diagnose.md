# codex 諮問 (diagnose): sern-b4-full-blocking

- **brief**: [`notes/reviews/briefs/2026-10-06-sern-b4-full-blocking.md`](../../notes/reviews/briefs/2026-10-06-sern-b4-full-blocking.md)
- **plan**: [`plans/active/tooling-sern-mesh-blocking.md`](../../plans/active/tooling-sern-mesh-blocking.md)
- **date**: 2026-10-06
- **commit**: `b29c6ecc` (feature/sern-design)
- **codex**: effort `xhigh`, 10.6 min, rc=0
- **結論**: **次は、全体の接続契約を先に定めた B4a として、MOC 接続模型の `t_te=0.005H` をメッシュだけで実証する。**
- **extra**: `plans/active/tooling-nozzle-sern-3d.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表です。**B4 は修正して着手可。ただし、接続模型の PASS を全体格子・runner・CFD の受入へ流用することは認めません。**

| 重大度 | 対象・採否 | 根拠と対案 |
|---|---|---|
| **Major** | 外部領域を最後に付け足す段取り：**修正採用** | 接続模型は側方領域の上を `vehicle` 壁で閉じています。一方、現行の全体模型は機体幅外の流体を上流から出口まで連続させています。[接続模型:199](/home/sano/work/forge-sern-design/case/46.sern_design/cad/hex_junction_model.py:199)、[全体模型:505](/home/sano/work/forge-sern-design/design/forge_design/meshing/mesh_sern3d.py:505)。**全体の固体・流体・共有面を先に定義し、上流・機体周囲・ベース後流・出口は一体で接続検証する。** |
| **Major** | ランプ／カウル関数の置換だけで MOC 対応：**却下** | 現実装は `YC=0`、`YCL=-TC` を前提とし、流れ方向の曲線も端点と中点の３点で再構成しています。[形状:119](/home/sano/work/forge-sern-design/case/46.sern_design/cad/hex_junction_model.py:119)、[曲線生成:338](/home/sano/work/forge-sern-design/case/46.sern_design/cad/hex_junction_model.py:338)。**MOC の輪郭、傾斜したカウルの壁法線、厚さ・テーパ、関連ブロックの写像を一緒に定義する。** |
| **Major** | `H1_END=t_te/5` と既存 `--scale` で解像条件を保証：**却下** | `march` は第一間隔の後から拡大し、後流の指定区間に上限を設けません。[間隔生成:85](/home/sano/work/forge-sern-design/case/46.sern_design/cad/hex_junction_model.py:85)。読み取り実行では、`t_te=0.5 mm`、第一間隔 `0.1 mm` でも、直後 `t_te` に掛かる区間で最大 **0.172 mm**。また `--scale` は明示した `H1_END` を変更せず、`scale=0.5` では設定上の `G` が **1.44** になります。[細分処理:547](/home/sano/work/forge-sern-design/case/46.sern_design/cad/hex_junction_model.py:547)。**後流区間全体の間隔上限と、端面を含む細分規則を実装し、出力実座標で検査する。** |
| **Major** | runner への接続をメッシャの呼び替えとして扱う：**却下** | 初期場は `NJ/nz/jm/k_sw` による節点番号の算術に依存します。[初期場:107](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern3d.py:107)。品質検査も終了コードだけで、ここには node 変換後の閉性検査がありません。[変換・検査:171](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern3d.py:171)。**領域情報・節点対応を返す backend 契約を設け、厳密な品質 PASS と閉性 PASS を必須にする。** |
| **Major** | 同名タグなら帳簿互換：**却下** | 模型の physID は自動割当です。[物理群:377](/home/sano/work/forge-sern-design/case/46.sern_design/cad/hex_junction_model.py:377)。runner は固定 ID を使い、新３タグを集計せず、出力欠落も黙って飛ばします。[力集計:285](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern3d.py:285)。**固定 ID、壁 BC、出力、力の帰属を共通の表から生成し、必要面の欠落をエラーにする。** |
| **Major** | 245 万節点の実績から 16 GB に載ると判断：**要再検証** | 旧模型だけで変換 RSS **9.0 GB**。全体領域も、今回の薄い後縁の解像も含みません。[規模実績:436](/home/sano/work/forge-sern-design/plans/active/tooling-sern-mesh-blocking.md:436)。**全ブロックの分割数から事前集計し、実際の変換経路でピークメモリを測る。** |
| **Minor** | 初版では断面隅フィレットを追加しない：**採用** | 接続模型の `PROFILE` は非ゼロ半径を含みます。流用するだけでは「追加なし」になりません。**y–z 断面の半径を明示的にゼロにし、既存の x–y 輪郭の `ramp_fillet` は別物として保持する。** |

採用する段階分けは以下です。**B4a～B4c は CFD なし**で進めます。

| 段階 | 内容 | 次段へ進むゲート |
|---|---|---|
| **B4a** | 全体の接続表を定義し、まず短い MOC 接続模型で有限厚後縁を検証 | 内壁・厚さ・端面位置が明確。下記 A/B で目標 `t_te` が合格。新しい間隔検査も合格 |
| **B4b** | 上流、ランプ全長、プルーム、機体上面・側面・ベース、外部境界をまとめて追加 | 全域の固体／流体分類、共有面、出口への接続が正しい。局所解像を保持し、遠方側を粗くした格子で共通ゲートに合格 |
| **B4c** | 生産用分布と３段階の細分列を確立 | 同一形状、端面を含む実際の細分、壁層・後流解像、メモリ予算が全水準で成立 |
| **B4d＋B5** | `blocks` backend、BC、restart、帳簿、`prepare_info` を統合 | タグと BC の完全対応、初期場の領域・組成基底、必要面の出力、入口積分、力の帰属が検証済み |

共通ゲートは、[plan §6](/home/sano/work/forge-sern-design/plans/active/tooling-sern-mesh-blocking.md:787) の条件をそのまま機械判定します。100% hex、欠落・重複ゼロ、内部面２所有者、外部面１所有者＋１タグ、float32 後の全頂点 Jacobian 正、双対体積正、全 CV 閉性 ≤ 1e−5、skew ≤ 0.90、規定の AR、壁別第一内部点距離・層数、継ぎ目を含む実測間隔比 ≤ 1.2です。`SOFT-PASS` と評価不能は不合格です。後流の `Δx ≤ t_te/5` と端面の壁解像は別々に検査します。

CFD は B4d＋B5 後の B6 で開始します。既決の３作動点、新形状の格子比較、同一実効設定の 20,000 step・10,000＋10,000 step の判定窓を維持してください。４係数の準定常性、全残差、床到達、局所 y₁⁺ を確認し、`m10_on` の低温点が新しいベースや側端へ移った場合も失敗です。途中段階の短い CFD を受入証拠にしません。

規模は長さ比から外挿せず、各ブロックの `nₓ nᵧ n_z` を合計し、節点は共有頂点・辺・面を一度だけ数えて求めます。旧実績の約 **3.7 kB/節点**を仮に使うと、200 万で約 7.4 GB、300 万で約 11.1 GB。これは変換器の概算で、常駐プロセスや余裕を含みません。現 runner は変換中も生成した配列を保持しています。[prepare:150](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern3d.py:150)

**予算超過時は、大容量 RAM 環境での生成・変換を第一選択**にします。不要配列の解放と検査配列の分割処理は併施できますが、変換器の大改修を B4 の前提にはしません。後縁・壁層・必要領域を削って載せる案は却下します。なお、固定分割数のまま下流の分布を緩めても、節点数は減りません。

タグと入力は次の契約に揃えます。

- 既存 [physID 1～18](/home/sano/work/forge-sern-design/design/forge_design/meshing/mesh_sern3d.py:26) を保持し、例えば `sidewall_end=19`、`cowl_side=20`、`cowl_base=21` を追加。新３面はノズル側の帳簿に含めます。
- 側壁厚みの外側に残るランプ高さの帯は、全体の固体定義に従って `vehicle`／内部面などへ分類します。現在の `ramp` のままでは、圧力は幅で分ける一方、摩擦は全面をノズルへ加算する不整合が生じます。[集計:301](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern3d.py:301)、[摩擦:331](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern3d.py:331)
- 入口の運動量・圧力項を実面積で積分し、理想推力の規格化も同じ入口定義に合わせます。現行の矩形積は残せません。[入口項:311](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern3d.py:311)
- `mesh3d.backend` は `legacy` 既定、`blocks` は選択式を維持。`mesh3d` 内で物理形状と格子設定を分離し、厚さ・テーパ・機体幅・端面位置と、第一間隔・層数・成長率を混ぜません。既存 `mesh` との優先順位を一度だけ解決し、未対応キーはエラーにします。
- 実効形状、単位、輪郭の識別情報、境界位置、分割数、各ゲートを `prepare_info` に保存します。現在、`L_up` や `ext_top` は `mesh` から読むため、単に `mesh3d` へ書いても効きません。[入力処理:133](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern3d.py:133)

結論: **次は、全体の接続契約を先に定めた B4a として、MOC 接続模型の `t_te=0.005H` をメッシュだけで実証する。**

第 1 仮説: **局所の19ブロック構成を維持し、輪郭写像と間隔生成を直せば、MOC＋有限厚後縁 `0.005H` を成立させられる。** 確度: **中**

根拠: 現模型は既に側壁跡・カウル跡を流体ブロックで埋める構成を持ちます。[ブロック定義:188](/home/sano/work/forge-sern-design/case/46.sern_design/cad/hex_junction_model.py:188)。一方、実証済みなのは二次ランプ・平坦カウルであり、MOC と目標後縁厚への適用は未確認です。

反証条件: 同じ MOC・分割方針で厚い後縁は全ゲートを通るのに、`0.005H` だけが反転、閉性、壁距離、間隔比のいずれかで落ちること。その場合、現在の写像をそのまま採用する案は棄却し、局所接続を見直します。

第 2・第 3 仮説:

- **第２：全体化の主な接続リスクは、機体幅外の流体を壁で閉じること。確度: 高。** 模型と全体模型で上面の所有関係が異なるためです。有限厚側壁は `W/2` の外側、現行機体は内側にあるので、側壁上端と機体下面の接合幅も明示が必要です。
- **第３：全体格子のメモリ増加は、局所の細分が共有辺の分割数を通じて外部領域へ伝わることに強く左右される。確度: 中。** 現コードは断面の同値類単位で分割数を共通化しています。[分割数:260](/home/sano/work/forge-sern-design/case/46.sern_design/cad/hex_junction_model.py:260)。全体構成での増加量は未測定です。

判別 A/B: **変更する入力は後縁厚 `t_te` だけ。**

- **A:** `t_te=0.02H`
- **B:** `t_te=0.005H`
- 共通: 実 MOC 輪郭、`L_sw=0.8H`、`L_cowl=1.2H`、下流端 `L_cowl+0.3H`、断面フィレットゼロ。内壁、上流板厚、テーパ区間、その他の分割設定は固定します。
- 小さい B の厚さを基準に、両者へ同じ厳しい間隔制約を適用します。後流 `0.02H` の区間全体で `Δx ≤ 0.001H` とし、各端面の第一内部点距離も同じ指定にします。
- **各１回のメッシュ生成・変換、CFD は０ step。** 見る量は全頂点 Jacobian、全 CV 閉性、タグ、skew／AR、壁別距離、実測間隔比、節点数、ピーク RSS です。

→ **両者 PASS:** 目標厚で現局所構成を使えると判断し、B4b へ進む。  
→ **A PASS・B FAIL:** 「同じ構成・写像で目標厚まで対応できる」という第１仮説を棄却。  
→ **A も FAIL:** 薄厚化だけを原因とする説明を棄却し、共通の MOC 写像・分布を修正する。全体化には進まない。

やらない方がよいこと:

- 旧端面間隔 `0.25 mm`、旧 `--scale`、模型の既存 PASS をそのまま流用する。
- `cowl_side` を `L_sw～L_cowl` だけに限定する。模型では側壁より下の板側面が上流区間にも露出しています。[タグ条件:208](/home/sano/work/forge-sern-design/case/46.sern_design/cad/hex_junction_model.py:208)
- 全体領域や帳簿が揃う前に、部分模型の CFD で生産受入を判断する。
- フィレット追加、ソルバ変更、物理形状の粗化を同時に行い、有限厚化の効果と混ぜる。

呼び出し側の前提への異議:

1. **「MOC 内壁を保持」の対象が曖昧です。** 旧メッシャのカウル実壁は中間線から `±tk/2` にあります。`cowl_xy` をそのまま新内壁にすれば、後縁以外も変わり得ます。[旧厚さ配置:208](/home/sano/work/forge-sern-design/design/forge_design/meshing/mesh_sern3d.py:208)。保持する輪郭を座標として固定し、旧実壁との差を記録してください。
2. **側方 `2.50H` は絶対位置です。** 現実装は `W/2+Z_ext`。外向き側壁厚を加えた位置からさらに `Z_ext` を足すと領域が変わります。[側方座標:191](/home/sano/work/forge-sern-design/design/forge_design/meshing/mesh_sern3d.py:191)。新形状でも領域独立性が成立するかは B6 の確認事項です。
3. **B1b 完了と、方式全体の受入完了は別です。** B0c の実測間隔検査、B2、B3 は残っています。また B1b の保存済み判定は `C_M` の評価幅 **5.23e−4 > 5e−4**、`VERDICT: 閾値超過あり → 記録して諮る` です。[判定記録](/home/sano/work/forge-sern-design/notes/investigations/2026-10-05-sern-b1b/B1B_VERDICT.txt)。変換器修正を再び疑う理由にはしませんが、全係数の同等性確認済みとも扱いません。

不足情報: 全体の機体半幅と側壁上端の接合形状、保持するカウル内壁の座標定義、厚さ方向・テーパ則、全体の分割数と実測 RSS が不足しています。参照された B1b run の生データはこの作業木になく、保存済み判定記録までの確認です。

ファイル変更・forge 実行はしていません。**plan 未反映**です。呼び出し側で `tooling-sern-mesh-blocking.md` の §4、§5.1 B4–B6、§6 に反映してください。
