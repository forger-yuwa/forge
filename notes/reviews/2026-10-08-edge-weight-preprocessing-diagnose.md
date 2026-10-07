# codex 諮問 (diagnose): edge-weight-preprocessing

- **brief**: [`notes/reviews/briefs/2026-10-08-edge-weight-preprocessing.md`](../../notes/reviews/briefs/2026-10-08-edge-weight-preprocessing.md)
- **plan**: [`plans/active/convection-zero-thickness-edge-reconstruction.md`](../../plans/active/convection-zero-thickness-edge-reconstruction.md)
- **date**: 2026-10-08
- **commit**: `3e38c41b` (feature/sern-design)
- **codex**: effort `high`, 3.0 min, rc=0
- **結論**: **前処理への移管を採用し、次は格子と `w` の結び付けを仕様化して、同数節点の別格子を誤って受理しない読込 A/B を先に通す。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0 / Major 4 / Minor 1）

| 重大度 | 論点 | 採否・根拠・対案 |
|---|---|---|
| Major | 前処理への移管と格子との整合 | **移管は採用。ただし現状の検査契約では不足。** 長さ・有限性・範囲・`w` 単独のハッシュでは、同数節点の別格子や節点番号の変更を検出できない。[plan:42](/home/sano/work/forge-sern-design/plans/active/convection-zero-thickness-edge-reconstruction.md:42)。**節点順序を保持した格子の意味内容のハッシュと `w` のハッシュを組にし、読込時に双方を再計算して照合する。** |
| Major | restart・段間継承 | **`VALUE` 外への配置は採用。読取元と来歴の確定方法を補う。** 格子は `meshFileName`、状態は `valueFileName` から別々に読む。[main.cpp:1627](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:1627)、[main.cpp:1748](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:1748)。`w` は必ず **`meshFileName` の `/AUX/…`** から読む。同一格子の継続では宛先の `w` を保持し、別格子では宛先で再生成する。各段の起動時に実際に読んだ格子・`w` の識別を確定して `stage_key` と評価来歴へ渡す。 |
| Major | 重みの適用位置と端点の保証 | **0〜1 の表現は採用。「乗算を入れれば旧介入と同じ」は却下。** 元の再構成の後に `contactBlend`・`lowMachThornber`・フォールバックがあり、既存の単面介入はその後にある。[SLAU:208](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:208)、[SLAU:279](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:279)。**既存介入位置で、補正済みの旧面速度を節点速度へ縮約する**契約に統一する。`w=1` は旧値をそのまま使用、`w=0` は節点値を直接代入、中間値だけ `u_i + w(u_f,old−u_i)` とする。 |
| Major | タグ検査・空集合・恒等フィールドの責任分担 | **前処理での空集合拒否は採用。汎用読込器での「ゼロ節点数 > 0」必須化は却下。** `w≡1` の受入試験と矛盾し、中間値だけのフィールドも拒否してしまう。[plan:49](/home/sano/work/forge-sern-design/plans/active/convection-zero-thickness-edge-reconstruction.md:49)、[plan:65](/home/sano/work/forge-sern-design/plans/active/convection-zero-thickness-edge-reconstruction.md:65)。端抽出ツールと SERN の prepare がタグ・E・発火数を検査し、ソルバは汎用フィールド契約を検査する。旧 `{tags, rings}` 形式が残る plan も `{field}` に統一する。 |
| Minor | ParaView での確認 | **目的は採用。HDF5 へ追加するだけでは未達。** 現行 XDMF は登録された変数を `VALUE/…` へ参照する。[output.cpp:247](/home/sano/work/forge-sern-design/solver_density_cuda/output/output.cpp:247)。前処理で `/AUX/w_recon_vel` を参照する **`Center='Node'` の XDMF** も生成し、元の節点番号で表示できることを確認する。 |

格子との結び付けは、次の契約まで具体化すべきです。

- `w` は ghost を含まない実節点順の一次元配列。型・長さ・有限性・範囲を検査し、不正値を clamp して受理しない。
- 格子署名には、少なくとも離散化・次元・節点順の座標・内部面の接続・境界の `physID` と元の境界面接続を含める。座標だけでは、座標が一致する別 ID の節点を区別できない。実際に SERN の継承コードにもこの取り違えの記録がある。[runner_sern.py:736](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern.py:736)。
- ハッシュの入力は、配列名・shape・規定した型と byte order・値を直列化したものとする。`VALUE`、生成時刻、HDF5 の圧縮・配置は除外する。**HDF5 ファイル全体の SHA256 を格子識別に使わない。** メモリ内の確認でも、同じ `w` を非圧縮／gzip で保存すると配列バイトは一致したが、ファイルサイズは 2448／4543 byte、ファイルハッシュは不一致だった。
- タグ名から `physID` への対応も記録する。現行 HDF5 の境界グループは数値 `physID` が正本で、名前だけの照合では足りない。[gmshReader.hpp:2673](/home/sano/work/forge-sern-design/solver_density_cuda/mesh/gmshReader.hpp:2673)。端抽出ツールは対応表を明示的に受け取り、対象面の存在と接続を検査する。
- `stage_key` は属性に書かれたハッシュを無検証で転記せず、読込時に検証した識別を使う。今回の読み取り実行では、現行 `stage_key()` は **OFF／ON で同一、フィールド名変更でも同一**だった。[stage_manifest.py:266](/home/sano/work/forge-sern-design/solver_density_cuda/tools/stage_manifest.py:266)。

中間値の意味は、**速度外挿量を縮める数値的な重み**です。物理モデルの係数ではありません。単純な再構成では実効リミッタが `wψ` になりますが、密度・圧力・組成との熱力学的整合、温度の正値性、冷却量の単調な改善は保証しません。また、滑らかな領域で `w<1` を固定すれば、線形場の厳密再構成も失います。初版の生産フィールドを 0/1 に限定し、中間値を使う生成則の受入れは別途必要です。周期は初版の拒否を維持してください。前処理へ移しただけで、周期 DOF 間の整合が保証されたことにはなりません。[periodicNode_d.cuh:14](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/periodicNode_d.cuh:14)。

§6.1 の物理的な判別条件と §6.2 の受入閾値は維持し、その前提試験として以下を追加します。

- 同じ格子の再生成・HDF5 再保存で、意味内容が同じなら同じ署名になること。節点番号・接続・対象境界の変更、`w` の改変、来歴不一致は検出すること。
- 同一格子 restart と cross-mesh 移植で、**宛先の `w` が保持され、元の `w` が移植されないこと**。prepare では品質確認用の cell 格子でなく、最終 node 変換後に生成すること。[runner_sern3d.py:166](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern3d.py:166)。
- OFF と `w≡1` の共通状態の面値・面流束がビット一致すること。`w=0` は指定側の速度が節点値とビット一致し、相手側・他変数・勾配配列を変更しないこと。中間値は規定式との比較を行う。
- 全段で実効フィールドの識別が追跡でき、変更時は区間が分離し、旧評価の無条件持ち越しを拒否すること。

結論: **前処理への移管を採用し、次は格子と `w` の結び付けを仕様化して、同数節点の別格子を誤って受理しない読込 A/B を先に通す。**

第 1 仮説: 現案で最も危険なのは、再構成式よりも、古い `w` を別の節点配置に適用しても正常入力として通ることである。確度: **高（設計上の穴として。実 run での発生は未確認）**  
  根拠: plan の検査項目には格子との再照合がない。[plan:42](/home/sano/work/forge-sern-design/plans/active/convection-zero-thickness-edge-reconstruction.md:42)。メモリ内の4節点例では、節点を入れ替えて古い `w` を保持すると、長さ・範囲・有限性・`w` ハッシュはすべて通り、正しい対応との差は2節点になった。  
  反証条件: 実際の読込器が格子の意味内容を再計算し、同数節点・同じ `w` ハッシュでも対応の異なる入力を、時間更新前に拒否すること。

第 2 仮説: 重みを元の勾配項にだけ掛ける実装では、後段の速度補正によって `w=0` の節点値保証が破れる。確度: **中**。補正の順序はコードで確認済みだが、新実装は未確認。

判別 A/B: **変更点は入力格子の節点番号付けだけ**。小さな端付き node 格子について、A は正しい格子＋生成済み `w`、B は座標・接続を整合して再番号付けした同じ格子＋更新せず残した同一の `/AUX`。読込検証だけ、時間更新 **0 回**。→ **結果A：正対応だけ受理し、取り違えを拒否**なら第1仮説の穴は塞がった／**結果B：取り違えも受理**なら整合契約は不十分で、物理 A/B へ進めない。両方拒否なら試験不成立。

やらない方がよいこと: `w` 単独のハッシュを格子整合の証明にすること、全1フィールドをソルバで一律拒否すること、中間重みやリング数を結果に合わせて調整すること。

呼び出し側の前提への異議: ブリーフの「変更の内容」は設計決定であり、実装済みの観測事実とは分ける必要がある。確認した HEAD `3e38c41b` には前処理ツールと読込実装がなく、ParaView 表示・周期整合・ビット同一は未検証である。

不足情報: 新実装の差分・バイナリ識別、`case/46.sern_design/run_1055_r7b_m10_A_c/` の入力と判定成果物。手元に対象 run がないため、冷却の解消・収束・定常性は認定していない。**ファイル変更なし、plan 未反映。**
