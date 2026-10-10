# codex レビュー: architecture-float-state-double-geometry (plan)

- **plan**: [`plans/active/architecture-float-state-double-geometry.md`](../../plans/active/architecture-float-state-double-geometry.md)
- **stage**: `plan`
- **date**: 2026-10-10
- **commit**: `2ac7f5a4` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `medium`, 1.8 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M3/m1
- **focus**: §4 (段 ①〜④・使う側ごとの係数・境界面ごとの壁の代表点・接続の生成) と §6 (V0〜V6 の事前登録) を、実装着手前の穴に絞って短く点検してほしい。特に V0 (状態を固定した演算の A/B) を forge でどう実現するか (状態・勾配・物性を固定して 3 つの腕で面ごとの粘性・拡散の流束を出す方法) が具体的に書けているか、段 ① で旧の係数と並べて比べる仕組みに抜けがないか
- **extra**: `notes/investigations/2026-10-10-geometry-precision-inventory.md`, `notes/reviews/2026-10-10-float-geometry-design-diagnose.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
目的と段①〜④の分割は妥当ですが、**段①だけで V0 を実行する具体的な仕組みが未設計**です。  
以下の順に修正し、V0 合格後に利用側の置換へ進むことを推奨します。

確認対象は `2ac7f5a4`。`run_0383*`・`run_0387_dqulp_B0`・`run_0183*` の実体は見つからず、速度・残差倍率は独立に実測確認できていません。ファイル変更・計算実行は行っていません。

1. **Major — V0 の三腕を動かす診断経路が段①に含まれていない。**

   **根拠:** [plan:88](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:88) は入力・配列追加のみ。一方、粘性は [viscousFlux_d.cu:380](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/viscousFlux_d.cu:380)、スカラー拡散は [scalarTransport_d.cu:129](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/scalarTransport_d.cu:129) で面流束を残差へ直接加算しています。既存の通常実行を三回回すだけでは、固定入力による面単位比較になりません。さらに `F1` は拡散直前に更新されます（[main.cpp:1967](/home/sano/work/forge-integ-1005/solver_density_cuda/main.cpp:1967)）。

   **対案:** 段①に診断専用 evaluator を明記してください。通常経路で一度準備した原始量・勾配・物性に加え、`F1`、`axisym_divU`、壁モデル入力、`fx`、半径重み適用後の `S/ss`、接続を保存し、三腕へ同じ入力を渡す。各実カーネルの流束式を共有し、**atomicAdd 前の応力・熱・k/ω 拡散を面 ID 付きで出力**、残差は別バッファへ集積する。参照腕も `S/ss` は `double(S32/ss32)` とし、元の double 面ベクトルへの交換を混ぜない。BC・EOS・勾配・乱流モデル・commit は再実行しない。診断旧腕と未変更カーネルの一致を先に確認する必要があります。

2. **Major — 「旧係数との並列比較」の精度契約と、段①／②の接続検査が曖昧。**

   **根拠:** [plan:60](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:60) は「double 差を float にしてから利用側で係数計算」、[plan:116](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:116) は「double 座標差から作った float の係数」と読め、丸め位置が確定していません。また block-DPLUR の旧経路は `ST` へ変換してから差を取ります（[timeIntegration_d.cu:952](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:952)）。単一の「旧 float 係数」では再現できません。接続変更は段②なのに、段①の V1 に接続一致も要求しています。

   **対案:** §4.2 の方式に統一し、利用側別に「差を取る型・係数演算型・ガード・対象面・向き」を表にする。段①では生産配列を変更せず、新旧の係数・分岐結果を別配列へ記録する。接続は段①で比較用の候補を別生成するか、合格ゲートを段②へ移す、と明記してください。壁の `(irep,y)` も候補なし・同率候補を含めて比較対象に加えるべきです。

3. **Major — closure の不整合を「記録」するだけでは、既存の自由流保持を壊せる。**

   **根拠:** 現行は半径重み適用後の実際の `sx/sy` を集積しています（[variables.cpp:637](/home/sano/work/forge-integ-1005/solver_density_cuda/variables.cpp:637)）。これは既存の [axisymmetric-freestream-hoop-gauge.md:32](/home/sano/work/forge-integ-1005/plans/active/axisymmetric-freestream-hoop-gauge.md:32) の設計そのものです。新 plan は生の double 幾何から closure を作り、差を記録するだけで、不合格条件がありません（[plan:126](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:126)）。

   **対案:** ソース用 closure は、**device に渡す最終面ベクトルを double で集積して最後に丸める**方式を基本にする。生の double 幾何からの closure は診断用として分離する。`P≠pRef` の一様静止場について、正規化残差と偽加速度の許容値を事前登録し、非劣化を段②の必須ゲートにしてください。

4. **Minor — V0 と V4 の判定に、まだ実装者が選べる余地がある。**

   **根拠:** [plan:118](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:118) の「誤差の比」はノルム・重み・壁層集合・ノイズ判別倍率が未定義。V4 は `drift` 指定だけですが、判定ツールには独立した `tail` と `osc` があり、既定はそれぞれ 0.4、0.10 です（[check_quasisteady.py:542](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/check_quasisteady.py:542)）。末尾二万 step の CSV を渡しても、既定ではその全窓を判定しません。

   **対案:** V0 の面流束／残差それぞれのノルム、固定した層マスク、ゼロ参照の扱い、判別不能時の停止条件を先に固定する。V4 は切り出した窓に `--tail 1` を使い、`--drift` と `--osc` を量別に明示する。V6 の基準 run・時刻／窓も、新版結果を見る前に固定してください。

**推奨は一つです。** 上記を §4〜§6 に反映し、まず段①を「三腕診断まで含む変更」として実装すること。double 入力から差を作る方向、利用側ごとの式の維持、境界面ごとの壁代表点、node 主体の検証選択は支持します。既存の median-dual 幾何改善は変換器側であり、本件の solver 内差分とは重複しません。V4 の VERDICT による結論制限も適切です。**本レビューは plan 未反映です。**

指摘数: Critical 0 / Major 3 / Minor 1
