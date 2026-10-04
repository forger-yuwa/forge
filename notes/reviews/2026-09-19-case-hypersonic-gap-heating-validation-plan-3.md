# codex レビュー: case-hypersonic-gap-heating-validation (plan)

- **plan**: [`plans/active/case-hypersonic-gap-heating-validation.md`](../../plans/active/case-hypersonic-gap-heating-validation.md)
- **stage**: `plan`
- **date**: 2026-09-19
- **commit**: `e7e4a0c1` (feature/sern-design)
- **codex**: effort `high`, 5.7 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M5/m1
- **focus**: 2 巡目 (NO-GO, C0/M8/m0) を全件採用した改訂稿。採用が実質的かを検証してほしい: §4.2b (T2-G0 成立ゲート), §4.4b (h と h_fp 相関による正規化), §4.5 (壁温接合部の形状定義と共有ノードの Ts 上書き回避), §4.7 (残余変化の推定と量別 IC 許容), §4.8 (G8 の誤差報告への降格と診断出力の別 plan 依存), §4.10 (事前登録と工学帯の根拠の限定), §4.11 (受け渡しの (A)/(B) 二分), G7b。まだ実装に進めない穴があれば指摘し、無ければ GO 判定の条件を明示してほしい。原報 PDF は papers/gap_heating/ にある
- **extra**: `notes/investigations/hypersonic-gap-cavity-heating-survey.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**

共通分母、T2 の二次元近似の限定、case/49 への誤差幅移植の撤回、G7b は実質的な改善です。  
ただし、試験条件の取り違えと収束区間判定器の欠陥が残り、「2巡目を全件解消」とは判断できません。  
ファイルは変更していません。以下は **plan 未反映**です。

目的と T1/T2 の構成は妥当です。`plans/README.md` と `plans/accepted/` を確認した範囲で、この検証を置き換える完了計画はありません。平面 node・SLAU・block-DPLUR も現在仕様に整合します。ソルバ共通部を変更しない本計画には、cell・周期・軸対称の回帰追加を要求しません。検証を node に限定する現行規則とも整合します。[検証規則:40](/home/sano/work/forge/procedures/verification/README.md:40)

1. **Major — `stage_manifest.py` は、今回区別すべき方程式・離散化の変更を見落とす**

   **根拠:** [plan:198](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:198) は `--segment` を収束ゲートにしています。しかし [stage_manifest.py:44](/home/sano/work/forge/solver_density_cuda/tools/stage_manifest.py:44) が探索するのは `turbulenceModel` で、現行設定の `turbulence.model` ではありません。`lowMachPrecond`、`prandtlLam` も判別対象にありません。

   実際の `stage_key()` と `segments()` をメモリ上で実行した結果は次のとおりでした。

   | 設定差分 | 同一キー判定 | 区間数 |
   |---|---:|---:|
   | `turbulence.model: none → sst` | `True` | 1 |
   | `lowMachPrecond: 0 → 2` | `True` | 1 |
   | `prandtlLam: 0.72 → 0.75` | `True` | 1 |

   したがって、manifest を生成するだけでは「同一方程式・BC・離散化の区間」という前提を保証できません。別設定の残差ピークが混入し、低下桁数を過大評価する経路が残ります。

   **対案:** 判定ツール側の計画に修正を先行依存として登録する。実効設定を YAML として解釈し、乱流モデル・輸送物性・RHS を変える前処理を区別させる。今回の三つの差分が別区間になり、CFL・反復数だけの変更は連結できることを確認してから使用してください。

2. **Major — T2 の「実測温度比」を一系列から全 Reynolds 数へ誤って一般化している**

   **根拠:** [plan:160](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:160) と残作業 #10 は、`1.00/1.18/1.34/1.52` を採用するとしています。しかし原報 Fig.8(a) の掃引角 0° の条件は以下です。[TH76 原報、p.25](https://ntrs.nasa.gov/api/citations/19760019344/downloads/19760019344.pdf#page=27)

   | \(Re'_\infty\,[\mathrm{m}^{-1}]\) | \(T_{surf}/T_{gap}\) |
   |---|---|
   | \(1.47\times10^6\) | 1.00 / 1.19 / 1.37 / **1.61** |
   | \(3.32\times10^6\) | 1.00 / 1.18 / 1.34 / **1.52** |
   | \(7.82\times10^6\) | 1.00 / 1.17 / 1.38 / **1.59** |

   中間 Reynolds 数の条件を他系列に適用すると、G12 の壁温感度と G13 の Reynolds 数感度が交絡します。前回レビューが Fig.5(b) の値を強調したことも、全系列への適用を正当化しません。

   **対案:** 条件表を図・パネル・系列単位で作り、各実測温度比と実測 Reynolds 数に対応した CFD を比較する。**G13 の一次判定は、共通条件が明確な温度比 1.00 の系列に固定する**ことを推奨します。非等温系列の Reynolds 数比較には、温度比補間の誤差を別途付けてください。

3. **Major — 事前登録と T2-G0 の依存順序が逆転し、ゲートもまだ定量化されていない**

   **根拠:** [plan:256](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:256) は CFD を見る前の事前登録を要求します。一方、[残作業表:300](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:300) は次の順です。

   `T0→T1` → `T2-G0` → `T2正規化` → `壁温接合部` → `事前登録`

   T2-G0 が使う分母・熱的条件を、そのゲートの後で確定する構成です。また [G0:332](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:332) の「整合」には、比較位置・分布誤差の尺度・数値許容がありません。

   加えて、TH76 の \(h_{fp}\) は温度比 1.0 の等温基準模型に由来します。加熱履歴を持つ smooth 計算にも一律に \(h/h_{fp}\approx1\) を要求すると、成立ゲート自体が別条件の比較になります。

   **対案:** 各系統を **条件・幾何・分母確定 → 許容の事前登録 → T0/T2-G0 → 本体計算** に並べ直す。G0 には参照模型の熱条件、比較断面、圧力・熱伝達・BL 分布それぞれの尺度と許容を明記する。未測定の \(k,\omega\) は実測再現量と区別し、仮定と感度範囲を固定してください。

4. **Major — 指数緩和の点推定を「残余変化の上限」と扱えず、G9 も旧基準のまま**

   **根拠:** [plan:205](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:205) は指数フィットの残余をゲートにしていますが、単一指数で表せない場合や、遅い成分を識別できない場合の扱いがありません。[G9:330](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:330) には、撤回したはずの「窓2倍」「相対2%かつ絶対 \(0.005q_{fp}\)」が残っています。

   合成系列
   \[
   v(s)=2-e^{-s/10}-e^{-s/10^6}
   \]
   を \(s=50\)〜100 の24点で評価すると、既存 `classify()` は **`STEADY`** を返しました。単一指数の最小二乗フィットでは、

   - 最大フィット誤差：\(7.24\times10^{-6}\)
   - 推定残余：\(4.06\times10^{-5}\)
   - 真の残余：\(0.999945\)

   でした。これは CFD run でも、IC 独立性を含む全ゲートへの反例でもありません。しかし、**よく合う指数フィットだけでは遅い緩和を排除できない**ことは示しています。[判定実装:278](/home/sano/work/forge/solver_density_cuda/tools/check_quasisteady.py:278)

   **対案:** フィット不能・時定数未識別・窓変更で漸近値が不安定な場合は「判定不能」とし、合格させない。物理時間による継続確認と IC 独立性を組み合わせ、識別できた残余だけを誤差予算と比較する。G9 は量別 `tolerances.json` と新基準への参照に置き換えてください。

5. **Major — 離散収支の「別 plan の先行依存」が実行可能な依存関係になっていない**

   **根拠:** [plan:229](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:229) は診断出力を先行依存としていますが、依存先の具体的な plan と解除条件がありません。しかも [§5:282](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:282) と残作業 #7 は離散収支を抽出器の仕事として残し、診断出力の起票は T2/T3 後の #16 です。

   技術的な制約は実在します。[nodeWallDirichlet_d.cu:85](/home/sano/work/forge/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:85) は壁エネルギー残差をゼロ化し、[output.cpp:55](/home/sano/work/forge/solver_density_cuda/output/output.cpp:55) は未登録量を `extraFields` に指定しても出力しません。

   G8 を物理収支の誤差報告に変えた点は正しいですが、それだけでは離散収支の取得経路は完成しません。

   **対案:** 診断出力の別 plan を具体的にリンクし、対象 CV・拘束前残差・境界を横切る数値流束・符号と単位・解除試験を定義する。その完了を離散収支評価の前に置き、§5 の「物理収支」と「離散収支」を明確に分離してください。

6. **Minor — 抽出器が必要とする物性と、計画された出力項目が一致しない**

   **根拠:** [plan:136](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:136) は \(c_p,\lambda\) を場から使う一方、[plan:186](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:186) は `output.level: 1` です。この出力には `thermCond` が含まれず、`cp` は現在の出力登録にもありません。[output.cpp:46](/home/sano/work/forge/solver_density_cuda/output/output.cpp:46)、[variables.hpp:223](/home/sano/work/forge/solver_density_cuda/variables.hpp:223)

   **対案:** `thermCond` を明示出力し、\(c_p\) は保存された \(T,Y\) と同一の物性 DB・基準から再計算する方式に固定する。抽出器の必須フィールド検査も仕様に加えてください。

§4.5 の共有ノードへの矛盾した `Ts` を禁止する方針は妥当です。ただし、実装時には `physID` の名称分割だけで済ませず、**各壁ノードに適用される温度が一意であること**を検査する必要があります。§4.10 の精度限界と統計量の区別、§4.11 の受け渡し (A)/(B)、G7b の量別格子評価は採用済みと評価します。

既存結果については、`case/48.flat_plate_cooled_m4/run_0011_Bplain_tw300_y3/` を再判定し、**本段単独の VERDICT は `NOT CONVERGED (stalled/plateau)`**、`rms_roe` 最終 \(1.42\times10^{-1}\)、低下0.1桁でした。既存 `CONVERGENCE_VERDICT.txt` と整合し、これを無条件の合格基準にしない修正は正しいです。[run 索引](/home/sano/work/forge/case/48.flat_plate_cooled_m4/README.md)

**推奨は、T1/T2 の構成を維持し、「条件台帳と判定ゲートの確定」を先行マイルストーンにすることです。** 優先順は、①区間判定器、②系列別実条件と G0・事前登録、③G9 と残余推定、④離散収支の依存関係、⑤抽出器の出力仕様です。これらを閉じてから、本体計算の実装へ進むことを勧めます。

指摘数: Critical 0 / Major 5 / Minor 1
