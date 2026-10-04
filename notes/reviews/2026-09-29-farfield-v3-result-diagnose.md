# codex 諮問 (diagnose): farfield-v3-result

- **brief**: [`notes/reviews/briefs/2026-09-29-farfield-v3-result.md`](../../notes/reviews/briefs/2026-09-29-farfield-v3-result.md)
- **plan**: [`plans/active/boundary-node-farfield-characteristic.md`](../../plans/active/boundary-node-farfield-characteristic.md)
- **date**: 2026-09-29
- **commit**: `c413bc1b` (feature/sern-design)
- **codex**: effort `high`, 4.6 min, rc=0
- **結論**: **V3の限定付き幅判定を保持し、生産切替の前に、2.50 H・farfieldを固定して初期場の履歴だけを変えるA/Bを一組行う。**
- **extra**: `case/46.sern_design/V3_EVAL.txt`, `case/46.sern_design/v3_farfield_setup.py`, `case/46.sern_design/v3_farfield_eval.py`, `case/46.sern_design/README.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：

| 対象・重大度 | 採否 | 根拠と対案 |
|---|---|---|
| **I1：2.50 H の幅判定** | **限定付き採用** | [事前規則](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:201)に対し、2.50 H と両広幅の全4量が `D ≤ ε`。[集計](/home/sano/work/forge-sern-design/case/46.sern_design/V3_EVAL.txt:12)の判定は整合する。記述は「**当該形状・m6_on・g3・固定リミッタ基準値・今回の初期場履歴で、試験した2.50／3.42／4.35 Hのうち最小の適合幅は2.50 H**」とする。無限遠との一致、初期場非依存、汎用境界の受入れは含めない。 |
| **Major：I2「反射をほぼ除去」** | **却下** | 測っているのは力係数の幅感度で、反射成分を分離していない。[境界実装](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/farfieldFlux_d.inc.cuh:98)は圧力・法線速度に加え、流入する組成・乱流量等も変える。代案は「**側方BC変更によって狭幅の力係数が広幅側へ移り、試験系列の幅感度が許容内になった。機序は未同定**」。 |
| **Major：V3通過をもって `accepted`** | **却下** | [ユーザ決定 §5.1 #3g](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:154)はV3への進行許可であり、V2未達の免除ではない。V2cの旧FAIL・新配置の判定不能、V2d-2の時間精度未達を保持し、現スコープでは `in_progress` を維持する。独立参照の収束不足は[§6](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:194)どおり境界合否とは別の未完了記録。結果レビューには提出できるが、完了認定の根拠は不足。 |
| **Major：旧格子差Gの転用** | **要再検証** | [旧g3/g4比較](/home/sano/work/forge-sern-design/plans/active/tooling-nozzle-sern-3d.md:1860)の `ΔC_L=0.00138`、`ΔC_M=−0.04049` は `side_far=slip`。farfieldへの変更量が格子によらず同じという証拠はない。生産の精度認定には、最終BC・同じリミッタ基準値で少なくともg3/g4のGを取り直し、R4dの `G+D` を評価する。g1/g3/g4列全体の結論を引き継ぐならg1も対象。 |
| **Major：評価器だけでは必要条件を保証しない** | **要再検証** | [`--pair`](/home/sano/work/forge-sern-design/case/46.sern_design/v3_farfield_eval.py:28)は前窓差を判定に含めない。読み取り専用の合成入力で、前窓差 `0.001 > 0.00005` でも全量D判定が真になることを確認した。今回の掲載値では窓条件を満たすため、これだけで結果は覆らない。対案は窓条件・標本被覆・有限性を必要条件に組み込み、各runの正式VERDICTと結び付けること。 |
| **Major：SERN帳簿のfarfield対応が未完了** | **修正を採用** | [`sern_momentum.py:58`](/home/sano/work/forge-sern-design/design/forge_design/metrics/sern_momentum.py:58)は境界節点の `ρ,u,P` から全境界の流束を組み、farfieldのHLLC流束を読まない。[plan §5](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:127)の実装事項が残る。診断面流束と圧力基準補正を使う経路を完成させる。これは今回の壁力比較を直ちに無効にする指摘ではないが、現帳簿でfarfieldの保存収支を認定できない。 |
| **Minor：C_M差の「εの45%」** | **訂正を採用** | [掲載値](/home/sano/work/forge-sern-design/case/46.sern_design/V3_EVAL.txt:5)から、平均差は `0.0016392–0.0016854`＝εの **32.8–33.7%**。振幅込みのDが `0.0022352／0.0021944`＝ **44.7／43.9%**。平均差とDを分けて記録する。 |

結論: **V3の限定付き幅判定を保持し、生産切替の前に、2.50 H・farfieldを固定して初期場の履歴だけを変えるA/Bを一組行う。**

第1仮説: 2.50 Hに残る `C_M` の平均差は、初期場履歴よりも境界位置・境界近傍の離散化への依存が主因である。 **確度: 中**

- 根拠: `case/46.sern_design/run_0996_ff_v3a_2p50_ff_cont20k/` と3.42／4.35 Hの平均差は約 `1.6–1.7e-3`。広幅同士の平均差は `4.62e-5`。狭幅を延長した後も差が残る。ただし、これは[集計原本](/home/sano/work/forge-sern-design/case/46.sern_design/V3_EVAL.txt:1)に基づき、初期場依存を除外した証拠ではない。
- 反証条件: 同一の2.50 H・farfield設定でも、広幅由来の初期場から始めた系列が、変動幅を超えて異なる平均値を維持する。

第2仮説: 初期場履歴によって異なる準定常状態に留まっている。**確度: 低、SERNでは未確認**。V2cの[共通初期場試験](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:148)は候補に残す理由になるが、SERNの多重定常解を証明しない。

判別A/B:

- **変えるのは初期保存量だけ。** 両者とも2.50 H・g3・farfield、同一バイナリ・BC・リミッタ基準値・CFL。
- Aは `run_0996` 最終場から同一格子restart。Bは `case/46.sern_design/run_0994_ff_v3b_3p42_ff/` 最終場を2.50 Hの共通領域へ制限して開始する。座標一致の双子節点を境界タグで区別し、保存量を直接コピーする。**既存の拡大用スクリプトはSRC全節点の対応を要求するため、そのまま縮小には使わない。**
- 新規run各20000 step、500 step出力。末尾10000 stepで全4量を比較し、前窓条件未達なら各20000 step延長。正式な準定常VERDICT、全残差の収束VERDICT、床・非有限・farfield置換回数を併記する。
- 追加診断の閾値を事前に `τ=0.2ε` とする。  
  **全量で `D_IC≤τ` なら**、今回の二履歴間の影響はτ以下として、第2仮説をその精度内で棄却する。  
  **いずれかで `|Δ平均|−a_A−a_B>τ` なら**、初期場依存を採用し、第1仮説だけで説明することを棄却する。  
  その中間、または準定常条件未達なら判定不能。V3本来の許容値は変更しない。

やらない方がよいこと: 反射ゼロ・無限遠一致・残差収束済みと表現すること、V2未達をV3で相殺すること、旧slipのGと新farfieldのDを足して生産精度を認定すること。

呼び出し側の前提への異議: 「V2cは試験形状側の問題」「V2d-2はfloat床」とする原因説明は未確定。また、2.50 Hの新旧slip比較で影響が小さくても、旧R4dの**全幅**における自動リミッタ基準値の交絡は除外できない。旧slipの必要幅3.42 Hは当時の設定系列での結果として残す。

不足情報: 手元には対象runディレクトリがなく、生成YAML、`force_history.csv`、全残差、場、`CONVERGENCE_VERDICT.txt`、準定常・品質判定原本を確認できなかった。したがって、I1は**掲載集計に対する規則適用として採用**し、独立した実測再認定は保留する。resultレビュー前にこれらと実行バイナリの識別情報を揃える必要がある。

ファイル変更・forge起動は行っていない。**plan未反映**。呼び出し側でfarfield plan §5.1／§6.1とSERN plan §5.1 R4dへ採否を記録する。
