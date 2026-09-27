# codex レビュー: boundary-node-farfield-characteristic (plan)

- **plan**: [`../forge-sern-design/plans/active/boundary-node-farfield-characteristic.md`](../../../forge-sern-design/plans/active/boundary-node-farfield-characteristic.md)
- **stage**: `plan`
- **date**: 2026-09-27
- **commit**: `e4ea8df8` (feature/gap-heating-precision)
- **codex**: effort `high`, 5.4 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M5/m2
- **extra**: `../forge-sern-design/methods/boundary.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**

専用境界流束とスカラー面値を同時に生成する改訂は妥当です。しかし、特性分岐に有限の不連続があり、TP 近似・保存収支・音響反射の検証にも欠落があります。実装前に設計を修正してください。

対象は `/home/sano/work/forge-sern-design`。ファイルは変更していません。`plans/README.md` と `accepted/` を確認した限り、既存の静圧出口計画とは目的が異なり、重複ではありません。node 限定も現行の検証方針と整合します。

実データでは、`case/46.sern_design/r4d_view/z2p50H/res_20000.h5` の最外面 **51,143 節点、Y0=0–0.128410、T=176.346–596.966 K、a²ρ/P=1.372255–1.404682** を再確認できました。ただし、`run_0986/0988/0989/0990` の残差履歴・判定原本はローカルに見つからず、既報の収束・準定常性や力係数差は独立に再認定していません。

1. **Major — 超音速への切替で境界流束が不連続になる**

   根拠: [plan:70](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:70) の「内部・自由流の両方が超音速なら全量内部、それ以外は亜音速式」。

   この式を float32 で評価しました。γ=1.4、ρᵢ=ρ∞=1、Pᵢ=P∞=1/γ、cᵢ=c∞=1、Uₙ,∞=2、面積=1 とすると、

   | Uₙ,ᵢ | 分岐 | Uₙ,b | P_b | 質量流束 |
   |---:|---|---:|---:|---:|
   | 0.999999 | 亜音速式 | 1.500000 | 0.341641 | 0.885735 |
   | 1.000001 | 全量内部 | 1.000001 | 0.714286 | 1.000001 |

   **入力差約 2×10⁻⁶ に対して圧力が約 2.09 倍**になります。丸め誤差ではなく分岐そのものの不連続です。V0u の「期待した分岐と有限値」では、この欠陥を合格させます。

   また、SU2 は構成状態をそのまま物理流束にせず、内部状態と構成状態を数値流束へ渡しています。したがって SU2 の構成式だけでは今回の `F(U_b)` を裏付けられません。[参照実装:5009](/home/sano/work/forge/.external/su2-src/SU2_CFD/src/solvers/CEulerSolver.cpp:5009)、[SU2 公式ソース](https://raw.githubusercontent.com/su2code/SU2/master/SU2_CFD/src/solvers/CEulerSolver.cpp)。

   **対案:** 内外の分類不一致・遷音速を、両側の波速に基づく一貫した Riemann 流束で扱う設計へ修正する。V0u は分岐確認に加え、音速近傍の連続性、流入・流出の極限、独立した一次元解との比較を必須にする。

2. **Major — V2b は frozen-γ 近似の誤差を測れない**

   根拠: [plan:60](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:60) は近似誤差を V2b で評価するとしていますが、[plan:135](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:135) の合格条件は収支・最終組成・有限性です。

   **誤った境界圧力・音速を返しても、その流束で保存量を更新すれば収支は閉じます。** 最終的に Y が自由流へ置換されることも、過渡の音響応答や TP 特性の正しさを証明しません。実測した γ の範囲だけから誤差を小さいとも判断できません。

   **対案:** 温度・組成差を持つ一次元問題について、同じ実物性を使う独立した参照解を設け、境界圧力・速度・波の振幅を比較する。許容値は実行前に固定し、保存性試験とは別の合否にする。CPG、単成分 TP、多成分 TP の対応範囲を試験表に明記する。

3. **Major — V2b の保存収支に SST ソースとエネルギーの保存対象が欠けている**

   根拠: [plan:136](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:136) は境界流束と体積積分の時間変化で k 収支を判定します。しかしコードは **(Pₖ−Dₖ)V** を加算し、Dₖ=β*ρkω です。[ransSource_d.cu:210](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/ransSource_d.cu:210)、[同:248](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/ransSource_d.cu:248)。

   一様流でも k・ω が正なら消滅項が残ります。さらに `sstEnergyIncludesK` 有効時も保存配列 `roe` は平均流エネルギーのままで、全エネルギーの積分対象は **`roe + roK`** です。[implementation.md:572](/home/sano/work/forge-sern-design/methods/turbulence/implementation.md:572)。

   **対案:** 各保存式について、蓄積・全境界流束・体積ソース・ピンによる補正を含む収支式を先に定義する。物理時間の積分方式に対応する時間重みと、相殺で消えない規格化分母も固定する。`sstEnergyIncludesK` の両設定、k と ω の双方を検証対象にする。

4. **Major — V2a の `slip` 対照は背景流と不整合で、反射率の基準にならない**

   根拠: [plan:132](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:132) はチャネル方向 M=0.3 の一様流に対し、端面 `slip` の反射率を約 1 と想定しています。

   その端面では、背景流の法線速度が非零なのに `slip` は質量流束を零にします。パルスを入れなくても背景流が変化し、振幅 10⁻³P∞ の音響反射と混ざります。また「block-DPLUR」だけでは物理時間積分を指定したことになりません。音響計算には `unsteady: 1` と、陰解法なら dual-time の設定が必要です。[recommended-settings.md:247](/home/sano/work/forge-sern-design/procedures/recommended-settings.md:247)。

   **対案:** M=0.3 の本試験は同じ背景流を維持する長領域参照と比較し、`slip` の反射検出対照は別の M=0 試験にする。物理時間刻み、パルス幅、観測窓、dual-time 内部反復の十分性を事前指定する。V2c も「高さ 2 倍」だけで参照扱いせず、参照側の反射が評価線へ届かない幾何条件を固定する。

5. **Major — `OPEN_KINDS` への追加だけでは、新しい境界の運動量収支を検算できない**

   根拠: [plan:107](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:107) が挙げる `sern_momentum.py` は、境界状態ではなく **owner の `VALUE/ro,U,P`** から流束を再構成します。[sern_momentum.py:58](/home/sano/work/forge-sern-design/design/forge_design/metrics/sern_momentum.py:58)。

   今回は内部状態 Uᵢ と構成状態 U_b が異なります。項目 1 の最初の例では、owner から求める質量流束は約 1、実際の計画流束は約 0.886 です。境界種別の登録ではこの差は直りません。既存ツール自身も厳密な離散収支ではないと明記しています。

   **対案:** `farfield_flux_d` が残差へ投入した面流束を、面 ID・評価時点付きで診断出力し、収支検証はその値を積分する。V2b の収支取得と SERN の境界帳簿を同じ仕組みにそろえ、§5 の実装ステップへ追加する。

6. **Minor — V3b が「境界条件の妥当性」と「2.50 H を採用できるか」を混同している**

   根拠: [plan:144](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:144) は 2.50→3.42 H の無感度を機能の合格条件にしています。

   正しく実装された遠方境界でも、近すぎる位置では有限振幅波・混合層・省略した粘性流束の影響が残り得ます。一方、三つの幅で同じ偏りを持つ可能性もあるため、無感度だけで無限遠との一致は示せません。

   **対案:** V1–V2 を境界機能の受入れ、V3 を配置・領域寸法の採否とする。V3b の結論は「試験した側方幅系列で許容内」に限定し、不合格なら幅を広げる。restart 方法、メッシュ品質、判定区間、run 索引更新も §5.1 に明示する。

7. **Minor — 改訂前の仕様が残り、対応範囲と退避処理の受入れ条件が曖昧**

   根拠: [plan:31](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:31) は依然「`bvar` を埋めればよい」とし、§4.0 と矛盾します。§2 と §4.3 の併用禁止一覧も一致していません。[plan:69](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:69) の異常状態・大膨張時の退避には、合格 run で許す回数がありません。

   **対案:** §2・§3・`methods/boundary.md` を改訂仕様へ統一し、対応外構成の起動拒否試験を追加する。退避は面 ID・理由・回数を記録し、検証の評価区間では発動零を要求する。有限値へ置換されたことを正常性の根拠にしない。

推奨は、**専用流束・同時刻スカラー面値という構造を維持し、計画を修正して再レビューすること**です。優先順は **①遷音速分岐、②TP の独立精度検証、③実流束による保存収支、④音響試験、⑤V3 と仕様の整理**。現行の ghostless A⁺ を近似前処理として残す判断自体は許容できますが、それで境界流束の不連続は解消しません。

本回答はレビュー提案であり、ユーザー指定に従い **plan 未反映**です。

指摘数: Critical 0 / Major 5 / Minor 2
