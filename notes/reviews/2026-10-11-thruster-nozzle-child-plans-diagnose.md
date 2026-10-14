# codex 諮問 (diagnose): thruster-nozzle-child-plans

- **brief**: [`notes/reviews/briefs/2026-10-11-thruster-nozzle-child-plans.md`](../../notes/reviews/briefs/2026-10-11-thruster-nozzle-child-plans.md)
- **plan**: [`plans/active/tooling-nozzle-bell-two-mode.md`](../../plans/active/tooling-nozzle-bell-two-mode.md)
- **date**: 2026-10-11
- **commit**: `a8eb7c53` (docs/thruster-nozzle-design)
- **codex**: effort `high`, 4.2 min, rc=0
- **結論**: **まず `outlet_statPress` の逆流状態が指定 `Tt` に応答するかを、固定内部状態からの境界カーネル A/B 一つで確認してください。**
- **extra**: `plans/active/tooling-nozzle-ambient-separation.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：**A＝NO-GO、B＝NO-GO。現行の §4・§6 のまま実装・比較へ進める案は受理しません。** ユーザが決めた TP・C2・2 モード・検証ケースは維持します。

以下はコードと設計文書のレビューです。禁止された run データ・PDF は読んでおらず、過去の数値と収束は再認証していません。共通指摘を含め、各 plan 最大 5 件です。

| 対象・重大度 | 採否・根拠 | 推奨する変更 |
|---|---|---|
| **共通 C1 — Critical：推力会計が閉じていない** | **却下。** A の [83 行](/home/sano/work/forge-thruster/plans/active/tooling-nozzle-ambient-separation.md:83)、B の [50 行](/home/sano/work/forge-thruster/plans/active/tooling-nozzle-bell-two-mode.md:50) は「入口運動量＋内壁力」としており、**入口圧力の面積積分が欠落**しています。既存の [sern_forces.py:88](/home/sano/work/forge-thruster/design/forge_design/metrics/sern_forces.py:88) は入口圧力項を含みます。また、出口の逆流は保存則による検算を弱める理由になりません。 | **同一の閉じた検査体積から両会計を定義する。** 入口は運動量・圧力・必要な粘性応力を含め、出口も逆流込みの符号付き積分で照合する。有限厚リップの底面力を含む総推力と、ノズル内部だけの流れ推力を区別する。付着流だけでなく剥離条件にも収支の受理基準を置く。 |
| **A1 — Major：外気 BC の説明が現行コードと違う** | **却下。** [A:64](/home/sano/work/forge-thruster/plans/active/tooling-nozzle-ambient-separation.md:64) の構成は双方向の静止大気境界ではありません。`inlet_Pressure` は流入を法線方向に限定し、流出時は速度をゼロにします（[boundaryCond_d.cu:1132](/home/sano/work/forge-thruster/solver_density_cuda/cuda_forge/boundaryCond_d.cu:1132)、[1155 行](/home/sano/work/forge-thruster/solver_density_cuda/cuda_forge/boundaryCond_d.cu:1155)）。`outlet_statPress` の逆流は内部エントロピー・組成を使い、指定 `Pt/Tt` から大気状態を作りません（[615 行](/home/sano/work/forge-thruster/solver_density_cuda/cuda_forge/boundaryCond_d.cu:615)、[655 行](/home/sano/work/forge-thruster/solver_density_cuda/cuda_forge/boundaryCond_d.cu:655)）。 | **軸対称の静止大気用に、流入時の温度・組成と流出時の内部状態を適切に扱う双方向開放 BC を別 plan で検証する。** `farfield` の拒否フラグだけを外さない。実物の有限厚リップを断熱壁にする案は採用するが、リップ外側の未定義の上流境界まで壁で塞がない。微小並行流への変更を単なる安定化として扱わない。 |
| **A2 — Major：主定義の置換を V1 では保証できない** | **却下。** [A:75](/home/sano/work/forge-thruster/plans/active/tooling-nozzle-ambient-separation.md:75) は MOC 参照曲線を forge の粘性全流れ解に置換しますが、[V1:116](/home/sano/work/forge-thruster/plans/active/tooling-nozzle-ambient-separation.md:116) が調べるのは X/R*＝1〜5 の壁圧です。交点付近での参照曲線の差や、最大勾配の抽出誤差を拘束しません。壁圧 3 % は位置誤差 2 % を保証しません。 | **実験と同じ MOC 全流れ壁圧を主参照に固定する。** forge 全流れ参照は副指標にする。ポート位置への再標本化、微分・平滑化、交点探索区間を事前登録し、抽出誤差を位置の単位で評価する。V1 の 3 % は上流場の整合確認に限定する。 |
| **A3 — Major：文献帯から誤差の原因は判定できない** | **却下。** [A:119](/home/sano/work/forge-thruster/plans/active/tooling-nozzle-ambient-separation.md:119) の帯内判定では、異なる実装・格子の誤差相殺を排除できません。さらに、採用予定レシピは `dilatationCorrection: 2, katoLaunder: 1`（[recommended-settings.md:136](/home/sano/work/forge-thruster/procedures/recommended-settings.md:136)）で、文献の「標準 SST」との方程式の同一性も未確認です。 | **TC1A の基準検証では、帯内でも同一格子・同一 BC・モデル項を照合した SU2 SST を必須にする。** 全設計候補で SU2 を回す必要はありません。文献帯は参考分布とし、帯内でも「モデルの偏り」「forge 固有誤差なし」と断定しない。 |
| **A4 — Major：U と非収束解の受理条件が不十分** | **却下。** [A:117](/home/sano/work/forge-thruster/plans/active/tooling-nozzle-ambient-separation.md:117) の 2 格子の差は誤差上限ではなく、複数の誤差の最大値も合成上限ではありません。[A:87](/home/sano/work/forge-thruster/plans/active/tooling-nozzle-ambient-separation.md:87) は定常反復の振動を物理的な周期現象として受理する余地がありますが、定常は局所擬似時間です（[solver-settings.md:37](/home/sano/work/forge-thruster/procedures/solver-settings.md:37)）。 | **3 水準の系統格子で収束傾向を確認し、格子・領域・反復・抽出の幅を保守的に合成する。** 0.132 R* は目標として保持できますが、実験不確かさではありません。定常計算の `OSCILLATING` を実験の平均±振幅へ対応させず、本 plan では判定不能とする。 |
| **B1 — Major：C2 の根探しが未定義で、既存関数も直接は使えない** | **要再検証。** [B:74](/home/sano/work/forge-thruster/plans/active/tooling-nozzle-bell-two-mode.md:74) は非負の「当てはめ残差」を小さくするだけで、根を取る符号付き関数と根の存在条件がありません。`joint_fit_wall` の始点拘束は任意の B の状態ではなく、`r′=x0/R, r″=1/R` に固定されています（[wall_axismach.py:1061](/home/sano/work/forge-thruster/design/forge_design/geometry/wall_axismach.py:1061)）。 | **親の根探し方針を維持し、g(κ_B)＝κ_Rao(B⁺;κ_B)−κ_B を明示する。** 各評価で B・D を解き直し、符号を挟めない候補は成立未確認とする。端点の r・r′・r″ を独立指定できる当てはめへ拡張する。κ＝r″/(1+r′²)³ᐟ² と r″ を混同しない。B を跨ぐ当てはめへの無条件置換は推奨しません。 |
| **B2 — Major：CPG 極限は TP 最適性の検証にならない** | **要再検証。** [B:73](/home/sano/work/forge-thruster/plans/active/tooling-nozzle-bell-two-mode.md:73) の「γ が陽に無ければ TP に使える」は未証明です。CPG 極限・保存則・50 形状との比較（[B:114](/home/sano/work/forge-thruster/plans/active/tooling-nozzle-bell-two-mode.md:114)）は、非定比熱条件での停留性を保証しません。`rao_planar.py` 自体も自前導出扱いです（[12 行](/home/sano/work/forge-thruster/design/forge_design/geometry/rao_planar.py:12)）。 | **凍結 TP の h(T)・s(T,p) を使って拘束付き変分を導出し、実際に cp が変化する条件で、拘束を保った壁摂動に対する推力の一次変分を検査する。** CPG 照合は回帰試験として残す。任意 50 点に勝ったことを最大推力の証明にしない。 |
| **B3 — Major：δ* 補正後の拘束回復が成立しない** | **却下。** [B:75](/home/sano/work/forge-thruster/plans/active/tooling-nozzle-bell-two-mode.md:75) の相似倍率では **ε＝(r_e/r_t)² は変化しません**。補正で変わった ε と L を相似倍率だけで同時に戻せません。また積分法は μ 以外にも温度回復の近似を持ちます（[deltastar_integral.py:24](/home/sano/work/forge-thruster/design/forge_design/feedback/deltastar_integral.py:24)）。 | **δ* 補正後の物理壁を、ε・L・スロート曲率を保持する拘束付き当てはめで確定する。** 両モードに共通の寸法決定を適用し、δ* の熱力学的近似と較正誤差も明記する。M の 1 変数と D の 4 変数という差自体は許容しますが、結論は「完成した設計手順の比較」とし、較正・探索を含む総費用を比較する。 |
| **B4 — Major：V2 が Euler にならない危険と、V3 の識別能不足** | **却下。** [B:118](/home/sano/work/forge-thruster/plans/active/tooling-nozzle-bell-two-mode.md:118) の `visc: 0`・`thermCond: 0`・slip だけでは不十分です。現行 runner は `viscMethod: 1`、SST が既定（[runner.py:90](/home/sano/work/forge-thruster/design/forge_design/evaluate/runner.py:90)、[103 行](/home/sano/work/forge-thruster/design/forge_design/evaluate/runner.py:103)）。さらに [B:120](/home/sano/work/forge-thruster/plans/active/tooling-nozzle-bell-two-mode.md:120) の各モード 1 形状の格子差＋停止閾値は、前線全体の誤差幅にはなりません。 | **Euler を `viscMethod: 0`・粘性/熱伝導ゼロ・`turbulence.model: none`・slip で明示する。** V2 の 0.2 % は粗い整合確認に限定する。V3 は共通の L で再評価し、差 Δη の格子・反復・補間・探索の不確かさを評価する。1.5e−3 を残すなら「工学的同等幅」とし、「差なし」は「その幅内で同等」へ変更する。 |

着手順の **評価器 TP 化 → D → M は採用**します。ただし、共通の推力定義・物理寸法・輸送物性を先に確定してください。M の導出成立前に比較キャンペーンを投入する根拠はありません。

結論: **まず `outlet_statPress` の逆流状態が指定 `Tt` に応答するかを、固定内部状態からの境界カーネル A/B 一つで確認してください。**

第 1 仮説: **現案 A は既存 BC の契約を誤認しており、下流逆流の大気温度を指定できていない。** 確度: **高**
  
根拠: [boundaryCond_d.cu:615](/home/sano/work/forge-thruster/solver_density_cuda/cuda_forge/boundaryCond_d.cu:615) は内部エントロピーから密度を作り、[660 行](/home/sano/work/forge-thruster/solver_density_cuda/cuda_forge/boundaryCond_d.cu:660) はその密度と静圧から温度を決めます。剥離位置への影響量は未確認です。
  
反証条件: 読んだソースに対応するビルドで、同じ内部状態・同じ `Ps` の逆流境界が、指定 `Tt` の変更に応じて大気由来のエンタルピーを変えること。

第 2・第 3 仮説: **追加しません。** SST の上流偏りや M/D の差の大きさは、現段階では検証済みの仮説に順位付けできません。

判別 A/B: **変える設定は下流 `Tt` だけ**。同じ非ゼロ逆流状態から、A＝T_a、B＝2T_a として境界カーネルを各 1 回評価する。見る量は ghost の T・ρ・全エンタルピー。**両者が不変なら「逆流 Tt が外気を規定する」を棄却／応答するなら今回読んだ経路だけで実行時挙動を説明できるという仮説を棄却**し、ビルド・呼出経路を確認する。時間発展は不要で、収束判定の対象ではありません。この診断結果を剥離位置への影響の証明には使いません。

やらない方がよいこと: **現行 BC のまま TC1A を回し、文献 SST の帯に入ったことを理由に外部領域・モデル・ソルバをまとめて受理すること。**

呼び出し側の前提への異議: 「現行評価器は粘性係数 1.8e−5 の定数」は誤りです。`viscMethod: 1` は空気の Sutherland です（[solver-settings.md:238](/home/sano/work/forge-thruster/procedures/solver-settings.md:238)）。また、case/27 の「本質的に非定常」「内部は定常」は今回確認できた範囲では README の解釈であり、VERDICT による根拠ではありません。Phase 2 の台地から今回の TP・ε≈50 の差を 0.1〜0.2 % と予測する根拠も不足しています。

不足情報: TC1A の参照 MOC 壁圧と抽出手順、実験不確かさ、リップ外形・外気組成・温度・乱流入口条件、比較する SST の具体的なモデル項、B の要求流量と CEA 全温・輸送物性、TP Rao の導出。マッハディスクの文献帯は原典未確認として扱います。

**ファイル変更・forge 起動なし。plan 未反映。反映先は A §4.2〜4.5・§6、B §4.1・4.4〜4.5・§6 です。**
