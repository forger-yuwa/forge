# codex レビュー: boundary-node-farfield-characteristic (plan)

- **plan**: [`../forge-sern-design/plans/active/boundary-node-farfield-characteristic.md`](../../../forge-sern-design/plans/active/boundary-node-farfield-characteristic.md)
- **stage**: `plan`
- **date**: 2026-09-28
- **commit**: `dffa3516` (feature/gap-heating-precision)
- **codex**: effort `high`, 5.2 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M3/m1
- **extra**: `../forge-sern-design/methods/boundary.md`, `../forge-sern-design/solver_density_cuda/tools/farfield_proto1d.py`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**

専用 HLLC・同時刻のスカラー輸送・面流束診断という構造は妥当です。
高温内部の音響反射改善も再現しましたが、正値性、超音速流出の独立性、音速分類の連続性に反例があります。
この3点をホスト段階で解消してから、本体実装へ進むべきです。

指定ファイル全文と参照コードを確認しました。以下は**添付された試作をメモリ上で実行した結果**で、forge 本体の検証結果ではありません。添付試作と隣接ワークツリーのファイルが一致することも確認しています。SERN の既存 run の収束・準定常性は再認定していません。

1. **Major — 正常な亜音速状態から負の外側圧力を構成する**

   根拠は [plan:77](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:77) の線形音響閉包と、[試作:143](/home/sano/work/forge-sern-design/solver_density_cuda/tools/farfield_proto1d.py:143) です。

   CPG、γ = 1.4、ρᵢ = ρ∞ = 1、Pᵢ = P∞ = 1/γ、したがって cᵢ = c∞ = 1 とし、法線速度を uᵢ = −0.9、u∞ = +0.9 にすると、

   `P_R = 1/γ − 0.9 = −0.185714286`

   となります。**内外の入力はどちらも正の密度・圧力を持つ亜音速状態**です。

   試作は密度計算の底だけをクランプするため、ρ_R ≈ 2.68×10⁻⁹ と負の P_R が組になり、音速が NaN、流束3成分もすべて NaN になります。HLL 退避は1回発動しますが、同じ不正状態・波速を使うため回復しません。

   [plan:97](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:97) の「非有限流束を0にする」は検出後の処置であり、境界閉包の成立性を保証しません。

   **対案:** 外側状態の正値性を保つ構成と、その適用範囲外での処理を先に定義してください。構成直後の ρ・P・T・EOS 妥当性を検査し、退避するなら**正の元状態から作る検証済みの保存的流束**を使います。V0h に上記反例と内外速度差の掃引を追加し、密度だけのクランプや流束ゼロ化で通さないことが必要です。

2. **Major — 内部が超音速流出でも、外気条件が流束を変える**

   [plan:82](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:82) は超音速流出の特別扱いを不要としています。しかし、**uᵢ−cᵢ > 0 でも、構成した外側状態の u_R−c_R が負なら Davis の S_L は負**です。

   上と同じ密度・圧力で内部速度を uᵢ = 1.1 に固定し、自由流速度だけを変更すると、添付試作は次を返しました。面積は1です。

   | u∞ | 質量流束 | 法線運動量流束 | エネルギー流束 |
   |---:|---:|---:|---:|
   | 1.1 | 1.100000 | 1.924286 | 3.415500 |
   | 0 | 0.889331 | 2.036991 | 2.981150 |

   後者では S_L = −0.534988。内部の物理流束に対する誤差は順に **−19.15%、+5.86%、−12.72%**です。float32・float64 の両方で再現し、HLL 退避は0回でした。

   [現在の T2:247](/home/sano/work/forge-sern-design/solver_density_cuda/tools/farfield_proto1d.py:247) は uᵢ = 2、u∞ = 3 の一点だけなので、この問題を検出できません。[残作業表:129](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:129) に掲げた「外側状態を変えても内部の物理流束に一致」という条件は満たしていません。

   **対案:** 内部特性がすべて外向きとなる極限で F = F(Uᵢ) を満たすよう、音響閉包と波の分類を修正してください。V0h は内部状態を固定して自由流の速度・圧力・密度を独立に掃引する試験へ広げます。局所 Mach による分岐を追加する場合も、音速通過時の連続性まで同時に確認する必要があります。

3. **Major — 自由流の超音速流入判定で有限の流束ジャンプが生じる**

   根拠は [plan:82](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:82) と [試作:146](/home/sano/work/forge-sern-design/solver_density_cuda/tools/farfield_proto1d.py:146) の `uinf <= -ainf` による外側状態の全置換です。

   γ = 1.4、ρᵢ = ρ∞ = 1、P∞ = 1/γ、Pᵢ = 1.2P∞、uᵢ = −0.3 として、u∞ を −1 の両側へ動かしました。

   | u∞ | 質量流束 | 法線運動量流束 | エネルギー流束 |
   |---:|---:|---:|---:|
   | −1−10⁻⁸ | −0.785513 | 1.884898 | −2.594795 |
   | −1+10⁻⁸ | −0.812740 | 1.670272 | −2.492519 |

   刻みを 10⁻⁴ → 10⁻⁶ → 10⁻⁸ と縮めても、運動量流束差は **0.215001 → 0.214631 → 0.214627** と残ります。丸めや有限差分の傾きではなく、分岐による不連続です。

   面ごとに分類を固定しても、自由流条件や面法線に対する不連続は残ります。[現在の T4:252](/home/sano/work/forge-sern-design/solver_density_cuda/tools/farfield_proto1d.py:252) は u∞ = 0 固定なので、この分岐を通りません。

   **対案:** 自由流 Mach の閾値だけで外側状態を全置換する方式を見直し、音速通過時にも流束が連続となる接続を定義してください。V0h/V0u に **Qₙ/a∞ = −1±ε** の掃引と面法線の連続回転を追加します。内部速度だけの掃引では不十分です。

4. **Minor — 採用方式とレビュー段階の記録が同期していない**

   [plan §5.1 #1g:129](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:129) は「自由流ゴースト＋HLLC が全項目合格→採用」のままです。[plans/README.md:28](/home/sano/work/forge-sern-design/plans/README.md:28) も固定自由流方式・plan-7 待ちを記載しています。現在の §4.2 は特性で構成する外側状態・plan-8 対象です。

   **対案:** 旧方式の結果は履歴として残し、現在の採用候補、通過した試験条件、未通過のゲートを明確に分けて同期してください。

推奨は一つです。**専用 HLLC・同時刻スカラー輸送の構造を維持し、正値性 → 超音速流出の独立性 → 音速分類の連続性の順に、CPG のホストゲートを修正・通過させてから TP、CUDA へ進めてください。**

高温内部の反射率は Δx = 5 / 2.5 mm で **0.14594 / 0.08189%**、退避0回と再現でき、前回指摘への改善は確認できました。目的は既存の静圧出口・ghostless 化と重複せず、node 限定、保存収支の検算、dual-time の刻み・反復数感度、V1–V2 後の V3、広幅対照による配置採否も妥当です。ただし、その改善を境界全体の成立性へ一般化するには上記ゲートが不足しています。

ファイルは変更していません。**plan 未反映**です。

指摘数: Critical 0 / Major 3 / Minor 1
