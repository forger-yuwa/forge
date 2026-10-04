# codex レビュー: boundary-node-farfield-characteristic (plan)

- **plan**: [`../forge-sern-design/plans/active/boundary-node-farfield-characteristic.md`](../../../forge-sern-design/plans/active/boundary-node-farfield-characteristic.md)
- **stage**: `plan`
- **date**: 2026-09-28
- **commit**: `39906fa6` (feature/gap-heating-precision)
- **codex**: effort `high`, 5.7 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M2/m1
- **extra**: `../forge-sern-design/methods/boundary.md`, `../forge-sern-design/solver_density_cuda/tools/farfield_proto1d.py`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定：NO-GO**

目的、node 境界専用流束、V1–V2 と SERN 配置評価 V3 の分離は妥当です。  
ただし、TRRS の正値性に再現可能な反例があり、状態混合の定義にも不整合が残っています。  
CUDA 実装前に、この２点をホスト試作で解消してください。

依頼ファイルに埋め込まれた plan・試作コードと、参照先ファイルが一致することを確認しました。以下の数値はホスト試作の再実行結果であり、forge の計算 run の結果ではありません。

1. **Major — TRRS は無条件に正値ではなく、試作の下限処理は float32 で破綻する**

   根拠：[plan:78](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:78)、[farfield_proto1d.py:168](/home/sano/work/forge-sern-design/solver_density_cuda/tools/farfield_proto1d.py:168)。

   TRRS の分子 B = cᵢ + cₚₒ − (γ−1)(u∞−uᵢ)/2 は、十分大きな離反速度で非正になります。「常に P_R > 0」は成立しません。

   γ = 1.4、ρᵢ = ρ∞ = 1、Pᵢ = P∞ = 1/γ、uᵢ = −0.9、u∞ = 10 では **B = −0.18**。現試作はこれを小さな正数に置換しますが、結果は次のとおりです。

   | 精度 | 構成した ρ_R | 構成した P_R | 流束・退避 |
   |---|---:|---:|---|
   | float64 | 約 1×10⁻⁶⁰ | 約 7.14×10⁻⁸⁵ | 有限、退避カウンタ 0 |
   | float32 | 0 | 0 | 全流束 NaN、退避カウンタ 1 |

   下限置換そのものは計数されていません。そのため、float64 の「退避 0」は、無補正で正値性を保った証拠になりません。plan の構成状態検査を実装すれば NaN を避けられますが、それは退避であり、現在の正値性の主張とは別です。

   **対案：** 真空・近真空の判定条件と処理を §4.2 に明記し、下限置換も必ず計数してください。無効な外側状態を HLLC/HLL に渡さず、正の入力状態を使う退避経路を定義する必要があります。V0h に float32、離反速度、密度比・圧力比の独立掃引を追加し、「有限」と「補正・退避なし」を別々に判定してください。

2. **Major — 混合する変数が未定義で、試作・多成分・SST の契約が一致しない**

   根拠：[plan:83](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:83)、[plan:102](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:102)、[farfield_proto1d.py:175](/home/sano/work/forge-sern-design/solver_density_cuda/tools/farfield_proto1d.py:175)。

   plan は `U_R ← (1−w)U_R + wU_i` と書きますが、試作が混合するのは **ρ・u・P の原始変数**です。U を保存量と解釈すると別の流束になります。

   内部 `(ρ,u,P) = (1, 0.95, 1/γ)`、自由流 `(10, −3, 10/γ)` では、w_q = 1、w_i = 0.5 となり、再計算結果は以下です。

   | 混合方法 | u_R | P_R | HLLC 質量流束 |
   |---|---:|---:|---:|
   | 試作の原始変数混合 | −1.025 | 3.92857 | −5.11729 |
   | 保存量混合 | −2.64091 | 5.34698 | −13.45933 |

   また、混合帯でも実際の質量流束は負になり得ます。このとき完全な状態を混合するなら、組成・k も混合後の値です。例えば原始変数混合で Yᵢ = 0.13、Y∞ = 0 なら Y_R = 0.065 ですが、§4.3 は流入種流束に Y∞ = 0 を使います。**エネルギー・EOS が参照する状態と、輸送する組成が分離します。** k についても同じ問題があります。

   **対案：** 試作に合わせて混合を原始変数で定義し、接線速度・全組成・k・ωまで処理を明文化してください。混合後の同一状態から EOS・音速・エネルギーとスカラー流束を作り、流入時はその `φ_R` を使う契約に統一することを推奨します。V0u の「常に自由流値」という条件も混合帯について更新し、上記の逆流例を TP・SST 両設定で検証してください。

3. **Minor — 現在仕様の一覧表が旧方式のまま**

   根拠：[methods/boundary.md:37](/home/sano/work/forge-sern-design/methods/boundary.md:37)。

   `farfield` の一覧表は依然として「境界節点の ρc で線形化」と説明しています。本文・plan の採用方式である TRRS と一致しません。

   **対案：** 上記２点の仕様確定後、一覧表も TRRS・混合変数・適用範囲に合わせて同期してください。

**推奨は、専用 HLLC 境界という構造を維持し、上記２点を修正したホストゲートを先に通すことです。**

Δx = 5 mm の反射率は、一様 M 0 で **0.501%**、M 0.3 で **0.152%**、高温内部で **0.155%** と、掲載値を再現できました。方式選定の根拠はあります。一方、この成功だけでは float32 と多成分・SST への移植を保証しません。

既存の静圧出口 plan と目的の重複はなく、node 限定・非対応構成の拒否、ghostless 対角近似、保存収支と非定常時間精度の検証方針は妥当です。修正順は **正値性・退避の定義 → 混合状態とスカラーの統一 → float32／TP／SST ホスト検証 → CUDA 実装**としてください。

ファイルは変更していません。指摘・推奨は **plan 未反映**です。

指摘数: Critical 0 / Major 2 / Minor 1
