# codex 諮問 (diagnose): contur-property-temperature

- **brief**: [`notes/reviews/briefs/2026-10-08-contur-property-temperature.md`](../../notes/reviews/briefs/2026-10-08-contur-property-temperature.md)
- **plan**: [`plans/active/tooling-nozzle-isothermal-wall-chain.md`](../../plans/active/tooling-nozzle-isothermal-wall-chain.md)
- **date**: 2026-10-08
- **commit**: `1dd09c45` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 6.0 min, rc=0
- **結論**: **係数と輸送物性を固定し、熱閉包だけを温度形からエンタルピー形へ替える CFD 0 step の A/B を先に行ってください。**
- **extra**: `plans/active/verification-m6-axis-wave-mesh-su2.md`, `design/forge_design/feedback/deltastar_integral.py`, `case/48.flat_plate_cooled_m4/README.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：**現案は修正付き採用。第1層の整合化を先行し、#8d の再較正と第2層の方式選定は保留**とします。Critical はありません。

| 重大度 | 対象・採否 | 根拠と対案 |
|---|---|---|
| **Major** | **2層分離：採用。ただし境界を修正** | 同じ組成・同じ `h(T)`・EOS・輸送物性を使うことは整合性の問題です。一方、**二次の `h(u/u_e)` 分布自体、回復係数 r、`a_crocco`、速度分布 N は閉包モデル**です。エンタルピー表現にすれば閉包まで一意になるわけではありません。`ρ/ρ_e = T_e/T` は、境界層内で圧力と組成が一定という仮定の下で成立します。現行の式は [deltastar_integral.py:49](/home/sano/work/forge-integ-1005/design/forge_design/feedback/deltastar_integral.py:49)。第1層には熱力学・輸送の整合、第2層にはこれらの閉包と適用範囲を置いてください。 |
| **Major** | **洗い出し a〜g：不足あり** | **運動量積分式の圧力勾配項にも定比熱近似が残っています。** [同ファイル:208](/home/sano/work/forge-integ-1005/design/forge_design/feedback/deltastar_integral.py:208) は `d ln u_e/dx` を一定 γ の式で `dM/dx` に変換しています。独立再計算では、case/45 の気体でこの変換係数の「現行／熱力学整合値」は M=3 で **0.96415**、M=6 で **0.98520**。対案は、項を `d ln ρ_e/dx + (2+H) d ln u_e/dx + d ln r_w/dx` で評価することです。等エントロピー・組成固定なら前二項は `(2+H−M²) d ln u_e/dx`。上流の `dMdx` 生成にも同じ近似があります（[同:116](/home/sano/work/forge-integ-1005/design/forge_design/feedback/deltastar_integral.py:116)）。上流では現行 RHS と代数的に相殺する部分があるため、片側だけ置換せず一貫して扱ってください。 |
| **Major** | **「NS の Pr_lam=0.72 だから r=0.72¹ᐟ³」：根拠を却下** | `viscMethod: 2` は μ と λ を混合輸送モデルから直接計算し、この経路で `prandtlLam` を使いません（[gasProperties_d.cu:107](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/gasProperties_d.cu:107)）。乱流熱伝導は別途 `λ + c_p μ_t/Pr_t` です（[viscousFlux_d.cu:259](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/viscousFlux_d.cu:259)）。独立参照でも、250 K の Pr は **0.75159**、1000 K は **0.75774**、1470 K は **0.74878**でした。**r=0.89628 は今回の診断では固定してよい**ですが、「NS と同じ分子 Pr」ではなく、検証対象の回復閉包として記録してください。混合気 Pr の立方根も、可変物性・圧力勾配下で自動的に正解になるわけではありません。 |
| **Major** | **μ を NS と統一：採用。欠落時の無条件 fallback は却下** | 独立参照で μ_mix/μ_air は 250 K **0.95565**、1470 K **1.08933**を再現しました。影響箇所は μ_e・μ_w に加え、N の Reynolds 数、入口 θ₀、`Re_theta_i` の床です。ただし `F_Rδ Re_θc = ρ_e u_e θ_c/μ_w` なので、**変換比の11%差がそのまま摩擦誤差になるわけではありません**（[deltastar_integral.py:163](/home/sano/work/forge-integ-1005/design/forge_design/feedback/deltastar_integral.py:163)）。TP の NS では `gas.transport` が必須という既存契約（[probdef.py:157](/home/sano/work/forge-integ-1005/design/forge_design/probdef.py:157)）に揃え、不足時は停止。CPG／明示された旧方式だけ Sutherland を使うのが妥当です。 |
| **Major** | **「T_aw 誤差は除外済み」：要再検証** | `tw_experiment()` は **T_w だけ**を観測値へ替え、温度分布内の T_aw は旧式のままです（[delta_contur_compare.py:231](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/delta_contur_compare.py:231)）。`taw_sensitivity()` も `T_e+r(T_t−T_e)` への置換で、エンタルピー分布への変更ではありません（[同:300](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/delta_contur_compare.py:300)）。したがって「壁温の単独置換では傾きが消えない」は観測ですが、**熱閉包全体を原因候補から外すことはできません**。 |
| **Major** | **θ の比較データ：修正が必要** | `integral_bl` の `theta` は既に θ/r_t です（[deltastar_integral.py:221](/home/sano/work/forge-integ-1005/design/forge_design/feedback/deltastar_integral.py:221)）。抽出スクリプトはさらに r_t で割って `theta_rt` に保存しています（[delta_contur_compare.py:99](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/delta_contur_compare.py:99)）。この列は **13.0437倍過大**です。現在の δ 比・C_f 比の集計では使われていないため、それらを撤回する理由にはなりませんが、H・N の診断に転用する前に直してください。 |
| **Major** | **第2層の交差検証：方向は採用、現データでの選定は保留** | case/48 は CPG の基準として有用ですが、断熱基準には壁温ドリフトが記録されています。case/44 の `run_0114` も熱負荷が変動中で、ω 残差が末尾上昇しています（[case/44 README:718](/home/sano/work/forge-integ-1005/case/44.windtunnel_nozzle/README.md:718)）。出口 M の静定だけでは熱閉包の検証になりません。**固定形状で、較正条件と検証条件を分離**し、δ_r・θ・C_f・T_w の時系列と抽出定義を揃えてください。`run_0115` の再設計成功は、壁温を跨ぐ閉包の予測成功とは別です。 |
| **Major** | **較正ガード：修正付き採用** | 壁温だけの照合では弱く、全 `prepare` を一律停止するのは強すぎます。**較正係数の流用と生産採用**を、熱条件・組成／物性DB・輸送モデル・熱閉包版・入口 θ₀・検証済みの作動範囲に結び付けて制限してください。未較正の初期壁による検証計算は可能にします。壁温変更後の「NS 1回」は再評価の開始であって、再較正の完了保証にはしません。現行の記録項目では不足しています（[deltastar_integral.py:228](/home/sano/work/forge-integ-1005/design/forge_design/feedback/deltastar_integral.py:228)）。 |
| **Major** | **壁温分布：配管は採用、物理的適用性は別途検証** | 単一の物理座標表から NS と CONTUR を生成する案は妥当です。ただし `Tw_table` は局所の T_w を代入するだけで、独立した熱境界層の発達履歴を持ちません（[deltastar_integral.py:132](/home/sano/work/forge-integ-1005/design/forge_design/feedback/deltastar_integral.py:132)）。急な冷却開始・加熱への切替まで妥当とは言えません。同条件判定は平均温度ではなく、**座標原点・単位・補間／外挿規則を含む区分関数の一致**で行います。物理 [m] を正本とし、CONTUR へ渡す時だけ x/r_t に変換。表の被覆不足を現在の `np.interp` の端値保持で黙認しないでください。 |
| **Minor** | **旧生産保存・新式で再較正・#8d 後回し：採用** | 式と係数を版で結び、旧壁・旧係数・生成条件を保存してください。ただし「NS と同じ Sutherland」と「旧 CONTUR と丸め差だけ」は両立しません。旧式は T₀=273.15、S=110.4（[deltastar.py:17](/home/sano/work/forge-integ-1005/design/forge_design/feedback/deltastar.py:17)）、NS は273、111です。1470 K の粘性差は **0.162%**。旧版再現と新モデルの CPG 極限試験を分けます。 |

第2層では、当面 **現行の van Driest II 型を基準として維持**することを推奨します。Eckert は単に μ の評価温度だけを替える方式として混ぜず、参照状態と摩擦相関を含む別の閉包として比較してください。参照エンタルピーの代表形は `h* = h_e + 0.5(h_w−h_e) + 0.22(h_aw−h_e)` です。[NASA TP-2914 §4.2](https://ntrs.nasa.gov/api/citations/19890017745/downloads/19890017745.pdf)

比較の事前基準は「最大誤差が小さい方」だけでは不十分です。較正係数と較正区間を凍結し、検証側で δ_r・θ・C_f の許容差を別々に設定してください。δ_r が零を跨ぐ上流では相対誤差を使わず絶対誤差にします。#8d の試験部トレンド **0.3%** は維持し、抽出感度・時間変動・格子差がそれより大きければ判定保留とします。

後段で追加する NS の最小候補は **case/45 の固定生産形状・300 K** です。既存断熱 run と同じメッシュが冷却時の壁解像を満たせるなら1本。再メッシュが必要なら、そのメッシュで断熱と300 Kの2本を揃えます。この対でも壁温分布の妥当性までは証明できません。

結論: **係数と輸送物性を固定し、熱閉包だけを温度形からエンタルピー形へ替える CFD 0 step の A/B を先に行ってください。**

第1仮説: **現行の熱閉包は、断熱から冷却壁への δ_r の変化率に1%以上の系統差を生む。** 確度: 中。  
根拠: `deltastar_integral.py:49,126`。同じ `GasSemiPerfect` で独立再計算すると、M=6、T_e=233.482 K、T_t=1600 K、r=0.89628 に対し、旧式 T_aw は **1663.033 K**、エンタルピー回復は **1472.684 K**です。ただし、これが δ_r に与える全効果は未測定です。  
反証条件: 下記 A/B で、試験部全域の冷却／断熱 δ_r 比の変化が1%未満となること。「1%以上の影響」という仮説を棄却します。熱力学的な不整合の存在を否定するものではありません。

第2仮説: **空気 Sutherland と混合輸送物性の差が壁温間の係数流用を妨げる。** 確度: 中。μ の差は再現しましたが、δ_r への寄与は未分離です。

第3仮説: **整合化後に残る厚さ・傾きの誤差は、N・温度分布・入口履歴の閉包にある。** 確度: 低。θ の比較修正と運動量式の整合化前には確定できません。

判別 A/B:  
A は現行の温度形、B は次のエンタルピー形とし、**変更因子は熱閉包の表現だけ**にします。

`h_aw = h_e + r(h_0−h_e)`  
`h(v) = h_w + a(h_aw−h_w)v + [h_e−a(h_aw−h_w)−h_w]v²`  
`v = u/u_e`、`T(v) = h⁻¹(h(v))`

両腕で r=0.72¹ᐟ³、a=1、k_f=k_N=1、現行 μ、入口 θ₀、縁条件、運動量式、求積・平滑化を固定。断熱／1000／600／300 Kについて入口から出口まで積分し、CFD は起動しません。CPG 極限の回帰と、積分精度による差が0.1%未満であることを前提に、試験部 [40,94] の δ_r・θ・H・C_f と各冷却／断熱比を比較します。

**B による300 K／断熱の δ_r 比の変化が最大1%以上なら第1仮説を支持、1%未満なら棄却。** NS への誤差が減るかは別欄に記録し、ここでは再較正しません。

やらない方がよいこと: **T_aw だけを T_t でクリップすること、旧 k_f の持ち越し、δ だけで4係数を合わせること、NS と近いという理由で分子 Pr を0.72に置き換えること。**

呼び出し側の前提への異議:  
「T_w を実測値に替えても傾きが残るので、T_aw は原因でない」は受け入れません。また、NS の摩擦に近いことは、H 側だけが原因という証明ではありません。

保存済み `extract.npz` の通算 **80000〜100000 step**、x/r_t=40・70・94 の δ 比と C_f 比を `check_quasisteady.py` で再判定した結果は **`OVERALL: ALL STEADY`** でした（全5点、drift／fluctuation許容0.1%）。一方、対象の

- `case/45.isobutane_m6_d155/run_0167_ns_n012_N2/`
- `case/45.isobutane_m6_d155/run_0179_ns_n012_N2_ext/`

は保存済み判定でともに **`NOT CONVERGED (stalled/plateau)`** です。ここで確認したのは保存された抽出量の静定までです。

不足情報: 元 run の HDF5・全残差系列・実効設定がこのワークスペースにないため、元場の NaN、T_w の時系列、接線壁応力による C_f、局所組成の一定性を再確認できません。case/44 の旧 run を現行輸送モデルの検証に使うには、当時の実効モデルとDBの照合も必要です。

ファイル変更・forge 起動はしていません。**plan 未反映**です。採用時の反映先は `tooling-nozzle-isothermal-wall-chain.md` §4.6・§5.1 と、`verification-m6-axis-wave-mesh-su2.md` §5.1 #8d です。
