# codex レビュー: boundary-cht-axisymmetric-graetz (plan)

- **plan**: [`plans/active/boundary-cht-axisymmetric-graetz.md`](../../plans/active/boundary-cht-axisymmetric-graetz.md)
- **stage**: `plan`
- **date**: 2026-09-27
- **commit**: `6d32c645` (feature/cht-axisym-graetz)
- **codex**: effort `high`, 6.5 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M4/m2
- **extra**: `case/63.graetz_cht/make_run.py`, `case/63.graetz_cht/gen_mesh.py`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
目的は妥当で、親計画の未検証範囲を補っています。  
ただし、入口温度の合格条件は設定した問題と矛盾し、格子収束・評価範囲・連成ゲートにも修正が必要です。

読み取り専用で確認した結果、`graetz_ref.py --selftest` は `VERDICT: PASS`、3格子の品質検査もすべて `VERDICT: PASS`（最大 AR 32、skewness 0）でした。有限 `Pe` の A/B は窓内最大 **0.11314 %** を再現しました。一方、確認時点の `case/63.graetz_cht/run_0001_dry_r16`、`run_0002_smoke_iso_r16`、`run_0003_smoke_iso0_r16` は、`check_convergence.py` がすべて **NOT CONVERGED**。これらから定常解の精度は判断していません。

`plans/README.md` と親計画の検証範囲に重複はありません。ソルバ変更を伴わない今回の検証を node・FP64 に限定することも妥当です。cell・FP32・周期境界への保証には拡張できません。

1. **Major — V-g1 の加熱開始断面「温度幅 ≤0.05 K」は、正しい基準解でも FAIL する。**

   **根拠:** [plan:155](/home/sano/work/forge-cht/plans/active/boundary-cht-axisymmetric-graetz.md:155) は `x=0` の断面温度一様性を要求します。しかし、その断面の壁は加熱壁であり、[温度ピン実装:80](/home/sano/work/forge-cht/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:80) は指定温度を壁ノードへ直接代入します。既存 `r32` メッシュにも `x=0, r=R` のノードが存在します。

   計画の `ellip` 条件（`nr=120, nx=2400, Pe=720`）を再計算すると、10 K 加熱時の `x=0` は軸で約 **300.000 K**、壁の第一内点で **306.779 K**、壁で **310.000 K**。温度幅は約 **10 K** です。壁点を除くだけでも解決しません。軸方向伝導による加熱開始点より上流への予熱は、[extended Graetz の一次文献](https://www.sciencedirect.com/science/article/pii/S0017931000001071)でも扱われています。

   **対案:** 一様入口温度の検査を実際の入口 `x=-L_up` に移す。`x=0` は有限 `Pe` 基準解の温度分布との比較に変更し、物理的な予熱と入口 BC の不整合を分離してください。

2. **Major — V-g3 は、モデル差を含む誤差の減少を「格子収束」の必要条件にしている。**

   **根拠:** [plan:160](/home/sano/work/forge-cht/plans/active/boundary-cht-axisymmetric-graetz.md:160) は古典 Graetz に対する誤差の単調減少を要求します。一方、[誤差予算:112](/home/sano/work/forge-cht/plans/active/boundary-cht-axisymmetric-graetz.md:112) は有限 `Pe` のモデル差を認めています。有限壁抵抗・圧縮性の差も残ります。

   したがって、`Nu_h−Nu_Graetz = δ_model + C h^p` です。例えば正規化したモデル差 `+0.001`、格子誤差 `−0.004, −0.001, −0.00025` は完全な2次収束ですが、基準解への絶対誤差は **0.3 % → 0 % → 0.075 %** となり、登録条件に落ちます。また、3格子から格子間差で得られる観測次数は1つで、「2組の次数の ±30 % 一致」は確認できません。

   **対案:** V-g2 の解析解との比較は維持し、V-g3 は共通位置での**格子間差**とその減少率を判定対象に変更する。3格子では観測次数を参考値に留め、漸近域確認済みとは扱わないでください。固定している固体半径8層の影響も、別途感度確認か誤差上限の根拠が必要です。

3. **Major — 窓全域の最大誤差を合否判定するのに、欠損と非定常を見逃す設計になっている。**

   **根拠:** [plan:95](/home/sano/work/forge-cht/plans/active/boundary-cht-axisymmetric-graetz.md:95) は `iface_ok=1` の節点だけを使用しますが、窓内の必要節点が欠けた場合の拒否条件がありません。実装では [conjugateWall.cpp:273](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:273) で診断値の設定を飛ばすため、NaN 検査だけでは欠損を検出できません。

   また、[V-g0:152](/home/sano/work/forge-cht/plans/active/boundary-cht-axisymmetric-graetz.md:152) の準定常検査は4点と総入熱ですが、V-g2・V-g4・V-g5 は窓内最大値です。観測点間の局所ドリフトや、総入熱で相殺する節点間振動は通過できます。

   **対案:** 窓内の期待節点集合をメッシュから固定し、欠損・重複・`iface_ok=0` があれば `REFUSED`。その全節点で分子・分母・`Nu` の準定常性と変動 ≤0.03 % を検査してください。負例には「最大誤差節点の欠落」と「4観測点を避けた局所ドリフト」を追加すべきです。

4. **Major — G-if / G-cons の事前登録が、実行可能な仕様まで閉じていない。**

   **根拠:** [plan:151](/home/sano/work/forge-cht/plans/active/boundary-cht-axisymmetric-graetz.md:151) の「親 plan の登録値の形」では数値が一意になりません。親の [G-if:199](/home/sano/work/forge-cht/plans/accepted/boundary-cht-axisymmetric-fem2d.md:199) は絶対許容 `0.46 W/m²`、確認時点の [case/63 テンプレート:47](/home/sano/work/forge-cht/case/63.graetz_cht/template/solverConfig_cht.yaml:47) は `1.0 W/m²` です。

   `G-cons ≤1e-9（相対）` にも分母・絶対床・絶対許容がありません。[check_cht_balance.py:139](/home/sano/work/forge-cht/solver_density_cuda/tools/check_cht_balance.py:139) は `--q-floor` を必須とし、相対・絶対条件の両方で判定します。`tol_solid=1e-9` の固体内部残差とは別の量です。

   **対案:** §6 に G-if の全パラメータと単位、G-cons の計算式・`q_floor`・絶対許容・対象出力を明記する。特に対照 run は熱量尺度が加熱 run と異なるため、同じ規格化を使う根拠も登録してください。テンプレートはその登録値から生成する形が適切です。

5. **Minor — `make_run.py --dT` と実際の固体外面温度の一致を保証していない。**

   **根拠:** [make_run.py:107](/home/sano/work/forge-cht/case/63.graetz_cht/make_run.py:107) は `Ts` と記録文を変更しますが、[同:117](/home/sano/work/forge-cht/case/63.graetz_cht/make_run.py:117) は固体ファイルをそのままコピーします。実際の外面温度は [solidFem2d.cpp:66](/home/sano/work/forge-cht/solver_density_cuda/conjugate/solidFem2d.cpp:66) が読む `ROBIN/TC` です。

   現在の [gen_solid.py:50](/home/sano/work/forge-cht/case/63.graetz_cht/gen_solid.py:50) は ΔT ごとの生成に対応していますが、異なる ΔT の固体を誤指定しても拒否されません。

   **対案:** run 作成前に `ROBIN/TC == T_IN + dT` を検査し、不一致は拒否する。記録にも固体から読み出した実値を残してください。

6. **Minor — 完了時の主張と §4.4 に、差し引き検証の限界が反映されていない。**

   **根拠:** [plan:89](/home/sano/work/forge-cht/plans/active/boundary-cht-axisymmetric-graetz.md:89) は依然として G2−G1 を「連成が入れた誤差」としていますが、V-g4 は有限壁抵抗の物理差を含むと訂正しています。また、[目的:21](/home/sano/work/forge-cht/plans/active/boundary-cht-axisymmetric-graetz.md:21) の「壁熱流束」の検証という表現は広すぎます。加熱・対照に共通の加算誤差は差し引きで消えるため、差分 `Nu` の合格だけでは絶対熱流束の精度を保証しません。

   **対案:** 完了時の主張を「登録した node・FP64 条件で、対照差し引き後の局所 `Nu` が許容差内」と明記し、§4.4 も V-g4 と統一してください。

**推奨:** この検証方針を維持し、実装前に **①入口検査の修正 → ②格子収束判定の分離 → ③全節点の評価・準定常ゲート → ④連成ゲートの数値登録** の順で plan を直してください。基準解と格子生成を作り直す必要はありません。

ファイルは変更していません。以上の提案は **plan 未反映** です。

指摘数: Critical 0 / Major 4 / Minor 2
