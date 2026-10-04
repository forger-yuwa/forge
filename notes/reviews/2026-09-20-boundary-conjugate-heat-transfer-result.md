# codex レビュー: boundary-conjugate-heat-transfer (result)

- **plan**: [`plans/active/boundary-conjugate-heat-transfer.md`](../../plans/active/boundary-conjugate-heat-transfer.md)
- **stage**: `result` (diff base `e3890d22`)
- **date**: 2026-09-20
- **commit**: `714e168f` (feature/sern-design)
- **codex**: effort `high`, 11.2 min, rc=0
- **判定**: **NO-GO**, 指摘 C1/M9/m1
- **focus**: 重点: (1) §4.2 の更新式と cht_loop の D_f 初期値。2026-09-20 に第一セル伝導 k_eff A/d1 から熱伝達係数 h*A へ変えた (31 倍過大で反復が止まっていた) が、安定性 (過小評価で発散しないか) の根拠が薄い。Df_safety=2 で足りるか、上界の条件を示せ。(2) §4.3 の界面熱量。現状 iface_q_compact を一次に使っているが、拘束反力込みの実効熱量 Q_f = ΣF^E - C は未接続 (依存 plan 待ち)。この状態で V5 の結論 (連成は正しい/誤差は遷移モデル) を出してよいか。(3) V5 の結論の妥当性。case/53 C3X と case/54 Mark II の README にある領域別の数字 (正圧面 +0.8 %/-0.2 K、負圧面層流域 +74.8 %/+42.8 K) と、内部条件を逆算した残差床 (7.8-9.9 K) の比較で「連成は正しい」と言えるか。逆算 (infer_internal_bc.py) が壁温を使っているので循環論法になっていないか検証せよ。(4) 残差は NOT CONVERGED (plateau) のまま check_quasisteady の STEADY で報告している。この扱いの是非。(5) accepted へ移してよいか、残すべき残作業は何か。読むべきパス: plans/active/boundary-conjugate-heat-transfer.md, case/53.c3x_vane_cht/README.md, case/54.markii_vane_cht/README.md, solver_density_cuda/tools/cht_loop.py, solver_density_cuda/tools/solid_fem2d.py, solver_density_cuda/conjugateWall.cpp, case/53.c3x_vane_cht/tools/{infer_internal_bc.py,compare_h.py,smooth_profile.py}

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**

保存的な界面熱量と合格ゲートが未実装で、反復処理にも再現可能な欠陥があります。  
README の領域別数値は再現できましたが、「連成は正しい／差の主因は遷移モデル」という結論は支持できません。`accepted` への移動は止めるべきです。

`e3890d22...714e168f` の diff、指定コード、保存済み run を確認しました。ファイル変更・新規 CFD 実行はしていません。固体単体試験は両方とも `VERDICT: PASS (all)`、V1 の専用検証も記載値を再現しました。ただし、以下の欠陥はその試験範囲外です。

1. **Critical — 実装が §4.3 の保存的界面契約を満たしていない**

   [conjugateWall.cpp:189](/home/sano/work/forge/solver_density_cuda/conjugateWall.cpp:189) はコンパクト差分を計算し、[cht_loop.py:268](/home/sano/work/forge/solver_density_cuda/tools/cht_loop.py:268) はそれに面積を掛けて固体へ渡しています。`iface_q_eff` の出力はなく、V5 の最終壁ダンプにも存在しません。ソルバ内連成も [conjugateWall.cpp:363](/home/sano/work/forge/solver_density_cuda/conjugateWall.cpp:363) の第一セル抵抗で更新しています。

   したがって、報告された「両側熱量の一致」は**採用したコンパクト差分と固体側の一致**です。拘束反力を含む流体のエネルギー授受が閉じた証拠ではありません。[verify_v1.py:54](/home/sano/work/forge/case/52.conjugate_slab/verify_v1.py:54) の `PASS` も、この限定された照合です。§4.3 が撤回した方式を合格根拠に戻しています。

   **対案:** 同一評価状態の \(Q_f=\sum F^E-C\) を接続し、dual-time／周期の解除試験、共有角の一回集計、G-cons を通す。それまでは V1・V4・V5 を保存的 CHT の合格実績として扱わない。

2. **Major — 棄却後に異なる状態の \(T\) と \(Q_f\) を組み合わせ、メリット関数の重みも混在させている**

   [solid_shell.py:423](/home/sano/work/forge/solver_density_cuda/tools/solid_shell.py:423) は棄却時に温度を `best` へ戻しますが、435 行の右辺は**棄却された温度で評価した `Qf`**のままです。さらに `Df` を倍増しても `best[0]` を再評価せず、次回は新しい重みの \(\Phi\) と古い重みの値を比較します。§4.2 の禁止事項そのものです。

   読み取り専用の再現計算で、\(A_s=1,\ b_s=0,\ Q(T)=3300-10T,\ T_0=301,\ D=1\) とすると、初回候補は 295.5 K。棄却後、301 K に戻しながら \(Q(295.5)=345\) を使い、315.6667 K を返しました。同じ状態の \(Q(301)=290\) を使う更新なら 297.3333 K です。

   初回にも同じ不整合があります。[cht_loop.py:215](/home/sano/work/forge/solver_density_cuda/tools/cht_loop.py:215) でテンプレートの分布を使う一方、245 行では一様 `T0` を作ります。実 run では次の組み合わせでした。

   | run | 初回 CFD に課した壁温 | ドライバが仮定した壁温 |
   |---|---:|---:|
   | `case/53.c3x_vane_cht/run_0004_cht/` | 512.61–611.82 K | 566 K 一様 |
   | `case/54.markii_vane_cht/run_0003_cht/` | 470.80–597.14 K | 537 K 一様 |

   **対案:** 温度・熱量・固体状態を同じ評価点の組として保存・復元する。重み変更時は基準メリットを再計算し、履歴を破棄する。初期温度は実際に課した分布から取得し、計画どおり line search と再試行上限を実装する。

3. **Major — `Df_safety=2` は上界でも安定保証でもない**

   [cht_loop.py:250](/home/sano/work/forge/solver_density_cuda/tools/cht_loop.py:250) の \(hA=|q|A/(T_g-T_0)\) は温度差に対する比であり、非対角応答 \(H=-\partial Q_f/\partial T\) の上界ではありません。

   線形・対称な応答で \(A_s\succ0,\ H\succeq0\) と仮定すると、素の反復の収束条件は

   \[
   \rho\!\left((A_s+D)^{-1}(D-H)\right)<1
   \quad\Longleftrightarrow\quad A_s+2D-H\succ0 .
   \]

   保守的な十分条件は \(D-H\succeq0\) です。例えば対角 \(D\) に対して
   \(D_{ii}\ge\sum_j|H_{ij}|\) なら保証できます。しかし `2*h*A` がこの条件を満たす証拠はありません。非対称・非線形・未収束 CFD には、この限定条件さえそのまま適用できません。

   計画の反例を再計算すると、

   \[
   H=\begin{pmatrix}1.1806&-1.0556\\-1.0556&1.1806\end{pmatrix},
   \quad A_s=0.1I,\quad hA=H\mathbf1=0.125\mathbf1
   \]

   に **安全率 2** を掛けても、反復固有値は **0.3571、−5.6749** で発散します。実履歴でも `Df_mean` は C3X で初期の128倍、Mark II で256倍になっています。

   **対案:** `hA` は初期推定として採用し、安定性は修正済みの受理・退避処理で担保する。低周波だけでなく交番温度摂動と CFD 緩和長依存を試験し、「安全率2で十分」という保証は撤回する。

4. **Major — `fem2d` の外部連成では局所 \(k_s(T)\) が解かれていない**

   [solid_fem2d.py:146](/home/sano/work/forge/solver_density_cuda/tools/solid_fem2d.py:146) は `self.u is None` のとき、全固体節点を界面平均温度にします。`self.u` を更新するのは `recover_interior()` ですが、`FixedPointDriver.advance()` はこれを呼びません。

   したがって V5 が使う経路は、実質的に**全域一様の \(k_s(\overline{T_w})\)** です。`solve()` の局所温度による Picard とは別の方程式を解いています。

   既存円環試験の形状に温度依存物性を与えて比較すると、全系 Picard は壁温 **494.8202 K**、外部ドライバは **492.7460 K**。**2.0742 K 違う状態で `converged=True`** でした。既存 FEM 試験は定数物性なので検出しません。

   **対案:** 各評価温度に対して内部温度と物性を自己整合させ、内部残差も判定する。温度依存円環で `driver` と全系求解を照合し、V5 と物性感度を再評価する。

5. **Major — 現在の収束ゲートは G-if を満たさず、熱量不釣合い100%でも合格できる**

   [solid_shell.py:419](/home/sano/work/forge/solver_density_cuda/tools/solid_shell.py:419) は残差を `max(|Qf|, |b|)` で割っています。これは計画の `max|Qf|` による規格化と異なり、背面温度を含む大きな `b` が判定を緩めます。局所面積あたりの絶対残差、固体内部残差、float32 量子化も判定していません。

   実コードで \(A_s=1000\) W/K、\(b=300000\) W、\(T=300.002\) K、\(Q_f=1\) W、\(D=10^9\) W/K を与えると、2 回で **`converged=True`**。表示 `res_rel=3.33e-6` に対し、物理的な \(|r|/|Q_f|\) は **約1＝100%** です。

   **対案:** §6 の G-if をそのまま実装し、絶対熱流束残差・相対残差・内部残差・温度更新・連続回数を独立に満たす場合だけ合格させる。`Df` 増大による微小更新を収束と認めない。

6. **Major — V5 は同定データへの再適合であり、独立した連成検証になっていない**

   [infer_internal_bc.py:154](/home/sano/work/forge/case/53.c3x_vane_cht/tools/infer_internal_bc.py:154) は実測 \(T_w\) を外部熱流束の生成に使い、[同:178](/home/sano/work/forge/case/53.c3x_vane_cht/tools/infer_internal_bc.py:178) は同じ \(T_w\) との差を最小化して孔の \(h_c\) を同定します。その条件で同じ壁温に近づくことを、[Mark II README:91](/home/sano/work/forge/case/54.markii_vane_cht/README.md:91) は連成の正しさの確認としています。**この独立検証の主張には循環があります。**

   数字自体は再現しました。Mark II 正圧面は \(h\) バイアス **+0.8246%** ですが RMS は **15.98%**。壁温バイアス **−0.1726 K** に対して RMS **12.04 K**、最大誤差 **32.53 K** です。平均の相殺を一致の根拠にはできません。また、全域の同定 RMS 7.78 K と一領域の CHT RMS 12.04 K は、同じ誤差尺度ではありません。

   7.8–9.9 K は**同定残差**であり、独立に評価した予測不確かさの下限ではありません。`Tc` と `hc` を再同定して差が小さいことも、縮退方向の感度を示すにとどまります。

   **対案:** V5 を「同定条件下での整合性評価」に格下げする。独立した冷却条件または未使用データによる検証を用意し、SU2 対照・遷移感度・格子感度を通してから誤差原因を判定する。現時点で遷移は有力仮説ですが、主因の確定はできません。

7. **Major — `STEADY` が残差収束と局所量の定常性の代用になっている**

   `check_convergence.py` を再実行した結果は次のとおりです。判定区間は各 run の現在の `residual_history.csv`、CHT は末尾の単一反復です。

   | run／区間 | 再判定 |
   |---|---|
   | `case/53.c3x_vane_cht/run_0003_smoothgeom/` | `NOT CONVERGED (stalled/plateau)` |
   | `case/54.markii_vane_cht/run_0002_shortexit/` | `NOT CONVERGED (stalled/plateau)` |
   | `case/53.c3x_vane_cht/run_0004_cht/it_024/` | `NOT CONVERGED (stalled/plateau)` |
   | `case/54.markii_vane_cht/run_0003_cht/it_021/` | `NOT CONVERGED (stalled/plateau)`、`rms_roK` は `RISING` |

   C3X の連成履歴末尾は `dTw=0.02170 K`、`res_rel=0.001203`、`converged=0`。既定温度許容 0.001 K を**約22倍**超えており、[README:201](/home/sano/work/forge/case/53.c3x_vane_cht/README.md:201) の「相対許容をわずかに超えた」という説明は不十分です。Mark II 末尾は `rejected=1`、`converged=0` です。

   積分量の `STEADY` も局所分布を保証しません。例として3%比較に対応する `drift=osc=0.006` で、保存された4時刻を `check_quasisteady.py` の判定関数に通すと、比較対象範囲内の結果は次でした。

   | 対象 | `STEADY` | `DRIFTING` | `TRANSIENT-UNSETTLED` | `OSCILLATING` |
   |---|---:|---:|---:|---:|
   | C3X 局所 \(q\)、402点 | 121 | 255 | 15 | 11 |
   | Mark II 局所 \(q\)、392点 | 276 | 109 | 7 | 0 |

   同じ条件で積分 `q_total` は C3X が `STEADY`、Mark II は `TRANSIENT-UNSETTLED` でした。

   **対案:** 現結果は「未収束の参考値」として報告する。報告対象の局所量にも事前登録した許容を適用し、流体残差・界面ゲート・準定常を別々に判定する。なお、両 CHT の全残差履歴と各反復最終場では NaN/Inf は0でしたが、これは収束の代替になりません。

8. **Major — 共有角の競合検査が連成中の競合を防げない**

   [conjugateWall.cpp:222](/home/sano/work/forge/solver_density_cuda/conjugateWall.cpp:222) は `wallProfile` と `interfaceDiag` が無い場合に検査を省略します。`conjugate` 自体は検査の起動条件ではありません。また244行は初期壁温が同じなら通します。

   連成壁と通常等温壁を同じ300 Kで起動すれば通過し、その後は [同:353](/home/sano/work/forge/solver_density_cuda/conjugateWall.cpp:353) が連成壁だけを更新します。共有 CV に異なる温度が入り、計画が防ぐはずの「後勝ち」が再発します。`conjugateGroup` と一意の温度 DOF も未実装です。初期化には、文書で拒否するとした周期・壁関数併用等の検査もありません。

   **対案:** 初期温度の一致ではなく、共有 CV の拘束所有者を検査する。対応する共有連成 DOF を実装するまでは競合構成を拒否し、未対応構成も明示的に拒否する。

9. **Major — 外部ループの「同一メッシュ index コピー」が実際は2次元最近傍転送**

   [cht_loop.py:208](/home/sano/work/forge/solver_density_cuda/tools/cht_loop.py:208) は `interp_field.py` を呼びます。同ツールの [centroids():31](/home/sano/work/forge/solver_density_cuda/tools/interp_field.py:31) は座標を `[:, :2]` に切り、最近傍探索で転送します。

   3D シェル連成では同じ \((x,y)\) の異なる \(z\) 節点を区別できず、温度・速度・乱流量を別の spanwise 位置から引き継ぎ得ます。現在の平面 V5 では露呈しません。

   **対案:** 同一メッシュの同一性を確認したうえで保存量を index コピーする。3D の spanwise 非一様場を使い、再開によって場が変わらないことを検証する。

10. **Major — cell の `wallProfile` が壁面重心でなく内部セル重心を参照する**

    [boundaryCond.cpp:331](/home/sano/work/forge/solver_density_cuda/boundaryCond.cpp:331) は、node でない場合に `msh.cells[icw].centCoords` を使います。計画 §4.5 と現在仕様が約束する「cell は面重心」と異なります。

    斜交セルや壁法線方向にも変化するプロファイルでは、壁面以外の位置の温度を境界値として与えます。

    **対案:** node のみノード座標を使い、cell は `msh.planes[ip].centCoords` に戻す。

11. **Minor — 文書と残作業表が現在状態に同期していない**

    [methods/index.md:22](/home/sano/work/forge/methods/index.md:22) は「仕様のみ確定・未実装」、[design/CAPABILITIES.md:58](/home/sano/work/forge/design/CAPABILITIES.md:58) は「選ばれたら専用 plan」のままです。§4.2／[methods/boundary.md:332](/home/sano/work/forge/methods/boundary.md:332) の初期推定・line search と実コードも一致しません。

    [残作業表:394](/home/sano/work/forge/plans/active/boundary-conjugate-heat-transfer.md:394) は #10 で `fem2d` を未完、#24 で完了としています。#16 のデータ整備も古い状態です。一方、V2・V3・V6、ソルバ内 `shell2d`、区間ハッシュは未完のままです。[stage_manifest.py:75](/home/sano/work/forge/solver_density_cuda/tools/stage_manifest.py:75) に外部入力ハッシュはありません。

    **対案:** 仕様・実装済み範囲・検証済み範囲を揃え、残作業に所有先と解除条件を付ける。壁ピン順序など既存 run にも効く変更と、別 plan 所管の limiter／設計ツール変更について、回帰証跡の参照先も明記する。

**推奨は、現在の plan を `active` に残し、保存的な連成の検証を完了させてから再レビューすることです。** 優先順は次に絞ります。

1. 実効界面熱量、同一状態の受け渡し、受理処理、局所物性、G-cons／G-if を修正する。
2. 欠陥を検出する単体試験と、V1・V2・V3・V4、共有角・周期・再開の検証を通す。
3. V5 の「連成の正しさを確認した」という結論を撤回し、独立検証と誤差要因の切り分けを行う。
4. V6・Phase 2 の残作業と文書を完了し、最終実装に対する回帰結果を揃える。

このレビューによる plan 更新は、read-only 指示に従い**未反映**です。

指摘数: Critical 1 / Major 9 / Minor 1
