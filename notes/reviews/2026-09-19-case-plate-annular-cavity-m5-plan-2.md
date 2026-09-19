# codex レビュー: case-plate-annular-cavity-m5 (plan)

- **plan**: [`plans/active/case-plate-annular-cavity-m5.md`](../../plans/active/case-plate-annular-cavity-m5.md)
- **stage**: `plan`
- **date**: 2026-09-19
- **commit**: `836d4d80` (feature/sern-design)
- **codex**: effort `high`, 7.2 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M9/m1
- **focus**: 前回 NO-GO (C2/M8/m1) を全件採用した改稿版。① 採用の実装が正しいか (3D restart 廃止/検査体積と q_w/面分類/変換時の壁タグ/準定常の列定義/EOS 決定ゲート) ② 新しく入れた 2D 前駆 + inletProfile による領域縮小 (§4.3) の妥当性と、外部 BL の法線解像 19 点で足りるかの判断基準 ③ URANS 非定常性チェック (§4.9) の設計 (物理 dt・時間窓・判定) ④ 偏心/すきま/マッハ数が可変というユーザ要件に対するチェーンの単一ソース化 ⑤ 残っている見落とし
- **extra**: `case/49.plate_annular_cavity_m5/cad/build_geom.py`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
前回の主要な誤りは修正されていますが、半割 URANS、境界層比較、偏心時の評価、準定常判定には実装前の修正が必要です。  
推奨は、下記のゲートを改稿してから評価ツールと Stage A に進むことです。現状のゲートで生産計算を承認することはできません。

`plans/README.md` と関連 `accepted/` を確認した範囲では、本目的を解決済みの計画との重複はありません。node・SLAU・低 Re SST・block-DPLUR の選択と、`case/48` を先行検証に使う方針は妥当です。共有コードを変更しない本件で、cell・周期・軸対称の検証を追加する必要はありません。

前回指摘への対応では、3D 最近傍 restart の廃止、最終壁タグでの変換、`pRef`、静止キャビティ IC、`cyl_top` の収支からの除外は適切です。`case/48.flat_plate_cooled_m4/run_0005_B_tw300_y3/res_wall_4_48000.h5` の `qwall/utau/ypls` は、今回も **各 1001 点すべてゼロ**と確認しました。ただし、勾配による熱流束評価と離散収支の関係は未解決です。

1. **Major — 半割 URANS を、半割仮定の検証開始条件にできない。**

   **根拠:** [plan:344](/home/sano/work/forge/plans/active/case-plate-annular-cavity-m5.md:344) は反対称モードが禁止されると認識しながら、354 行では「反対称的な挙動を示唆したら」全周計算を行います。`y=0` の `slip` は法線流速を除去するため、その制約で排除したモードを検出できません。[boundaryCond_d.cu:68](/home/sano/work/forge/solver_density_cuda/cuda_forge/boundaryCond_d.cu:68)

   また、`dT_up/dT_dn` は流れ方向の上流／下流差であり、半割が禁止する左右非対称の指標ではありません。

   **対案:** **Stage A 全周 URANS を無条件の必須ゲート**にする。小さな左右非対称擾乱を与え、左右の温度・開口交換流量と反対称成分の成長／減衰を測る。全周でも完全対称 IC のまま無擾乱で回すだけでは、検査として弱いです。

2. **Major — URANS の設定と時間精度・観測窓が、非定常性の否定に不足する。**

   **根拠:** [plan:347](/home/sano/work/forge/plans/active/case-plate-annular-cavity-m5.md:347) は物理 CFL、`nSubIterDualTime=10`、音響往復時間の約 10 倍を指定するだけです。ところが実装は **`unsteady:1` と `deltaT.control:0` が必要**で、定常設定の `control:1` を継承すると例外終了します。[main.cpp:1724](/home/sano/work/forge/solver_density_cuda/main.cpp:1724)

   内部反復 10 回の根拠は別ケースの限定条件であり、現行手順も **`nSub` 倍増比較**を要求しています。[recommended-settings.md:184](/home/sano/work/forge/procedures/recommended-settings.md:184)  
   また、音響時間は深部の交換・熱緩和時間を保証しません。粗い Stage A で振動が減衰しても、空間・時間離散化による減衰を排除できません。

   **対案:** 実行可能な URANS 設定を別 YAML として定義し、`dt` 半減、`nSub` 倍増、開口せん断層の細化を必須にする。平均・振幅・卓越周波数の変化を比較し、時間・内部反復誤差は結論量の **2 % 以下**を目標にする。2.5 ms は最初の確認時刻とし、観測窓を倍増して深部温度・入熱・交換量の統計が安定するまで延長する。時間刻みの細化による確認は、[NASA の検証手順](https://www.grc.nasa.gov/www/wind/valid/tutorial/spatconv.html)とも整合します。

3. **Major — 指定した `gen --table` 経路では、前駆解の `Ps` が失われる。**

   **根拠:** [plan:195](/home/sano/work/forge/plans/active/case-plate-annular-cavity-m5.md:195) の原始変数表を渡す場合、[gen_inlet_profile.py:198](/home/sano/work/forge/solver_density_cuda/tools/gen_inlet_profile.py:198) は `Ps` を直接出力列から除外します。その後、`Tt/M` による換算を行わない経路では `Ps` を出力しません。[同:276](/home/sano/work/forge/solver_density_cuda/tools/gen_inlet_profile.py:276)

   実関数をメモリ内で実行した結果:

   ```text
   入力: z ro Ux Ps k omega
   出力: z k omega ro Ux
   ```

   このままでは圧力分布が BC の一様値に戻ります。さらに計画の抽出量には、発達境界層の法線速度がありません。

   **対案:** case 側で **`z ro Ux Uy Uz Ps k omega` を直接書き出す**。2D の法線速度を 3D の `Uz` に写し、スパン速度を `Uy=0` とする。生成列・単位・`applied:` ログを検査し、入口の実際の `U,T,k,omega` と前駆解を照合する。亜音速部の閉包も CPG と TP で異なるため、「上書きするから成立する」だけでは検証になりません。[boundaryCond_d.cu:817](/home/sano/work/forge/solver_density_cuda/cuda_forge/boundaryCond_d.cu:817)

4. **Major — 「50 mm 下流でも入口の厚さから 5 % 以内」は、BL 解像度の検査になっていない。19 点の採否も判断できない。**

   **根拠:** [plan:199](/home/sano/work/forge/plans/active/case-plate-annular-cavity-m5.md:199) は、流入後の物理的な境界層発達を比較基準に含めていません。

   既存の断熱平板 `case/48.flat_plate_cooled_m4/run_0004_A_ad_y3_cont/` を再評価すると:

   | 量 | x≈0.30 m | x≈0.35 m | 下流増加 |
   |---|---:|---:|---:|
   | `delta*` | 2.8439 mm | 3.1819 mm | 11.88 % |
   | `theta` | 0.27668 mm | 0.31467 mm | 13.73 % |

   両断面・両量の全スナップショットを `check_quasisteady.py` の判定関数に渡した **VERDICT は `STEADY`**。一方、run 全体の `check_convergence.py` は **`NOT CONVERGED (stalled/plateau)`**です。これは M5 の予測値ではありませんが、下流発達をゼロと扱えない実例です。

   さらに「19 点」は `plate_in` の概算で、入口側の `plate` は 2.5 mm サイズです。[plan:228](/home/sano/work/forge/plans/active/case-plate-annular-cavity-m5.md:228) 面サイズだけでは tet 領域の法線配置も保証されません。

   **対案:** 前駆 station を `x₀` とし、3D の 50 mm 下流は **前駆の `x₀+50 mm`**と比較する。まずキャビティを塞いだ同じ 3D メッシュ方式で入口移送・BL 発達を検証する。厚さだけでなく `U,T,k,omega` の分布と壁摩擦を比較し、外層・prism→tet 遷移を細化した差を確認する。**19 点を先験的に合格とはしないが、試験格子として禁止する根拠もありません。採否はこの比較で決めるべきです。**

5. **Major — 偏心に対する単一ソース化が、評価領域まで届いていない。**

   **根拠:** [plan:327](/home/sano/work/forge/plans/active/case-plate-annular-cavity-m5.md:327) の `Ri < r < Ro` は同心時だけ正しい定義です。偏心時の正しい断面は

   \[
   x^2+y^2<R_o^2,\qquad (x-x_{\rm off})^2+y^2>R_i^2
   \]

   です。既定半径・`x_off=1 mm` では、同心マスクは真の半環流体断面を **22.498 mm²、約 12.06 %**取りこぼし、同面積の固体側を誤って含めます。面積総和が同じでも誤りを検出できません。

   `conditions.json` と `geom_config.json` の関係、すきま中央の定義、`delta/gap` の基準幅も未確定です。

   **対案:** 条件と形状を束ねる最上位入力を一つ定義し、`setup.py` が解決済み manifest を生成する。CAD・メッシュ・IC・BC・評価の全段がそれを読む。CV マスク、偏心円との交点による測線、周方向平均の重み、測点の深さ依存、`delta/gap_nom` 等の定義を共通関数にする。M・すきま変更時の適用範囲と再検証条件も manifest に含める。

6. **Major — 壁温ピンを持つ node 離散化で、勾配入熱と離散エネルギー収支の同一視が残る。**

   **根拠:** [plan:334](/home/sano/work/forge/plans/active/case-plate-annular-cavity-m5.md:334) の片側二次差分は、物理壁熱流束の推定方法として妥当です。しかし forge は壁ノードのエネルギー残差をゼロ化します。[nodeWallDirichlet_d.cu:85](/home/sano/work/forge/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:85)  
   ソルバの壁流束も、後処理の三点差分とは別の勾配を使用します。[viscousFlux_d.cu:521](/home/sano/work/forge/solver_density_cuda/cuda_forge/viscousFlux_d.cu:521)

   したがって、物理壁の勾配積分が、温度拘束による実効的なエネルギー授受と有限格子上で一致する保証はありません。特にリップでは、開口をまたぐ dual CV と等温／断熱共有ノードの扱いが必要です。

   **対案:** 評価ツールに **物理的な壁熱流束評価と離散保存検査を別々に実装**する。非拘束 CV 群と壁隣接界面の流束で離散収支を検算し、勾配入熱との差が細化で減ることを確認する。平板に加え、角・曲面を持つ伝導問題で法線、重み、共有ノード処理を検証する。開口の `h0` は保存された値と採用 EOS の datum を使う。

7. **Major — 準定常判定がゼロ近傍で破綻し、URANS の合格条件も互いに矛盾する。**

   **根拠:** [check_quasisteady.py:261](/home/sano/work/forge/solver_density_cuda/tools/check_quasisteady.py:261) は `max(abs(mean),1e-30)` で正規化します。[plan:315](/home/sano/work/forge/plans/active/case-plate-annular-cavity-m5.md:315) の `mdot_net` をそのまま全列 `STEADY` の対象にすると、ゼロへの収束を拒否します。

   実関数による合成系列の再現では、**±1e−12 のゼロ平均系列も、1e−10→1e−14 に減衰する系列も `DRIFTING`**でした。`dT≈0`、入熱≈0、侵入深さゼロでも同種の問題があります。

   また §4.9 は「定常解との差 10 %」、§6.4 は「振幅 10 %」、§8 は「全量 `STEADY`」で、同じ合否条件ではありません。[plan:351](/home/sano/work/forge/plans/active/case-plate-annular-cavity-m5.md:351)、[plan:442](/home/sano/work/forge/plans/active/case-plate-annular-cavity-m5.md:442)

   **対案:** `mdot_net` は交換流量で正規化した不釣合いと絶対許容で判定する。温度・入熱・侵入深さにも、事前に固定した絶対許容と非ゼロの基準尺度を設ける。既存ツールには適切に尺度化した検査列を渡し、生値も併記する。定常解と周期統計の合格条件を分け、URANS では平均差・振幅・窓間 drift をすべて判定する。

8. **Major — Stage A の「EOS 差 3 % 未満」で CPG を確定するには、誤差と評価量の範囲が不足する。**

   **根拠:** [plan:367](/home/sano/work/forge/plans/active/case-plate-annular-cavity-m5.md:367) は三つの量だけで EOS を確定しますが、数値側には格子差 5 %、収支差 5 %、温度 drift の例として 5 K が許されています。これでは **3 % 未満という比較結果の確かさ**が示せません。

   また深部温度が壁温近傍なら、`dT_mid` と `zpen` の差が小さくても、開口温度分布や局所熱流束の差は残り得ます。成果物は温度場全体です。

   **対案:** Stage A の EOS 選定を暫定判定にする。CPG を採用する条件には、数値不確かさを含む比較、`dT_mouth`・温度プロファイル・壁別入熱を追加する。境界近傍の判定なら生産解像度で再比較する。前駆解・入口 CSV も EOS ごとに整合させ、同じ作動条件・BL 条件で比較したことを記録する。**誤差込みで 3 % 未満を示せない場合は TP を採用**する。

9. **Major — 追加メッシュ品質ゲートが、まだ実行可能な定量仕様になっていない。**

   **根拠:** [plan:242](/home/sano/work/forge/plans/active/case-plate-annular-cavity-m5.md:242) の sliver 指標には下限値がありません。`V/Lmax³` は正常な薄い prism でも小さくなるため、全要素共通の閾値では壁層と sliver を区別できません。

   また変換器の閉性ログは、**全体の最大面積**で正規化します。[gmshReader.hpp:2082](/home/sano/work/forge/solver_density_cuda/mesh/gmshReader.hpp:2082) 粗い外部セルと微小壁セルが混在する本ケースでは、これだけで微小 CV の局所閉性を保証できません。

   **対案:** tet/prism 別に指標・数値閾値・不合格条件を確定する。prism は局所 Jacobian の符号も検査し、AR 緩和は対象セルだけに適用する。閉性は CV ごとの面積尺度でも判定する。層あり率の対象には、等温の `cyl_top` を含む **壁 6 グループ**を明記する。

10. **Minor — 参照と幾何数値に修正漏れがある。**

    **根拠:** [plan:9](/home/sano/work/forge/plans/active/case-plate-annular-cavity-m5.md:9) の `methods/turbulence.md` は存在しません。[plans/README.md:26](/home/sano/work/forge/plans/README.md:26) も旧版の「高度確定・CPG」のままです。  
    また第一層 10 µm・伸長率 1.25・14 層の等比和は **0.869495 mm**、最終層は **0.181899 mm**で、記載値と一致しません。

    **対案:** 現行の `methods/turbulence/theory.md`・`implementation.md` を参照し、索引を同期する。層厚・派生幾何量は manifest から生成して転記誤差をなくす。

**推奨は、評価可能な Stage A を先に成立させる順序へ改稿することです。** 実装前の修正優先順は次のとおりです。

1. 単一入力・偏心 CV・入口 CSV の仕様を確定する。
2. 熱収支、ゼロ近傍の準定常判定、メッシュ品質の検査仕様を確定する。
3. 前駆→3D 平板の移送検証と、全周 URANS の時間・格子・観測窓検証を必須化する。
4. それらの誤差評価を通した EOS 決定を、Stage B 移行条件にする。

ファイル変更・新規 CFD 計算投入は行っていません。本回答はレビュー提案であり、**plan 未反映**です。

指摘数: Critical 0 / Major 9 / Minor 1
