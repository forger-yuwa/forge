# 非平衡凝縮ソースの θ 律速を残差から外し、定常解を擬似時間刻みに依存させない

## メタ

- **area**: `condensation / time_integration`
- **status**: `in_progress`
- **related_docs**:
  - `methods/condensation.md` (実装 §4 ソース項 / §4c θ 律速の dt 依存)
  - `methods/time_integration/` (scalar-DPLUR の対角 `src_jac`)
- **related_plans**: [condensation-followups.md](condensation-followups.md) F-cf7 (起票元), [condensation-equilibrium-eos.md](../accepted/condensation-equilibrium-eos.md) (緩和形の θ)
- **created**: `2026-09-15`
- **owner**: Claude (session 011spCYH)

## 1. 目的

非平衡凝縮ソース (`condensationSourceKernels_d.cuh`) の律速係数 θ は「1 擬似ステップあたりの Δg ≤ `dg_max`=5e-3、
潜熱 ΔT ≤ `dT_max`=1 K、蒸気枯渇 Δg ≤ 0.9(Y_w−g)」を `Δg = S_g·Δτ_loc/ρ` で評価し、**その θ を定常残差のソース
$S_{Q_0..Q_2}, S_g$ に掛けている**。Δτ_loc (局所擬似時間刻み) に比例するため、定常局所時間刻み (`unsteady 0`) の
収束解が `cfl_pseudo` とセル体積に依存する。case/44 (6 m ノズル, Δτ_loc 2.6e-5 s) では ΔT/step 3.4 K → θ 0.25 (内側) /
0.55 (壁 2 ノード: Δτ_loc が半分) で、成長が 1/4 に絞られ「壁第一層だけ液相が速く増える」偽の壁異常と凝縮完了の
3–4 $r_t$ 遅れを生んだ (cfl 2/1/0.5 で θ 0.25/0.5/1、出口 g 0.437/0.574/0.584 %; `case/44 run_0127/0130/0131`)。
完了時: 収束解が Δτ に依存せず (cfl 2 と 0.5 で g 場が一致)、起動時の安定性 (onset の急な潜熱放出) は従来どおり保たれる。

## 2. スコープ

- **やる** (定常・非平衡 `condEquilibrium 0` に限定): 凝縮・**蒸発**の両ソースを「状態だけで決まる瞬間速度」に書き直し Δτ を残差から完全に外す (§4.2)。
  θ は**未クリップの 4 モーメント候補増分をまとめて縮める更新クランプ**に限定する (§4.2-3)。ヤコビアン `sj_g`/`sj_Q1` の θ 倍も外す。
  診断 `condLim_<s>` は更新クランプ係数 (蒸発分岐でも書く)。double / float (`condFloat 0/1`) の両カーネルを同時に直す。
- **やる**: 「状態を固定して Δτ だけ変えても残差 (凝縮・蒸発・消滅) が変わらない」単体試験を必須にする (§6)。
- **やらない (スコープ外として明記)**: 核生成率・成長則・物性の変更。平衡緩和形 `condEquilibrium 1` (§4.3: 旧モデル互換のため据え置き、輸送との釣り合いが
  Δτ 依存なのは既知の制約として記録)。**dual-time / RK 陽解法の時間精度保証** (凝縮モーメントには BDF 物理時間項が無い = 別課題 F-cf8; 新経路の初回適用は
  定常 point-implicit 更新に限定し、RK 経路は `condLimiterMode 0` を既定のまま残す)。

## 3. 関連 docs と前提

- [methods/condensation.md](../../methods/condensation.md) 実装 §4「ソース項」「安定化」: θ の定義と `dg_max`/`dT_max`
  (ハードコード `condensationSource_d.cu` L33–35)。§4c (本 plan で追加) に dt 依存の記録。
- [methods/time_integration](../../methods/time_integration/) scalar-DPLUR: 対角 $D_\phi = V/\Delta\tau + V\,\mathrm{src\_jac} + \mathrm{transport\_diag}$
  (`scalarTransport_d.cu` L282–286)。`src_jac_g` は `condensationTransport_d_wrapper` で 0 クリア → `condensationSource_d_wrapper` が
  `sj_g = −θ ∂S_g/∂T · ∂T/∂(ρg)` を書き、後段の陰的スカラー更新が使う (`main.cpp` L1205–1206 の順)。
- 潜熱はエネルギー方程式のソースではなく二相 EOS の温度反転 (`dependentVariables_d.cu`) で T に現れる。したがって「ΔT ≤ 1 K/step」は
  ρg の更新量にだけ掛ければよく、roe には触らない。
- Wysłouzil (case/16, mm ノズル) は Δτ_loc が小さく ΔT/step ≤ 0.86 K で θ≡1 → 過去の検証は無影響 (本 plan の回帰基準)。

## 4. 設計方針

### 4.1 何が間違っていたか

定常擬似時間反復は $V\,\Delta(\rho\phi)/\Delta\tau = -R(\phi)$ の固定点 $R(\phi^\ast)=0$ を求める。ソースに Δτ 依存の係数
$\theta(\Delta\tau)$ を掛けると $R$ 自体が Δτ の関数になり、固定点 $\phi^\ast(\Delta\tau)$ が動く。「1 step の変化量を抑える」意図の
安全弁を残差側に置いたのが原因。

### 4.2 修正: 残差は Δτ を含まない瞬間速度、θ は未クリップ候補増分のクランプ

1. **凝縮側の残差**: $S_{Q_0..Q_2}, S_g$ をそのまま `res_*` に加える (θ 倍を削除)。蒸気枯渇 ($Y_w-g\le0$ → $S=0$) と $J$ 上限 (`Jmax`) は
   Δτ を含まない物理条件なので残差側に残す。
2. **蒸発側の残差** (codex M2): 現行 `cond_evap_source` は $\lambda=1+\dot r\,\Delta\tau/r_{30}$, $S_g=\rho g(\lambda^3-1)/\Delta\tau$ で、
   制限が非作動でも $S_g=\rho g(3a+3a^2\Delta\tau+a^3\Delta\tau^2)$, $a=\dot r/r_{30}$ と Δτ が残る。これを**瞬間速度形** (成長側と同じ一様 $\dot r$ のモーメント形;
   $\dot r<0$ を $r_{30}$ で評価) $S_g=4\pi\rho_l Q_2\,\dot r$, $S_{Q_1}=Q_0\dot r$, $S_{Q_2}=2Q_1\dot r$, $S_{Q_0}=0$ に書き換える (実装 `cond_evap_source_rate{,_f}`;
   初回実装の自己相似形 $aq_1,2aq_2,3a\rho g$ は monodisperse でのみ一致するので codex result M4 で本式に揃えた。多分散の値検証は単体 (g))。
   完全蒸発による数密度の消滅は従来どおり実現可能性クランプ ($r_{30}<2r_{min}$, $Q_0=0$) が確定する。半径半減 ($\lambda\ge\tfrac12$)・Δg・ΔT の制限は更新クランプに移す。
   **ソース積分の体積は `volumePartial_d`** (node 周期 seam の合併体積二重計上の既存欠陥; codex result M6。case/09 1 step 試験で旧 2.0/8.0 → 1.000)。
   **新経路は `condEquilibrium 0` のみ** (平衡形は従来更新; codex result M5)。
3. **ヤコビアン**: `sj_g`, `sj_Q1` の θ 倍を外す (負帰還はそのまま point-implicit の対角へ)。蒸発側も瞬間速度形の $\partial S/\partial(\rho\phi)$ に合わせる。
4. **更新クランプ** (新設 `condensation_update_limiter_d`, codex M1/M6 反映):
   - 位置: 凝縮モーメントは NS の block-DPLUR sweep 内ではなく**流れ・化学種更新後の独立した point-implicit 更新** (`main.cpp` の `condensationTransport`→
     `condensationSource`→ scalar 更新)。その更新で **floor/硬クランプを掛ける前の候補増分** $\delta_\phi=\Delta(\rho\phi)$ (4 本) を取り出し、
     クランプ後に確定 → primitive 同期、の順に固定する。
   - 評価状態 (codex result M3 で確定): $\Delta g=\delta_g/\rho^{new}$ = **更新済みの流れ $\rho^{new}$ を固定したモーメント修正 (相変化分)**。
     移流による密度変化分 $g^{old}\Delta\rho/\rho$ は潜熱を伴わず流れの CFL が抑えるので数えない (初稿の $(\delta_g-g^{old}\Delta\rho)/\rho^{new}$ は撤回)。
     $T$, $c_p$ は前回の従属変数評価 (流れ更新前) を使う: 安全弁の見積りであり固定点 ($\delta\to0$) には影響しない。潜熱 ΔT は EOS と同じ有効比熱
     $c_{v,\rm eff}=c_v+g(R_w-dL/dT)$ で $\Delta T=\Delta g\,L/c_{v,\rm eff}$。密度が変わっても $\delta_g=0$ なら無作用 (単体 (i))。
   - 係数: $\Delta g>0$ のとき $\theta_u=\min(1,\ dg_{max}/\Delta g,\ dT_{max}/\Delta T,\ \mathrm{avail}/\Delta g)$、$\Delta g<0$ (蒸発) のとき
     $\theta_u=\min(1,\ dg_{max}/|\Delta g|,\ dT_{max}/|\Delta T|,\ (1-\lambda_{min}^3)\,g^{old}/|\Delta g|)$。**$\theta_u$ は常に $>0$** (停止穴を作らない:
     `avail=0` は残差側で $S=0$ なので候補増分が 0)。4 本の増分を同じ $\theta_u$ で縮め、その後に既存の非負 floor / $\rho g\le\rho Y_w$ / $0.99\rho$ を
     「補正量を計測しながら」適用する (補正量の総和を診断に残し、収束時 0 を確認)。
   - 消滅: クランプ後 $r_{30}<r_{min}$ かつ $g<g_{min}$ のセルは 4 モーメントを 0 (現行の消滅処理を移設)。
5. **診断**: `condLim_<s>` = $\theta_u$ (更新クランプ)、`condClampCorr_<s>` = このステップの**全**硬クランプ (更新 floor + 実現可能性 $g\le Y_w$/$0.99$ + 液滴消滅) による
   $|\Delta\rho g|/\rho$ の累積、`condClampCorrQ_<s>` = $Q_0..Q_2$ の最大相対補正 (codex result M2)。従来の θ (残差律速) は消える。
6. **設定**: `condDgMaxStep` / `condDTmaxStep` (既定 5e-3 / 1 K), `condLimiterMode: 0` (旧: 残差に θ) / `1` (新)。**既定は定常 point-implicit 経路で 1、
   RK 陽解法・dual-time では 0 のまま** (§2)。回帰後に旧経路を削除する時期は followups で決める。

### 4.3 平衡緩和形 (`condEquilibrium 1`) は据え置き — 理由は互換性とスコープ分離 (codex M3/m1 で訂正)

$S_g=\alpha\rho(g_{eq}-g)/\Delta\tau_{loc}$ は輸送残差と釣り合うので、定常条件は $R_{\rm transport}+V\alpha\theta\rho(g_{eq}-g)/\Delta\tau=0$ で
一般には $g\ne g_{eq}$ — **固定点も Δτ 依存** (「θ は接近速度だけ」は輸送の無い局所緩和にしか成り立たない)。本 plan では旧モデル互換のため触らず、
この制約を [methods/condensation.md](../../methods/condensation.md) に既知の制約として記録する。平衡凝縮を選ぶ場合の推奨は EOS 拘束形
`condEquilibrium 2` (代数拘束、Δτ 非依存)。なお `condEquilibrium` の**設定既定値は 0 (非平衡)** であり 2 は既定ではない。

### 4.3b dual-time / RK は初回スコープ外 (codex M4/M6)

現行の BDF 物理時間項と時間レベルのシフトは平均流 5 変数と `roK/roOmega` のみで、`rog`/`roQ0..2` には物理時間項が無い (サブ反復ごとに
`N/M` を現在値へコピーして定常と同じ point-implicit 更新)。したがって dual-time で凝縮モーメントが物理時間で積分される保証は現状無く、
本 plan の新旧一致を時間精度の保証には使わない。followups に **F-cf8: 凝縮モーメントの dual-time 物理時間項** を登録する。RK 陽解法は
未制限残差の累積バッファを持つため段ごとの制限では不十分で、専用試験なしに新経路を既定適用しない。

### 4.4 安定性の見立て

θ を外した残差は onset 直後に大きくなるが、(i) implicit の対角に `sj_g` (潜熱負帰還) が入る、(ii) 更新クランプで 1 step の
Δg / ΔT は従来と同じ上限に抑えられる、ので**起動時の安定性は従来と同等**のはず。RK 陽解法・dual-time は本 plan の対象外で旧経路に自動降格する (§4.3b)。
懸念は「クランプが常時効くセル (θ_u≪1 が収束後も残る)」= 残差が消えないのに更新が止まる状態。$\theta_u>0$ を保証しても収束が
遅くなるだけで固定点は変わらないが、§6 では **未加工残差** (`res_rog` 等の RMS) と `condClampCorr` の収束時ゼロ化を合否に含め、`condLim≈1` だけを固定点の証拠にしない。

## 5. 実装ステップ

1. `methods/condensation.md` 実装 §4 に「4c. θ 律速の dt 依存と更新クランプ化」を追記 (本 plan の §4.1–4.2 の要約 + 発見経緯)。
2. `cuda_forge/condensationSourceKernels_d.cuh` (double / float): θ 倍を残差・ヤコビアンから外す。`avail ≤ 0` と `Jmax` は残す。
   **蒸発 `cond_evap_source{,_f}` を瞬間速度形に書き換え** (§4.2-2)、消滅と半径半減制限を更新クランプへ移す。`diagLim` は更新クランプ側で書く。
3. `cuda_forge/condensationSource_d.cu` + `scalarTransport_d.cu` の凝縮モーメント更新: **候補増分 (floor 前) を取り出す経路**を作り、新 kernel
   `condensation_update_limiter_d` で 4 本同率クランプ → 硬クランプ (補正量計測) → 確定 → primitive 同期。前段値は既存 `...N`/`...M` バッファ。
4. `main.cpp`: 定常 point-implicit のモーメント更新直後に update limiter を呼ぶ。RK / dual-time 経路は `condLimiterMode 0` のまま (§4.3b)。
5. `input/solverConfig.{hpp,cpp}`: `condDgMaxStep`, `condDTmaxStep`, `condLimiterMode` (既定 1)。
6. 単体テスト `tests/unit/test_cond_update_limiter.cu` / `test_cond_source_dt_invariance.cu`: (a) **状態固定で Δτ を 1e-7〜1e-3 に振っても凝縮・蒸発・消滅の
   残差が不変** (double/float)、(b) Δg 上限内で θ_u=1・増分不変、(c) 上限超で 4 本同率縮小、(d) 収束状態 (Δ=0) で無作用、(e) 輸送とソースが非ゼロで釣り合う
   1 セル移流モデルで固定点が Δτ 非依存、(f) 蒸気枯渇・負増分・密度/組成変化・float32 極小増分。
7. 回帰 (§6) → docs 同期 → `condLimiterMode 0` 経路の削除は次 plan (followups) に送る。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 |
| --- | --- | --- |
| 1 | ~~codex plan レビュー~~ | 2026-09-15 実施 (§6.1)。M1–M6/m1 を全て採用し §2/§4/§6 に反映済み。**実装着手可** (ユーザ確認後) |
| 2 | ~~methods §4c 追記~~ | 実装後の仕様に書き換え済み (2026-09-15) |
| 3 | ~~実装 (ステップ 2–5)~~ | 済 (2026-09-15): `condLimiterMode`/`condDgMaxStep`/`condDTmaxStep`、`cond_evap_source_rate{,_f}`、`cond_moment_update_limited_d` (`condensationUpdateLimiter_d.cuh`)、RK/dual-time は起動時に 0 へ自動降格、診断 `condClampCorr_<s>` |
| 3b | followups F-cf8 登録 | 凝縮モーメントの dual-time 物理時間項 (BDF 残差・対角・時間レベルシフト) — 本 plan の外 |
| 4 | ~~単体テスト~~ | 済: `tests/unit/test_cond_limiter_steady.cu` (Δτ 不変性 double/float ビット一致 + 更新クランプ 8 ケース ALL PASS)、`test_cond_float_device.cu` 両モード PASS |
| 5 | 回帰 run (result-2 で再開) | 1 巡目 (2026-09-15): case/44 cfl 0.5/2/4 × condFloat 0/1 一致 (`run_0132`–`0136`); Arthur N2 cell/node・空気 cell (`case/34 run_0105`–`0107`) 参照とノイズ床内 PASS・onset 同一; Wysłouzil (`case/16 run_0470`) 報告量同一・場の差 ≤1e-4 |
| 6 | codex result レビュー (3 回目) → accepted | 2 回目 NO-GO (M3/m2) を反映中: `run_0147` (真の継続) / `0148`・`0149` (cell condFloat 0) / `0150` (cell mode 0 床) と単体 (h)(j) 拡張。 1 回目 NO-GO の M1–M7/m1 は反映済み (§6.1)。残: Wys の mode 1/0 の step 時間差 (同時実行時 7.4 vs 4.4 ms/step → 単独計測は §9 に記載)、F-cf9 (旧 mode 0 削除) は followups |

## 6. 検証

- **単体 / ビルド**: `cmake --build build` (double + float 経路)、§5-6 のテスト (Δτ 不変性を含む)。
- **収束・定常の扱い (codex M5, result-2 M1 で明文化)**: 比較する run は `check_convergence.py` **PASS** を原則とする。PASS が取れない case では、各モーメント残差が「3 桁以上低下」**または**「同一条件の旧経路 (mode 0) run と同じ床 (plateau 値が ±30 %)」に到達していることを残差ゲートとし、報告は「定常解の一致」でなく **「プラトー上の最終場一致 + 凝縮固有量 series STEADY + 補正 0」** と表現する (床が case 固有 = 修正と無関係であることを mode 0 との比較で示す)。case/44 のような warm 床 plateau
  (rms_roe 床 0.42 が dry 一様 run と同値) で PASS が取れない場合は、(i) 未加工残差 `res_rog`/`res_roQ*` の RMS が 3 桁以上低下、
  (ii) 凝縮固有量 (出口 g 平均 / g max / onset x / 壁 M) の時系列を CSV 化して `check_quasisteady.py --series-csv` で
  **ドリフト・変動 0.2 % 以内の STEADY**、(iii) `condClampCorr` 総和が収束時 0、を全て満たすことを PASS 相当とし、そう報告する。
  終了は step 数でなくこのゲートで決める。既存 `run_0131` は**参考** (最終場の 1 点で θ<1、`NOT CONVERGED` plateau) であり正解データではない。
- **差のノルム**: 凝縮領域 $\Omega_c$ = 両 run の $g>10^{-4}$ の和集合。$\|g_a-g_b\|_{L^1(\Omega_c)}/\|g_b\|_{L^1(\Omega_c)}\le1\%$、
  出口 g 平均 (質量流束重み, $x=x_{max}-2r_t$) の相対差 ≤1 %、onset x (壁流線, $g>10^{-3}Y_w$) の差 ≤0.1 $r_t$、$\Omega_c$ の $|\Delta T|\le0.5$ K、$|\Delta M|\le5\times10^{-3}$。
- **検証ケース**:
  1. **Δτ 非依存** (主目的): case/44 `run_0127` プロトコル (入口 Tt 分布, node Euler TP) を新バイナリで cfl_pseudo 2 と 0.5 で回し、上のノルムで一致。
     node と cell、`condFloat 0/1` の 4 組合せ (cell はメッシュを cell 変換して同条件)。
  2. **回帰 (無影響)**: case/16 Wysłouzil `run_0335` プロトコル (θ≡1 だった; 収束 PASS が取れる) を新バイナリで再実行し、壁圧・onset・g が
     同一バイナリ 3 反復のノイズ床以内。case/34 Arthur (N2 空気, CPG carrier, 蒸発分岐あり) で onset ±0.02 in 以内 — **蒸発の瞬間速度形化の影響は
     ここで初めて出るので、差が出たら量と符号を記録し採否を判断** (蒸発が効く Wysłouzil 出口側も同様)。
  3. **起動安定性**: case/44 一様 IC (`run_0092` の mid 段場) から cfl 2 で入口 Tt 分布を与えて起動、NaN 0、`condLim` 最小値の推移
     (起動時 <1 → 収束時 ≈1)。cfl 4 の挙動も記録 (合否外)。
  4. **dual-time / RK**: 合否に含めない (§4.3b)。`condLimiterMode 0` のままで従来と同一 (ビット一致) を確認するだけ。
- **判定基準**: 上のノルム・ノイズ床・NaN 0・収束ゲート。`check_convergence.py` / `check_quasisteady.py` の VERDICT を添付。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| result (2 回目) | `2026-09-16` | [2026-09-16-condensation-source-limiter-steady-result.md](../../notes/reviews/2026-09-16-condensation-source-limiter-steady-result.md) | **NO-GO**, C0/M3/m2 | **全採用 (2026-09-16 反映)**: M1 (残差ゲート未達: node roQ0 2.7 桁, cell 1.1–2.3 桁) → §6 のゲートを「3 桁低下、または旧経路 (mode 0) と同じ床 (plateau ±30 %) に到達」に明文化し、cell は mode 0 の同条件 run (`run_0150`) で床を実測、報告は「プラトー上の最終場一致 + series STEADY」と表現; M2 (`run_0146` が roY を引き継がず液相が初期クランプで消えていた: `interp_field.py` が roY* を新設しない) → `interp_field.py` 修正 (roY* 新設)、`run_0146` は「流れを引き継ぎ化学種・凝縮場を再発達させた 48000 step」と訂正、修正後の真の継続 `run_0147` で不変を確認; M3 (cell × condFloat 0 なし、固定点試験がクランプを通さず残差/補正を判定しない、corQ 未検証) → `run_0148`/`0149` (cell condFloat 0, cfl 2/0.5)、単体 (h) に `cond_realizability_clamp_d` を通し未加工残差 <1e-6・補正 0 を判定、(j) Q 負値 floor の corQ=1 を検証; m1 (旧経路の補正診断が累積) → 旧経路でもステージ 0 でリセット (更新 floor は未記録と明記); m2 (run_0131 正本の残記述、F-cf7 状態、Wys の正式ツール判定) → README/methods/followups 同期、`wys_cond_series.py` + `check_quasisteady --series-csv` の VERDICT を run に保存 (lim1b/lim0b ALL STEADY) |
| result (1 回目) | `2026-09-15` | [2026-09-15-condensation-source-limiter-steady-result.md](../../notes/reviews/2026-09-15-condensation-source-limiter-steady-result.md) | **NO-GO**, C0/M7/m1 | **全採用 (2026-09-16 反映)**: M1 (収束/定常ゲート未達・Wys 出力 2 枚・Arthur 反復基準) → 凝縮固有量の時系列 CSV (`cond_series_csv.py`) + `check_quasisteady --series-csv` 0.2 % を全 run に、Wys は outStepInterval 4000 で再実行し `--series` STEADY、Arthur は新バイナリ反復 3 本の広がりを床に; M2 (補正記録が ρg 負値のみ) → 実現可能性クランプ (g≤Y_w/0.99・負値・消滅) も `condClampCorr`/`condClampCorrQ` に記録; M3 (Δg 評価状態) → §4.2-4 を「更新済み密度固定の相変化分」と確定 (密度変化分は数えない理由を記載, 単体 (i)); M4 (蒸発式) → plan の一様 ṙ 形に実装を揃え多分散の値検証 (単体 (g)); M5 (eq=1 に新クランプ) → `condEquilibrium 0` に限定し eq=1 の lim1/lim0 回帰 (`run_0143`–`0145`); M6 (周期 seam 二重計上) → `volumePartial_d` に修正し case/09 1 step 試験で旧 2.0/8.0 → 1.000; M7 (単体欠落) → 1 セル固定点 (h)・密度変化 (i)・cell 離散化 (`run_0141`/`0142`)・condFloat 0 cfl 0.5 (`run_0140`); m1 (文書同期) → §4.4/§7/README/recommended-settings 修正, F-cf9 登録 |
| plan | `2026-09-15` | [2026-09-15-condensation-source-limiter-steady-plan.md](../../notes/reviews/2026-09-15-condensation-source-limiter-steady-plan.md) | GO-with-changes, C0/M6/m1 | **全採用**: M1 (停止穴・未クリップ候補増分・補正量計測) → §4.2-4; M2 (蒸発の Δτ 依存 → 瞬間速度形) → §4.2-2, §5.1 #3; M3 (eq=1 の理由訂正) → §4.3; M4 (dual-time 物理時間項無し → スコープ外, F-cf8) → §4.3b, §5.1 #3b; M5 (収束 PASS/series ゲート・ノルム定義・run_0131 は参考) → §6; M6 (更新位置・評価状態・RK/node/cell/condFloat) → §4.2-4, §6; m1 (eq 既定 0) → §4.3 |

## 7. 影響範囲

- `cuda_forge/condensationSourceKernels_d.cuh`, `condensationSource_d.cu`, `main.cpp`, `input/solverConfig.{hpp,cpp}`, `tests/unit/`
- 凝縮 ON の定常 run は成長が速くなる (大型ノズルで顕著)。case/44 va3 の入口 Tt 分布 run の旧 `run_0131` (cfl 0.5, θ≡1) は比較参考 (正本は新バイナリの `run_0137` 系)。
  Wysłouzil / Arthur は無影響 (θ≡1)。
- docs: `methods/condensation.md` §4 (安定化の記述を更新クランプに書き換え), `methods/index.md` は変更なし、`procedures/recommended-settings.md` 凝縮節に
  「`condLim` が全域 ≈1 を確認」を追記。

## 8. 完了条件

- [ ] `methods/condensation.md` §4c 追記済み
- [ ] 実装・§6 の検証 1–4 を満たす
- [ ] codex レビュー 2 回 (`plan` / `result`) を §6.1 に記録し、Critical / Major の採否を §5.1 に反映済み
- [ ] `status: done`、§9 に変更ログ
- [ ] `plans/active/` → `plans/accepted/` へ移動、`plans/README.md` 同期

## 9. 変更ログ

- `2026-09-15` — 初稿 (followups F-cf7 から独立 plan 化)。case/44 cfl A/B (run_0127/0130/0131) で dt 依存を確定。
- `2026-09-15` — **実装** (commit 後述): 残差から θ を撤去 (蒸気枯渇→0 のみ残す)、蒸発を瞬間速度形 `cond_evap_source_rate{,_f}` に、更新クランプ kernel `cond_moment_update_limited_d` (未クリップ 4 本候補増分を同率縮小、θ_u≥1e-12、c_v,eff で ΔT)、config 3 キー、RK/dual-time 自動降格、診断 `condClampCorr`。単体: Δτ 1e-7/1e-3 で double/float ともビット一致 (H2O 962 状態・N2 300 状態; 旧 mode 0 は 525/169 セルで差)、更新クランプ 8 ケース、float/double 一致テスト両モード PASS。CFD (case/44 入口 Tt 分布, node Euler TP, 24000 step from 一様 IC): **mode 1 cfl 2 (`run_0132`) は旧 θ≡1 の cfl 0.5 (`run_0131`) と g L1 1.6e-5 / |ΔT| 0.007 K / |ΔM| 6e-5 で一致**、cfl 4 (`run_0134`) も cfl 2 と g L1 3.9e-5 / |ΔT| 0.010 K で一致・起動 NaN 0; 収束時 `condLim` は凝縮域で 1.000 (θ<1 が残る 62 セルは g≤2e-23 の数値塵, S 0.2–0.35, 総液相の 7e-25)、`condClampCorr` 0。旧 mode 0 cfl 2 (`run_0127`) との差は g L1 50 % (絞られていた分)。
- `2026-09-16` — **codex result レビュー (NO-GO, M7/m1) を全採用した 2 巡目**: 蒸発を一様 ṙ 形に (`cond_evap_source_rate{,_f}`: $q_0\dot r, 2q_1\dot r, 4\pi\rho_l q_2\dot r$; 多分散で自己相似形と 8 % 差), 新経路を `condEquilibrium 0` に限定, ソース積分を `volumePartial_d` に (node 周期 seam の二重計上を修正: case/09 `run_0060`/`0061` 1 step 試験 旧 2.0 (面)/8.0 (角) → 新 1.000), 補正診断を全クランプに拡張 (`condClampCorr`=|Δρg|/ρ 累積, `condClampCorrQ`=Q の最大相対補正)。単体: (g) 多分散蒸発の値 (float 差 7e-7), (h) 1 セル輸送+ソース固定点が dt 1e-6/1e-5/1e-4 で 3.8e-5 一致, (i) 密度半減で δ=0 なら無作用; 全 PASS。CFD (全て新バイナリ, 凝縮固有量 series 0.2 % で STEADY): node cfl 2 (`run_0137`) vs cfl 4 (`0138`) g L1 4.2e-5 / |ΔT| 0.007 K, vs cfl 0.5 warm (`0139`) 1.9e-3 / 0.21 K / |ΔM| 2.0e-3, 1 巡目 (`0132`) との差 2.1e-5 (蒸発形の効果は小)、`condLim` 凝縮域 1.000、`condClampCorr`/`Q` 凝縮域 0 (Q の floor は g≤6e-24 の乾き塵 222 セルのみ)。cell (`run_0141` cfl 2 STEADY; `0142`+`0146` cfl 0.5 72000 step STEADY): cell cfl 0.5 vs 2 で g L1 6.8e-4 / |ΔT| 0.12 K / |ΔM| 9.6e-4 (24000 step 時点は未収束で 2.8 K)。condFloat 0 cfl 0.5 (`run_0140`) は float cfl 0.5 と 1.8e-5 で一致・STEADY。eq=1 (`run_0143` lim1 / `0144` lim0 / `0145` lim1 反復): lim1−lim0 (ro 1.4e-6, g 5.7e-5) = 反復ノイズ (1.8e-6, 6.2e-5) 以内。Arthur (`lim1b`): N2 cell 3 反復とも 28/28 PASS (反復広がり ro 2.4–4.6e-4)・N2 node 28/28・空気 cell 28/28・onset 同一。Wysłouzil (`lim1b`/`lim0b`, 26 枚): 報告量同一, `--series` STEADY, rms_rog 5.1 桁, lim1b−lim0b g 1.2e-5, `condClampCorr` 凝縮域 0。
- `2026-09-16` — 蒸発の一様 ṙ 形の実現可能性ガード: `test_cond_float_device` が mode 1 で塵状態 (g=1e-13, Q0=1e19, S=0.999) の evap-edge ヤコビアン差 (|Δsj|·dt 7e3) で FAIL → 原因は不整合モーメントで $S_g=4\pi\rho_l q_2\dot r$ が液相の 1e6 倍/s・sj が 1e15 に発散すること。べき平均不等式による上限 ($q_{1e}=\min(q_1,q_0r_{30})$, $q_{2e}=\min(q_2,q_0r_{30}^2)$) と $r_{30}<2r_{min}$ の消滅待ち ($S=0$) を `cond_evap_source_rate{,_f}` に入れ、g 摂動幅を $10^{-3}\max(g,10^{-9})$ に。両単体 ALL PASS。**この最終バイナリで検証 run を全て取り直す** (`run_0151`–`0161`, Arthur `lim1c`, Wys `lim1c`/`lim0c`, case/09 `run_0062`)。
- `2026-09-16` — **codex result レビュー 2 回目 (NO-GO, M3/m2) を全採用**: `interp_field.py` が roY* を新設しない欠陥 (cell 継続 `run_0146` で液相が初期クランプで消えていた) を修正、旧経路の補正診断をステージ 0 でリセット、単体 (h) を実現可能性クランプ込みの全経路 + 未加工残差/補正判定に、(j) Q floor の corQ 検証を追加、§6 の残差ゲートを床基準で明文化。追加 run: `run_0147` (cell cfl 0.5 真の継続), `run_0148`/`0149` (cell condFloat 0, cfl 2/0.5), `run_0150` (cell cfl 2 mode 0 = 床の参照)。Wys の正式 `check_quasisteady --series-csv` (`wys_cond_series.py`) は lim1b/lim0b とも ALL STEADY。
- `2026-09-16` — step 時間 (Wysłouzil 2D node SST 凝縮, 4000 step 単独計測, RTX 3060): mode 1 = 3.28 / 3.30 / 3.27 ms/step、mode 0 = 3.09 / 3.10 / 3.12 ms/step → 更新クランプの追加コストは約 +6 % (48000 step 同時実行時の 7.4 vs 4.4 は GPU 共有の影響)。
- `2026-09-15` — **検証完了 (1 巡目; codex result で不足指摘)**: (§6-1) case/44 cfl 0.5 (`run_0133`+`0136`, 72000 step) vs cfl 2 (`run_0132`): g L1 2.1e-3 / |ΔT| 0.23 K / |ΔM| 2.2e-3 (ゲート 1 % / 0.5 K / 5e-3 内、なお緩やかに接近中)、cfl 4 (`run_0134`) 3.9e-5、condFloat 0 (`run_0135`) 1.8e-5; 報告量 (onset 15.92 r_t, 出口 g 0.584 %, M 4.050) は全て同一。(§6-2) Arthur N2 node (`case/34 run_0106`) lim1/lim0 28/28 PASS; N2 cell (`run_0105`) lim1 27/28 (ro 6.2e-4 = 床 2.1 倍) だが同一バイナリ反復 `lim1_r2` が同じセルで 5.6e-4 差 = cell atomicAdd 反復ノイズ (`lim1_r2` 28/28 PASS), onset 2.13 in = 参照; 空気 cell (`run_0107`, 蒸発分岐あり) lim1/lim0 28/28 PASS, onset 2.21 in = 参照。Wysłouzil (`case/16 run_0470`) lim1/lim0 とも報告量すべて run_0335 と同一 (onset 22.51 mm, 壁偏差 +2.46/−4.5/+5.1/+5.4 %); **lim1 vs lim0 は ro 1.9e-6 / g 9.9e-6** (バイナリ世代差 lim0 vs run_0335 の g 9.7e-5 より小) = 蒸発の瞬間速度形化を含め無影響。(§6-3) cfl 4 起動 NaN 0。(§6-4) RK/dual-time は自動降格で旧経路のまま。
- `2026-09-15` — codex plan レビュー (GO-with-changes, M6/m1) を全採用: 蒸発ソースも瞬間速度形へ、更新クランプは未クリップ候補増分に、dual-time/RK はスコープ外 (F-cf8)、検証ゲートを収束 PASS/series STEADY + ノルム定義に強化。

## 10. 未確定事項

- ~~更新クランプの呼び出し粒度 (RK 段ごと vs 段の最後だけ)~~ 決着 (2026-09-15, §4.3b): RK は初回スコープ外、定常 point-implicit のみ。
- `condLimiterMode 0` (旧経路) を残す期間。回帰が揃えば即削除でもよい。
