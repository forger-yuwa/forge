# 空気 (N2/O2) 自体の凝縮: CPG carrier 形 (N2 選択凝縮 + O2 キャリア)・N2 潜熱/飽和圧の低温整合・Arthur / Daum–Gyarmathy 検証

## メタ

- **area**: `condensation`
- **status**: `in_progress`  <!-- codex plan 1 回目 NO-GO (C1/M8/m1) → v2 → 2 回目 GO-with-changes (M5/m2) 全件採用 → v3 で実装 -->
- **related_docs**:
  - `methods/condensation.md` (§8 N2 物性「潜熱と飽和圧の低温整合」節・新設「空気: CPG carrier 形」節)
  - [notes/investigations/condensation-carrier-kantrowitz-air-survey.md](../../notes/investigations/condensation-carrier-kantrowitz-air-survey.md) §3 (onset データ、O2 物性、forge との差分)
- **related_plans**:
  - [condensation-nonequilibrium.md](../accepted/condensation-nonequilibrium.md) (N2 モデルの導入元、case/34 Fig.2 検証)
  - [condensation-kantrowitz-gamma-twophase-sonic.md](condensation-kantrowitz-gamma-twophase-sonic.md) (二相音速; pure/CPG は `n2_latent` の $L'>0$ で対象外 → 本 plan で潜熱は直すが CPG 二相音速の適用は後続 §5.1)
- **created**: `2026-09-12`
- **owner**: `CFD Dev`
- **branch**: `feature/condensation-air` (worktree `../forge-cond`)

## 1. 目的

極超音速ノズルの強膨張で**試験ガスの空気そのもの**が凝縮する解析を forge で行えるようにし、Arthur (case/34) と Daum & Gyarmathy (1968) の onset データで検証する。

1. **空気 = CPG carrier 形 (N2 選択凝縮 + 非凝縮 O2 キャリア)**: Daum & Gyarmathy の「低圧では空気は純 N2 として振る舞い N2 の自発核生成が onset を決める」に従い、
   気相は空気 (CPG, $R_{air}$, $c_{p,air}$)、凝縮するのは N2 だけ (質量分率 $Y_w=Y_{N_2}$ 固定)、O2 (+Ar) は凝縮しないキャリアとする。
   $g$ は**総混合物に対する液 N2 の質量分率**で、EOS・核生成・成長・蒸発・枯渇上限・流束・実現可能性の全てで同じ定義を使う (codex C1)。
   O2/N2 理想溶液の露点線 (混合液モデル) は本 plan から外し後続に置く。
2. **N2 潜熱と飽和圧の低温整合**: `n2_latent` (Lin 式 26) は ≈61 K で $L'>0$ となり、55 K 未満で $L'>c_{p,v}$ となって液比熱 $c_l=c_{p,v}-L'$ が負になる。70 K 未満を $c_l$ 一定の線形外挿にし、
   **同じ $L(T)$ で飽和圧の低温 Clausius–Clapeyron 外挿も再構成**する (現在の `n2_psat` は `n2_latent(50)` を C–C に直接使うので、潜熱だけ変えると飽和圧が動く; codex M2)。
   旧物性一式は `condN2LatentLowT: 0` で回帰用に保持。
3. **slip 境界の二相整合**: `slip` の ghost は $T=p/(\rho R)$, $\rho e=p/(\gamma-1)+\rho e_k$ で再構成しており、二相 ($p=\rho R_{eff}T$) では内部と矛盾する。
   内部の熱力学状態 ($T$, $c$, $\rho e$; 反射で $|u|$ 不変なので $\rho e$ も不変) を保持する形に改める (codex M4)。Arthur は壁・対称面・面外が全て slip。
4. **検証**: (a) case/34 Arthur N2 の回帰 (旧物性 = 現行) と潜熱/飽和圧整合の効果 — 比較対象は **Arthur 実験記号** (Lin 計算曲線とは別に報告; codex M6)、
   (b) 同ノズルの空気版 (CPG carrier) の onset $(p,T)$ を Daum & Gyarmathy の **$\dot P$ を揃えた理論 onset 線** と比較、(c) node/cell 両離散化 (codex M9)、
   (d) 報告量の時系列定常判定を専用ツールで (codex M7)。

## 2. スコープ

- **やる**:
  - `condVaporMassFraction` (>0 で CPG carrier 形: 凝縮種の質量分率を固定値で与える。既定 −1=従来 pure)。CPG 分岐に carrier 版の EOS
    ($e=(c_{v,air}+gR_w)T-gL$, $p=\rho T(R_{air}-gR_w)$, $p_v=\rho(Y_w-g)R_wT$)、source kernel (carrier=1, $Y_w$ 定数)、SLAU CPG 二相補正
    ($h_{2\phi}=c_{p,air}T-gL$; 単相 $h=\frac{\gamma}{\gamma-1}\frac{p}{\rho}$ からの補正 $g(c_{p,air}\frac{R_w}{R_{air}}T-L)$)、実現可能性 $g\le Y_w$、蒸発の枯渇上限。
  - `condN2LatentLowT` (既定 1): `n2_latent` 70 K 未満線形化 + `n2_psat` の C–C 外挿を新 $L(T)$ の積分で再構成 (接続点で値・微分連続)。0 で旧一式。
  - `slip` ghost の状態保持 (CPG/TP 共通: `T[ig]=T[ic]`, `sonic[ig]=sonic[ic]`, `roe[ig]=roe[ic]` (反射で $e_k$ 不変), `Ht` は $(\rho e+p)/\rho$)。dry 回帰で反復ノイズ以内を確認。
  - 空気の定数は **二成分 N2/O2 = 0.79/0.21 mol** から生成: $M$=28.850 g/mol, $R_{air}$=288.19, $c_{p,air}$=1008.7 (γ 1.4), $Y_{N_2}$=0.7671 (codex m1)。
  - 検証 run (§6)、`onset_analysis.py --series` (時系列定常判定)、node メッシュ、case/34 README、methods、plan §9。
- **やらない**:
  - CPG への二相 frozen 音速 (`condSonicModel`) の拡張: CPG 固有系が config γ・χ=0 を仮定し `gasProperties` が γ を上書きするため、整合した frozen Jacobian にするには
    一般 EOS 経路への切替が要る (codex M3)。**本 plan では pure/CPG の音速は旧式のまま**とし、後続 (§5.1) で「γ_2φ + 実 Ht + 一般 EOS 固有系 + N2/空気での流束 FD 照合」を一括で行う。
  - O2/N2 理想溶液の露点線 (混合液モデル)、O2 と N2 の 2 成分凝縮 (`nCondSpecies 2`)。
  - TP (`thermalMethod 2`) の空気凝縮: NASA-9 は 200 K 未満を線形外挿し、`DEPVAR_TMIN`=50 K で Arthur の 27 K を表現できない。CPG のまま。
  - Iland 補正の空気向け再較正 (データ無し。空気 ≈ N2 の前提で N2 の係数を使う)。
  - `outlet_statPress` の超音速判定 (法線 Mach でなく $|u|$ の Mach) の変更: Arthur の出口は軸方向流で $u_n\approx|u|$。全スナップショットで出口面の $u_n/c$ を監視するに留める。

## 3. 関連 docs と前提

- N2 物性は methods §8 (Lin 2014: Jacobsen $p_{sat}$ (50 K 以上) + C–C 外挿、Nowak $\rho_l$、Stansfield σ、式 26 潜熱、Iland 補正)。液 N2 の実測 $c_{p,l}$ ≈ 2.0 kJ/kg/K (63–77 K)。
  低温側 (30–60 K) の $c_l$ は測定が無く、$c_l$ 一定は**正の熱容量を保証する閉包**で物性精度は保証しない → $c_l$ の感度 (1.5 / 2.0 / 2.5 kJ/kg/K) を単体と run で見る (codex M2)。
- 現行の CPG pure 経路: `cond_vapor_state(carrier=0)` は全気相を凝縮可能と数える ($p_v=p$, $\rho_v=(1-g)\rho$)。空気に流用すると O2 も凝縮してしまう → carrier 形が必要。
- onset データと $\dot P$ の出典はノート §3。**理論 onset 線 (Grossir Fig. 4b, 300 dpi 再読み, ±1 K)**: `case/34.arthur_n2_nozzle/daum_gyarmathy_theory_onset_Pdot20000_n2.csv`
  ($\dot P$=20000 /s: 100 Pa 28.5 K / 1 kPa 38.5 K / 10 kPa 49 K; 1000 /s: 32 / 42 / 52 K; 飽和線 44.5 / 52 / 61.5 K)。初稿の CSV (1 kPa で 46 K) は誤読で、本文の「~40 K」とも
  食い違っていた (codex M5) → 差し替え済み。最小実験 onset 曲線 (Fig. 4a 破線, 1 kPa 43 K) は Grossir の説明では「達成された最大過冷却の限界」で、
  $\dot P$ の大きい conical/wedge 群の下端。
- 参照 run (`run_0008_ref_n2`, 現行物性, `onset_analysis.py`): onset (Δp/p_dry>1 %) x=2.37 in, $p_{dry}$=678 Pa, $T_{dry}$=37.9 K, $\dot P$=1.73×10⁴ /s。
  再読み後の理論線 ($\dot P$=20000, 678 Pa) は 36.8 K → **forge は +1.1 K (理論線よりやや浅い過冷却)**、最小実験曲線 (~41 K) より 3 K 深い。
  ただしこの run は 0/4000/8000 の 3 スナップショットしか無く定常性未確認 (codex M7) → 本 plan の run は 1000 step 毎に保存し `--series` で判定する。
- Arthur 実験記号 (`arthur_fig2_digitized.csv` の `exp` 列): cond/dry = 1.133 (3 in) / 1.250 (4 in) / 1.500 (5 in)。Lin 計算曲線は 1.200 / 1.308 / 1.447 (別物; codex M6)。
- 空気 dry (`run_0013_air_dry`) は密度だけ換算し入口速度を N2 のままにしていたので $M_{in}$=1.068, $P_0$=862 kPa になっていた (codex M8) → 速度も $\sqrt{R_{air}/R_{N_2}}$ 倍にして取り直す (実 config: $u$ 325.12 m/s, $\rho$ 6.137 kg/m³ = M 1.05 スロート条件を空気で再計算)。

## 4. 設計方針

### 4.1 CPG carrier 形 (空気)

$Y_w$ (=$Y_{N_2}$) を config 定数とし、既存 carrier 経路 (TP の H2O–N2) と同じ定義で CPG に載せる:

$$
p_v=\rho\,(Y_w-g)\,R_wT,\qquad p=\rho T\,(R_{air}-gR_w),\qquad e=(c_{v,air}+gR_w)\,T-g\,L(T),\qquad 0\le g\le Y_w
$$

- 温度反転は `cond_T_from_e_cpg` を $R\to R_w$ で使う ($a=c_v+gR_w$; pure は $R_w=R$)。$p$ は $R_{eff}=R_{air}-gR_w$。
  **反転は括弧付き Newton + 二分法退避に改める** (codex v2 M2): 既存の 30 回 Newton は物性クランプ (臨界直下・45 K 床) をまたいで往復すると収束せずに温度を返し
  (例: g=0.75, T=122 K で 99.2 K を返し e が −28 kJ/kg ずれる)、dependentVariables がその温度で `roe` を上書きして保存量を壊す。$e(T,g)$ は $T$ に単調なので
  [T_min, T_max] で括弧を作り、残差 $|e(T)-e_{in}|<10^{-8}|e_{in}|+1$ J/kg を成功条件として返す。**失敗時は保存量を上書きしない** (前ステップの T を保持しフラグを診断に出す)。
  pure N2 も同じ関数を通る (収束している既存 run は結果不変: 単体で旧 Newton 収束点と 1e-10 一致を確認)。
- source kernel: `carrier=1`, `Yw` 定数 (`roY_w=nullptr`)、`Rw=R_{N_2}`、`cprops=N2`。核生成の $R$, $M$, $\rho_v$ は N2 の値 (Iland 補正込み)、成長は Goodheart。成長・蒸発・二温度が呼ぶ $k_{gas}$ は
  `CondSpeciesProps.gasKgasModel` (0=N2, 1=空気 Sutherland $\mu_0$ 1.716e-5 @273 K, C 111 × $c_p$/Pr) で `cond_kgas` がディスパッチし、CPG carrier 空気では kernel が 1 を設定 (codex v2 m2)。
- SLAU CPG 二相の面エンタルピーは**面状態で一貫して構成**する (codex v2 M3; 既存は再構成 $p_f/\rho_f$ の単相項にセル温度 $T_{cell}$ の補正を足す混在で、
  $g$=0.1・$T_{cell}$=40 K・面二相温度 45 K で 1 kJ/kg ずれる): $g_f$ はセル値 (1 次) のまま、
  $T_f=p_f/[\rho_f(R_{air}-g_fR_w)]$, $h_f=c_{p,air}T_f-g_fL(T_f)$ (pure は $R_{air}=R_w$)。TP carrier の面温度修正 (2026-08) と同じ流儀。
- 実現可能性 `cond_realizability_clamp_d`: $g_{max}=Y_w\rho$ (定数)。蒸発・枯渇上限も $Y_w-g$。
- **初版の受付範囲** (codex v2 M4): CPG (`thermalMethod 0`) × `solver: SLAU` × `condModel 0` (N2) × `nCondSpecies 1` × `condEquilibrium 0` × `condKantrowitz ≤ 1`
  × 検証する境界 (`inlet_uniformVelocity` / `outlet_statPress` 超音速 / `slip`)。Roe/KEEP (pure 補正 $g(c_pT-L)$ と $p/[\rho(\gamma-1)]$ の内部エネルギーのまま)、
  `condKantrowitz 2/3` (CPG では `carrierSum=0` で O2 の冷却が入らない) は config で拒否し後続に明記。$Y_w\in(0,1]$、$R_{air}-Y_wR_w>0$ を検査。
  block-DPLUR は CPG 固有系 (config γ, χ=0) のままで**近似 Jacobian** (整合した frozen Jacobian ではない; 収束経路にのみ効く)。
- 空気の定数 (二成分 0.79/0.21 mol): $M_{air}$=28.850e-3, $R_{air}$=288.19, $c_{p,air}$=1008.7, $Y_{N_2}$=0.7671, $Y_{O_2}$=0.2329。config は `physProp: {cp: 1008.7, gamma: 1.4}` + `condensation: {condModel: 0, condVaporMassFraction: 0.7671}`。

### 4.2 N2 潜熱と飽和圧の低温整合 (`condN2LatentLowT: 1`)

- $T\ge T_a$=70 K: 現行多項式。$T<T_a$: $L(T)=L(T_a)+(c_{p,v}-c_l)(T-T_a)$, $c_l$ 既定 2000 J/kg/K → $L'=-961$ J/kg/K (45 K で 233 kJ/kg, 30 K で 248; 旧 188 / 187)。
  接続は **C0 (値のみ連続)**: 多項式側の $L'(70)$=−1072 と線形側 −961 は不連続 (codex v2 M1; 微分連続は要求しない)。片側微分がともに負 = $c_l>0$ を検査。
- 飽和圧: 50 K 以上は Jacobsen のまま。50 K 未満の C–C 外挿を新 $L(T)$ で積分 ($L$ が $T$ に線形なので閉形式:
  $\ln\frac{p_{sat}(T)}{p_{sat}(T_s)}=\frac{1}{R}\Big[\frac{L_a-c'T_a}{T_s}-\frac{L_a-c'T_a}{T}+c'\ln\frac{T}{T_s}\Big]$, $c'=c_{p,v}-c_l$, $T_s$=50 K, $L_a=L(T_a)$)。
  接続 (50 K) も C0 で、$(\ln p_{sat})'$ は 0.308 (低温側) vs 0.279 (Jacobsen 側) と不連続。**38 K の $p_{sat}$ は旧 (一定 $L_{old}(50)$=204 kJ/kg の C–C) の 0.518 倍**
  (閉形式の値; 初稿の 0.59 は新潜熱を一定値で外挿した数)。厳密な C–C 整合は $T<50$ K の外挿域だけの主張。単調性 ($p_{sat}$ が T に単調増) を検査。
- キー (codex v2 m2): `condN2LatentLowT` (1: 新潜熱), `condN2PsatLowT` (1: 新飽和圧; 診断用に 0 で旧 C–C のまま = 「潜熱だけ新」の R2), `condN2LiquidCp` ($c_l$, 既定 2000)。
  `CondSpeciesProps` に載せ、EOS・流束・核生成・成長・蒸発・実現可能性の全経路が同じ構造体を読む。各 run は解決済み設定を README に記録。
- $c_l$ 感度: 1500 / 2000 / 2500 J/kg/K で onset の移動を run で見る (物性精度の不確かさ幅)。

### 4.3 slip の状態保持

CPG/TP とも: ghost は `roe[ig]=roe[ic]` (反射で $|u|$ 不変なので $e_k$ 同一), `T[ig]=T[ic]`, `sonic[ig]=sonic[ic]`, `Ht[ig]=(roe[ig]+P[ig])/ro[ig]`。
bvar は接線速度 (法線成分を除去) なので `roeb = roe[ic] − ½ρU_n²` (codex v2 m1)、`Tsb=T[ic]`。dry では旧式 ($p/(\gamma-1)+\rho e_k$) と丸め差のみ →
凝縮 OFF 回帰 (node/cell) で反復ノイズ以内を確認。二相では ghost が内部の $g$ と整合する。非零法線速度を与えた境界単体試験を追加。

### 4.4 onset の定義と比較

中心線の $p/p_{dry}-1>1\%$ となる最初の点を onset とし、その点の dry 状態 $(p,T)$、$\dot P=-(u/p)dp/dx$ を出す (副: $g>10^{-4}$)。Daum & Gyarmathy との比較は
forge の $\dot P$ (≈1.8×10⁴ /s) に対応する理論線 ($\dot P$=20000) を基準にし、最小実験曲線 (41 K) との差は観測。**±3 K はモデル間比較の基準で、空気モデルの
「実験検証合格」とは分ける** (単一条件で選抜しない; 条件を増やすのは後続)。

## 5. 実装ステップ

1. `condensationProperties_d.cuh`: `CondSpeciesProps.latentLowT`, `n2_latent`/`n2_psat` の新旧切替 (接続点連続)、`air_kgas`。
2. `input/solverConfig.hpp/.cpp`: `condN2LatentLowT` (1), `condVaporMassFraction` (−1)。CPG carrier × `condEquilibrium≠0` / `nCondSpecies≠1` / `thermalMethod 2` は拒否。
3. `dependentVariables_d.cu` CPG 分岐: carrier 版 EOS ($R_w$, $R_{eff}$, $g\le Y_w$)。`condensationEOS_d.cuh`: `cond_T_from_e_cpg` の $R_w$ 引数分離。
4. `condensationSource_d.cu`: CPG carrier (`carrier=1`, `Yw` 定数)。`condensationTransport_d.cu`: 実現可能性 $g_{max}=Y_w\rho$。
5. `convectiveFlux_slau_d.inc.cuh` CPG 二相補正の carrier 形。`convectiveFlux_d.cu` の `CondArgs` に $Y_w$, $R_w$。
6. `boundaryCond_d.cu` `slip_d`: 状態保持。
7. `tests/unit/test_cond_air.cpp`: (a) 新 $L$/$p_{sat}$ の接続連続 (70 K, 50 K で値・微分)、$c_l>0$ 全域 (25–125 K)、(b) CPG carrier EOS 往復 ($e\to T\to e$ 1e-10, $g$ 0–0.99·$Y_w$, T 25–125 K)、
   (c) pure ($Y_w$=1, $R_w=R$) が既存 pure 式と一致、(d) SLAU 補正項が $h_{2\phi}$ と一致、(e) 旧一式 (`latentLowT 0`) が旧関数値と bitwise。
8. `case/34.arthur_n2_nozzle/onset_analysis.py --series` (codex v2 M5 で改修): **壁圧比は壁セル列** (中心線ではなく) から、dry/cond は config (`condensation`) で決め、
   判定対象に onset $x$/$T$/$p$・$\dot P$・壁 $p/p_{dry}$ (3/4/5 in)・$g_{exit}$・**出口境界隣接セルの $u_n/c$ 最小 (全スナップショットで >1 を独立の合否)** を入れる。
   閾値 (plan で固定; onset は x ビン ≈0.013 in の離散位置なので 1 ビン分の交差ゆらぎ Δp/p≈Ṗ Δx/u≈0.9 %, ΔT≈0.25 K を許容): onset 0.02 in / $T_{on}$ 0.4 K /
   $p_{on}$ 2 % / $\dot P$ 10 % / 壁圧比 0.2 % / $g_{exit}$ 0.5 %。末尾窓 4 枚 (最低 3 枚)。
   `check_quasisteady.py` への統合は後続 (dry 参照場を要する case 固有量のため; 1000 step 間隔の保存は「保存時刻での」判定であることを明記)。
9. node メッシュ (`mesh/arthur_nozzle.geo` → node 変換) + `check_mesh_quality`。
10. run (§6)、README、methods、plan §9。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 |
| --- | --- | --- |
| 1 | ~~codex plan レビュー 2 回目~~ | 決着 (2026-09-12, §6.1): GO-with-changes M5/m2 全件採用 (v3) |
| 2 | ~~実装 (§5 1–7)~~ | 決着 (2026-09-13, §9): 実装・full rebuild・単体 `test_cond_air` ALL PASS |
| 3 | ~~ツール・メッシュ (§5 8–9)~~ | 決着 (2026-09-13): `--series` (壁セル列・config dry・出口 u_n/c)、node 平面メッシュ QC PASS |
| 4 | ~~回帰 (§6 R)~~ | 決着 (2026-09-13, §9): slip 変更は dry で cell 1e-3 / node 9e-6 (ノイズ水準)、旧物性再現 onset 同一 |
| 5 | ~~空気 (§6 E)~~ | 決着 (2026-09-13, §9): 空気 onset は理論線 +1.7 K、N2 との差 0.7 K |
| 6 | codex result レビュー | 1 回目 NO-GO (§6.1) → 全件採用。2 回目 NO-GO (2026-09-13, M3/m2, §6.1) → 全件採用・反映済 (§9 2026-09-13 ②)。**3 回目待ち** → GO で `status: done` |
| 6b | 境界試験を実装経由で・流束収支 (result ① 要求, 未実施) | `test_cond_air` (e) の slip 試験はテスト内の代数式比較のまま。slip ghost/bvar を実装 (boundaryCond_d の kernel) 経由で通す試験と、node/cell の質量・エネルギー流束収支 (入口−出口, 凝縮 run で h0 保存) を追加する (codex result ② m5) |
| 6c | `check_quasisteady.py` 統合 (後続) | onset/壁圧比の判定を `--quantity` に追加 (dry 参照場を要する case 固有量; 保存時刻での判定であることを明記)。現状は `onset_analysis.py --series` が代替 (codex result ② m5) |
| 7 | (後続) CPG 二相音速 | γ_2φ + 実 Ht + 一般 EOS 固有系を CPG 二相へ、N2/空気で流束 FD 照合 (double/float32) |
| 8 | (後続) 混合液モデル | O2/N2 理想溶液露点線、2 成分凝縮 |
| 9 | (後続) 条件の拡張 | Longshot 級 (M 10–14, $\dot P$ 小) と Daum の複数条件で理論線・実験点との比較を増やす |
| 10 | (後続) Arthur 3–4 in の ~9 % 過大の**原因切り分け** | 今回の物性変更 (潜熱/飽和圧の低温整合, c_l ±500) では偏差を解消しなかった (物性一般の除外までは立証していない; codex result ② m5)。候補: 核生成 $J$ (CNT×Iland)、成長 $\dot r$ (Goodheart, $\alpha$)、壁圧の抽出位置 (壁セル列 vs 実験の壁静圧孔)、Arthur 記号の読み取り。切り分けは $J$/$\dot r$ 各 ×0.5/×2 の感度と Fig.2 の再デジタイズから (codex 2026-09-13 m1: 「レート側」と断定しない) |

後続 (#6b, #6c, #7–#10) の正本は [condensation-followups.md](condensation-followups.md) §5.1 (codex result ② の要求)。

## 6. 検証

- **単体**: §5 (7)。
- **検証ケース**: `case/34.arthur_n2_nozzle/` (dry 収束場から restart, Euler, `convMethod 1`, cfl_pseudo 1, **12000 step, 1000 step 毎保存**)。全 run で NaN 0、
  `check_convergence.py` VERDICT (2 次のリミットサイクル plateau は既知) と `onset_analysis.py --series` VERDICT を貼る。

  | # | run | 離散化 | 内容 |
  | --- | --- | --- | --- |
  | R0 | `run_0014_n2_ref_cell` / `run_0020_n2_ref_node` | cell / node | 現行物性 (`condN2LatentLowT 0`) の N2 参照 (run_0006 プロトコル; node は新メッシュ, dry も取る) |
  | R1 | `run_0015_dry_slip_cell` / `run_0021_dry_slip_node` | cell / node | 凝縮 OFF、slip 変更前 (branch 先端バイナリ) vs 後 → 反復ノイズ以内 |
  | R2 | `run_0016_n2_latent_only` | cell | 潜熱だけ新・飽和圧旧 (診断) |
  | R3 | `run_0017_n2_new` / `run_0022_n2_new_node` | cell / node | 潜熱+飽和圧新 (最終形) |
  | R4 | `run_0018_n2_cl15` / `run_0019_n2_cl25` | cell | $c_l$ 1.5 / 2.5 kJ/kg/K 感度 |
  | E2 | `run_0023_air_dry` (cell) / `run_0025_air_dry_node` | cell / node | 空気 dry (入口 $u$ 325.12 m/s, $\rho$ 6.137 kg/m³ = `bcondConfig.yaml` 実値, $P_0$=844 kPa, $T_0$=290 K) |
  | E1 | `run_0024_air_cpgcarrier` / `run_0026_air_cpgcarrier_node` | cell / node | 空気 CPG carrier 形 (新物性) |
  | E1' / R3' | `run_0027_air_cpgcarrier_v2` (+反復 `run_0030`) / `run_0028_air_cpgcarrier_node_v2`, `run_0029_n2_new_v2` / `run_0032_n2_new_node_v2` | cell / node | result レビュー①反映後バイナリでの再取得 + 同一バイナリ反復 (ノイズ床) |
  | N1 | `run_0031_dry_slip_cell_oldbin_rep` | cell | 旧バイナリ dry の同一バイナリ反復 (R1 判定のノイズ床) |
  | N2 | `run_0033_air_cpgcarrier_v2_rep2` | cell | E1 v2 の 3 回目反復 (床を 3 run の全ペア最大で定義; codex result ② M2) |
  | E1'' | `run_0035_air_cpgcarrier_final` | cell | result ② 反映後 (反転の非有限入力拒否・面 R_eff 床撤去) の最終バイナリ; 標準空気では両変更は不活性 → run_0027 と床内 |
  | C1 | `run_0034_cfg_default_check` | cell (20 step) | 凝縮セクション有効で `condKantrowitz`/`condKantrowitzGammaMode` を省略 → 起動ログの実効値 0/0 を確認 (読込試験) |

- **判定基準 (合否)**:
  1. 単体 PASS。
  0. **場差の判定方法** (codex result ② M2): 同一バイナリ・同一 config の反復を 3 run 以上取り、変数別に全ペアの max|Δ|/max|ref| の最大をノイズ床 JSON にする (`diff_res.py --dump` → 変数別 max)。判定は `diff_res.py REF NEW --tolfile 床.json --factor 2` の **exit code** (表示丸めで判断しない)。床に無い変数は判定しない。
  2. R1: slip 変更で凝縮 OFF の場差が同一バイナリ反復ノイズ (case/34 は ~1e-3 [README]) 以内。R0 が旧結果 (run_0006: cond/dry 1.21/1.33/1.43 vs Lin 曲線) を再現。
  3. 全 run NaN 0、`--series` STEADY (未達なら延長; plateau は「未収束の準定常比較」と明記)。出口面 $u_n/c>1$ を全スナップショットで確認。
  4. R3/E1: onset $T$ が $\dot P$ 対応理論線から **±3 K** (モデル間比較の基準)。R3 と E1 の onset 差 ≲2 K (「空気 ≈ N2」の再現)。
  5. R3 の Arthur 実験記号 (1.133/1.250/1.500) との cond/dry 比の差を報告 (合否にせず観測; 旧 Lin 曲線比も併記)。悪化・改善の方向を §9 に書く。
  6. node/cell の onset 差 ≲1 K・cond/dry 比差 ≲2 % (共有実装の整合)。
- **観測項目**: onset $x$/$T$/$p$/$\dot P$、$g_{exit}$、cond/dry 比、$c_l$ 感度の onset 移動、E1 と R3 の差。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| result | `2026-09-13` | 2 回目 [`notes/reviews/2026-09-13-condensation-air-result-2.md`](../../notes/reviews/2026-09-13-condensation-air-result-2.md) | NO-GO, C0/M3/m2 | **全件採用** (2026-09-13 ②, §9)。M1 (e=±Inf が tol=Inf で反転成功に化ける) → 入口で非有限 e/g/T_guess・非正熱容量を拒否し、成功条件にも T,G の有限性を要求。`test_cond_air` (i) に ±Inf/NaN/到達不能 e の 7 例 (全て ok=false・有限 T)。M2 (U_x 1.684e-4 > 1 ペア床×2=1.665e-4 で diff_res FAIL、「ちょうど ok」は誤記) → 判定未達を §9 に明記し、3 回目反復 run_0033 で床を 3 反復の全ペア最大に再定義 (U_x 床 1.34e-4) → exit 0 で PASS。方法を §6 に明記。M3 (面エンタルピーの R_eff<1 床が受付範囲 (γ1.4, cp 1038.67, Y_w=g=0.999 → R_eff 0.297) で EOS と食い違う) → 床を撤去し正で有限な R_eff をそのまま使う、非正/非有限は乾き面へ退避。`test_cond_air` (j) に R_eff<1 の 6 状態 (h_f=e+p/ρ 2e-16)。m4 (文書) → methods の CPG 音速を $\sqrt{\gamma R_{air}T}$ に、README run_0023 の入口 325.12/6.137、c_l 感度は同符号 (+0.58/−0.72 K)、`condKantrowitz` の省略時既定 0 (Wysłouzil 参照設定は 1 明示) と診断名 `condTheta_0`/`condLim_0` に訂正。m5 (§5.1 の抜け) → 6b 境界試験の実装経由化・流束収支、6c `check_quasisteady` 統合を残作業表に戻し、#10 を「今回の物性変更では解消しなかった」に限定 |
| result | `2026-09-13` | 1 回目 [`notes/reviews/2026-09-13-condensation-air-result.md`](../../notes/reviews/2026-09-13-condensation-air-result.md) | NO-GO, C0/M5/m2 | **全件採用** (2026-09-13, §9)。M1 (T 反転失敗でも T/P/sonic を更新) → 失敗セルは原始量・roe とも前ステップ値を保持し `g_condTinvFail` を集計して警告 (再取得 run は全て 0 件)。M2 (消滅判定が全圧) → `cond_clamp_vapor_pressure` で source と同じ N2 分圧に統一、全圧/分圧で判定が分かれる状態を `test_cond_air` (h) に追加。M3 (境界種別の受付未実装・Yw NaN) → `main.cpp` で inlet_uniformVelocity/outlet_statPress/slip 以外を拒否、`isfinite` 検査。M4 (node 中心線・壁抽出) → 中心線 |y|<1e-9 ノード、壁は x ビン毎の y 最大ノード (全 slip の node メッシュは wall_dist が無い) で再計算 (onset 表は変わらず)。M5 (cell ノイズ根拠・float32 試験) → 同一バイナリ反復 run_0027/0030 (凝縮) と run_0015o/0031 (旧 dry) で床を実測し `--tolfile` 2 倍で判定、SLAU 面エンタルピーを `cond_face_h_cpg` に抽出して float 入力の試験 (g) を追加。m1 (「レート側」断定) → §5.1 #10 を原因切り分けに書き換え。m2 (文書の古い値) → methods「計画中」/0.59、README 旧注記、入口 325.12 m/s・6.137 kg/m³、c_l 符号の条件 (L'>c_pv)、procedures の新キーを同期 |
| plan | `2026-09-12` | 2 回目 (v2) [`notes/reviews/2026-09-12-condensation-air-plan-2.md`](../../notes/reviews/2026-09-12-condensation-air-plan-2.md) | GO-with-changes, C0/M5/m2 | **全件採用 → v3**。M1 (接続の微分連続は成立しない・0.59 倍は誤値) → C0 接続に改め 0.518 倍に訂正、片側微分と単調性を検査 (§4.2)。M2 (Newton が非収束のまま T を返し roe を上書き) → 括弧付き Newton+二分法・成功フラグ・失敗時は保存量不変 (§4.1)。M3 (SLAU の面温度混在) → 面状態で $T_f$, $h_f$ を一貫構成 (§4.1)。M4 (受付範囲) → SLAU/CPG/N2/単一種/Kantrowitz≤1 に限定し他は拒否、近似 Jacobian と明記、$Y_w$ 範囲検査 (§4.1)。M5 (`--series` の壁圧が中心線・dry 推定・出口判定) → 壁セル抽出・config で dry・$u_n/c$ を独立合否・閾値固定 (§5 8)。m1 (slip bvar の運動エネルギー) → $\rho E_b=\rho E_i-\tfrac12\rho U_n^2$ (§4.3)。m2 (R2/R4 の設定経路・kgas 伝播) → `condN2PsatLowT`, `condN2LiquidCp`, `cond_kgas` ディスパッチ (§4.1, §4.2) |
| plan | `2026-09-12` | 1 回目 [`notes/reviews/2026-09-12-condensation-air-plan.md`](../../notes/reviews/2026-09-12-condensation-air-plan.md) | NO-GO, C1/M8/m1 | **全件採用 → v2**。C1 (pure 擬似種は N2 選択凝縮と EOS が両立しない) → CPG carrier 形 ($Y_w$ 定数, 全経路で同じ $g$) に置換、露点線は後続 (§1, §4.1)。M2 (潜熱修正は飽和圧も変える) → 新 $L$ の積分で C–C を再構成、2 段 A/B と $c_l$ 感度 (§4.2)。M3 (CPG `sonic` のみは Jacobian と不整合) → CPG 二相音速を本 plan から外し後続 #7 (§2)。M4 (slip ghost の二相不整合) → 状態保持に修正、出口は $u_n/c$ 監視 (§4.3)。M5 (理論線 CSV と本文の食い違い) → 300 dpi 再読みで CSV 差し替え (1 kPa: 38.5 K)、±3 K はモデル間比較基準に限定 (§3, §4.4)。M6 (Lin 曲線との比較を実験一致と呼んでいた) → 実験記号と直接比較、Lin 曲線は別報告 (§3, §6)。M7 (準定常判定不能) → 1000 step 毎保存 + `onset_analysis.py --series` (§5 8, §6)。M8 (空気 dry の入口速度未換算) → 324.48 m/s に修正して取り直し (§3, §6 E2)。M9 (cell のみ・掃引不足) → node/cell 両方、単体掃引 g≤0.99·Y_w / T 25–125 K、EOS 往復 (§5 7, §6)。m1 (露点 78.8 K・空気定数) → 82.2 K に訂正、二成分定数に統一 (§4.1) |

## 7. 影響範囲

- `condN2LatentLowT` 既定 1 は**凝縮 ON の pure N2 run (case/34)** の 50 K 未満の飽和圧と 70 K 未満の潜熱を変える (旧は 0)。H2O・TP carrier には無関係。
- `slip` の状態保持は全ケースの slip 境界に効く (dry は丸め差のみ; 回帰 R1 で確認)。
- `condVaporMassFraction` は指定しなければ無影響。CPG 二相音速は変えない (pure/CPG は旧式のまま)。

## 8. 完了条件

- [x] methods 更新 (潜熱/飽和圧の整合・空気 CPG carrier 形・slip)
- [x] 実装・単体・run (§6; 2026-09-13)
- [ ] codex plan (2 回目) / result
- [ ] `status: done`、accepted、README 同期

## 9. 変更ログ

- `2026-09-12` — 初稿。ユーザ要望「空気そのものの凝縮の検証 (case/34 + 文献)」。文献調査 (ノート §3) に基づき設計。
- `2026-09-13` — **実装・検証完了** (case/34 README「空気凝縮 (CPG carrier 形) と N2 低温物性の整合」節、図 `compare_air_n2_wall.png`):
  - 実装: `CondSpeciesProps` に `latentLowT/psatLowT/liquidCp/gasKgasModel`、`n2_latent_ex`/`n2_psat_ex` (C0 接続, 閉形式 C–C)、`cond_kgas` ディスパッチ、
    `CondPropOpts` を source/transport/dependentVariables/SLAU (`CondArgs.cprops`) に値渡し、`cond_T_from_e_cpg` を括弧付き Newton+二分法 (成功フラグ, 失敗時 roe 不変)、
    CPG carrier 形 (`condVaporMassFraction`; EOS $p=\rho T(R_{air}-gR_w)$, $e=(c_v+gR_w)T-gL$, $g\le Y_w$, source `carrier=1`, SLAU 面状態 $T_f=p_f/(\rho_fR_{eff})$, $h_f=c_pT_f-g_fL(T_f)$)、
    `slip` ghost の状態保持 (`roe[ig]=roe[ic]`, bvar $\rho E_b=\rho E_i-\tfrac12\rho U_n^2$)、受付範囲の config 検査 (SLAU/CPG/N2/単一種/非平衡/Kantrowitz≤1, $R_{air}-Y_wR_w>0$)。
  - 単体 `tests/unit/test_cond_air.cpp` ALL PASS: L/p_sat の C0 接続と単調性、$c_l>0$ (25–125 K)、38 K の p_sat 比 0.5176、旧一式 bitwise、EOS 往復 420/420 (悪い初期推定込み、codex の反例 g=0.75/122 K 収束)、
    面状態 h_f=e+p/ρ (旧混在形は −1.1 kJ/kg ずれていた)、slip エネルギー、空気定数、kgas 切替。既存 `test_cond_sonic`/`test_cond_evaporation`/`test_cond_equilibrium_eos` も PASS。
  - 回帰: slip 変更は dry で cell ρ 3e-4 / U_y 2e-3 (case/34 の反復ノイズ ~1e-3 水準), node ≤8.6e-6。旧物性 (run_0014) は旧バイナリ run_0008 と onset 同一 (2.367 in, 37.85 K)、場差 ≤1e-3。
  - **結果 (12000 step, 1000 毎保存, 全 run NaN 0・`--series` STEADY・出口 $u_n/c$>1; 残差は既知のリミットサイクル plateau)**:

    | run | 物性 | onset [in] / p [Pa] / T [K] | 理論線 ($\dot P$=20000) との差 | cond/dry @3/4/5 in | g_exit |
    | --- | --- | --- | --- | --- | --- |
    | run_0014 (cell) / run_0020 (node) | N2 旧 | 2.367 / 678 / 37.85 ; 2.336 / 680 / 37.89 | +1.0 K | 1.241/1.357/1.451 ; 1.245/1.364/1.456 | 0.081 |
    | run_0016 (cell) | 潜熱のみ新 | 2.367 / 678 / 37.85 | +1.0 K | 1.158/1.345/1.446 | 0.059 |
    | run_0017 (cell) / run_0022 (node) | N2 新 (c_l 2000) | 2.112 / 800 / 39.67 ; 2.102 / 793 / 39.59 | **+2.1 K** | 1.242/1.371/1.476 ; 1.248/1.379/1.481 | 0.065 |
    | run_0018 / run_0019 (cell) | c_l 1500 / 2500 | 38.95 / 40.25 K | +1.7 / +2.6 K | 1.2375 / 1.246 @3 in | 0.068 / 0.062 |
    | run_0024 (cell) / run_0026 (node) | 空気 CPG carrier | 2.207 / 750 / 38.94 ; 2.197 / 744 / 38.88 | **+1.7 K** | 1.214/1.351/1.455 ; 1.218/1.357/1.459 | 0.059 |

    合否: 単体 PASS ✓; R1 ノイズ水準 ✓; NaN 0・STEADY・出口超音速 ✓; **N2 新/空気とも $\dot P$ 対応理論線から ±3 K 以内 (+2.1 / +1.7 K)** ✓; **N2 と空気の onset 差 0.7 K (≲2 K)** ✓
    (空気は $S=y_{N_2}p/p_{sat}$ で S が小さい分だけ遅い = Daum & Gyarmathy の「空気 ≈ N2」); node/cell: onset ΔT ≤0.08 K, Δx ≤0.03 in, cond/dry ≤0.6 % ✓。
    Arthur 実験記号 (1.133/1.250/1.500) との差: 旧 +9.5/+8.6/−3.2 %、新 +9.6/+9.7/−1.6 %、空気 +7.1/+8.1/−3.0 % (観測; Lin 計算曲線との差は旧 +3.4/+3.8/+0.3 %)。
    **潜熱整合の効果**: 潜熱だけ新では onset 不変で g_exit 0.081→0.059 (放出熱 +25 % で凝縮量減; 3 in の cond/dry 1.241→1.158 と実験 1.133 に近づく)、飽和圧も新にすると
    S が 1/0.52 倍で onset が 0.25 in 上流 (+1.8 K) に移り 3 in は 1.242 に戻る。$c_l$ 1500 / 2000 / 2500 J/kg/K で onset 38.95 / 39.67 / 40.25 K (**同符号**: $c_l$ +500 で +0.58 K, −500 で −0.72 K; $c_l$ 大 → 低温側の $L$ が大きく $p_{sat}$ が下がり S 増)。3–4 in の実験に対する ~9 % 過大は物性修正前後で同程度で、
    N2 モデル (CNT×Iland, Goodheart) の較正問題として残る (§5.1 後続)。
- `2026-09-13` ② — codex result 2 回目 NO-GO (M3/m2, §6.1) を全件採用。コード: `cond_T_from_e_cpg` の非有限入力拒否と有限性を含む成功条件、`cond_face_h_cpg` の R_eff 床撤去 (正で有限ならそのまま、非正/非有限は乾き面)、起動ログに凝縮キーの実効値を出力。単体 `test_cond_air` (i)(j) 追加 ALL PASS (test_kwc/sonic/evap/eq も PASS)。**M2 の訂正**: E1 の場差判定は 1 ペア床で FAIL していた → 3 反復 (run_0027/0030/0033) の全ペア最大を床にして PASS (方法は §6)。最終バイナリ `run_0035_air_cpgcarrier_final` は run_0027 と 2 倍床内 (exit 0)、onset 表同一 (2.207 in / 38.94 K)、警告 0 件。読込試験 `run_0034_cfg_default_check`: 凝縮有効・キー省略で `condKantrowitz=0 condKantrowitzGammaMode=0` をログで確認。文書同期 (m4) と §5.1 6b/6c/#10 (m5)。
- `2026-09-13` — codex result 1 回目 NO-GO (M5/m2, §6.1) を全件採用。コード: T 反転失敗セルの原始量凍結 + `g_condTinvFail` 警告、消滅判定の N2 分圧化 (`cond_clamp_vapor_pressure`)、
  `cond_face_h_cpg` 抽出、境界種別・Yw 有限性の受付検査。ツール: node の中心線/壁抽出 (y=0 ノード / x ビン毎 y 最大)。単体 `test_cond_air` に (g) float32 面エンタルピー・(h) 分圧判定を追加 ALL PASS。
  **再取得 (最終バイナリ)**: E1 cell `run_0027_air_cpgcarrier_v2` は onset 2.207 in / 38.94 K・cond/dry 1.2134/1.3514/1.4547 で run_0024 と同一、node `run_0028` も run_0026 と同一 (場差 ≤1e-5);
  R3 `run_0029_n2_new_v2` は onset 2.131 in / 39.52 K (run_0017 2.112 / 39.67; 1 % 閾値交差が 1 セル動く = cell の onset ノイズ ±0.02 in / ±0.15 K)、壁比・g_exit は同一。
  **ノイズ床の実測** (codex M5): 同一バイナリ反復 `run_0027` vs `run_0030` (凝縮空気) ρ 4.2e-4 / U_y 2.0e-3 / g 7.6e-4 / Q0 9e-3 (`noise_cell_air_12000.json`)、`run_0015o` vs `run_0031` (旧バイナリ dry)
  ρ 2.7e-4 / U_y 2.8e-3 (`noise_cell_dry_oldbin_12000.json`)。この 2 倍に対し slip 変更 (dry) は全変数 ok。E1 (run_0024 vs run_0027) は 1 ペア床では U_x 1.684e-4 > 2×8.33e-5=1.665e-4 で `diff_res.py` **FAIL (exit 1)** だった (この時点の「ちょうど ok」は表示丸めの誤記; codex result ② M2)。
  → 3 回目反復 run_0033 を追加し、床を **3 反復の全ペア最大 (変数別)** に定義し直した (U_x 8.3e-5 / 1.34e-4 / 1.26e-4 → 床 1.34e-4): 2 倍床 2.7e-4 に対し run_0024 vs run_0027 は全変数 ok (exit 0)、run_0017 vs run_0029 も exit 0 (§9 2026-09-13 ②)。
  T 反転失敗の警告は全 run 0 件。§3 の理論線比較・§9 の結論 (+2.1 / +1.7 K, N2–空気 0.7 K) は変わらない。
- `2026-09-12` — codex plan 2 回目 GO-with-changes (M5/m2) を全件採用して v3 (§6.1): C0 接続と 0.518 倍、括弧付き反転、SLAU 面状態の一貫構成、受付範囲の限定、
  `--series` 改修、slip bvar、診断キー。`status: in_progress`。node 平面メッシュ (`mesh_node/`, 24000 双対 CV, cell 版 QC PASS AR 21.9 / skew 0.07) と
  node dry 参照 `run_0021_dry_slip_node` (PASS, 出口 M 6.93) と空気 dry 修正版 `run_0023_air_dry` (入口 325.12 m/s, ρ 6.137) を取得済み。
- `2026-09-12` — codex plan 1 回目 NO-GO (C1/M8/m1) を全件採用して v2 (§6.1)。空気を pure 擬似種 → CPG carrier 形 (N2 選択凝縮 + O2 キャリア) に変更、
  潜熱と飽和圧の低温外挿を連動、slip の状態保持、CPG 二相音速を後続へ、理論線 CSV を再読みで差し替え (forge N2 参照は $\dot P$=20000 線に対し +1.1 K)、
  Arthur は実験記号と直接比較、空気 dry の入口速度換算、node/cell 両方、`--series` 判定ツール。2 回目レビューは codex 使用上限 (2026-09-13 01:24 解除) 待ち。

## 10. 未確定事項

- ~~空気の核生成に N2 の Iland 補正をそのまま使うか (Daum & Gyarmathy の「空気≈N2」に依拠)。~~ 決着 (2026-09-13, §4.1/§9): 初版は N2 の Iland 補正をそのまま使う (空気向け較正データが無く、E1 の onset が理論線 +1.7 K・N2 との差 0.7 K で「空気 ≈ N2」と整合)。空気固有の再較正は Longshot 級条件の比較 (#9) で判断。
- ~~過飽和度の基準 (N2 分圧 vs 露点線)~~ 決着 (2026-09-12, codex C1): 初版は N2 分圧基準 (CPG carrier 形)。露点線 (混合液) は後続 #8。
- 低温 (30–60 K) の液 N2 比熱 $c_l$: 実測が無く 2.0 kJ/kg/K は閉包。感度 R4 で幅を出す。
