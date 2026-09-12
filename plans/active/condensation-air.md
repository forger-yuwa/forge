# 空気 (N2/O2) 自体の凝縮: N2 潜熱フィットの整合修正・空気擬似種・Arthur / Daum–Gyarmathy 検証

## メタ

- **area**: `condensation`
- **status**: `draft`  <!-- plan 段 codex レビュー待ち -->
- **related_docs**:
  - `methods/condensation.md` (§8 N2 物性「潜熱」節・新設「空気擬似種」節・§5 二相音速の pure 適用)
  - [notes/investigations/condensation-carrier-kantrowitz-air-survey.md](../../notes/investigations/condensation-carrier-kantrowitz-air-survey.md) §3 (onset データ、O2 物性、forge との差分)
- **related_plans**:
  - [condensation-nonequilibrium.md](../accepted/condensation-nonequilibrium.md) (N2 モデルの導入元、case/34 Fig.2 検証)
  - [condensation-kantrowitz-gamma-twophase-sonic.md](condensation-kantrowitz-gamma-twophase-sonic.md) (二相音速; pure N2 は `n2_latent` の $L'>0$ で対象外にした → 本 plan で解消)
- **created**: `2026-09-12`
- **owner**: `CFD Dev`
- **branch**: `feature/condensation-air` (worktree `../forge-cond`)

## 1. 目的

極超音速ノズルの強膨張で**試験ガスの空気そのもの**が凝縮する解析を forge で行えるようにし、検証する。

1. **N2 潜熱フィットの熱力学整合**: `n2_latent` (Lin 2014 式 26, 4 次多項式) は 60 K 未満で $L'=dL/dT>0$ となり液の比熱 $c_l=c_{p,v}-L'$ が負になる
   (45 K で −3100 J/kg/K)。Arthur (case/34) の最低温 28–35 K はこの域。低温側を $c_l$ 一定の線形外挿に置き換え、$c_{v,2\phi}>0$ を保証する。
2. **二相 frozen 音速を pure (CPG) にも適用**: 1. により `cond_twophase_sonic` の前提 ($c_{v,2\phi}>0$) が満たせるので、CPG 分岐にも実装し
   `resolveCondSonicModel` の自動 ON を pure N2/空気 + 検証済み境界に広げる。
3. **空気擬似種 `condModel 2`**: Daum & Gyarmathy (1968) の「低圧では空気は純 N2 として振る舞う」に基づき、**核生成・成長は N2 の式**を使い、
   **飽和線は O2/N2 理想溶液の露点線** (Hansen & Nothwang 1952 式 A1–A2)、$\rho_l$・σ・$L$・$R$・$M$ は空気 (N2 0.79 / O2 0.21 mol) の混合値にした擬似種。
4. **検証**: (a) case/34 Arthur N2 (Fig. 2 壁圧 cond/dry 比) の回帰 (潜熱修正前後・二相音速 on/off)、(b) 同ノズルの空気版で onset (p, T) を
   Daum & Gyarmathy の最小 onset 曲線 (ノート §3.2) と比較、(c) 膨張率 $\dot P$ を出して Fig. 4b の理論 onset 線と対応づける。

## 2. スコープ

- **やる**:
  - `n2_latent` の低温整合 (`condN2LatentLowT` 1=新 (既定) / 0=旧多項式のまま: A/B 用)。
  - `cond_twophase_sonic` の CPG 分岐適用 + resolver の pure 対応 (`condModel` 0/2、CPG、検証済み境界 → 自動 1)。単体 sweep で N2/空気の $c_{v,2\phi}>0$ を全有効域で確認。
  - `condProps_AIR()` (M 28.9647e-3, R 287.05, cv 717.6, cp 1004.6, Tc 132.6 (擬似), Tt 63.15) と物性関数 `air_psat` (露点線), `air_rho_cond`, `air_sigma`, `air_latent`, `air_kgas`。
    O2 側: Antoine (NIST, 54–154 K) + 54 K 未満 C–C 外挿、$\rho_l$ 線形 (1141 @90.19 K, 1306 @54.36 K)、σ $=\sigma_0(1-T/154.58)^{1.25}$ ($\sigma_0$ を 13.2 mN/m @90.19 K に合わせる)、
    $L$ = 213 kJ/kg @90.19 K + $(c_{p,v}-c_l)(T-90.19)$ ($c_l$ 1.7 kJ/kg/K)。混合: 露点線は $p_{dew}=1/(y_{N_2}/p^{sat}_{N_2}+y_{O_2}/p^{sat}_{O_2})$、
    液相組成 $x_{O_2}=y_{O_2}p_{dew}/p^{sat}_{O_2}$ で $\rho_l$・σ・$L$ を $x$ (モル) / 質量分率で混合。
  - 核生成: CNT × Iland 補正 (N2 と同じ; 空気≈N2 の前提)。成長: Goodheart、$k_{gas}$ は空気 Sutherland。
  - 検証 run (§6)、case/34 README、methods、plan §9。
- **やらない**:
  - O2 と N2 を別々に凝縮させる 2 成分凝縮 (`nCondSpecies 2`) — 理想溶液の単一擬似種で十分か、まず onset で判断。
  - TP (`thermalMethod 2`) の pure 凝縮 — NASA-9 が 200 K 未満で無効なので CPG のまま。
  - Iland 補正の空気向け再較正 (データ無し)。Daum & Gyarmathy の $\dot P$ 相関の再現は「対応づけ」まで。

## 3. 関連 docs と前提

- N2 物性は methods §8 (Lin 2014: Jacobsen $p_{sat}$, Nowak $\rho_l$, Stansfield σ, 式 26 潜熱, Iland 補正)。液 N2 の実測 $c_{p,l}$ ≈ 2.0 kJ/kg/K (63–77 K)。
- onset データと O2 物性の出典はノート §3。Arthur 条件 ($P_0$ 844 kPa, $T_0$ 290 K) の onset 近傍は $p\approx1$ kPa, $T\approx43$ K で
  Daum & Gyarmathy の 1 kPa 点 (onset 43 K, 過冷却 ~8 K) に相当。
- Arthur nozzle の $\dot P=-(1/p)dp/dt=-(u/p)\,dp/dx$ は dry 場から算出。**実測 (run_0008_ref_n2, `onset_analysis.py`, 2026-09-12): onset 付近で
  $\dot P\approx1.7$–$1.8\times10^4$ /s** (Fig. 4b の最上位クラス 16000–24000 /s、conical 群)。参照 run の onset は Δp/p_dry>1 % で x=2.37 in,
  $p_{dry}$=678 Pa, $T_{dry}$=37.9 K (g>1e-4 では 2.27 in, 722 Pa, 38.5 K)。Daum & Gyarmathy **最小**実験 onset 曲線は 700 Pa で 41.3 K なので
  forge は 3 K 深く過冷却するが、Fig. 4b の $\dot P$=20000 /s 理論 onset 線 (700 Pa で ~40 K) には近い。**比較は $\dot P$ を揃えて行う**
  (最小曲線は $\dot P$ の小さい風洞も含む包絡なので、$\dot P$=2×10⁴ の Arthur には理論線の方が対応する)。

## 4. 設計方針

### 4.1 `n2_latent` の低温整合

$T\ge T_a$ (=70 K) は現行多項式。$T<T_a$ は $L(T)=L(T_a)+(c_{p,v}-c_l)(T-T_a)$、$c_l$=2000 J/kg/K → $L'=-961$ J/kg/K で単調 (45 K で 233 kJ/kg; 旧 188)。
45 K の物性クランプ (`cond_clamp_Tprop`) は他物性に残すが $L$ は線形を続ける。**影響**: Arthur の 30–45 K 域で潜熱 +5〜+25 %
→ Fig. 2 の cond/dry 比が動く可能性 → A/B で定量化し、実験 (Arthur 壁圧) との一致が悪化しないか見る。

### 4.2 pure への二相音速

CPG 二相分岐: $R_{eff}=(1-g)R$, $c_p^{全蒸気}=c_p$、$L'$ は 4.1 の式。`gam_array` は CPG では未使用 (Jacobian は config γ) なので `sonic` のみ。
resolver: `condModel∈{0,2} && thermalMethod==0 && condGasSpecies<0 && condEquilibrium==0 && 検証済み境界` → 1。Arthur の bcond は
`inlet_uniformVelocity` (超音速入口: ρ,U,Ps 全量固定、g=0 なので二相の影響なし) / `outlet_statPress` (出口は M≈6.9 の超音速で全量外挿、
亜音速セルは無い) / `slip`。前 2 種は carrier 用の verified リストに無いので、**pure CPG 用のリストに `inlet_uniformVelocity` を追加**し、
`outlet_statPress` は「出口が全面超音速なら ghost は内部コピー」という条件付きで許可する (超音速判定は実行時にしかできないので、
resolver は warn を出したうえで 1 にし、run 後に出口面の亜音速セル数 0 を確認する手順を §6 に置く)。

### 4.3 空気擬似種

`CondSpeciesProps.model = COND_MODEL_AIR (2)`。ディスパッチ関数 (`cond_psat/rho_cond/latent/sigma/kgas`) に分岐追加。核生成の Iland 補正は N2 と同じ係数。
$k_{gas}$: 空気 Sutherland ($\mu_0$ 1.716e-5, 273 K, C 111) × $c_p$/Pr。

**飽和線の定義 (2026-09-12 追記, 起票後の数値確認)**: 理想溶液露点線 (N2 Jacobsen + O2 Antoine) を Python で組むと 1 atm の露点 82.2 K (実在空気 81.7 K と整合) だが、
**全圧域で N2 飽和線より 4.3–4.8 K 高温側**になる。一方 Daum & Gyarmathy (Grossir Fig. 4a) の「air saturation line」は N2 線より ~1 K 高温側で、
onset データも air と N2 が重なる (「低圧では空気は純 N2 として振る舞う」)。露点線 (O2 濃厚な最初の液滴) を過飽和度の基準にすると、N2 核生成で空気の
onset が N2 より ~4 K 早くなり実験と矛盾する。したがって過飽和度の定義を 2 つ用意して比較する:
- `condAirSatModel: 0` (既定候補): **N2 分圧基準** $S=y_{N_2}p/p^{sat}_{N_2}(T)$ (O2 は核生成に関与しないキャリア扱い、凝縮相物性は N2)。Daum & Gyarmathy の経験則に直接対応。
- `condAirSatModel: 1`: 理想溶液露点線 $S=p/p_{dew}(T)$、凝縮相は $x_{O_2}$ 混合物性 (Hansen & Nothwang の枠組み)。
どちらが onset データに合うかを §6 E1 で判定し、既定を決める (§10)。

### 4.4 onset の定義 (Daum & Gyarmathy との比較)

実験は静圧の等エントロピーからの逸脱で onset を検出するので、forge でも **中心線の $p/p_{dry}-1$ が 1 % を超える最初の点**を onset とし、
その点の $(p_{dry}, T_{dry})$ を Fig. 4a の onset 曲線と比べる。副として $g>10^{-4}$ の点も出す。

## 5. 実装ステップ

1. `condensationProperties_d.cuh`: `n2_latent` 低温線形化 (キー連動)、`condProps_AIR`、`air_*` 物性、ディスパッチ、`cond_kgas` ディスパッチ (N2/空気)。
2. `input/solverConfig.hpp/.cpp`: `condN2LatentLowT` (既定 1)、`condModel` 2 の許容。`condSonicResolve.hpp`: pure CPG 条件。
3. `dependentVariables_d.cu`: CPG 二相分岐に `cond_twophase_sonic`。
4. `tests/unit/test_cond_air.cpp`: 露点線 (N2/O2 純成分極限、1 atm で ~78.8 K)、物性の連続性、$c_{v,2\phi}>0$ sweep (N2/空気, 30–120 K, g≤0.3)、旧多項式との切替。
5. case/34 run (§6)、`postproc.py` に onset 抽出と $\dot P$ を追加、README、methods、plan §9。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 |
| --- | --- | --- |
| 1 | codex plan レビュー | 1 回目は 2026-09-12 に codex 使用上限で中断 (log のみ; 17:26 以降に再実行)。**再実行して GO を得るまで実装に入らない** |
| 2 | 実装 (§5 1–4) | worktree、full rebuild |
| 3 | case/34 回帰 (§6 R) | 現行バイナリ参照 → 潜熱新旧 A/B → 二相音速 on/off |
| 4 | 空気 run と onset 比較 (§6 E) | Arthur 条件の空気、$\dot P$ |
| 5 | codex result レビュー | → `status: done` |
| 6 | (後続) 2 成分凝縮 (N2/O2 別) | 露点線近似で不足なら |
| 7 | (後続) Longshot 級 (M 10–14, $\dot P$ 小) 条件 | 輪郭ノズルで過冷却が小さい側の検証 |

## 6. 検証

- **単体**: §5 (4)。
- **検証ケース**: `case/34.arthur_n2_nozzle/` (run_0006 のプロトコル: 収束 dry `restart_dry.h5` から restart, CPG N2, Euler, `convMethod 1`, cfl_pseudo 1, 8000 step)。

  | # | run | 内容 |
  | --- | --- | --- |
  | R0 | `run_0008_ref_n2` | 現行バイナリ (ブランチ先端, 変更前) の run_0006 再現 |
  | R1 | `run_0009_n2_latent_old` | 新バイナリ + `condN2LatentLowT: 0` + `condSonicModel: 0` → R0 と一致 (回帰) |
  | R2 | `run_0010_n2_latent_new` | 潜熱整合 (既定) + `condSonicModel: 0` → 潜熱修正の単独効果 (Fig. 2 比) |
  | R3 | `run_0011_n2_new_sonic` | 潜熱整合 + 二相音速 (既定) → 音速の単独効果 |
  | E1 | `run_0012_air_ref` | `condModel 2` (空気), 同ノズル・同 $P_0,T_0$ (dry は N2 と γ 同一なので `restart_dry.h5` を流用可: R=287 との差は密度スケールのみ → 空気用 dry を別途取る) |
  | E2 | `run_0013_air_dry` | 空気 dry (onset 判定と $\dot P$ の基準) |

- **判定基準 (合否)**:
  1. 単体 PASS。R1 ≡ R0 (cell run なので反復ノイズ ~1e-3 [case/34 README] 以内)。
  2. 全 run NaN 0、`check_convergence` VERDICT を貼る (case/34 は 2 次でリミットサイクル ~3e-4 プラトーが既知 → 準定常判定 `--series` 相当を壁圧比で行い STEADY を要求)。
  3. R2: Fig. 2 の cond/dry 比 (3/4/5 in) が実験 ±5 % 以内を維持 (現状 1–2 %)。悪化したら潜熱修正の妥当性を再検討 (合否は「±5 % 以内」)。
  4. R3: 音速 on/off の場差を報告 (観測)。$c$ の変化は g~0.1 で −5 % 級の見込み。
  5. E1/R2: onset $(p,T)$ を Daum & Gyarmathy の Fig. 4 と **$\dot P$ を揃えて**比較する: forge の $\dot P$ (≈1.8×10⁴ /s) に対応する理論 onset 線
     ($\dot P$=20000 /s: 700 Pa で ~40 K, Grossir Fig. 4b 読み) から **±3 K** を合否帯とし、最小実験曲線 (41.3 K @700 Pa) との差は観測に置く
     (参照 N2 run は理論線より ~1.5 K 深く、最小曲線より 3 K 深い)。N2 (R2) と空気 (E1) の onset 差 ≲2 K (Daum & Gyarmathy の「air ≈ N2」) は観測。
- **観測項目**: onset 位置 [in]、$T_{onset}$、$p_{onset}$、$\dot P$、$g_{exit}$、cond/dry 比、$c$ の変化。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |

## 7. 影響範囲

- `n2_latent` の既定変更は **凝縮 ON の pure N2 run (case/34)** の 60 K 未満の潜熱を変える (旧は `condN2LatentLowT: 0`)。H2O・carrier には無関係。
- pure CPG の二相音速自動 ON は case/34 系のみ。凝縮 OFF・TP carrier はビット不変。
- 新規 `condModel 2` は選ばなければ無影響。

## 8. 完了条件

- [ ] methods 更新 (潜熱・空気擬似種・pure 音速)
- [ ] 実装・単体・run (§6)
- [ ] codex plan / result
- [ ] `status: done`、accepted、README 同期

## 9. 変更ログ

- `2026-09-12` — 初稿。ユーザ要望「空気そのものの凝縮の検証 (case/34 + 文献)」。文献調査 (ノート §3) に基づき設計。

## 10. 未確定事項

- 空気の核生成に N2 の Iland 補正をそのまま使うか (Daum & Gyarmathy の「空気≈N2」に依拠)。
- 露点線 (理想溶液) と実在混合物 (Lemmon 2000 air EOS) の差 (~1 K) を許容するか。
- 過飽和度の基準 (N2 分圧 vs 露点線): §4.3 追記のとおり露点線は N2 線より 4.5 K 高温側で Daum & Gyarmathy の「air ≈ N2」と合わない可能性 → E1 で判定。
