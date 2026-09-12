# carrier 中の非等温核生成補正 (Feder 形) と H2O 表面張力の小半径妥当性

## メタ

- **area**: `condensation`
- **status**: `in_progress`  <!-- plan 段 codex GO-with-changes (M6/m3) 全件採用 → 実装 -->
- **related_docs**:
  - `methods/condensation.md` (§2 核生成「Kantrowitz 非等温補正」・「表面張力の妥当性」節 — 本 plan 起票時に更新)
  - [notes/investigations/condensation-carrier-kantrowitz-air-survey.md](../../notes/investigations/condensation-carrier-kantrowitz-air-survey.md) (文献調査: Feder 一般形、Wysłouzil 条件の見積り、過冷却水 σ の実測状況、Tolman 長)
- **related_plans**:
  - [condensation-kantrowitz-gamma-twophase-sonic.md](condensation-kantrowitz-gamma-twophase-sonic.md) (前段: 純蒸気形の γ_v 修正。本 plan はその carrier 拡張)
  - [condensation-nonequilibrium.md](../accepted/condensation-nonequilibrium.md) (親)
- **created**: `2026-09-12`
- **owner**: `CFD Dev`
- **branch**: `feature/condensation-air` (worktree `../forge-cond`)

## 1. 目的

1. **carrier 中の非等温補正**: 現行 `condKantrowitz: 1` は純蒸気形 (蒸気分子との衝突だけがクラスタの潜熱を持ち去る) で、H2O–N2 では
   N2 分子の衝突による冷却を無視するため $\theta$ を 40 倍過大に見積もる (Wysłouzil 230 K: $\theta$ 167 vs Feder carrier 形 4.1;
   $J/J_{iso}$ 0.006 vs 0.20)。Feder et al. (1966) の一般形 ($b^2$ にキャリア衝突項) を `condKantrowitz: 2/3` として実装し、
   Wysłouzil 2D で等温 (0) / 純蒸気 (1) / carrier (2, 3) の onset・壁圧を比較する。
2. **H2O 表面張力の妥当性 (ユーザ質問「50 nm 級でも成り立つか」)**: 文献調査の結論 (成長液滴 50 nm では平面 σ で十分、臨界核 ~1 nm は
   capillarity 近似の限界、過冷却域の IAPWS 外挿は 242 K まで実測で支持) を methods に記し、σ の不確かさが onset に与える感度を
   `condSigmaScale` (σ 倍率) の A/B で定量化する。Tolman 補正は実装しない (δ の符号・大きさが未確定)。

完了時: `condKantrowitz` 0/1/2/3 が選べ、Wysłouzil 1.0 kPa で 4 者の onset (mm) と壁 p/p₀ 偏差の表、mode 3 基準の σ ±3 % 局所感度が README と plan にある。
**位置づけ: mode 3 = 物理モデル (Feder/Wedekind 式 8–10 に対応)、mode 2 = 表面仕事項を落とした比較用、mode 1 = 旧結果再現用 (Kantrowitz 原形)、mode 0 = 等温**。
グローバル既定値は **0 (等温) のまま**変えない (Wysłouzil 参照 config が 1 を明示しているだけ)。2/3 を推奨値にするかは結果を見て別途決める。

## 2. スコープ

- **やる**:
  - `cond_nucleation` に Feder 形 $\theta=q^2/b^2$ を追加: mode 2 = $q=m_vL-k_BT/2$ (Kantrowitz と同じ $q$)、mode 3 = $q=m_vL-k_BT(\tfrac12+\ln S)$ (表面仕事項込み)。
    衝突項は**種 DB (NASA-9 $c_v(T)$, $M_i$) と種質量分率から種別に集計**する (混合 $c_p$ からの引き算はしない; codex M1)。
  - 本体評価と差分評価 (`src_jac` 用の温度・モーメント摂動、`cond_source_vector` 経由) の**全呼び出しで同じモデル**を使う (codex M2)。
  - pure ($N_c=0$): mode 2 は Feder 純蒸気形 $(b_L-\tfrac12)^2/(\tilde c_{v,v}+\tfrac12)$ (mode 1 と 2 % 差)、mode 3 は表面項が残る。各式の解析値を単体で個別に検査。
  - `condSigmaScale` (既定 1.0): `CondSpeciesProps` に載せ `cond_sigma` (`condensationProperties_d.cuh`) に掛ける倍率。核生成・Kelvin・蒸発に一貫。
    「一定倍率による局所感度」であり、温度依存・曲率依存の誤差モデルや信頼区間ではない (codex M5)。
  - Wysłouzil 2D node (run_0335 プロトコル, IC=run_0213 dry, 48000 step) で `condKantrowitz` 0/1/2/3 と **mode 3 基準**の σ ×0.97/×1.03、cell の mode 1/3 対照。
  - methods §2 更新 (Feder 形の式・見積り・表面張力節)、README run 一覧、本 plan §9。
- **やらない**:
  - pV 仕事 (Wedekind の pressure effect): Wysłouzil で $v_lp_c/(k_BT\ln S)\sim3\times10^{-5}$ で無視 (ノート §1.4)。
  - Tolman 補正 (`σ(r)`) の実装。感度試験で代替。
  - 既定値の変更 (別途決定)。分圧スイープ (0.5 / 0.26 kPa) は dry 場の再取得が要るので後続 (§5.1)。
  - 一般の多成分擬似種 (MIXDRY のような合成種) の厳密扱い: 衝突項は種 DB の各エントリ (MIXDRY は実質 N2 の物性) で集計し、擬似種の内部組成までは分解しない (限界として明記)。
  - N2 pure への carrier 形 (キャリア無しなので mode 1 と同値。自動で同値になる)。

## 3. 関連 docs と前提

- Feder 一般形の導出・記号はノート §1。純蒸気極限で forge の mode 1 と 2 % 差 ($(b_L-\tfrac12)^2$ vs $b_L(b_L-\tfrac12)$) — mode 2 は
  Feder の $(b_L-\tfrac12)^2$ 形を採り、mode 1 は Kantrowitz 原形のまま残す (旧結果の再現性)。多原子・複数キャリアの一般形は Horsch et al.
  (CO2–air 核生成) 式 14 の種別和と同形 (codex 参照)。
- 参照 run (run_0335) で中心線が最初に $g>10^{-3}$ となる節点は **T=213.8 K, ln S=5.03** (codex 実測) — 230 K・ln S=3.4 の手計算点より低温・高過飽和。
  単体掃引はこの域 (T 200–260 K, ln S 2–6) を覆う。
- 表面仕事項: 臨界核で $\gamma\,\partial A/\partial n=2\sigma v_l/r_*=k_BT\ln S$ (Kelvin–Thomson)。$\ln S\approx3.4$ (Wysłouzil onset) で $q$ が 14 % 減。
- 現行の kernel 入力: `cp_cell`/`Rmix_cell` (全蒸気混合)、`roY_w` (総水)、`g`。キャリアの物性はこれらから逆算する (下記)。
- 表面張力: `h2o_sigma` = IAPWS R1-76 形の過冷却外挿。実測 (Hrubý 2014, Vinš 2015/2020) は 241.8 K まで外挿と一致 (ノート §2)。

## 4. 設計方針

### 4.1 Feder carrier 形 (`condKantrowitz: 2, 3`)

$$
\frac{J}{J_{iso}}=\frac{1}{1+\theta},\quad \theta=\frac{q^2}{b^2},\qquad
q=\begin{cases} m_vL(T)-\tfrac12k_BT & (\text{mode 2})\\ m_vL(T)-k_BT\big(\tfrac12+\ln S\big) & (\text{mode 3})\end{cases},\qquad
b^2=k_B^2T^2\Big[\tilde c_{v,v}+\tfrac12+\frac{N_c}{N_v}\sqrt{\frac{m_v}{m_c}}\big(\tilde c_{v,c}+\tfrac12\big)\Big]
$$

- 実装形 (codex M1 の対案; 蒸気枯渇で有限): $\hat q=q/(k_BT)$、$a_v=(Y_w-g)/M_v$ として

  $$
  \theta=\frac{a_v\,\hat q^2}{a_v(\tilde c_{v,v}+\tfrac12)+\sum_{i\ne v}\frac{Y_i}{M_i}\sqrt{\frac{M_v}{M_i}}\big(\tilde c_{v,i}+\tfrac12\big)}
  $$

  $\tilde c_{v,i}=c_{v,i}(T)M_i/R_u$ は **種 DB の NASA-9** から $c_{p,i}(T)-R_i$ で評価 (蒸気 H2O も同じ DB: 230 K で 3.0)。混合 $c_p$ からの引き算は
  純蒸気極限で破綻するので使わない。$Y_w-g\to0$ で $\theta\to0$ (等温) に連続に落ち、$Y_c\to0$ で Feder 純蒸気形に落ちる。
- kernel は `SpeciesThermo* sp, nSpecies, roY[]` を受け取り (dependentVariables と同じ `thermo_species_device_ptr()`)、種別和をその場で計算する。
  CPG (pure N2, `sp==nullptr`) は和が空で純蒸気形。
- $q$ は式のまま $q^2$ を使う ($q<0$ は表面仕事が潜熱を上回る $\ln S>b_L-\tfrac12$ の域で、Feder 形の適用外。**黙って等温に落とさず**、単体で適用域
  $\ln S<b_L-\tfrac12$ を確認し、run では診断 (`condTheta_<s>` 出力) で $q<0$ セルが無いことを見る)。
- 伝播 (codex M2): `cond_nucleation` の新引数は `cond_source_vector` とその全呼び出し元 (本体 + `src_jac` の温度・モーメント摂動) に同じ値で渡す。
  温度摂動時の $\tilde c_{v,i}(T)$ は**摂動前 T で凍結** (摂動幅 0.1 K で $c_v$ の変化は 1e-5、差分係数の一貫性を優先)。
- 診断出力 `condTheta_<s>` (θ) と `condLim_<s>` (ソース律速係数 θ_lim) を追加し、核生成域で律速が結果を支配していないか見る (codex M6)。
- mode 0/1 は既存経路のまま (ビット不変)。

### 4.2 σ 倍率 (`condSigmaScale`)

`CondSpeciesProps.sigmaScale` を `cond_sigma()` (`condensationProperties_d.cuh`) の戻り値に掛ける (核生成・Kelvin 項・蒸発の全経路に一貫して効く)。既定 1.0 で
ビット同一 (乗算 ×1.0 は IEEE で恒等)。感度: $\Delta\ln J\approx-3(\Delta G^*/k_BT)\,\Delta\sigma/\sigma$、$\Delta G^*/k_BT\sim50$–70 → ±3 % で $J$ が 10²–10³ 倍動く。
これは「一定倍率による局所感度」で、σ の温度外挿誤差 (Vinš 2020 は −20 °C 未満で IAPWS 外挿からの偏差を報告) や曲率依存 (Tolman) の代替ではない。
Tolman 補正を今回入れないのは、対象温度・臨界核サイズで採用する曲率モデルと係数の検証が不足しているため (異常が否定されたからではない; codex M3)。

### 4.3 Wysłouzil 条件の見込み (観測項目)

局所 $\theta$ (230 K, ln S 3.4, コード物性; codex 再計算): mode 1 167.4 / mode 2 4.05 / mode 3 2.98。$J/J_{iso}$: 0.006 / 0.20 / 0.25。
参照 run の実際の onset 状態 (213.8 K, ln S 5.0) では別値になるので単体掃引で覆う。等温 (mode 0) は先行研究で onset 過早 (x≈2.0 cm)、mode 1 は 22.5 mm。
**結合 run の onset 序列は観測項目** (核生成率の局所序列 $J_0\ge J_3\ge J_2\ge J_1$ とは別物: 成長・潜熱・律速が絡む; codex M4)。逆転したら
局所状態・成長・律速係数を調べる。実験の onset 帯 (壁 p/p₀ @21 mm −4.5 %) に近づくかも観測。

## 5. 実装ステップ

1. `cuda_forge/condensationProperties_d.cuh`: `CondSpeciesProps.sigmaScale`、`cond_sigma` に倍率。
2. `cuda_forge/condensationSource_d.cuh`: `cond_nucleation` に mode 2/3 (Feder θ; 種別和は呼び出し側で作った `carrierSum`=$\sum_{i\ne v}\frac{Y_i}{M_i}\sqrt{M_v/M_i}(\tilde c_{v,i}+\tfrac12)$ と $a_v$, $\tilde c_{v,v}$ を受ける)、
   `cond_source_vector` に同じ引数を通す。
3. `cuda_forge/condensationSource_d.cu`: kernel に `sp, nSpecies, roY` を渡し、セルごとに $a_v$・`carrierSum`・$\tilde c_{v,v}$ を種 DB から計算 (T 凍結)。本体と `src_jac` 摂動の全呼び出しに同じ値。
   診断 `condTheta_<s>`, `condLim_<s>` を登録・出力。
4. `input/solverConfig.hpp/.cpp`: `condKantrowitz` の許容値 0–3、`condSigmaScale` (既定 1.0)。struct 変更 → full rebuild。
5. `tests/unit/test_cond_kantrowitz_carrier.cu` (nvcc, host+device): (a) 純蒸気極限の解析値 (mode 2: $(b_L-\tfrac12)^2/(\tilde c+\tfrac12)$, mode 3: 同+表面項) と 1e-12、
   (b) Wysłouzil 条件 (230 K, ln S 3.4) で θ = 167.4 / 4.05 / 2.98 ± 1 %、(c) 掃引 T 200–260 K × ln S 2–6 × $Y_w$ 0.005–0.05 × g/Y_w 0–0.99 で θ 有限・単調 (Yv→0 で θ→0, Yc→0 で純蒸気形)、
   $q<0$ が掃引域に無い、(d) mode 0/1 と `sigmaScale` 1.0 で旧関数値とビット同一、(e) **同じ入力を device kernel で評価し CPU double と一致** (double 経路 1e-12; float 入力からの経路 1e-6)、
   (f) キー省略時の既定 (0, 1.0)。
6. run (§6)、README、methods、plan §9。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 |
| --- | --- | --- |
| 1 | ~~codex plan レビュー~~ | 決着 (2026-09-12, §6.1): GO-with-changes M6/m3 全件採用 |
| 2 | ~~実装 (§5 1–5)~~ | 決着 (2026-09-12, §9): 実装・単体 ALL PASS |
| 3 | ~~Wysłouzil run (§6)~~ | 決着 (2026-09-12, §9): node 6 run + 凝縮 OFF + cell mode 1/3 対照 完了 |
| 4 | codex result レビュー | → `status: done` |
| 5 | (後続) 分圧スイープ | 0.5 / 0.26 kPa の dry node 場を作ってから mode 1/2/3 の分圧応答 |
| 6 | (後続) 既定値の決定 | 2/3 を既定にするかはスイープ結果と実験一致で判断 (plan に記録) |

## 6. 検証

- **単体**: §5 (4)。Feder 純蒸気極限の解析値と 1e-12 一致、Wysłouzil 条件の $\theta$ がノート §1.3 の手計算 (4.1 / 3.0) と 5 % 以内。
- **検証ケース**: `case/16.nozzle_wys/` 2D node SST 凝縮、run_0335 プロトコル (IC=run_0213 dry 場 index コピー, `outflow`, cfl 2, TP MIXDRY+H2O, HK 成長, `condSonicModel` 自動=1, 48000 step, 出力 level 2)。

  | run | `condKantrowitz` | σ 倍率 | 目的 |
  | --- | --- | --- | --- |
  | `run_0350_kw1_ref` | 1 | 1.0 | 参照 (=run_0335 の再現、新バイナリ) |
  | `run_0351_kw0_iso` | 0 | 1.0 | 等温 CNT (下限) |
  | `run_0352_kw2_feder` | 2 | 1.0 | Feder carrier |
  | `run_0353_kw3_feder_surf` | 3 | 1.0 | Feder carrier + 表面仕事 |
  | `run_0354_kw3_sig097` / `run_0355_kw3_sig103` | 3 | 0.97 / 1.03 | σ の一定倍率局所感度 (mode 3 基準; codex M5) |
  | `run_0356_dry_regress` | — | — | 凝縮 OFF (run_0213 継続 1000 step) 新旧バイナリ回帰 |
  | `run_0357_cell_kw1` / `run_0358_cell_kw3` | 1 / 3 | 1.0 | cell 対照 (run_0341 プロトコル, dry 場から 48000; 非ゲート・未収束準定常比較) |

- **判定基準 (合否)**:
  1. 単体全 PASS (§5 (5) a–f)。特に局所序列 $J_0\ge J_3\ge J_2\ge J_1$ を掃引域全点で確認 (これが序列の合否; codex M4)。
  2. 回帰: `run_0350` (mode 1, σ 1.0) が main 側 run_0335 と `diff_res.py --tolfile noise_node_48000.json --factor 2` で PASS; `run_0356` 凝縮 OFF が新旧バイナリで反復ノイズ以内。
  3. 全 node run NaN 0、`check_convergence.py` PASS (凝縮列込み)、`compare_condfix.py --series` STEADY (未達なら延長、届かなければ「未収束」と明記)。cell は非ゲート。
  4. 核生成域で `condLim_<s>` (律速係数) が 1 に近い (律速が結果を支配していない)。支配していれば cfl を下げた対照を追加。
  5. 実 CUDA 照合: run の `condTheta_<s>` を res の (T, P, ρ, Y, g) から host 式で再計算し一致 (≤1e-5)。
  6. メッシュ品質 `MESH_QUALITY.txt` PASS (継承)、IC は同一メッシュ index コピー、`check_quasisteady.py pmax,machmax` も併記。
- **観測項目 (合否にしない)**: 結合 run の onset [mm] とその序列、壁 p/p₀ の実験偏差 (平均・21/42/52 mm)、$g_{exit}$、σ ±3 % の onset 移動 (符号・量)、cell/node の一致。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan | `2026-09-12` | [`notes/reviews/2026-09-12-condensation-kantrowitz-carrier-plan.md`](../../notes/reviews/2026-09-12-condensation-kantrowitz-carrier-plan.md) | GO-with-changes, C0/M6/m3 | **全件採用**。M1 (混合 cp からの引き算は純蒸気極限で破綻) → 種 DB から種別和、$a_v$ 形で Yv→0 有限 (§4.1)。M2 (src_jac 経路) → `cond_source_vector` と全呼び出しに同一モデル、T 凍結 (§4.1, §5)。M3 (Vinš 2020 の誤読・Tolman 引用) → methods/ノート訂正、Tolman 見送りの根拠を「検証不足」に (§4.2)。M4 (onset 序列) → 局所 J 序列を合否、結合 onset は観測 (§6)。M5 (σ 感度の基準) → mode 3 基準 (§6)。M6 (検証不足) → 掃引・CUDA 照合・cell・凝縮 OFF・律速診断・品質/IC/定常 (§5, §6)。m1 (mode の位置づけ, q<0) → §1/§4.1。m2 (既定 0) → §1。m3 (1 nm=139 分子, 障壁 +17–23 kT) → ノート §2 訂正 |

## 7. 影響範囲

- 触るファイル: `condensationProperties_d.cuh` (σ 倍率), `condensationSource_d.cu(h)`, `variables` (診断 2 本), `input/solverConfig.hpp/.cpp`, `tests/unit/test_cond_kantrowitz_carrier.cu`, `methods/condensation.md`, case/16 README。
- グローバル既定 (`condKantrowitz` 0, `condSigmaScale` 1.0) と mode 1 の経路はビット同一。既存 run への影響なし。

## 8. 完了条件

- [x] methods §2 更新 (Feder 形・表面張力の妥当性・結果)
- [x] 実装・単体・run (§6; 2026-09-12)
- [ ] codex レビュー plan / result
- [ ] `status: done`、accepted へ移動、README 同期

## 9. 変更ログ

- `2026-09-12` — 初稿。文献調査 (ノート) に基づき Feder carrier 形 (mode 2/3) と σ 感度キーを設計。
- `2026-09-12` — **実装・node 検証完了** (case/16 README「carrier 中の非等温核生成補正」節、図 `compare_kantrowitz_carrier.png`, 表 `compare_kantrowitz_carrier.txt`):
  - 実装: `cond_kantrowitz_theta` (mode 1 = 旧演算順でビット不変, 2/3 = Feder $a_v$ 形)、`CondNucCarrier` を種 DB から kernel 内で集計し本体・src_jac 摂動の全呼び出しに同一値、
    `CondSpeciesProps.sigmaScale` → `cond_sigma`、診断 `condTheta_<s>`/`condLim_<s>`、キー `condKantrowitz` 0–3 (範囲検査)・`condSigmaScale`。
  - 単体 `tests/unit/test_cond_kantrowitz_carrier.cu` ALL PASS: 純蒸気極限の解析値 1e-12、Wysłouzil 条件 θ1/θ2/θ3 = 167.4/4.03/2.97 (codex 再計算 167.4/4.05/2.98 と 0.5 %)、
    掃引 400 状態で θ 有限・局所序列 $J_0\ge J_3\ge J_2\ge J_1$・$\hat q>0$ (最小 14.6)・Yv→0 で θ→0・Yc→0 で純蒸気形、mode 0/1 と σ 1.0 は旧関数とビット同一、device==host 1e-12。
  - 回帰: run_0350 (mode 1) は main 側 run_0335 と 29 変数 `diff_res --tolfile noise_node_48000.json --factor 2` PASS; 凝縮 OFF run_0356 は旧バイナリ run_0339 と反復ノイズ以内。
  - 実カーネル照合: `verify_theta.py` で res の `condTheta_0` が host 式と 1.2e-4 (θ>0.1) / 絶対 4e-4 以内 (float 保存) → PASS (mode 1/2/3)。律速係数 `condLim_0` は核生成域で 1.000。
  - **結果 (node 48000 step, 全 run PASS+STEADY; σ 0.97 のみ未収束)**:

    | run | mode / σ | onset [mm] | 壁 p/p₀ 偏差 @21 / 42 / 52 mm [%] |
    | --- | --- | --- | --- |
    | run_0350 | 1 (Kantrowitz 純蒸気) | 22.51 | −4.5 / +5.1 / +5.4 |
    | run_0351 | 0 (等温) | 12.28 | +21.6 / +2.6 / +2.1 |
    | run_0352 | 2 (Feder carrier) | 14.75 | +15.0 / +3.7 / +3.0 |
    | run_0353 | 3 (Feder carrier + 表面仕事) | 14.28 | +17.5 / +3.5 / +2.9 |
    | run_0354 | 3, σ×0.97 (未収束・過渡) | 12.09 | +21.4 / +2.5 / +2.1 |
    | run_0355 | 3, σ×1.03 | 16.81 | +6.2 / +4.4 / +3.7 |

    cell 対照 (run_0357/0358, dry 場から 48000, 既知の床で plateau=未収束準定常比較, 非ゲート): mode 1 は main run_0341 と報告量同一・場差は cell 床の ≤1.2 倍、
    **mode 3 の onset 14.28 mm は node と同一**。
    結合 run の onset 序列は局所 θ の序列 (0 < 3 < 2 < 1) と一致 (観測)。**carrier 形 (θ 3–4) は等温 (θ=0) 側に寄り、Wysłouzil 実験の onset 帯 (壁圧 @21 mm で
    forge mode 1 が −4.5 %) から +15〜18 % 早い側へ外れる**。σ ±3 % で onset は ∓2.3 mm 動き、carrier 形と純蒸気形の差 (8 mm) は σ 換算で ~+10 % に相当。
    **解釈**: キャリア冷却を物理どおり入れると J は等温 CNT に近づき、水の CNT が低温で J を桁で過大評価する既知の傾向 (Wölk & Strey 2001) がそのまま出る。
    純蒸気 Kantrowitz (mode 1) の実験との「良い一致」は、この J 過大を非等温抑制で偶然補償していたと見るのが整合的。したがって物理モデルとしては
    mode 3 + 別途の J 較正 (CNT 補正・σ の低温補正) が筋で、mode 1 は経験的較正として残す。既定は 0 のまま (§10)。
- `2026-09-12` — codex plan レビュー GO-with-changes (M6/m3) を全件採用 (§6.1)。衝突項を種 DB から種別和で評価する $a_v$ 形に変更、src_jac 経路への伝播、
  σ 感度を mode 3 基準に、局所 J 序列を合否に、掃引/CUDA 照合/cell/凝縮 OFF/律速診断を追加、文献要約 (Vinš 2020, Wilhelmsen 2015) を訂正。`status: in_progress`。

## 10. 未確定事項

- mode 2/3 を既定にするか: 結果 (§9) は「mode 3 は物理的に正しい方向だが CNT の J 過大が露出し実験より 8 mm 早い」なので、**既定化は J 較正 (Wölk–Strey 型補正
  か σ 低温補正) と組にする必要がある**。分圧スイープ (0.5/0.26 kPa) で mode 1 の補償が偶然かを見るのが次の判断材料。
- Tolman 補正 (δ<0, MD) を入れるか — 実測の裏付けが無いので保留。
