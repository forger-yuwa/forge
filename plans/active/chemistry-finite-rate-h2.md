# 有限速度化学 (H₂ 燃焼・ノズル化学非平衡) の導入

## メタ

- **area**: `thermophysics / chemistry`
- **status**: `in_progress` (Phase 0・1 完了 2026-09-04、Phase 2 未着手)
- **related_docs**:
  - `methods/chemistry.md` (現在仕様: 理論・実装方針)
  - `methods/thermophysics.md` (多成分 TP gas、sensible datum、種輸送・陰解法結合)
  - `notes/investigations/chemistry-finite-rate-h2-survey.md` (文献調査・方針の根拠)
- **related_plans**:
  - `plans/accepted/thermophysics-multicomponent-tpgas.md` (前提。§10 に反応は将来拡張と明記)
  - `plans/accepted/convection-multispecies-contact-pressure.md`
- **created**: `2026-09-04`
- **owner**: `CFD Dev`

## 1. 目的

forge の多成分 TP gas に化学反応ソース項を加え、(B) ノズル膨張中のラジカル再結合・凍結 (試験部の残留 OH/H/O/NO と $T$, $M$ への影響)、
(A) H₂ 燃焼加熱器 (低マッハ・RANS) を解けるようにする。完了時には frozen / equilibrium / finite-rate の三者比較が
ノズル設計チェーン (`methods/design/`) の確認 CFD メニューとして使える状態になる。

## 2. スコープ

- **やる**: 反応機構 (Cantera YAML サブセット) 読込、Arrhenius・三体・逆反応 ($K_c$)、反応熱の sensible datum 整合、解析 Jacobian、
  種ブロック point-implicit (定常) と Strang 分離 (非定常)、0-D / Q1D ノズル / Burrows–Kurkov 検証、No-TCI → PaSR。
- **やらない**: fully coupled $(5+n_s)$ block-DPLUR、EDC / flamelet、熱的非平衡 (多温度)、炭化水素機構 (CH₄ 加熱器は後続)、
  燃焼器の着火過渡 (定常は高温平衡 IC から始める)。

## 3. 関連 docs と前提

- 理論・実装方針: `methods/chemistry.md` (本 plan と同時に起草)。
- 既存資産: `speciesTransport_d` (種輸送・拡散・`src_jac_Y{s}` 配管・scalar-DPLUR sweep・案C EOS-JVP)、
  `condensationSource_d` / `ransSource_d` (point-implicit ソース項の雛型)、`thermo_d` (NASA-9, `THERMO_HD` host/device 共通)。
- 既知の制約: 多成分 TP 陰解法の安定 `cfl_pseudo` は 1–2 (thermophysics plan §10)。反応で剛性が上乗せされる。

## 4. 設計方針 (差分のみ。本文は `methods/chemistry.md`)

1. **datum**: `thermoHrefTemp: 298.15` を維持し、$\dot Q=-\sum_s c_s\dot\omega_s$ ($c_s=h^{\rm abs}_s(T_{\rm ref})$) をエネルギー式に陽注入。
   $K_c$ は焼き込み前の絶対 $a_7$ (`a7_abs`) で評価。
2. **剛性**: 既存 `src_jac_Y{s}` (対角) を $n_s\times n_s$ 密ブロック `src_jac_Y{s,k}` に拡張し `species_dplur_sweep_d` でセル毎 LU。
   流れ側は案C の EOS-JVP に反応寄与を載せ、$(5,5)$ 対角に $\partial\dot Q/\partial(\rho e_t)$。
3. **機構**: Jachimowski 1988 (9 種 20 反応 / 13 種 33 反応、falloff 無し) から。Burke 2012 (Troe) は Phase 2 で falloff 実装後。
4. **精度**: 濃度・Jacobian は double、$\ln k_f$/$\ln K_c$ は FP32。
5. **参照解**: Cantera (`.venv-chem`, CPU) と CEA (`.venv-cea`)。forge と同一 YAML・同一 NASA-9 を使う。

## 5. 実装ステップ

| Phase | 内容 | 主要ファイル | 状態 |
| --- | --- | --- | --- |
| **0** 前提整備 | CEA 平衡/凍結スクリーニング、`thermo.inp→species_db.yaml` ツール、Jachimowski YAML 化 + Cantera 検証、方針文書 | `tools/cea_thermo_to_species_db.py`, `tools/mechanisms/*.yaml`, `notes/investigations/scripts/cea_kinetics_screen/`, `methods/chemistry.md` | **done** (§9) |
| **0b** BC | 指定組成 `Yb` 配線: `inlet_Pressure` は M5/M7 で配線済 (case/44 `run_0091` が `Y0/Y1` を使用中、thermophysics plan §10 の記述は古い)。残りは `inlet_Pressure_dir` のみ (燃焼室出口組成は `inlet_Pressure` で与えられるので低優先) | `cuda_forge/boundaryCond_d.cu` | `inlet_Pressure` done / `_dir` todo |
| **1** ソース項 | `chemistry_d.{cuh,cu}` (機構読込・SI 換算・device 定数)、`chemistrySource_d` ($\dot\omega_s$, $\dot Q$, 解析 Jacobian)、`h_datum` 保持、config キー、ホスト 0-D テスト `tools/test_chemistry.cpp`、GPU 0-D 着火 (case/35 run_0049) | `cuda_forge/chemistry_d.*`, `chemistry_mech_io.hpp`, `thermo_d.*`, `input/solverConfig.*`, `main.cpp`, `periodicNode_d.cu` | **done** (§9) |
| **2** 陰解法結合 | 種ブロック point-implicit (LU)、案C への反応寄与、$(5,5)$ 反応熱 Jacobian、Strang 分離 (陽解法/dual-time)、falloff (Troe) | `speciesTransport_d.cu`, `species_eos_coupling_d.cuh`, `timeIntegration_d.cu` | todo |
| **3** 燃焼器 RANS | 低マッハ前処理 + SST + No-TCI、PaSR、13 種 (NO) | `ransSource_d` 連携, `chemistry_d` | todo |
| **4** 応用 | 加熱器 → ノズル → 試験部 end-to-end、設計チェーン確認メニューへ | `design/`, `methods/design/` | todo |

## 6. 検証

- **単体 / ビルド**: `tools/test_chemistry.cpp` (0-D 定積着火遅れ・定圧平衡到達 vs Cantera 同一 YAML、Jacobian 有限差分照合 $<10^{-6}$ 相対)。
- **検証ケース**:
  1. Q1D ノズル再結合: 新設 `case/46.nozzle_h2o2_kinetics` (H₂/O₂ 平衡入口、軸対称 node)。参照は Cantera PFR (面積則)。
     frozen / equilibrium / finite-rate の $T$, $Y_{OH}$, 凍結位置。
  2. Burrows–Kurkov (M2.44 vitiated air + H₂ 壁噴射、NASA GRC 検証アーカイブ)。
  3. 回帰: `case/28.cutler_coaxial_jet` (非反応多成分) と `case/44` (TP 単一擬似種) がビット不変 (`chemistry.enabled: 0`)。
- **判定基準**: 0-D 着火遅れ $\pm2$ %、平衡温度 $\pm$ 数 K; Q1D ノズルの $T$ 差 $<1$ %、$Y_{OH}$ 凍結値 $\pm10$ %;
  反応 ON で `cfl_pseudo` 上限が非反応と同等 (剛性が陰化で吸収されていること); 全残差 `check_convergence.py` PASS。

## 7. 影響範囲

- 新規: `cuda_forge/chemistry_d.*`, `tools/test_chemistry.cpp`, `tools/mechanisms/`, `tools/cea_thermo_to_species_db.py`。
- 変更: `speciesTransport_d.cu` (ブロック Jacobian), `thermo_d.*` (`a7_abs`), `species_eos_coupling_d.cuh`, `solverConfig.*`, `main.cpp`,
  `boundaryCond_d.cu` (`Yb`)。既定 (`chemistry.enabled: 0`) では全経路ビット不変。
- ドキュメント: `methods/chemistry.md`, `methods/index.md`, `methods/thermophysics.md` (datum 節に反応時の注記), `procedures/solver-settings.md` (config キー)。

## 8. 完了条件

- [x] `methods/chemistry.md` 起草 (理論・実装方針)
- [x] Phase 1 実装・検証 (0-D 着火)
- [ ] Phase 2–3 実装・検証 (§6)
- [ ] `status: done` + §9 変更ログ
- [ ] `plans/active/` → `plans/accepted/` へ移動、`plans/README.md` 同期

## 9. 変更ログ

- `2026-09-04` — 初稿。Phase 0 完了:
  - **CEA スクリーニング** (`notes/investigations/scripts/cea_kinetics_screen/run_screen.py`, 結果 `screen_table.md`, $p_c$ 11.39 bar,
    出口 $A/A^*$ 15.07 (M≈4.2) と 53.2 (M≈6)): 燃焼室温度 $T_c$ ごとの平衡 vs 凍結 (nfz=1 燃焼室凍結 / nfz=2 スロート凍結) の出口差。
    | $T_c$ [K] | 組成 | $X_{OH,c}$ | $X_{NO,c}$ | ΔM (M4.2, eq−fz1) | ΔT (M4.2) [K] | ΔM (M6) |
    |---|---|---|---|---|---|---|
    | 1161 (case/44 va3) | va3 | 0 | 1.3e-4 | 0.000 | 0.1 | 0 (凝縮を除く) |
    | 1595 | va3 | 9e-5 | 1.8e-3 | −0.003 (0.07 %) | +2 | −0.004 |
    | 1979 | va3 | 8.9e-4 | 6.7e-3 | −0.013 (0.3 %) | +10 | −0.019 |
    | 2422 | va3 | 5.1e-3 | 1.8e-2 | −0.043 (1.1 %) | +41 | −0.067 |
    | 2049 (H₂ 加熱器, H₂O 26 %) | B40 | 2.6e-3 | 7.0e-3 | −0.017 (0.4 %) | +14 | −0.025 |
    → **$T_c\le1600$ K では化学非平衡は無視可 (0.1 % 未満)、$T_c\approx2000$ K で 0.3–0.4 %、2400 K で 1 %** (設計ゲート 0.5 % Md と同桁)。
    有限速度解は eq と fz2 の間に入る。case/44 (1161 K) と case/45 (1600 K) は現状の frozen 設計で十分、
    本機能が効くのは $T_c\gtrsim2000$ K の加熱器案件。M6 出口の A1161/B10 eq 列は CEA が H₂O(cr) 凝縮を含めるため T が不連続 (forge は凝縮を別モデルで扱う)。
  - **熱力学 DB ツール** `tools/cea_thermo_to_species_db.py`: CEA `thermo.inp` → `species_db.yaml`。内蔵 N2 係数と相対差 0、
    Cantera と $c_p/R$ 4 桁一致 ($h/RT$ の OH/HO₂ 2–5 % 差は生成エンタルピーの出典差)。生成物 `tools/mechanisms/species_db_h2air_cea.yaml`。
  - **機構ファイル** `tools/mechanisms/h2air_jachimowski1988_13sp33r.yaml` / `h2o2_jachimowski1988_9sp20r.yaml` (TP-2791 Table I、熱力学は CEA NASA-9)。
    Cantera 定積着火遅れ (量論 H₂-air, `notes/investigations/scripts/chem_ignition_check.py`):
    | $T_0$ [K] / p [atm] | Jach13 | Jach9 | Cantera h2o2 (GRI 系) |
    |---|---|---|---|
    | 1000 / 1 | 152 µs | 152 | 305 |
    | 1200 / 1 | 32.2 | 32.2 | 44.2 |
    | 1500 / 1 | 10.1 | 10.1 | 12.8 |
    | 2000 / 1 | 3.5 | 3.5 | 3.6 |
    Jach13 = Jach9 (N 化学は着火に無関係) で YAML 化の整合を確認。Jachimowski は低温側で GRI 系より 1.4–2 倍速い (Slack データに合わせた 1988 校正、既知傾向)。
  - **YAML の罠**: 種名 `NO` / `N` は YAML 1.1 (PyYAML) で真偽値に化ける → DB 生成はキーを引用符付きにした。yaml-cpp (forge) 側も同様の扱いを Phase 1 で確認する。
- `2026-09-04` — **Phase 1 完了** (反応ソース項・反応熱・解析 Jacobian・対角 point-implicit)。
  - 実装: `cuda_forge/chemistry_d.cuh` (`ReactionTable`, `chem_source`: Arrhenius・三体・$K_c$・$\dot Q$・$n_s\times n_s$ 解析 Jacobian・$\partial\omega/\partial T$)、
    `chemistry_mech_io.hpp` (Cantera YAML サブセット読込)、`chemistry_d.cu` (device ソース項、`res_roY`/`res_roe`/`src_jac_Y` 加算、`chemQdot`/`chemTau` 出力)、
    `SpeciesThermo::h_datum` (sensible datum の除去分を保持 → $\dot Q=-\sum c_s\dot\omega_s$、$K_c$ は絶対 $H$)、`physProp.chemistry` キー、
    `registerSpecies(n, chemistry)`。既定 (`enabled: 0`) では全経路ビット不変 (構造体変更のため full rebuild)。
  - ホスト検証 (`tools/test_chemistry.cpp`, 参照 `tools/chem_reference_cantera.py`): Jacobian 有限差分照合 $\partial\omega/\partial\rho Y$ 1.8e-8・$\partial\omega/\partial T$ 1.3e-9 (PASS)、
    $\sum\dot\omega_s\sim10^{-13}$ (質量保存)。0-D 定積 BDF1 反応器 (量論 H₂-air 1200 K 1 atm) の着火遅れは刻みを締めると Cantera へ収束
    (29.4 → 31.1 → 32.0 µs; Cantera 32.2)、平衡 T 2945.2 vs 2943.9 K。13 種 33 反応 (1500 K 2 atm) も 4.7 vs 5.0 µs・3065 vs 3064 K。
    絶対 datum と sensible datum で同一の平衡 T (反応熱経路の検証)。
  - **GPU 検証** `case/35.uniform_periodic_box/run_0049_node_h2_ignition` (node 8³ 全面 periodic = 定積反応器, explicit RK3 dt 5e-9 s, 10000 step 94 s):
    着火遅れ 32.0 µs (Cantera 32.2)、平衡 T 2948.4 K (Cantera 2943.9, +0.15 %)、全ノード一様 (T 差 0.004 K, |u| 1e-4 m/s)。
    図 `h2_ignition_vs_cantera.png`。
  - **発見・修正したバグ (化学以前からの欠落)**: node 周期の `periodicNodeGather` が化学種残差 `res_roY{s}` を合算していなかった
    (流れ 5 変数と k/ω のみ)。seam ノードでは化学種の移流・反応が部分体積分 (面 1/2・辺 1/4・角 1/8) しか効かず、エネルギーだけ合算されて
    不整合 (run_0049 初回: seam が着火せず圧力波で場が乱れた)。`res_roY{s}`・凝縮モーメント残差の gather と `roY{s}` の root→member ミラーを追加。
    体積ソースは `volumePartial_d` を使う (bodyForce/ransSource と同じ規約)。**既往の node 周期×多成分 run (case/28 等は非周期なので影響なし) を要確認**。
  - 残: Phase 2 (種ブロック point-implicit、案C 結合、falloff、Strang)、`chemQdot` の陰解法 $(5,5)$ Jacobian、FP32 化。

