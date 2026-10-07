"""case/66 回帰ハーネスの構成表 (plan architecture-solver-host-memory §6 の構成表を具体化したもの)。

prepare_inputs.py (ローカル: 元 run から入力を複製して inputs/<入力名>/ を作る) と
run_matrix.py (AWS: 入力から run_* を作って回す) の両方がこのファイルを読む。正本はここ。

用語:
  入力 (INPUTS)  … inputs/<名前>/ に置くテンプレート。元 run から複製し、起動に要る最小の修正だけを加える。
  構成 (CONFIGS) … 入力 + 環境変数 + 種類 (forge / 変換器)。同じ入力を env 違いで複数の構成が共有する。

元入力は別セッションのワークツリー (/home/sano/work/forge/case/...) にある。**読むだけで書かない**。
"""

SRC_ROOT = "/home/sano/work/forge/case"

# ---- 設定ファイルの修正 (複製側だけ。元には触らない) ------------------------------------------
# (正規表現, 置換, 期待一致数)。置換後に PyYAML で読んで期待値を検査する (prepare_inputs.py)。


def steps(n, out=None):
    """nStepOuter と time.outStepInterval を書き換える (初期出力 res_0 と最終 res_n の 2 回に絞る)。"""
    out = n if out is None else out
    return [
        (r"nStepOuter:\s*\d+", f"nStepOuter: {n}", 1),
        (r"(?m)^(\s*)outStepInterval:\s*\d+", rf"\g<1>outStepInterval: {out}", 1),
    ]


# 旧キー体系 (LESorRANS/LESmodel/RANSmodel) は現行 solver が起動時に拒否する (recommended-settings §9)。
OLD_TURB_SST = [
    (r"(?m)^(\s*)LESorRANS:\s*2\b.*$", r'\g<1>model: "sst"   # 複製時に旧キー LESorRANS 2 / RANSmodel 1 から移行', 1),
    (r"(?m)^\s*LESmodel:\s*\d+.*\n", "", 1),
    (r"(?m)^\s*RANSmodel:\s*\d+.*\n", "", 1),
]
OLD_TURB_WALE = [
    (r"(?m)^(\s*)LESorRANS:\s*1\b.*$", r'\g<1>model: "wale"   # 複製時に旧キー LESorRANS 1 / LESmodel 1 から移行', 1),
    (r"(?m)^\s*LESmodel:\s*1\b.*\n", "", 1),
]
OLD_TURB_NONE = [
    (r"(?m)^(\s*)LESorRANS:\s*0\b.*$", r'\g<1>model: "none"   # 複製時に旧キー LESorRANS 0 から移行', 1),
    (r"(?m)^\s*LESmodel:\s*\d+.*\n", "", 1),
]
# SST 壁関数は使わない (2026-09-20 ユーザ方針、recommended-settings §2)。複製時に低 Re (0) へ。
WALL_SST0 = [(r"wallTreatmentSST:\s*1\b", "wallTreatmentSST: 0   # 複製時に 1 (壁関数, 使用禁止) から変更", 1)]

# 凝縮 ON では凝縮する気体 H2O は組込み (液相と同じ CEA 基準) でなければならず、speciesDBFile の H2O は拒否される
# (現行 solver の [speciesDB] 規約)。複製側の species_db.yaml から H2O を外す (起動ログの移行手順 (1))。
DROP_H2O = [(r'(?ms)^"H2O":.*\Z', "", 1)]

PROBE_C36 = """# case/66 で追加した点 probe (plan §6: case/36 に 2〜3 点を足す)。outStepInterval は 20 step ごと
outStepInterval: 20
outStepStart: 0
points:
  p0:
    x: 0.10
    y: 0.0
    z: 0.0
  p1:
    x: 0.35
    y: 0.005
    z: 0.0
  p2:
    x: 0.60
    y: -0.01
    z: 0.0
surfaces:
"""

# ---- 入力テンプレート ----------------------------------------------------------------------
# src     : 元 run (SRC_ROOT 相対)
# copy    : 複製するファイル (src 相対。{"to": 名前} で改名)。run には「copy の宛先 + 生成物」を写す
# seed    : 同一メッシュ restart の種 (src 相対の res を seed_src.h5 として置き、AWS で restart_field.py を掛けて dst へ)
# ckpt    : dual-time の再開に使う checkpoint (base r1 の ckpt100 の res_100.h5 を AWS で置く。run_matrix.py set-ckpt)
# edits   : {ファイル: [(正規表現, 置換, 期待数)]}
# expect  : 修正後の solverConfig.yaml で検査する値 (ドット区切りのキー: 値)
# files   : 上書き生成するファイル {名前: 内容}
# note    : README に写す説明 (直した点)
INPUTS = {
    # ---- 2D node 標準 (表の A・G・H・probe・壁出力) ----
    "c36node": dict(
        src="36.passive_pseudoshock_control/run_sym_H_2up_node",
        copy=["solverConfig.yaml", "bcondConfig.yaml", "passive_solid.h5"],
        seed=dict(src="res_80000.h5", dst="passive_solid.h5"),
        edits={"solverConfig.yaml": steps(200) + OLD_TURB_SST + WALL_SST0},
        files={"probe.yaml": PROBE_C36},
        expect={"turbulence.model": "sst", "turbulence.wallTreatmentSST": 0, "time.last.nStepOuter": 200},
        note="旧キー LESorRANS/LESmodel/RANSmodel → `model: sst`、`wallTreatmentSST` 1→0、probe 3 点を追加 (20 step ごと)",
    ),
    # ---- dual-time checkpoint (表の F・I・J・K・L・B・D・N) ----
    "c09ckpt100": dict(
        src="09.Taylor-Green/run_0160_passiveG_fct_ckpt100",
        copy=["solverConfig.yaml", "bcondConfig.yaml", "probe.yaml", "species_db.yaml",
              "Taylor-Green.h5", "TG_stepxi_gaussY_seam.h5"],
        edits={"solverConfig.yaml": steps(100)},
        expect={"time.dualTime": 1, "output.level": 2},
        note="無修正 (100 step、res_100 に /CHECKPOINT)",
    ),
    "c09cont200": dict(
        src="09.Taylor-Green/run_0162_passiveG_fct_cont200",
        copy=["solverConfig.yaml", "bcondConfig.yaml", "probe.yaml", "species_db.yaml",
              "Taylor-Green.h5", "TG_stepxi_gaussY_seam.h5"],
        edits={"solverConfig.yaml": steps(200, 100)},
        expect={"time.dualTime": 1},
        note="無修正 (連続 200 step、res_100/res_200)",
    ),
    "c09restart100": dict(
        src="09.Taylor-Green/run_0164_passiveG_fct_restart100_fixed",
        copy=["solverConfig.yaml", "bcondConfig.yaml", "probe.yaml", "species_db.yaml", "Taylor-Green.h5"],
        ckpt=dict(cfg="c09ckpt100", file="res_100.h5", dst="ic_ckpt.h5"),
        edits={"solverConfig.yaml": steps(100)},
        expect={"mesh.valueFileName": "ic_ckpt.h5"},
        note="無修正。`ic_ckpt.h5` = base r1 の c09ckpt100 の res_100.h5 (全ビルド・全反復で同じファイル)",
    ),
    # ---- 軸対称・多成分・凝縮 ----
    "c44steady": dict(
        src="44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX",
        copy=["solverConfig.yaml", "bcondConfig.yaml", "probe.yaml", "species_db.yaml", "nozzle.h5"],
        seed=dict(src="res_24000.h5", dst="nozzle.h5", force_species=True,
                  precreate=["rog_0", "roQ2_0", "roQ1_0", "roQ0_0"]),
        edits={"solverConfig.yaml": steps(200), "species_db.yaml": DROP_H2O},
        expect={"condensation.condensation": 1, "output.level": 2},
        note="無修正 (定常 active run、res_24000 から restart)。凝縮モーメント rog_0/roQ*_0 は変換直後の nozzle.h5 に無いので 0 で作ってから"
             "写す。種の記録が無い場なので restart_field は --force-species、solver は FORGE_ALLOW_UNVERIFIED_SPECIES=1。"
             "species_db.yaml から H2O を外す (凝縮 ON は組込み H2O が必須、現行 solver が拒否する)",
    ),
    "c44dual_ckpt100": dict(
        src="44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float",
        copy=["solverConfig.yaml", "bcondConfig.yaml", "probe.yaml", "species_db.yaml", "nozzle.h5", "ic.h5",
              "inlet_profile_1.csv"],
        edits={"solverConfig.yaml": steps(100), "species_db.yaml": DROP_H2O},
        expect={"time.dualTime": 1, "output.level": 2},
        note="species_db.yaml から H2O を外す (凝縮 ON は組込み H2O が必須)。run_0376 型 (plan §6) のうち推奨設定に合う run_0468 (run_0376 は dual-time で implicitRelax 0.7・SFR 2 なのに convMethod 0 で check_solver_config が FAIL)。200 step を 100 + 100 に分割 (前半)",
    ),
    "c44dual_restart100": dict(
        src="44.vitiated_air_wt/run_0468_sweep_cflp12_nsub20_float",
        copy=["solverConfig.yaml", "bcondConfig.yaml", "probe.yaml", "species_db.yaml", "nozzle.h5",
              "inlet_profile_1.csv"],
        ckpt=dict(cfg="c44dual_ckpt100", file="res_100.h5", dst="ic_ckpt.h5"),
        edits={"solverConfig.yaml": steps(100) + [(r"valueFileName:\s*ic\.h5", "valueFileName: ic_ckpt.h5", 1)],
               "species_db.yaml": DROP_H2O},
        expect={"mesh.valueFileName": "ic_ckpt.h5"},
        note="species_db.yaml から H2O を外す。後半。`valueFileName` を ic.h5 → ic_ckpt.h5 (= base r1 の c44dual_ckpt100 の res_100.h5)",
    ),
    # ---- 共役伝熱 ----
    "c52cht": dict(
        src="52.conjugate_slab/run_0007_fxhalf",
        copy=["solverConfig.yaml", "bcondConfig.yaml", "probe.yaml", "mesh.h5", "wall_profile_3.csv"],
        seed=dict(src="res_8000.h5", dst="mesh.h5"),
        edits={"solverConfig.yaml": steps(200)},
        expect={"output.interfaceDiag": 1},
        note="無修正 (res_8000 から restart。CHT の内部状態は res に無いので界面は再初期化)",
    ),
    # ---- cell モード ----
    "c36cell": dict(
        src="36.passive_pseudoshock_control/run_sym_H_2up_cell",
        copy=["solverConfig.yaml", "bcondConfig.yaml", "probe.yaml", "passive_solid.h5"],
        seed=dict(src="res_60000.h5", dst="passive_solid.h5"),
        edits={"solverConfig.yaml": steps(200) + OLD_TURB_SST + WALL_SST0},
        expect={"turbulence.model": "sst", "turbulence.wallTreatmentSST": 0},
        note="旧キー → `model: sst`、`wallTreatmentSST` 1→0",
    ),
    "c20cell_rk3": dict(
        src="20.naca_ml/001.test/run_slau",
        copy=["solverConfig.yaml", "bcondConfig.yaml", "probe.yaml", "naca.h5"],
        edits={"solverConfig.yaml": steps(200)},
        expect={"time.timeIntegration": 3},
        note="無修正 (cell RK3、変換時の初期場から)",
    ),
    "c20cell_dual": dict(
        src="20.naca_ml/001.test/run_case04_les_unsteady_dualtime",
        copy=["solverConfig.yaml", "bcondConfig.yaml", "probe.yaml", "naca.h5"],
        edits={"solverConfig.yaml": steps(100) + OLD_TURB_WALE},
        expect={"turbulence.model": "wale", "time.dualTime": 1},
        note="旧キー LESorRANS 1 / LESmodel 1 → `model: wale`",
    ),
    "c20cell_impl": dict(
        src="20.naca_ml/001.test/run_slau_20260511_003420_implicit_diag_cfl2_smoke",
        copy=["solverConfig.yaml", "bcondConfig.yaml", "probe.yaml", "naca.h5"],
        edits={"solverConfig.yaml": steps(200) + OLD_TURB_NONE},
        expect={"turbulence.model": "none", "time.timeIntegration": 11},
        note="旧キー LESorRANS 0 → `model: none` (cell 陰解法。FORGE_IMPLICIT_DIAG_CSV の cell 版)",
    ),
    # ---- 遷移モデル ----
    "c57lm": dict(
        src="57.transition_flat_plate/run_0014_t3a_lm_unitcheck",
        copy=["solverConfig.yaml", "bcondConfig.yaml", "probe.yaml", "mesh.h5"],
        edits={"solverConfig.yaml": steps(200)},
        expect={"turbulence.transition": "lm2009", "output.level": 2},
        note="無修正 (mesh.h5 は run_0011 の場を種にしたもの = roGamma を読む経路)",
    ),
    "c57lm_fromsst": dict(
        src="57.transition_flat_plate/run_0013_t3b_lm",
        copy=["solverConfig.yaml", "bcondConfig.yaml", "probe.yaml",
              {"from": "../run_0012_t3b_sst/mesh.h5", "to": "mesh.h5"}],
        seed=dict(src="../run_0012_t3b_sst/res_20000.h5", dst="mesh.h5"),
        edits={"solverConfig.yaml": steps(200)},
        expect={"turbulence.transition": "lm2009"},
        note="run_0013 の LM 設定 + run_0012 (SST だけ) の mesh.h5 に res_20000 を写した場 (roGamma が無い = 初期化する経路)",
    ),
    # ---- line-implicit と extraFields ----
    "c56lineimp": dict(
        src="56.gap_tp1187/run_0019_lineimplicit",
        copy=["solverConfig.yaml", "bcondConfig.yaml", "probe.yaml", "mesh.h5"],
        seed=dict(src="res_100000.h5", dst="mesh.h5"),
        edits={"solverConfig.yaml": steps(200)},
        expect={"time.deltaT.lineImplicit": 1},
        note="無修正 (res_100000 から restart)",
    ),
    "c56extra": dict(
        src="56.gap_tp1187/run_0027_s6_f32_b",
        copy=["solverConfig.yaml", "bcondConfig.yaml", "probe.yaml", "mesh.h5"],
        edits={"solverConfig.yaml": steps(200)},
        expect={"time.deltaT.qAccumulatorFP64": 1},
        note="無修正 (extraFields に res_ro 等。元は 1 step → 200 step)",
    ),
    "c48absorb": dict(
        src="48.flat_plate_cooled_m4/run_0903_absorb2",
        copy=["solverConfig.yaml", "bcondConfig.yaml", "probe.yaml", "mesh.h5"],
        edits={"solverConfig.yaml": steps(200)},
        expect={"turbulence.wallTreatmentSST": 0},
        note="無修正 (extraFields に roN.. と dq_block_old_*。元は 1 step → 200 step)",
    ),
    # ---- 環境変数で登録が変わる診断 ----
    "c26optin": dict(
        src="26.flat_plate_sst/run_0086_optin_base",
        copy=["solverConfig.yaml", "bcondConfig.yaml", "probe.yaml", "flat_plate_ny52_planar.h5"],
        edits={"solverConfig.yaml": steps(200) + WALL_SST0},
        expect={"turbulence.wallTreatmentSST": 0},
        note="`wallTreatmentSST` 1→0 (env による名前の登録・除去は壁処理に依らない: variables.cpp:301-356)",
    ),
    # ---- 変換器 (variables.cpp を共有する) ----
    "v36node": dict(
        src="36.passive_pseudoshock_control",
        copy=[{"from": "mesh/passive_solid.msh", "to": "mesh.msh"}],
        derive=dict(input="c36node", files=["solverConfig.yaml", "bcondConfig.yaml"]),
        note="c36node の config (node) で mesh/passive_solid.msh を変換",
    ),
    "v36cell": dict(
        src="36.passive_pseudoshock_control",
        copy=[{"from": "mesh/passive_solid.msh", "to": "mesh.msh"}],
        derive=dict(input="c36cell", files=["solverConfig.yaml", "bcondConfig.yaml"]),
        note="c36cell の config (cell) で mesh/passive_solid.msh を変換",
    ),
    "v09": dict(
        src="09.Taylor-Green",
        copy=[{"from": "mesh/Taylor-Green.msh", "to": "mesh.msh"}],
        derive=dict(input="c09ckpt100", files=["solverConfig.yaml", "bcondConfig.yaml", "species_db.yaml"]),
        note="c09ckpt100 の config (node 周期・2 種・トレーサ) で mesh/Taylor-Green.msh を変換",
    ),
    "v52": dict(
        src="52.conjugate_slab",
        copy=[{"from": "mesh/slab.msh", "to": "mesh.msh"},
              {"from": "mesh/solverConfig.yaml", "to": "solverConfig.yaml"},
              {"from": "mesh/bcondConfig.yaml", "to": "bcondConfig.yaml"}],
        note="case/52 mesh/ の変換 config のまま",
    ),
    "v44": dict(
        src="44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX",
        copy=[{"from": "nozzle.msh", "to": "mesh.msh"}],
        derive=dict(input="c44steady", files=["solverConfig.yaml", "bcondConfig.yaml", "species_db.yaml"]),
        note="c44steady の config (軸対称・2 種・凝縮) で nozzle.msh を変換",
    ),
}

# ---- 構成 (run の単位) ---------------------------------------------------------------------
# kind: forge (run_case.sh) / convert (変換器)。N は入力の nStepOuter。checks は README の「確認する成果物」。
UNVERIFIED = {"FORGE_ALLOW_UNVERIFIED_SPECIES": "1"}
CONFIGS = {
    "c36node":            dict(input="c36node", env={}, group="2D node 標準"),
    "c36node_impdiag":    dict(input="c36node", env={"FORGE_IMPLICIT_DIAG_CSV": "diag.csv"}, group="2D node + 陰解法診断 CSV (表 M)"),
    "c36node_psidual":    dict(input="c36node", env={"FORGE_DIAG_PSI_DUALEVAL": "2,4"}, group="2D node + ψ 二重評価 (R2 の pdeSize)"),
    "c09ckpt100":         dict(input="c09ckpt100", env={}, group="dual-time 100 step (checkpoint 書出し)"),
    "c09restart100":      dict(input="c09restart100", env={}, group="dual-time 再開 100 step (checkpoint 復元)"),
    "c09cont200":         dict(input="c09cont200", env={}, group="dual-time 連続 200 step"),
    "c09ckpt100_outres":  dict(input="c09ckpt100", env={"FORGE_OUT_RESIDUALS": "1", "FORGE_RESID_SNAP": "0"}, group="dual-time + 残差出力 (表 K・L)"),
    "c09ckpt100_rawdiag": dict(input="c09ckpt100", env={"FORGE_SPECIES_RAW_DIAG": "1"}, group="dual-time + roYraw (表 I)"),
    "c09ckpt100_pindiag": dict(input="c09ckpt100", env={"FORGE_PIN_DIAG": "1"}, group="dual-time + ピン診断 (表 N)"),
    # c44: 元の場に化学種の記録 (species_hash) が無いので solver の照合を env で許可する (host c の名前集合には効かない: 監査 §1 の env 一覧)
    "c44steady":          dict(input="c44steady", env=UNVERIFIED, group="軸対称・多成分・凝縮 (定常)"),
    "c44dual_ckpt100":    dict(input="c44dual_ckpt100", env=UNVERIFIED, group="軸対称・凝縮 dual-time 前半"),
    "c44dual_restart100": dict(input="c44dual_restart100", env=UNVERIFIED, group="軸対称・凝縮 dual-time 後半 (再開)"),
    "c44dual_pindiag":    dict(input="c44dual_ckpt100", env=dict(UNVERIFIED, FORGE_PIN_DIAG="1"), group="軸対称・凝縮 dual-time + ピン診断 (入口あり)"),
    "c52cht":             dict(input="c52cht", env={}, group="共役伝熱"),
    "c36cell":            dict(input="c36cell", env={}, group="cell 定常 SST"),
    "c20cell_rk3":        dict(input="c20cell_rk3", env={}, group="cell RK3"),
    "c20cell_dual":       dict(input="c20cell_dual", env={}, group="cell LES dual-time"),
    "c20cell_impdiag":    dict(input="c20cell_impl", env={"FORGE_IMPLICIT_DIAG_CSV": "diag.csv"}, group="cell 陰解法 + 診断 CSV"),
    "c57lm":              dict(input="c57lm", env={}, group="遷移 (roGamma を読む)"),
    "c57lm_fromsst":      dict(input="c57lm_fromsst", env={}, group="遷移 (SST の場から初期化)"),
    "c56lineimp":         dict(input="c56lineimp", env={}, group="line-implicit + extraFields"),
    "c56extra":           dict(input="c56extra", env={}, group="extraFields (res_* を含む) + FP64 アキュムレータ"),
    "c48absorb":          dict(input="c48absorb", env={}, group="extraFields (roN・dq_block_old_*)"),
    "c26optin":           dict(input="c26optin", env={}, group="env 診断なし (名前が消えて警告)"),
    "c26optin_env":       dict(input="c26optin", env={"FORGE_WI_FORCE_DIAG": "1", "FORGE_WF_CLOSURE_DIAG": "1",
                                                      "FORGE_OMEGA_BUDGET": "1", "FORGE_WF_REP_DIAG": "1"},
                               group="env 診断 4 種 ON (表 I)"),
    "v36node":            dict(input="v36node", env={}, kind="convert", group="変換器 node"),
    "v36cell":            dict(input="v36cell", env={}, kind="convert", group="変換器 cell"),
    "v09":                dict(input="v09", env={}, kind="convert", group="変換器 node 周期・種"),
    "v52":                dict(input="v52", env={}, kind="convert", group="変換器 node (CHT 用スラブ)"),
    "v44":                dict(input="v44", env={}, kind="convert", group="変換器 軸対称・種・凝縮"),
}
