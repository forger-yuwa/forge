"""冷却壁の腕の擬似 CFL の試行と、質量の収支の監査の準備 (plan tooling-nozzle-isothermal-wall-chain §5.1 #27)。
cold_pair.py の prep_ext (延長) と同じ手順 (restart_field でビット一致・FP64 の型のまま、段なし) で、変えてよい設定を
nStepOuter・cfl・cfl_pseudo・outStepInterval と output.extraFields に限る (implicitRelax は変えない — codex 2026-10-08。例外は --relax を明示した §6.13 の run)。
走行中の run が cold_pair.py を使っているので、cold_pair.py は書き換えずにここで包む。

usage (AWS の case dir、別バイナリは COLD_ALT_BINARY=<キー> を前に付ける):
  python3 cold_cfl.py prep <src_run> <run> --steps N --cfl C [--out 5000] [--extra res_ro,volume] [--limiter-ref-from <run>]
    --limiter-ref-from: リミッタの基準値 (limiterRoRef・limiterPRef・limiterARef) を指定した run の forge_run.log の値に固定する。
      既定 (自動) では開始場から決まるので、restart した run は親と別の作用素になる (forge の警告; run_0182 → run_0183 で a_ref が 1.6 % 違った)。
    --line dir|only: 壁法線のライン陰解法 (lineImplicit 1)。dir は方向別の擬似 dt (lineDtDirectional 1) も足す (§5.1 #27 の試行)。
    --isp 0|1: time.deltaT.implicitSolvePrecision (陰解法の行列の組立て・解法の精度。既定 0 = float は FP64 のビルドでも float)。
  python3 cold_pair.py run <run>       (投入は既存の run_one; FORGE_DUMP_MASSFLUX は投入側の環境変数で渡す)
  python3 cold_cfl.py run <run>        (同じ run_one を、COLD_ALT_BINARY の登録を効かせて呼ぶ)
"""
import argparse
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cold_pair as CP  # noqa: E402
NS = CP.NS

# 別バイナリの登録 (plan time_integration-implicit-thermal-jacobian §6 の検証用)。環境変数 COLD_ALT_BINARY にキーを書くと、
# cold_pair.binary_record の照合 (FORGE_SHA・FP64_TREE) をこの登録値に差し替える。変換器は従来の FP64 のもの (CONV_SHA) のまま。
ALT_BINARIES = {
    # キー: (forge の sha256, ソースの作業ツリー)
    # implicitThermalJacobian の検証用: commit 1ad0b9bb + typedef double (座標の stod は HEAD に入っている)、2026-10-09 AWS でビルド
    "thermjac_fp64": ("9ffc4d1efcec6ca4799f4026ef418931cf99bd27ad61f29d264cbfc4639018b7", "~/forge-thermjac-fp64"),
    # 06d1b149 (ビット 4 を追加) + typedef double、2026-10-09 に同じ作業ツリーで再ビルド (上のバイナリは上書きされて残っていない)
    "thermjac5_fp64": ("985aca0f2ebba912851042ec9abc8b7deab9707b20f19f757641dc9b3521215e", "~/forge-thermjac-fp64"),
    # b4771052 (lineDtDirectionalCap を追加) + typedef double、2026-10-09 に同じ作業ツリーで再ビルド
    "thermjac_cap_fp64": ("35e498b14b5f3cdaa09b754f7bbaa6455f54631fd2edf1ad4929964e8828d4bc", "~/forge-thermjac-fp64"),
    # 4d394a71 (ライン上の節点の対角を storeLU の sweep 以外で組まない、plan time_integration-line-implicit-speed 案 A) + typedef double、2026-10-09 AWS で新しい作業ツリーにビルド
    "linespeed_fp64": ("d8b06ebcb91cfc3151cfcedf441b3e79f20f4c3702aa89e287801cb1c4f1d10b", "~/forge-linespeed-fp64"),
    # 5ab83056 (lineViscCoupling 2 = 薄層の粘性・熱伝導の Jacobian、案 A を含む) + typedef double、2026-10-09 AWS でビルド
    # (上のバイナリは 2026-10-09 に下の lineB で上書きされて残っていない)
    "linevisc_fp64": ("6631a87eb1279b2fa4eff22080378ac12db4937890ac4470cf0d9185daf55653", "~/forge-linevisc-fp64"),
    # 64ed2cd6 (ライン内で並列にした Thomas・新旧の比較モード・ライン行列の書き出し、値 2 と案 A を含む) + typedef double、~/forge-linevisc-fp64 を上書きしてビルド
    # (上の lineB は下の lineB2 で上書きされて残っていない)
    "lineB_fp64": ("7a7e9eb9150104babc5b1c9877e15eb826eaa46fde8f2fe96dce74b65ad61f6d", "~/forge-linevisc-fp64"),
    # 874ae90a (並列 Thomas の代入を行の分担に、行のポインタを 1 回だけ選ぶ、書き出しのフラグの修正) + typedef double、~/forge-linevisc-fp64 を上書き
    # (上の lineB2 は下の lineB3 で上書きされて残っていない)
    "lineB2_fp64": ("f3765961a601cff3f63e5ed86e0f0cf84d62138f504ba4a03dc317330ea45fad", "~/forge-linevisc-fp64"),
    # 3764852d (lu5 の行の入れ替えを静的な添字にしてレジスタに置く、既定は 1 ライン 1 スレッド、並列は FORGE_LINE_PAR=1) + typedef double
    # (上の lineB3 は下の lineD で上書きされて残っていない)
    "lineB3_fp64": ("9a094160caa9b3a68cbad20cd7b035a529f8f31a3de34fcf82590039040f6d3f", "~/forge-linevisc-fp64"),
    # d5538001 (lu5 を元に戻す。既定 = 案 A + 値 2 の経路 + 並列 Thomas は FORGE_LINE_PAR=1 + 比較・書き出しのデバッグ) + typedef double
    # (上の lineD は下の lineE で上書きされて残っていない)
    "lineD_fp64": ("9ef80d5f148f73b1ee330a141a6f5fde31737a996ea0b5611e749e79b9442e8d", "~/forge-linevisc-fp64"),
    # 6d49738e (診断用の lineViscCoupling 3 を追加) + typedef double、~/forge-linevisc-fp64 を上書き
    # (上の lineE は下の lineF で上書きされて残っていない)
    "lineE_fp64": ("8e989dc688b7ccef70c444aecf06fb9ec406fa0a7f4d2ad0dc8158cef7851178", "~/forge-linevisc-fp64"),
    # 4fdc2c0e (診断のマスク FORGE_LVC_TERMS) + typedef double、~/forge-linevisc-fp64 を上書き
    # (上の lineF は下の lineG で上書きされて残っていない)
    "lineF_fp64": ("89ea94385435c16babde906465b58a54748d0e4b2d50b5edc661b46903ae19b0", "~/forge-linevisc-fp64"),
    # 846727be (マスクのビット 8 = 熱伝導の K の密度の列を外す) + typedef double、~/forge-linevisc-fp64 を上書き
    "lineG_fp64": ("c0b649052fdcf6092c7c73317563d77b570df60aca31e1bf709cf47586000c8c", "~/forge-linevisc-fp64"),
    # 5197e00e (診断の FORGE_DIAG_FACE_H_DOUBLE = TP の面エンタルピーを double で、plan time_integration-line-viscous-jacobian §6.15) + typedef double、
    # ~/forge-linevisc-fp64 を上書き (上の lineG は残っていない。切替なしの残差が lineG と同じことは run_0323_jph_a_q0 で確かめる)
    "lineH_fp64": ("561813b3564473420824ccca302b795bf9c344f69776f6381a43cb23e0ee6456", "~/forge-linevisc-fp64"),
    # 270f1d75 (opt-in の逆行列の保存 FORGE_LINE_INV と比較の経路の非有限・全ラインの η、plan time_integration-line-implicit-speed §6.2) + typedef double、~/forge-linevisc-fp64 を上書き (lineH は残っていない)
    "lineI_fp64": ("b467c3d7a1b32b595ac0ed8cc5db33f434329ef39865d727abfd11ef835b7661", "~/forge-linevisc-fp64"),
    # d1d0de7c (opt-in の float の Thomas FORGE_LINE_F32=1/2 と比較の経路の非有限の検査、plan time_integration-line-implicit-speed §6.4) + typedef double、~/forge-linevisc-fp64 を上書き (lineI は残っていない)
    "lineJ_fp64": ("3d045221b8677ca108af365a4012010c3e8f0e65e5d9e3ef53827d9cef963c0f", "~/forge-linevisc-fp64"),
    # 00da938b (opt-in の Thomas の配列の並べ替え FORGE_LINE_LAYOUT=1 とビット列の比較・因子の比較、plan time_integration-line-implicit-speed §6.7) + typedef double、~/forge-linevisc-fp64 を上書き (lineJ は残っていない)
    "lineK_fp64": ("c62eaf5a3910f9b1573ebbb1e740727d3c81b0a356a9352ae167ed67e93eec27", "~/forge-linevisc-fp64"),
    # f9be0c4f (opt-in の並べ替え + 連鎖の短縮 FORGE_LINE_LAYOUT=2、plan time_integration-line-implicit-speed §6.10) + typedef double、~/forge-linevisc-fp64 を上書き (lineK は残っていない)
    "lineL_fp64": ("e10a195dbbec9bd4b67650348985d68d6a5354e1f6884ab5e73f11974f085e3d", "~/forge-linevisc-fp64"),
}


def _use_alt_binary():
    import os
    k = os.environ.get("COLD_ALT_BINARY", "")
    if not k:
        return None
    if k not in ALT_BINARIES:
        raise SystemExit(f"COLD_ALT_BINARY={k} は登録されていない ({sorted(ALT_BINARIES)}) — 止める")
    sha, tree = ALT_BINARIES[k]
    CP.FORGE_SHA = sha
    CP.FP64_TREE = Path(tree).expanduser()
    return k


def limiter_refs(run: Path) -> dict:
    """forge_run.log の最後の limiterRoRef・limiterPRef・limiterARef (forge が貼るよう促す行) を読む。"""
    import re
    ref = {}
    for line in (run / "forge_run.log").read_text(errors="replace").splitlines():
        m = re.match(r"^\[limiter\]\s+(limiterRoRef|limiterPRef|limiterARef):\s*([0-9.eE+-]+)\s*$", line)
        if m:
            ref[m.group(1)] = m.group(2)
    if set(ref) != {"limiterRoRef", "limiterPRef", "limiterARef"}:
        raise SystemExit(f"{run.name}/forge_run.log にリミッタの基準値の 3 行がそろっていない: {ref}")
    return ref


def prep(src: Path, run: Path, steps: int, cfl: float, out_int: int, extra: list[str], ref_from: Path | None = None,
         line: str = "", isp: int | None = None, inner: int | None = None, conv: int | None = None, itj: int | None = None,
         cap: float | None = None, lvc: int | None = None, field_from: Path | None = None, relax: float | None = None) -> dict:
    NS.check_dry_env(False)
    binrec = CP.binary_record()
    if not NS.RUN_RE.match(run.name) or run.exists():
        raise SystemExit(f"{run} の名前が不正か既にある — 止める")
    srec = NS.jload(src / CP.RECORD)
    rs = NS.res_files(src)
    if not rs:
        raise SystemExit(f"{src} に res が無い")
    src_h5 = rs[-1]
    if field_from is not None:              # 設定は src、場は別の run の最終の res (同じ格子。切り戻し試験用)
        fr = NS.res_files(field_from)
        if not fr:
            raise SystemExit(f"{field_from} に res が無い")
        src_h5 = fr[-1]
    ys = NS.yaml_strict()
    ptext = (src / "solverConfig.yaml").read_text()
    ctext = ys.replace_scalars(ptext, {NS.NSTEP: str(int(steps)), NS.CFL: repr(float(cfl)), NS.CFLP: repr(float(cfl)),
                                       NS.OUTINT: str(int(out_int))})
    pcfg = ys.load(ptext)
    if extra:
        if "output" in pcfg:
            raise SystemExit("親の設定に output がある — extraFields の足し方を決めていないので止める")
        ctext = ctext.rstrip("\n") + "\noutput: {level: 1, extraFields: [" + ", ".join(extra) + "]}\n"
    refs = limiter_refs(ref_from) if ref_from is not None else {}
    if refs:
        if any(k in pcfg.get("space", {}) for k in refs):
            raise SystemExit("親の設定に既にリミッタの基準値がある — 止める")
        if ctext.count("space: {") != 1:
            raise SystemExit("space が 1 行のフロー形式でない — 基準値の足し方を決めていないので止める")
        ctext = ctext.replace("space: {", "space: {" + ", ".join(f"{k}: {v}" for k, v in refs.items()) + ", ")
    # line: "" = なし、"dir" = lineImplicit + lineDtDirectional、"only" = lineImplicit だけ (方向別の擬似 dt なし)
    line_keys = {"dir": {"lineImplicit": 1, "lineDtDirectional": 1}, "only": {"lineImplicit": 1},
                 "dirvisc": {"lineImplicit": 1, "lineDtDirectional": 1, "lineViscCoupling": 1}, "": {}}[line]
    if line_keys:                           # 壁法線のライン陰解法 + 方向別の擬似 dt (procedures/solver-settings.md「lineImplicit」)
        dt = pcfg["time"]["deltaT"]
        if int(dt.get("blockDPLUR", 0)) != 1 or int(pcfg["time"].get("timeIntegration", 0)) != 11 or int(dt.get("lowMachPrecond", 0)) >= 2:
            raise SystemExit("lineImplicit は timeIntegration 11 + blockDPLUR 1 + lowMachPrecond < 2 専用 — 止める")
        if any(k in dt for k in line_keys) or ctext.count("deltaT: {") != 1:
            raise SystemExit("deltaT に既にライン陰解法のキーがあるか、deltaT が 1 つのフロー形式でない — 止める")
        ctext = ctext.replace("deltaT: {", "deltaT: {" + ", ".join(f"{k}: {v}" for k, v in line_keys.items()) + ", ")
    if isp is not None:                     # 陰解法の行列の組立て・解法の精度 (0 = float、1 = double; 既定 0 は FP64 のビルドでも float)
        if "implicitSolvePrecision" in pcfg["time"]["deltaT"] or ctext.count("deltaT: {") != 1:
            raise SystemExit("deltaT に既に implicitSolvePrecision があるか、deltaT が 1 つのフロー形式でない — 止める")
        ctext = ctext.replace("deltaT: {", f"deltaT: {{implicitSolvePrecision: {int(isp)}, ")
    if itj is not None:                     # implicitThermalJacobian (plan time_integration-implicit-thermal-jacobian、別バイナリ)
        if "implicitThermalJacobian" in pcfg["time"]["deltaT"] or ctext.count("deltaT: {") != 1:
            raise SystemExit("deltaT に既に implicitThermalJacobian があるか、deltaT が 1 つのフロー形式でない — 止める")
        ctext = ctext.replace("deltaT: {", f"deltaT: {{implicitThermalJacobian: {int(itj)}, ")
    if cap is not None:                     # 方向別 dt の伸びの上限 R (lineDtDirectionalCap、line dir が要る)
        if not line_keys.get("lineDtDirectional"):
            raise SystemExit("--cap は --line dir と一緒に使う — 止める")
        ctext = ctext.replace("deltaT: {", f"deltaT: {{lineDtDirectionalCap: {float(cap)!r}, ")
    if lvc is not None:                     # lineViscCoupling (2 = 薄層の粘性 Jacobian、plan time_integration-line-viscous-jacobian、line が要る)
        if not line_keys.get("lineImplicit") or "lineViscCoupling" in pcfg["time"]["deltaT"]:
            raise SystemExit("--lvc は --line と一緒に使う (親に lineViscCoupling があっても止める) — 止める")
        ctext = ctext.replace("deltaT: {", f"deltaT: {{lineViscCoupling: {int(lvc)}, ")
    if relax is not None:                   # implicitRelax を明示して変える (plan time_integration-line-implicit-speed §5.1 #19・§6.13 の事前登録の run だけ)
        if ctext.count("deltaT: {") != 1 or "implicitRelax" not in pcfg["time"]["deltaT"]:
            raise SystemExit("deltaT が 1 つのフロー形式でないか親に implicitRelax が無い — 止める")
        import re as _re
        ctext, nsub = _re.subn(r"implicitRelax:\s*[0-9.eE+-]+", f"implicitRelax: {float(relax)!r}", ctext)
        if nsub != 1: raise SystemExit("implicitRelax の書き換えが 1 か所でない — 止める")
    one = {}
    if inner is not None:
        one[("time", "nStepInner")] = str(int(inner))
    if conv is not None:
        one[NS.CONVP] = str(int(conv))
    if one:                                  # 単因子の切り分け用 (nStepInner・convMethod)
        ctext = ys.replace_scalars(ctext, one)
    allowed = ({NS.NSTEP, NS.CFL, NS.CFLP, NS.OUTINT} | set(one) | ({("output",)} if extra else set())
               | {("space", k) for k in refs} | {("time", "deltaT", k) for k in line_keys}
               | ({("time", "deltaT", "implicitSolvePrecision")} if isp is not None else set())
               | ({("time", "deltaT", "implicitThermalJacobian")} if itj is not None else set())
               | ({("time", "deltaT", "lineDtDirectionalCap")} if cap is not None else set())
               | ({("time", "deltaT", "lineViscCoupling")} if lvc is not None else set())
               | ({("time", "deltaT", "implicitRelax")} if relax is not None else set()))
    diff = set(NS.MK.diff_paths(pcfg, ys.load(ctext)))
    if not diff <= allowed:
        raise SystemExit(f"許していない設定の差がある: {sorted(diff - allowed)} — 止める")
    if relax is None and ys.load(ctext)["time"]["deltaT"]["implicitRelax"] != pcfg["time"]["deltaT"]["implicitRelax"]:
        raise SystemExit("implicitRelax が変わった — 止める")
    if relax is not None and float(ys.load(ctext)["time"]["deltaT"]["implicitRelax"]) != float(relax):
        raise SystemExit("implicitRelax が指定どおりでない — 止める")
    run.mkdir(parents=True)
    for fn in NS.EXT_COPY + ("wall_repr.json", "bcondConfig.yaml", "species_meta.yaml"):
        if (src / fn).is_file():
            shutil.copy2(src / fn, run / fn)
    for p in sorted(src.glob("resolved_species_*.yaml")):
        shutil.copy2(p, run / p.name)
    (run / "solverConfig.yaml").write_text(ctext)
    cmd = [sys.executable, str(NS.TOOLS / "restart_field.py"), str(src_h5), str(run / "nozzle.h5"), "--dst-run", str(run), "--keep-src-dtype"]
    r = subprocess.run(cmd, capture_output=True, text=True, env=NS.runner()._ENV)
    (run / "restart_field.log").write_text(r.stdout + r.stderr)
    if r.returncode != 0 or "ビット一致" not in (r.stdout + r.stderr):
        print((r.stdout + r.stderr)[-3000:])
        raise SystemExit(f"restart_field がビット一致を確認していない (rc {r.returncode}) — 止める")
    info = NS.jload(run / "prepare_info.json")
    info.update(stages={"stages": "none", "ramp": None, "ramp_steps": 1000}, extends=src.name, restart_from=f"{src_h5.parent.name}/{src_h5.name}")
    NS.jdump(run / "prepare_info.json", info)
    keep = ("plan", "kind", "problem", "problem_sha256", "delta_r_csv", "euler_ref", "implicit_relax", "mesh_checks",
            "geometry_vs_production", "wall_thermal", "bcond_wall")
    rec = {**{k: srec[k] for k in keep if k in srec},
           "tool": "cold_cfl.py prep", "plan_item": "§5.1 #27", "created": NS.now(), "git_head": NS.git_head(), "binary": binrec,
           "stages": "none", "parent": src.name, "field_from": field_from.name if field_from is not None else None, "parent_res": src_h5.name, "parent_res_sha256": NS.sha256_file(src_h5),
           "ext_steps": int(steps), "cfl_main": float(cfl), "cfl_parent": srec.get("cfl_main"), "out_interval": int(out_int),
           "extra_fields": extra, "limiter_ref_from": ref_from.name if ref_from is not None else None, "limiter_refs": refs, "line_keys": line_keys, "implicit_solve_precision": isp, "n_step_inner": inner, "conv_method": conv, "implicit_thermal_jacobian": itj, "line_dt_directional_cap": cap, "line_visc_coupling": lvc, "implicit_relax_override": relax,
           "config_diff": sorted("/".join(p) for p in diff),
           "restart_field_tail": (r.stdout + r.stderr).strip().splitlines()[-1:], "nozzle_sha256_after_prep": NS.sha256_file(run / "nozzle.h5")}
    NS.jdump(run / CP.RECORD, rec)
    print(f"[cold_cfl prep] {run.name} ← {src.name}/{src_h5.name}: cfl {cfl}・{steps} step・出力 {out_int} ごと・extra {extra}; "
          f"{rec['restart_field_tail']}")
    return rec


if __name__ == "__main__":
    alt = _use_alt_binary()
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    p = sp.add_parser("run"); p.add_argument("run")
    p = sp.add_parser("prep"); p.add_argument("src"); p.add_argument("run")
    p.add_argument("--steps", type=int, required=True); p.add_argument("--cfl", type=float, required=True)
    p.add_argument("--out", type=int, default=5000); p.add_argument("--extra", default="")
    p.add_argument("--limiter-ref-from", default=None)
    p.add_argument("--line", choices=("dir", "only", "dirvisc"), default="",
                   help="dir = lineImplicit 1 + lineDtDirectional 1、only = lineImplicit 1 だけ、dirvisc = dir + lineViscCoupling 1")
    p.add_argument("--isp", type=int, choices=(0, 1), default=None, help="time.deltaT.implicitSolvePrecision を書く")
    p.add_argument("--inner", type=int, default=None, help="time.nStepInner を変える")
    p.add_argument("--conv", type=int, default=None, help="space.convMethod を変える (0 = 1 次)")
    p.add_argument("--itj", type=int, default=None, help="time.deltaT.implicitThermalJacobian を書く (別バイナリ COLD_ALT_BINARY が要る)")
    p.add_argument("--cap", type=float, default=None, help="time.deltaT.lineDtDirectionalCap を書く (別バイナリが要る)")
    p.add_argument("--field-from", default=None, help="場だけをこの run の最終の res から取る (設定は src、切り戻し試験用)")
    p.add_argument("--lvc", type=int, default=None, help="time.deltaT.lineViscCoupling を書く (2 = 薄層の粘性 Jacobian、別バイナリが要る)")
    p.add_argument("--relax", type=float, default=None, help="time.deltaT.implicitRelax を書き換える (§6.13 の事前登録の run だけ。既定は親のまま)")
    a = ap.parse_args()
    if a.cmd == "run":                      # cold_pair.run_one を (別バイナリの登録を効かせて) 呼ぶ
        sys.exit(CP.run_one(HERE / a.run))
    prep(HERE / a.src, HERE / a.run, a.steps, a.cfl, a.out, [s for s in a.extra.split(",") if s],
         HERE / a.limiter_ref_from if a.limiter_ref_from else None, a.line, a.isp, a.inner, a.conv, a.itj, a.cap, a.lvc, HERE / a.field_from if a.field_from else None, a.relax)
