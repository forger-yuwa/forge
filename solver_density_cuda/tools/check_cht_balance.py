#!/usr/bin/env python3
r"""CHT の**収支ゲート G-cons** — 界面で受け渡した熱量が閉じているかを判定する。

plan `boundary-conjugate-heat-transfer.md` §4.3 / §6 の実体化。同一状態について

$$\varepsilon = \Big|\sum_i Q_{f,i} - Q_{\rm solid}\Big|,\qquad
  \text{分母} = \max\Big(\sum_i |Q_{f,i}|,\ Q_{\rm floor}\Big)$$

を作り、**$\varepsilon/\text{分母} \le$ `--tol-rel` かつ $\varepsilon \le$ `--tol-abs`** で合格とする。
**正味量で割らない**のは、正負が相殺する構成で分母が消えるため。$Q_{\rm floor}$ は
**ケースごとに計算前に登録**する絶対床 (`--q-floor`)。

流体側 $\sum Q_f$ は壁ダンプの `iface_q_eff` × 面積 (host・double)、
固体側 $Q_{\rm solid}$ は孔 Robin から出て行く熱 $\sum_k h_{c,k}\int(T-T_{c,k})\,ds$。

**`iface_q_eff` が NaN の節点があれば不合格**にする (適用範囲外の構成をそのまま通さない。
dual-time / 周期壁ノード / 軸対称は `methods/boundary.md` 参照)。

usage:
  python3 solver_density_cuda/tools/check_cht_balance.py <cht_run_dir> \
      --solid <solid.json> --solid-mode fem2d --phys-id 5 --phys-name wall \
      --q-floor 100.0 [--tol-rel 0.005] [--tol-abs 50.0]
"""
import argparse
import glob
import json
import sys
from pathlib import Path

import numpy as np
import h5py

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from solid_fem2d import Fem2DOperator          # noqa: E402
from cht_loop import build_fem2d, read_wall_dump, latest_wall_dump   # noqa: E402


def _dump_area(coords, faces):
    """壁ダンプの可視化要素から節点あたりの面積 (2D なら長さ) を集中させる。"""
    n = len(coords)
    a = np.zeros(n)
    for f in faces:
        idx = [int(k) for k in f if int(k) >= 0]
        if len(idx) == 2:
            L = float(np.linalg.norm(coords[idx[1]] - coords[idx[0]]))
            a[idx[0]] += 0.5*L; a[idx[1]] += 0.5*L
        elif len(idx) >= 3:
            p = coords[idx]
            ar = 0.0
            for k in range(1, len(idx)-1):
                ar += 0.5*float(np.linalg.norm(np.cross(p[k]-p[0], p[k+1]-p[0])))
            for k in idx:
                a[k] += ar/len(idx)
    return a


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--solid", default=None,
                    help="fem2d のとき必須 (固体 JSON)。local1d では run の solverConfig.yaml を読む")
    ap.add_argument("--solid-mode", default="fem2d", choices=["fem2d", "local1d"],
                    help="fem2d: `cht_loop` の外部連成 (固体 JSON)。"
                         "local1d: **ソルバ内連成** (`conjugate:` 節)。"
                         "固体側は $\\sum_i (T_{w,i}-T_b)/R_{\\rm tot}\\,A_i$ で作る")
    ap.add_argument("--phys-id", type=int, required=True)
    ap.add_argument("--phys-name", default="wall")
    ap.add_argument("--flux", default="q_eff",
                    help="収支に使う界面熱量 (既定 q_eff = 保存形。他は診断)")
    ap.add_argument("--q-floor", type=float, required=True,
                    help="規格化の**絶対床** [W or W/m]。ケースごとに計算前に登録する")
    ap.add_argument("--tol-rel", type=float, default=0.005)
    ap.add_argument("--tol-abs", type=float, default=None,
                    help="絶対許容 [W or W/m]。既定は q_floor")
    a = ap.parse_args()
    tol_abs = a.tol_abs if a.tol_abs is not None else a.q_floor

    run = Path(a.run_dir).resolve()
    its = sorted(glob.glob(str(run / "it_*")))
    src = Path(its[-1]) if its else run
    dump = latest_wall_dump(src, a.phys_name, a.phys_id)
    coords, faces, vals = read_wall_dump(dump)

    key = "iface_" + a.flux
    if key not in vals:
        sys.exit(f"{dump}: {key} が無い (output.interfaceDiag: 1 が要る)")
    q = np.asarray(vals[key], float)
    # 積分済み荷重 (ソルバが出す `iface_Qf_eff`)。`--flux q_eff` のときだけ正本にする。
    Qdirect = (np.asarray(vals["iface_Qf_eff"], float)
               if (a.flux == "q_eff" and "iface_Qf_eff" in vals) else None)
    n_nan = int(np.sum(~np.isfinite(q)))

    if a.solid_mode == "local1d":
        # ---- ソルバ内連成 (`conjugate: {mode: local1d, ...}`) ----
        # 固体は節点ごとに独立な直列抵抗 $R_{\rm tot}=t/k_s+R_{\rm back}$ なので、
        # 固体が持ち去る熱は $\sum_i (T_{w,i}-T_b)/R_{\rm tot}\cdot A_i$。
        # 界面の面積は壁ダンプの三角形/線分から作らず、**forge と同じ bplane 面積**が要るので
        # `iface_q_eff` と `Q` の比から逆算せず、mesh.h5 の面積を使う (下で読む)。
        import yaml
        cfgp = src / "solverConfig.yaml"
        cfg = yaml.safe_load(cfgp.read_text())
        cj = cfg.get("conjugate", {})
        t_s = float(cj["thickness"]); k_s = float(cj["k_solid"])
        Tb = float(cj["T_b"]); back = str(cj.get("back", "isothermal"))
        Rb = 0.0 if back == "isothermal" else float(cj.get("R_back", 0.0))
        Rtot = t_s/k_s + Rb
        Tw = np.asarray(vals["Ts"], float)
        with h5py.File(src / "mesh.h5", "r") as mf:
            pass
        # 面積: 壁ダンプの線分長 (2D) / 三角形面積 (3D) から作る
        area = _dump_area(coords, faces)
        Qf = q * area
        sum_Qf = float(np.sum(Qf[np.isfinite(Qf)]))
        sum_absQf = float(np.sum(np.abs(Qf[np.isfinite(Qf)])))
        Q_solid = float(np.sum((Tw - Tb)/Rtot * area))
        eps = abs(sum_Qf - Q_solid); denom = max(sum_absQf, a.q_floor); rel = eps/denom
        ok = (rel <= a.tol_rel) and (eps <= tol_abs) and (n_nan == 0)
        print(f"=== G-cons: {run.name}  ({dump.name}, flux={key}, local1d) ===")
        print(f"  固体: R_tot = t/k_s + R_back = {t_s}/{k_s} + {Rb} = {Rtot:.6f} m2K/W, T_b = {Tb} K")
        print(f"  流体側  sum Q_f     = {sum_Qf:14.6f} W   (sum |Q_f| = {sum_absQf:.6f})")
        print(f"  固体側  Q_solid     = {Q_solid:14.6f} W")
        print(f"  不釣合い eps        = {eps:14.6f} W")
        print(f"  規格化  eps/denom   = {100*rel:14.6f} %     (denom = max(sum|Q_f|, q_floor={a.q_floor}))")
        print(f"  NaN 節点            = {n_nan} / {len(q)}")
        print(f"  許容                : rel <= {100*a.tol_rel:.3f} % かつ abs <= {tol_abs}")
        print(f"VERDICT: {'PASS' if ok else 'FAIL'}")
        return 0 if ok else 1

    if a.solid is None:
        sys.exit("--solid-mode fem2d には --solid が要る")
    spec = json.loads(Path(a.solid).read_text())
    op, perm = build_fem2d(spec, coords)
    # **積分済み荷重があればそれを使う** (codex result 2 巡目 M4)。面積で割って集中辺長を
    # 掛け直すと角で +30 % 歪む。ソルバ内連成と同じ契約にする。
    if Qdirect is not None:
        Qf = Qdirect[perm]                      # [W/m]
        area = op.area
    else:
        qs = q[perm]
        area = op.area                          # 界面の集中長さ [m] (単位奥行き)
        Qf = qs * area                          # [W/m]
    # **実際に集計する配列**の非有限を検査する (codex result 4 巡目 M3)。
    # `iface_q_eff` だけ見ていたので、`iface_Qf_eff` の 1 点を NaN にしても「NaN 0/480」で PASS した。
    Qf = np.asarray(Qf, float)
    n_nan_load = int(np.count_nonzero(~np.isfinite(Qf)))
    if n_nan_load > 0:
        print(f"=== G-cons: {run.name} -> REFUSED ===")
        print(f"  集計に使う荷重に非有限が {n_nan_load} / {Qf.size} 点ある "
              f"(flux={key}{', 積分済み荷重' if Qdirect is not None else ''})")
        print("\nVERDICT: REFUSED")
        return 1
    sum_Qf = float(np.sum(Qf))
    sum_absQf = float(np.sum(np.abs(Qf)))

    # 固体側: 孔 Robin が持ち去る熱。
    # **ソルバ内連成では保存された固体状態 `conjugate_state_<physID>.h5` を使う** (3 巡目 M2)。
    # 初期入力の `wall_profile_*.csv` から復元すると、**比較相手が最新の流体荷重でなく初期条件**になる
    # (実測: C3X で 43430.6 W/m vs 保存状態からの 43473.5 W/m)。
    # **ソルバ内連成では、照合した固体 h5 からそのまま集計する** (codex result 4 巡目 M2)。
    # `--solid` の JSON / npz は**ハッシュ照合の対象外**なので、そこから作った作用素で計算すると
    # 孔の `h` を 10 % 変えるだけで固体放熱量が 43474 → 47821 W/m に変わっても拒否されない。
    # **ソルバ内連成 (`conjugate.mode: fem2d`) かどうかを run の config で判別する** (5 巡目 M1)。
    # チェックポイントの有無で経路を選ぶと、**無いときに初期壁温から復元した別状態で合格**してしまう
    # (実測: チェックポイントを隠すと 43474.1 → 43430.6 W/m・0.0984 % で PASS した)。
    insolver, solid_cfg = False, None
    cfgp = src / "solverConfig.yaml"
    if cfgp.exists():
        try:
            import yaml
            cj = (yaml.safe_load(cfgp.read_text()) or {}).get("conjugate")
            if isinstance(cj, dict) and str(cj.get("mode", "")) == "fem2d":
                insolver = True
                solid_cfg = cj.get("solid")
        except Exception:
            pass

    state = src / f"conjugate_state_{a.phys_id}.h5"
    solid_src = "conjugate_state"
    if insolver and not state.exists():
        print(f"=== G-cons: {run.name} -> REFUSED ===")
        print(f"  mode: fem2d なのに {state.name} が無い。"
              "初期壁温から復元した固体で代替しない (別状態を合格にしてしまう)")
        print("\nVERDICT: REFUSED")
        return 1
    if state.exists():
        with h5py.File(state, "r") as fs:
            u = np.asarray(fs["SOLID/T"][:], float)
            st_sha = str(fs.attrs.get("content_sha1", ""))
            st_step = int(fs.attrs["step"]) if "step" in fs.attrs else None
        # **必須属性が欠けていたら合格にしない** (5 巡目 M2: 欠落時に照合を飛ばしていた)
        if st_step is None or st_step < 0:
            print(f"=== G-cons: {run.name} -> REFUSED ===")
            print(f"  {state.name} に有効な step 属性が無い (時刻一致を確認できない)")
            print("\nVERDICT: REFUSED")
            return 1
        if not st_sha:
            print(f"=== G-cons: {run.name} -> REFUSED ===")
            print(f"  {state.name} に content_sha1 が無い (固体の同一性を確認できない)")
            print("\nVERDICT: REFUSED")
            return 1
        if len(u) != op.N:
            sys.exit(f"{state}: SOLID/T が {len(u)} 点 (固体は {op.N} 点)")
        # **並べ替えを戻す**。`SOLID/T` は固体 h5 の RCM 並べ替え後の順、`op` は npz の元の順。
        # ここを取り違えると孔の持ち去りが 66740 W/m (正 43473) のように桁違いに出る。
        # **固体ファイルは config の `conjugate.solid` から解決する** (固定名 `solid.h5` に頼らない)
        solid_h5 = (src / str(solid_cfg)) if solid_cfg else (src / "solid.h5")
        if not solid_h5.exists():
            print(f"=== G-cons: {run.name} -> REFUSED ===")
            print(f"  固体 HDF5 {solid_h5} が無い (conjugate.solid = {solid_cfg})")
            print("\nVERDICT: REFUSED")
            return 1
        with h5py.File(solid_h5, "r") as fh:
            perm_h5 = np.asarray(fh["MESH/PERM"][:], int)     # new[i] = old index
            mesh_sha = str(fh.attrs.get("content_sha1", ""))
        if not mesh_sha:
            print(f"=== G-cons: {run.name} -> REFUSED ===")
            print(f"  {solid_h5.name} に content_sha1 が無い (古い変換器で作った)")
            print("\nVERDICT: REFUSED")
            return 1
        if mesh_sha != st_sha:
            print(f"=== G-cons: {run.name} -> REFUSED ===")
            print(f"  content_sha1 が違う: 状態 {st_sha[:16]} / 固体 {mesh_sha[:16]} (別の固体)")
            print("\nVERDICT: REFUSED")
            return 1
        u_raw = u.copy()                     # h5 (RCM) 順
        u_old = np.empty_like(u)
        u_old[perm_h5] = u
        u = u_old
        # 孔 Robin の持ち去りを**固体 h5 の幾何と条件から直接**計算する (JSON を通さない)。
        with h5py.File(solid_h5, "r") as fh:
            xy = np.asarray(fh["MESH/COORD"][:], float)
            re_ = np.asarray(fh["ROBIN/EDGES"][:], int)
            rh = np.asarray(fh["ROBIN/H"][:], float)
            rtc = np.asarray(fh["ROBIN/TC"][:], float)
        u_h5 = np.asarray(u_raw, float)          # h5 (RCM) 順のまま使う
        L = np.hypot(xy[re_[:, 1], 0] - xy[re_[:, 0], 0], xy[re_[:, 1], 1] - xy[re_[:, 0], 1])
        Q_solid_direct = float(np.sum(rh * L * (0.5 * (u_h5[re_[:, 0]] + u_h5[re_[:, 1]]) - rtc)))
        solid_src = f"conjugate_state (step {st_step}) + solid.h5 の Robin 辺"
        direct = True
        # **時刻の一致を要求する** (4 巡目 M1: 壁が 960 step 古い組合せを PASS にしていた)。
        with h5py.File(dump, "r") as fw:
            w_abs = int(fw.attrs["step_abs"]) if "step_abs" in fw.attrs else None
        if w_abs is None:
            print(f"=== G-cons: {run.name} -> REFUSED ===")
            print(f"  壁ダンプ {dump.name} に step_abs 属性が無い (古い run)。"
                  "固体チェックポイントとの時刻一致を確認できないので合格にしない")
            print("\nVERDICT: REFUSED")
            return 1
        if w_abs != st_step:
            print(f"=== G-cons: {run.name} -> REFUSED ===")
            print(f"  時刻が違う: 壁ダンプ step_abs={w_abs} / 固体チェックポイント step={st_step}")
            print("\nVERDICT: REFUSED")
            return 1
    else:
        prof = sorted(glob.glob(str(src / f"wall_profile_{a.phys_id}.csv")))
        if not prof:
            sys.exit(f"{src}: wall_profile_{a.phys_id}.csv も conjugate_state も無い")
        d = np.loadtxt(prof[-1], skiprows=1)
        Tw_file = d[:, 3]
        # 壁プロファイルは固体界面節点の順 (cht_loop が op.coords で書いている)
        Tw = np.asarray(Tw_file if len(Tw_file) == op.n else Tw_file[perm], float)
        # `recover_interior` は `assemble` が張る Schur の中間量を使うので、必ず先に組む
        op.assemble(Tw); op.recover_interior(Tw)
        op.assemble(Tw)                      # 局所 k_s(T) で組み直してから最終復元
        u = op.recover_interior(Tw)
        solid_src = "wall_profile (外部ループ想定)"
        direct = False
    K = parts = None
    if not direct:
        K, parts = op.parts(T_ref=float(np.mean(u[op.iface])),
                        groups=[[(int(e[0]), int(e[1])) for e in np.load(spec["mesh_npz"])[k]]
                                for k in sorted([k for k in np.load(spec["mesh_npz"]).files
                                                 if k.startswith("hole")],
                                                key=lambda z: int(z[4:]))])
    if direct:
        Q_solid = Q_solid_direct
        K = parts = None
    else:
      Q_solid = 0.0
      for k, hole in enumerate(spec["holes"]):
        M, v = parts[k]
        h, Tc = float(hole["h"]), float(hole["T_c"])
        Q_solid += h * (float(np.asarray(M.dot(u)).sum()) - Tc * float(np.asarray(v).sum()))

    eps = abs(sum_Qf - Q_solid)
    denom = max(sum_absQf, a.q_floor)
    rel = eps / denom
    ok = (rel <= a.tol_rel) and (eps <= tol_abs) and (n_nan == 0)

    print(f"=== G-cons: {run.name}  ({dump.name}, flux={key}"
          f"{', 積分済み荷重' if Qdirect is not None else ''}) ===")
    print(f"  固体状態の出所      : {solid_src}")
    print(f"  流体側  sum Q_f     = {sum_Qf:14.4f} W/m   (sum |Q_f| = {sum_absQf:.4f})")
    print(f"  固体側  Q_solid     = {Q_solid:14.4f} W/m   (孔 Robin の持ち去り)")
    print(f"  不釣合い eps        = {eps:14.4f} W/m")
    print(f"  規格化  eps/denom   = {100*rel:14.4f} %     (denom = max(sum|Q_f|, q_floor={a.q_floor}))")
    print(f"  NaN 節点            = {n_nan} / {len(q)}"
          + ("   <-- 適用範囲外の構成 (methods/boundary.md)" if n_nan else ""))
    print(f"  許容                : rel <= {100*a.tol_rel:.3f} % かつ abs <= {tol_abs} W/m")
    print(f"VERDICT: {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
