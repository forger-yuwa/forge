#!/usr/bin/env python3
"""F1 の読み出し試験 (plan boundary-node-periodic-gradient-fix §4.2 / codex m3、§5.1 #5a)。

`sstF1` は出力変数に無い (`variables.hpp` output_cellValNames) ので、輸送が読む F1 を残差から判定する。
TGV 32^3 (三重周期・壁なし)、U=0、ρ=1、粘性 0 (μ=0 → 拡散係数は σ·μ_t だけ)、k = 1 + 0.3 sin x cos y + 0.2 cos z、ω = 64 k。

- wall_dist ≡ ∞ (h5 の 0 をソルバが 1e30 で埋める): arg1 = 0 → **計算値 F1 = 0 (厳密)**。σ_k = 1.0、σ_ω = 0.856。
- wall_dist ≡ 1e-6: arg1 ≫ 1 → 計算値 F1 = 1 (厳密)。σ は k-ω 側 0.85 / 0.5 と同じ。
sstSigmaBlend 0 (σ = 0.85/0.5 固定) と 1 (F1 ブレンド) の対を 2 step 回し、res_1 (1 回目の残差組立)・res_2 (2 回目) の
res_roK / res_roOmega を比べる。源項の F1 は ransSource が自分で計算するので対の両方で同じ (差は輸送の σ だけ)。

判定:
- 対 ∞ (計算値 F1=0): 輸送が初期値 1 (旧 buildScalarDescs の上書き) を読むなら blend 1 = blend 0 と同じ σ で残差は一致する。
  計算値 0 を読むなら差が出て、ω = 64k で μ_t が両式共通なので Δres_ω/Δres_k = 64·(0.856−0.5)/(1.0−0.85) = 151.89 になる
  (拡散は場について線形、64 は 2 の冪なので離散でも同じ比)。1・2 回目とも中央値がこの比の ±1 % に入ること。
- 対 1e-6 (計算値 F1=1): 両者は同じ σ → 差は atomicAdd の順序の床 (同一設定の再実行との差) と同程度。
"""
import os
import sys

import h5py
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gharness as G   # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
SCR = "/tmp/claude-1000/-home-sano-work-forge/4b0c8643-66fb-4fde-8878-7c8a5104c060/scratchpad"
SRC = os.path.join(SCR, "lsqseam_m1", "Taylor-Green.h5")
BC = os.path.join(HERE, "bcondConfig.yaml")
RATIO = 64.0 * (0.856 - 0.5) / (1.0 - 0.85)


def make(tag, blend, wd):
    run = os.path.join(SCR, f"g_harness_f1_{tag}")
    cfg = G.base_solver_cfg(sst=True, nstep=2)
    cfg["turbulence"]["sstSigmaBlend"] = blend
    h5 = G.make_run(run, SRC, cfg, BC)
    msh = G.Mesh(h5, os.path.join(run, "bcondConfig.yaml"))
    x = msh.xyz.astype(np.float64)[msh.root]            # 周期像は root の座標で評価 (group 内で同値)
    k = 1.0 + 0.3 * np.sin(x[:, 0]) * np.cos(x[:, 1]) + 0.2 * np.cos(x[:, 2])
    k = k.astype(np.float32).astype(np.float64)
    om = 64.0 * k                                        # 2 の冪倍なので f32 でも ω = 64k が厳密
    z = 0.0 * k
    f = G.prim_state_fields(msh, 1.0, z, z, z, 101325.0, k=k, om=om)
    f["wall_dist"] = np.full(msh.nCells, wd)
    G.bake(h5, f)
    G.run_forge(run)
    return run


def res(run, step):
    with h5py.File(os.path.join(run, f"res_{step}.h5"), "r") as f:
        return f["VALUE/res_roK"][()].astype(np.float64), f["VALUE/res_roOmega"][()].astype(np.float64)


def main():
    runs = {}
    for tag, blend, wd in (("inf_b0", 0, 0.0), ("inf_b1", 1, 0.0), ("inf_b1_repeat", 1, 0.0),
                           ("tiny_b0", 0, 1e-6), ("tiny_b1", 1, 1e-6)):
        runs[tag] = make(tag, blend, wd)
    out = ["# F1 の読み出し試験 (plan boundary-node-periodic-gradient-fix §4.2、codex m3)",
           f"harness revision: {G.git_rev()}", G.provenance(runs["inf_b1"]).rstrip(), "",
           f"期待比 Δres_ω/Δres_k = 64·(0.856−0.5)/(1.0−0.85) = {RATIO:.3f} (計算値 F1=0 を読んだとき)", "",
           "| 残差組立 | 対 | max|Δres_k| / max|res_k| | max|Δres_ω| / max|res_ω| | Δres_ω/Δres_k 中央値 [p5, p95] | 床 (再実行との差) / max|res_k| | 判定 |",
           "| --- | --- | --- | --- | --- | --- | --- |"]
    ok = True
    for step in (1, 2):
        k0, w0 = res(runs["inf_b0"], step); k1, w1 = res(runs["inf_b1"], step); kr, wr = res(runs["inf_b1_repeat"], step)
        dk, dw = k1 - k0, w1 - w0
        floor = np.abs(kr - k1).max() / np.abs(k1).max()
        sel = np.abs(dk) > 1e-3 * np.abs(dk).max()
        r = dw[sel] / dk[sel]
        med, p5, p95 = np.median(r), np.percentile(r, 5), np.percentile(r, 95)
        good = np.abs(dk).max() / np.abs(k0).max() > 100 * floor and abs(med / RATIO - 1) <= 0.01
        ok &= good
        out.append(f"| {step} 回目 | wall_dist ∞ (計算値 F1=0): blend1 − blend0 | {np.abs(dk).max() / np.abs(k0).max():.3e} | "
                   f"{np.abs(dw).max() / np.abs(w0).max():.3e} | {med:.3f} [{p5:.3f}, {p95:.3f}] | {floor:.1e} | {'ok' if good else 'NG'} |")
        t0k, t0w = res(runs["tiny_b0"], step); t1k, t1w = res(runs["tiny_b1"], step)
        e = np.abs(t1k - t0k).max() / np.abs(t0k).max()
        ew = np.abs(t1w - t0w).max() / np.abs(t0w).max()
        good2 = e <= 10 * max(floor, 1e-7) and ew <= 10 * max(np.abs(wr - w1).max() / np.abs(w1).max(), 1e-7)
        ok &= good2
        out.append(f"| {step} 回目 | wall_dist 1e-6 (計算値 F1=1): blend1 − blend0 | {e:.1e} | {ew:.1e} | - | {floor:.1e} | "
                   f"{'ok' if good2 else 'NG'} (床の 10 倍以内) |")
    out.append("")
    out.append("run: " + ", ".join(f"{k} = {os.path.basename(v)}" for k, v in runs.items()))
    out.append(f"VERDICT F1 読み出し (1・2 回目の残差組立で輸送が計算値 F1 を読む): {'PASS' if ok else 'FAIL'}")
    txt = "\n".join(out) + "\n"
    p = os.path.join(HERE, "G_f1_read.txt")
    open(p, "w").write(txt)
    print(txt)
    print("->", p)


if __name__ == "__main__":
    main()
