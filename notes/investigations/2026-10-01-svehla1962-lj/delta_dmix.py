"""SERN 代表 m6_on ノズル入口組成での混合平均拡散係数 D_i,mix の legacy_v1 → 既定 [gri30, svehla1962] の相対変化。
ソルバと同じ補数形 D_i,mix = (Σ_{j≠i} X_j) / Σ_{j≠i} X_j/D_ij (`thermo_Dmix_species_f`)。SERN R9 の予測値用 (#14-L2, diagnostician 2026-10-02)。
usage: python3 delta_dmix.py (delta_table.py と同じ場所で)"""
import delta_table as DT

TS = (300.0, 1000.0, 2000.0, 2500.0)


def dmix(species, Y, lj_source, T):
    props = [DT.transport_lj(s, lj_source) for s in species]
    mw = [p[0] for p in props]
    n = [y / m for y, m in zip(Y, mw)]
    X = [v / sum(n) for v in n]
    out = []
    for i in range(len(species)):
        num = sum(X[j] for j in range(len(species)) if j != i)
        den = sum(X[j] / DT.dbin(props[i], props[j], T) for j in range(len(species)) if j != i)
        out.append(num / den)
    return out


def main():
    sp, Y = DT.SERN, DT.SERN_Y
    print("| 種 | X (入口) | " + " | ".join(f"ΔD_mix/D @{int(T)} K" for T in TS) + " |")
    print("| --- | --- | " + " | ".join("---" for _ in TS) + " |")
    old = {T: dmix(sp, Y, DT.OLD, T) for T in TS}
    new = {T: dmix(sp, Y, DT.NEW, T) for T in TS}
    mw = [DT.transport_lj(s, DT.OLD)[0] for s in sp]
    n = [y / m for y, m in zip(Y, mw)]
    X = [v / sum(n) for v in n]
    for i, s in enumerate(sp):
        print(f"| {s} | {X[i]:.3g} | " + " | ".join(f"{new[T][i] / old[T][i] - 1:+.2%}" for T in TS) + " |")


if __name__ == "__main__":
    main()
