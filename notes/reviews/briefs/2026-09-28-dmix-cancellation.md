# 諮問ブリーフ: 気相組成の拡散係数試験 G2 の FAIL — 許容の見直しか式の補数形化か (2026-09-28)

作業ツリー `/home/sano/work/forge-species` (`feature/species-transport`, HEAD 28897a88)。plan: `plans/active/condensation-two-phase-transport.md` §5.1 #3・§6 物性。
AGENTS.md エスカレーション 3 (事前に書いた許容との比較が FAIL) と 6 (直すなら `cuda_forge/thermo_d.cuh` の数値の振る舞いを変える)。

## 問い (1 つ)

§6 の「拡散係数も float 格納 ≤1e-5 で独立 double 参照と一致」が 1 点で FAIL した。**(A) 式を補数形 `(Σ_{j≠i} X_j)/Σ_{j≠i}(X_j/D_ij)` に変えて許容を維持する か、(B) 許容を「1−X_i ≥ 1e-2 で ≤1e-5、それ以外は ε₃₂/(1−X_i) の数倍」に書き直すか**。推奨を 1 つと、選んだ場合の合格条件を示してほしい。

## 観測事実

- 試験 `solver_density_cuda/tests/unit/test_transport_gas_phase.py` G2: 独立 double 参照との差 最大 2.93e-4 (合格 ≤1e-5)。判別 (総組成で参照した値との差) 0.31 で組成は正しく入っている。1−X_s ≥ 1e-2 の点では 6.97e-7。外れるのはほぼ純粋な点 (X_MIXDRY = 0.99990) だけで、誤差は 1.65 × 6e-8/(1−X_s) = 桁落ちの限界どおり (implementer 報告)。
- 式: `cuda_forge/thermo_d.cuh:566-577` `thermo_Dmix_species_f`: `(1.0f − X[i])/denom`, `denom = Σ_{j≠i} X_j/D_ij`。float。凝縮以外の全 TP 多成分 NS run が使う既存経路。
- 同じ試験の μ・λ (G1) は float 2.3e-7 で PASS、g=0 のビット一致 (G0) も PASS。
- ほぼ純粋な点の拡散流束は `ρD∇Y` で、その種の ∇Y は小さい (支配種)。ただし種の分子流束の補正 (Σ J = 0) に使われる。

## 当方の見立て (棄却してよい)

- (A) は数学的に同値で、分子と分母が同じ小さな和になるので桁落ちが無い。変わるのは X_i→1 の点の丸めだけで、既存 run の結果は丸め程度しか動かない見込み。ただし既存 TP 多成分 run のビット一致は崩れる (g=0 のビット一致試験 G0 は μ・λ のみで D を含まない)。
- (B) は結果を見てから合格条件を作ることになる (AGENTS.md が禁じる形)。

## 禁止事項 (厳守)

- ファイルを変更しない。**`*.log`, `residual_history.csv`, `res_*.h5`, `*.vtu`, `plans/README.md`, `trans.inp`・`thermo.inp` 全体を読まない**。下の範囲以外のファイル読みをしない (grep は可)。
- 推奨は 1 つに絞る。根拠は `ファイル:行` か本ブリーフで示す。

## 読んでよいもの (パスは `/home/sano/work/forge-species/` 基準)

- `sed -n '530,580p' solver_density_cuda/cuda_forge/thermo_d.cuh`
- `sed -n '180,320p' solver_density_cuda/cuda_forge/speciesTransport_d.cu`
- `solver_density_cuda/tests/unit/test_transport_gas_phase.py`
- `sed -n '100,140p' plans/active/condensation-two-phase-transport.md`
