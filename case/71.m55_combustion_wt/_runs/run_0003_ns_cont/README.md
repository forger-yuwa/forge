# run_0003_ns_cont — M5.5 の NS の延長を全部 double で (AWS A)

plan [`tooling-nozzle-card-local-mesh`](../../../../plans/active/tooling-nozzle-card-local-mesh.md) §5.1 #3d (ユーザ決定 2026-10-15
「どうせやるなら全部doubleかな」「thermoだけとは限らないよね」)。延長の条件は codex 諮問
[`2026-10-15-m55-ns-unsettled-2-diagnose.md`](../../../../notes/reviews/2026-10-15-m55-ns-unsettled-2-diagnose.md) の「延長へ進む場合」。
本体は AWS A の `~/forge-cardmesh/case/71.m55_combustion_wt/run_0003_ns_cont` (保存場・格子はそちらだけ)。このディレクトリは台本と記録の写し。

## 作り方

- 台本 [`cont_0003.py`](cont_0003.py) (`prepare` → `run`)。A の `~/forge-cardmesh` は f85e4f1b (git_dirty 0)。
- 出発点: `run_0002_ns/res_183000.h5` (sha256 `c29044a8…`) の保存量を、run_0002_ns の `nozzle.h5` の複製へ
  `restart_field.py --keep-src-dtype` で index コピー → **VERDICT: OK (9 量を移した、SRC とビット一致)**
  (ro・roUx・roUy・roUz・roe・roK・roOmega・roY0・roY1、`wall_dist` は格子のもの; 同一メッシュの検査 MESH/COORD 一致 441,600)。
  run_0002_ns の段の restart では roY が float32 の入力に丸められていた (相対 3〜6e−8) が、この run では float64 のまま。
- 前段は回さない: `runner_axismach.run_staged_ns(stages="none", stop_rule=True)` (本段だけ、停止規則の見張りごと)。
- バイナリ: FP64 のビルド `~/forge-integ-dflow2-fp64/solver_density_cuda/build/forge` (sha256 `3a9e347b…`、run_0002_ns と同じ)。

## run_0002_ns の本段との差 (これ以外は同じ; `prepare` が YAML を読み直して検査)

| 項目 | run_0002_ns | この run | 効き先 (forge_run.log の起動の行) |
| --- | --- | --- | --- |
| `physProp.thermoFloat` | 1 (既定: float Newton + double 研磨) | **0** | `'thermoFloat' in 'physProp': 0 (double Newton)` |
| `time.deltaT.implicitSolvePrecision` | 0 (既定: block-DPLUR の solve を float) | **1** (条件 timeIntegration 11・blockDPLUR 1・lowMachPrecond 0 を確認) | `'implicitSolvePrecision' in 'time.deltaT': 1` |
| 環境変数 `FORGE_TRANSPORT_TABLE` | 無し (輸送物性を float の表で引く) | **0** (double の `transport_mix_Y`) | `[thermo_d] transport tables disabled (FORGE_TRANSPORT_TABLE=0): double evaluation` |
| 環境変数 `FORGE_DIAG_FACE_H_DOUBLE` | 無し (SLAU の面の全エンタルピーを float の写しで) | **1** (double の `thermo_h_mix`) | `[DIAG] FORGE_DIAG_FACE_H_DOUBLE: …` |
| `time.outStepInterval` | 1,000 | **2,000** | — |
| 停止規則 | 既定の要求 (explicit false)、euler_ref はディレクトリだけ、帯は本段の res_1000 から | **明示** (explicit true)・窓 2 万・保持 4 万・上限 200,000 (= 追加分)・出力 2,000、euler_ref は `run_0001_euler/res_12000.h5` (sha256 `93cb88b0…`) に固定、帯はこの run の最初の出力から | `stop_rule.json` の settings |

環境変数は `RUN_PROVENANCE.txt` に記録される (`run_case.sh` の記録の対象に 2 つを足した、f85e4f1b)。

**残る float** (この run でも float のまま): 化学種の拡散の D_s・h_s・X の float の写し (`speciesTransport_d.cu:240-356`)、
リミッタの長さの `sqrtf` (`limiter_d.cu:398`)、XMF の型の宣言 (表示だけ)。

判定区間はこの run の本段だけ (設定が違うので run_0002_ns の残差とは連結しない; restart の直後の跳ねから測ることになる)。
θ_r はこの run の帯で取る別の系列。追加 20 万 step は予算の上限で、到達の予測ではない (2 万 step ごとに再評価)。

## 入口の面の実効の流入とエネルギー流束 (forge 0 step、run_0002_ns の res_183000)

台本 [`inlet_face_check.py`](inlet_face_check.py)、原文 [`inlet_face_check.json`](inlet_face_check.json)・[`inlet_face_check.stdout.txt`](inlet_face_check.stdout.txt)。
`inlet_Pressure_d` (境界の状態) と `convectiveFlux_boundary_d` (境界の流束) の式を写し、節点の状態は A = 出力の原始量、B = 保存量から
thermoFloat 0 の経路で作り直したもの。入口の面 92 (列 i = 0 の全節点、外向きの法線は全て −x)。違反の目安 (事前に置いたもの):
v1 逆流 (Un_c > 0)・v2 流入の面で |Ht_R − h0_in| > 1e−6 h0_in・v3 合計のエネルギー流束と ṁ·h0_in の差 > 1e−6・v4 mach ≥ 1。

| 量 | A | B |
| --- | --- | --- |
| 流入の面 / 逆流 / U = 0 | 91 / 0 / 1 (壁との角の節点、壁の Dirichlet で U = 0 → 逆流の分岐で ṁ = 0) | 同じ |
| 流入の面の \|Ht_R − h0_in\| / h0_in の最大 | 9.7e−14 | 9.7e−14 |
| 面の流入の合計 ṁ (2π 倍) | 121.572701 kg/s | 121.572699 kg/s |
| 合計のエネルギー流束 − ṁ·h0_in (相対) | −9.5e−14 | −9.5e−14 |
| 節点の値の ṁ_in (停止規則と同じ式) / 面の ṁ | 121.573953 / 0.99999 | 121.573949 / 0.99999 |
| mach_new の最大 | 0.0658 | 0.0658 |
| 違反 v1〜v4 | なし | なし |

T0 > Tt + 1 K の入口の 9 点 (壁から 19〜27 点、r 0.29988〜0.29997 m): 全て流入の面 (Un_c −44.7〜−44.8 m/s、mach 0.0653〜0.0655)。
境界の状態は Ts_b 1,249.25 K・Ps_b 3,490,397〜3,490,447 Pa (節点の静温 1,250.39〜1,252.62 K・静圧 3,490,377.6〜3,490,381.4 Pa)、
面から入るエンタルピーは h0_in ちょうど (差 0)。9 点の面の流入の割合は合わせて 0.07 %。つまり入口の面の流束は Tt を超える
エンタルピーを入れていない (節点の静温の超過がどこから来たかはこの検査の範囲外)。
壁の近く (壁から 1〜5 点) の流入の面は節点が遅く冷たい (T 392〜588 K、mach 0.006〜0.026) ので境界の状態が淀みに近く、
Ps_b が節点の静圧より最大 +0.27 % 高く Ts_b ≈ 1,250 K になる。これらの面の面積の割合は 1.4e−5、流入の割合は 3.5e−6
(W で重み付けした Ps_b − P の平均は相対 3.4e−7)。

## 投入 (2026-10-15 00:24:13 UTC、implementer)

```
cd ~/forge-cardmesh/design
export FORGE_BIN=/home/ubuntu/forge-integ-dflow2-fp64/solver_density_cuda/build/forge FORGE_CUDA_BLOCKSIZE=128 \
       FORGE_TRANSPORT_TABLE=0 FORGE_DIAG_FACE_H_DOUBLE=1 LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial
nohup python3 ../case/71.m55_combustion_wt/_runs/run_0003_ns_cont/cont_0003.py run \
      --run-dir ../case/71.m55_combustion_wt/run_0003_ns_cont \
      > ../case/71.m55_combustion_wt/run_0003_ns_cont.driver.log 2>&1 < /dev/null &
```

- PID (A): ドライバ (python、停止規則の見張りを含む) 3320、`run_case.sh` 3329、forge 3377。
- 起動の記録: `RUN_PROVENANCE.txt` に `FORGE_CUDA_BLOCKSIZE=128`・`FORGE_TRANSPORT_TABLE=0`・`FORGE_DIAG_FACE_H_DOUBLE=1`・`FORGE_BIN`、
  forge の sha256 `3a9e347b…`、git_head f85e4f1b (dirty 0)。`forge_run.head300.log` に 4 つの切り替えの起動の行。
  `MAIN_STAGE_STARTED` は stages none・前段なし・本段の設定 sha256 `29d9390b…`・初期の場 `nozzle.h5` sha256 `b883877f…`。
- 序盤の NaN: 残差 (先頭 43,437 行 = step 0〜14,478、内反復を含む全 rms 列) に非有限 0、`res_nan_*` 無し、見張りの行
  (res_2000〜10000) の非有限 0。step 0 の rms_ro 1.70e−6 (run_0002_ns の最後 1.56e−6)、step 14,478 で 6.45e−7。
- 1 step の時間: 26.94 ms/step (起動から 342 s、step 12,649; run_0002_ns の本段は平均 14.4 ms/step (2,648 s / 184,000 step) なので約 1.9 倍)。
  上限 200,000 step まで回ると forge は約 90 分 (5,390 s)。停止規則の見張りは最初の出力 res_2000 で帯を決めるのに 226 s かかった。
  出力ごとの δ の抽出 (記録だけ) は 1 本数分なので、forge の終了の後も最後の 2 窓ぶんの抽出を待つ (上限 2 時間)。
- 容量 (A の `/`): 投入の前の空き 16,356,311,040 B → 格子の複製の後 16,049,217,536 B。この run の見込みは全場 100 枚 + res_0
  (1 枚 100,605,412 B) で約 10.2 GB、壁の出力 1,170,416 B × 100、出口の出力、残差履歴 (約 0.12 GB)・forge の log を足して約 10.4 GB
  → 終わりの空きは約 5.6 GB (他のセッションが使わない場合)。
- 停止規則の帯: res_2000 (この run の最初の出力) × `run_0001_euler/res_12000.h5` (sha256 照合済み)、`stop_rule_band/`。
