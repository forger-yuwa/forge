# 諮問ブリーフ: 凝縮二相拡散 #4e (定常専用初版カーネル) の受入判定 (2026-10-02)

作業ツリー `/home/sano/work/forge-species` (HEAD d9fd68dd)。plan: `plans/active/condensation-two-phase-transport.md` §5.1 #4e。設計メモ `notes/investigations/2026-10-02-twophase-diffusion-kernel-design.md` §14 (実装の選択 §14.1、事前固定の受入基準 §14.2、実測、#1b 手順 §14.4)。前回諮問 `notes/reviews/2026-10-02-twophase-diffusion-4c-diagnose.md`。

## 問い

1. 受入 1 の **3 セル B・1000 更新の累積保存 3.14e-6 > 1e-6 (FAIL)** をどう扱うか。実装者の診断: 試験側の BE 解法から擬似時間項 V/Δτ を除くと同じ 6ε 床で 2.83e-7・上限 0 → 累積は試験側の BE 解法・停止則で決まり (#4a §6.2 と同じ現象)、面流束・更新写像は原因でない。初版は dual-time を拒否し、カーネルに物理時間の更新を積み重ねる経路はない。#4b の判断「相変化込みの累積保存を初版の必須ゲートにしない」は輸送単独の S3 型にも及ぶか。及ばないなら直し方 (試験の解法を直すか、カーネルに何か要るか) を 1 つ。
2. 蒸気残差の監視 `rms_roYv` は、格納状態から毎反復組み直した float32 残差の差であり、#4c で決めた「独立 float64 再評価」ではない。CFD の収束判定は rms の数桁低下で行うので丸め床より十分上、という前提でよいか。よくないなら最小の直し方。
3. 実装上の選択 (メモ §14.1) の妥当性: vapour/liquid/Q は点対角前処理 (DPLUR に載せない)、モーメントに φ_N δρ を足さない (再正規化係数で密度変化を二重計上するため)、θ_b 撤去で θ ≥ 0、定 Schmidt で ρ_g D = μ/Sc、Σz 正規化。
4. 以上で #4e を受け入れて、CFD の A/B (#1b、AWS) に進んでよいか。#1b の交絡 (B は前処理も変わる; 固定点は前処理に依らないので両者 PASS/STEADY なら場の比較は有効、収束速度は比較しない) の扱い。

## 実測 (要約; 詳細はメモ §14)

- 受入 2: T1 面流束 (3 種、D 比 1/3、15600 面) vs 同じ float 入力の double 参照で最大 0.175×許容。T3 の 1D 問題を GPU 版の拡散・更新に置換: DT_MAX 1 K で GPU 69/host 69 反復、0.01 K で 911/911 (θ<1 が 9759 セル·反復)、出口 g・T が表示桁まで一致、GPU 版も `accept()` PASS。
- 受入 3: キー OFF の新バイナリは旧と G0 4 構成バイト一致・test_transport_gpu 69 PASS。run_0482 入力の 2 step で res_0 の全 120 データセットがビット一致 (res_2 は旧同士と同じ 83 データセットが atomicAdd で不一致)。
- 受入 4: dual-time + キー ON で rc 1、res 書かれず。
- 受入 1: 構造 0.004/0.041×許容、面恒等式 ≤0.132、非負 2 セル min ρv 0、エネルギー接線 6.9e-4、3 セル 1 更新 4.92e-7 PASS、**1000 更新 3.14e-6 FAIL**。

## 禁止事項 (厳守)

- ファイルを変更しない。`*.log`, `residual_history.csv`, `res_*.h5`, `*.vtu`, `plans/README.md` を読まない。推奨は問いごとに 1 つ。根拠は `ファイル:行` か本ブリーフ。

## 読んでよいもの

- 設計メモ §14、plan、前回諮問、`solver_density_cuda/cuda_forge/twoPhaseDiffusion_d.cuh`、`speciesTransport_d.cu:2015-2164`、`condensationTransport_d.cu:670-837`、`main.cpp:2139-2176`、`solver_density_cuda/tests/unit/test_twophase_kernel.cu`
