# 諮問ブリーフ: 凝縮二相拡散 #4a の結果と #4 カーネル着手の可否 (2026-10-02)

作業ツリー `/home/sano/work/forge-species`。plan: `plans/active/condensation-two-phase-transport.md` §5.1 #4・#4a。前回諮問: `notes/reviews/2026-10-02-twophase-diffusion-kernel-diagnose.md`。設計メモ: `notes/investigations/2026-10-02-twophase-diffusion-kernel-design.md` (§6 結果、§6.1 B の定義)。

**問い**: (Q1) §6.1 の B の limiter・commit の定義で #4 のカーネル実装に進んでよいか。(Q2) float32 の累積保存 ≤1e-6・上限到達なしを、
停止則で満たすのか (6ε は窓が狭い)、残差の持ち越し (離散式の変更) を採るのか、許容を見直すのか。(Q3) 補正項の種間結合で緩和なし点対角が収束しない問題を、
緩和で済ませるか、前処理に入れるか。

**観測事実** (ホストのみ; HEAD `9fd131fc` + 未 commit のハーネス変更、`python3 solver_density_cuda/tests/unit/test_twophase_diffusion_harness.py`、PASS 42 / FAIL 2):
- S7: A (P1) ΔρY_w −3.3333302e-3 (f32) / −3.3333333e-3 (f64)。B は全増分・補正が厳密に 0 (f32/f64)。
- S3 float32 1000 更新: 床 16/8/6/4ε で 9.98e-7 / 6.41e-7 / 3.18e-7 / 2.83e-7 (上限到達なし)、持ち越し 7.9e-8。
- S5 float32 1000 更新: 床 16ε 6.84e-6、8ε 1.25e-6、6ε 5.69e-7、4ε 5.00e-7 (上限到達 113 回 = FAIL)、16ε + 持ち越し 1.53e-8 (上限 0)。
  残差の和を丸めの大きさまで下げる停止 (`sum_q` = 4/2/1) は 5.25e-6 / 5.38e-6 / 1.60e-6 (q=2,1 は上限到達あり)。
- S1 (dx 1 mm): 入力量子化 3.082e-10 (許容 3.4e-9)、演算誤差 5.3e-10 (許容 9.2e-9)。旧案は許容の 18 倍。S1a (D_k 固定) も同様に PASS。
- 反例 (両セル g=0, D 比 100, 風上 z): 点対角 1 回 ρv_L +3.2442e-4。緩和なし反復は 5000 回で停滞 (max|r_k| 2.9 から下がらない、2 周期)。
  緩和 0.5 で 232–247 回で収束し Newton の BE 解 (ρv = [1.195754e-3, 8.804246e-3]) と一致。

**期待値と出典**: plan §5.1 #4a の合格条件 (S7: A −0.003333・B 全 0、S3/S5 float32: 累積 ≤1e-6・上限到達なし)。diagnose 2026-10-02 の採否表。

**実施済みの操作**: 本メモ §6 (ハーネスの変更点・B の定義・実測)。カーネル・plan・methods は未変更 (methods §7c の「未確定」は diagnose の裁定をまだ反映していない)。

**仮説**:
- H1: B の limiter・commit は固定点を保つ (S7)。非固定点の総水分の増減は前処理の性質で、固定点の保存は S3/S5 の float64 で成立している。
- H2: float32 の累積保存は、停止の床では格納の丸めの偏りを抑えられない。持ち越しは丸めの偏りを次ステップで打ち消すので効く (ただし 1 例・閉じた箱の拡散単独での結果)。
- H3: 反例の非収束は補正項 $-z_c\Sigma j^0$ の種間結合が点対角に入らないことによる (緩和で収束、解は存在)。D 比の大きい組成 (H2 を含む燃焼ガス) で実用上どれだけ効くかは未測定。

**呼び出し側の前提で疑わしいもの**: 「解き切った」= 上限到達なしを、反復回数に上限のある擬似時間 (定常) とサブ反復数固定の dual-time にどう写すか。
持ち越しを採ると dual-time の BDF 履歴 (FCT の流束形) とどう整合させるか。

## 禁止事項 (厳守)

- ファイルを変更しない。`*.log`, `residual_history.csv`, `res_*.h5`, `*.vtu`, `plans/README.md` を読まない。推奨は問いごとに 1 つ。根拠は `ファイル:行` か本ブリーフ。

## 読んでよいもの

- `notes/investigations/2026-10-02-twophase-diffusion-kernel-design.md` 全体、plan 全体、前回諮問記録、`solver_density_cuda/tests/unit/test_twophase_diffusion_harness.py`、メモが挙げた `ファイル:行` の範囲
