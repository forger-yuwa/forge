# 諮問ブリーフ: B0 の停滞の原因の切り分け (続き) — 右辺は変わらないと分かった後の、反復の側の診断

- **plan**: `plans/active/tooling-nozzle-core-grid.md` (§4.33・§4.34)
- **date**: 2026-10-11
- **前回の諮問**: `notes/reviews/2026-10-11-core-grid-b0-cause-diagnose.md` — 「S1 (同じ保存量から上限 50/1 で更新前の残差ベクトルを比べる) を先に」
- **諮問の条件**: AGENTS.md「エスカレーション条件」4 (登録に無い run の前)

## 禁止事項 (厳守)

- ファイルを変更しない (read-only サンドボックスで動いている)。
- **`*.log`, `residual_history.csv`, `res_*.h5`, `*.vtu`, `plans/README.md` を読まない** (巨大)。
- 下に列挙した `sed -n 'A,Bp' <file>` 以外のファイル読みをしない。grep は可。
- 両論併記で逃げず、**推奨は 1 つに絞る**。根拠は `ファイル:行` か下に示した数値で示す。

## 読んでよいもの

```
sed -n '2318,2345p' solver_density_cuda/main.cpp            (組立 → 記録 → setDT → 更新の順)
grep -rn 'dt_local' --include=*.cu --include=*.cpp solver_density_cuda | grep -v tests | cut -c1-120
grep -rn 'lineDtDirectionalCap\|lineDtCap' --include=*.cu --include=*.cpp solver_density_cuda | grep -v 'tests\|input/'
sed -n '495,510p' solver_density_cuda/main.cpp              (rms_dq_* は常に 0 で記録している)
```

## 1. 観測事実

- **S1 はコードで答えが出た (実行していない)**: `lineDtDirectionalCap`・`lineDtDirectional` は `setDT_d.cu` の中だけで使われ (366〜371 行で setDT に渡すだけ)、`dt_local` を読むのは更新・左辺の側 (`update_d.cu`・`update.cpp`・`timeIntegration_d.cu`) と、この case で無効な凝縮・遷移・化学種の輸送の更新だけ。残差の組立 (`assembleResidual`) は `dt_local` を読まない。組立は setDT の前に行われる (main.cpp 2322〜2350)。AWS の実行バイナリ (`~/forge-fgeom7-fp64`) のソースでも `dt_local` を読むファイルの組は同じ、上限は `setDT_d.cu` の 5 か所だけ。したがって、この case では上限は右辺を変えない (前回の S1 の A/B は step 0 の残差が自明に同じになり、情報を持たない)。
- 残差の記録の `rms_dq_*` は常に 0 (main.cpp 503)。更新量の大きさは手元のデータに無い。
- 止まった状態の残差 (最後の 4000 step の中央値): B0 (`run_0608`・`run_0615`) rms_ro 2.1e-5・roe 13・roOmega 114、変動 (std/mean) 0.09〜0.20。上限 1 (`run_0616`) rms_ro 2.8e-7・roe 0.30・**roOmega 1.5e3**。point (`run_0614`) rms_ro 6.8e-8・roe 0.063・roOmega 130。Grefx の B0 (`run_0610`) rms_ro 1.2e-7・roOmega 42。
- 出力できる場: 保存量・原始量 (level 1)、`res_ro`・`res_roUx`・`res_roUy`・`res_roe`・`res_roK`・`res_roOmega`・`dt_local`・`limiter_*`・`omg_prod`・`omg_dest` など。`outStepStart`・`outStepInterval` で毎 step の書き出しもできる。

## 2. 期待値と出典

- 右辺が同じなら、B0 が R ≠ 0 のまま量が STEADY になる理由は反復の側にある。候補: (a1) 行き過ぎによる有界の振動 (増幅率 < −1 の モードが step ごとに符号を変えて揺れ続ける。平均がずれる)、(a2) 更新の切り詰め (正値の下限・EOS の床・k/ω の下限・SST の点陰解法のクリップ) で R ≠ 0 なのに更新が止まる節点がある (止まる位置が dt に依る)、(a3) 遅いモード (減衰率が 1 に近い)。

## 3. 判断してほしいこと

推奨を 1 つに絞り、事前の読み方まで示してほしい。候補 (縛られなくてよい):

- D1: `run_0608` の res_60000 から上限 50 (B0) と上限 1 で各 2003 step、`outStepStart 2000`・`outStepInterval 1` で 2000〜2003 の場を書き出す (`extraFields: [res_ro, res_roUx, res_roUy, res_roe, res_roK, res_roOmega, dt_local]`)。毎 step の更新 ΔU_k = U_{k} − U_{k−1} (k = 2001〜2003) を節点ごとに出し、(i) ΔU の符号の反転 (隣り合う ΔU の内積・相関が負の節点の割合と場所)、(ii) |R| が大きいのに |ΔU| が小さい節点 (切り詰めの疑い)、(iii) 壁からの節点番号 j・x・変数ごとの分布、を見る。再開の跳ねを避けるため 2000 step 回してから書き出す。数分で済む。
- 前回の推奨の S1 は上記のとおりコードで答えが出たので回さない、でよいか。
