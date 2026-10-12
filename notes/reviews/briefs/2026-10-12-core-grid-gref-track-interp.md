# 諮問ブリーフ: Gref の波束の追跡の結果の解釈と、仕組みについて今言えること

- **plan**: `plans/active/tooling-nozzle-core-grid.md` (§4.62〜§4.71、とくに §4.66・§4.68・§4.71)
- **date**: 2026-10-12
- **諮問の条件**: AGENTS.md「エスカレーション条件」3 (判定不能)・7 (解釈の確定の前)
- **ユーザの依頼**: 「なぜ流れ方向に細かいとライン陰解法の反復が遅くなるのか、という仕組みです。つきとめてくれ」

## 禁止事項 (厳守)

- ファイルを変更しない (read-only サンドボックスで動いている)。
- **`*.log`, `residual_history.csv`, `res_*.h5`, `*.vtu`, `*.f64`, `plans/README.md` を読まない** (巨大)。
- 下に列挙した `sed -n 'A,Bp' <file>` 以外のファイル読みをしない。grep は可。
- 両論併記で逃げず、**推奨は 1 つに絞る**。根拠は `ファイル:行` か下に示した数値で示す。

## 読んでよいもの

```
grep -n '^### 4.6[4-9]\|^### 4.7[01]' plans/active/tooling-nozzle-core-grid.md   (その節を sed で 40 行ずつ)
python3 -c "import json;j=json.load(open('case/45.isobutane_m6_d155/_band_ab/core_grid/twin_track_gref.json'));print({k:j[k] for k in ('center_x','xcorr_max','qualification','result','center_cols_per_step','packet_dt_over_dx')});print({k:{kk:v[kk] for kk in ('lambda','R2')} for k,v in j['follow'].items()}, j['fixed']['lambda'])"
python3 -c "import json;j=json.load(open('case/45.isobutane_m6_d155/_band_ab/core_grid/twin_track.json'));print({k:j[k] for k in ('center_x','qualification','result')});print({k:{kk:v[kk] for kk in ('lambda','R2')} for k,v in j['follow'].items()}, j['fixed']['lambda'])"
python3 -c "import json;[print(r['step'],round(r['island_D_T'],5),r['max_x']) for r in sorted((json.loads(l) for l in open('case/45.isobutane_m6_d155/_band_ab/core_grid/twin_main_cap50.jsonl')),key=lambda r:r['step'])]"
python3 -c "import json;[print(r['step'],round(r['island_D_T'],5),r['max_x']) for r in sorted((json.loads(l) for l in open('case/45.isobutane_m6_d155/_band_ab/core_grid/twin_main_grefc.jsonl')),key=lambda r:r['step'])]"
```

## 1. 事実 (§4.71)

- Gref の回し直し (§4.62 と同じ起点・摂動、3,000 step): 判定は判定不能 (最後の点 step 2750 だけ max|ΔT/T| 1.057e-3 で線形の資格を外れた。2500 までは ≤ 4.6e-4)。追従窓 λ 3.63e-3 /step (R² 1.000、窓幅で ±3 %)、固定窓 3.89e-3。
- 線形区間の波束の中心の移動: Gref 1.5e-3 r_t/step (0.0665 列/step)、Grefc 1.76e-3 r_t/step (0.037 列/step)。物理の距離ではほぼ同じ。
- 帯 (壁から 70〜80 節点目) の dt_local/Δx は両格子の同じ x で一致 → Δτ ∝ Δx。Grefc の擬似 dt は同じ場所で Gref の約 2 倍。
- 島の D_T は step 3000 まで両格子でほぼ同じ (0.23 対 0.21 K)。その後 Gref は 4500 まで育ち 11 K、Grefc は 3500〜3750 で 0.44〜0.49 K で頭打ち。最大点の位置では x −6〜−4.5 を Gref 約 1.2e-3、Grefc 約 1.8e-3 r_t/step で進む (非線形の区間)。
- 残差のプラトーの比 Grefc/Gref (同じ起点からの B0 の継続): ρ 0.028、ρE 0.016 (§4.37)。上限 1 (Gref) は増幅なし (§4.64)。

## 2. 判断してほしいこと

1. §4.70 の判定不能 (最後の点の線形の資格) の扱い。値は結果 A の条件 (λ ≥ 1e-3・R² ≥ 0.95・前後半とも正) を満たしているように見えるが、窓の変更は事後になる。
2. 以上から「流れ方向に細かいと止まる」について、ユーザに今言えること・言えないことを、短く (3〜5 文) 整理してほしい。私の読み: 増幅の速さ (1 step あたり) と線形区間の波束の物理の速さは両格子でほぼ同じ。違いは step 3000 以降 (x −6〜−4.5) の振る舞いで、Gref だけが育ち続ける。「列数に比例して滞在が延びる」はそのままでは成り立たない。
3. 調査をここで区切るか、次の 1 手があるなら 1 つ (事前の読み方つき)。
