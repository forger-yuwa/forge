# codex 諮問 (diagnose): te-gradient-hypothesis

- **brief**: [`notes/reviews/briefs/2026-10-08-te-gradient-hypothesis.md`](../../notes/reviews/briefs/2026-10-08-te-gradient-hypothesis.md)
- **plan**: [`plans/active/convection-zero-thickness-edge-reconstruction.md`](../../plans/active/convection-zero-thickness-edge-reconstruction.md)
- **date**: 2026-10-08
- **commit**: `6370dc0c` (feature/sern-design)
- **codex**: effort `xhigh`, 6.9 min, rc=0
- **結論**: **次は、冷却中の共通状態を固定し、実際のGPU勾配と同じ近傍・重み・打切りを使う独立倍精度LSQを照合する、局所作用素のA/Bを一つ行う。**
- **extra**: `plans/active/tooling-nozzle-sern-3d.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表

| 重大度 | 論点・採否 | 根拠と対案 |
|---|---|---|
| **Major** | 勾配の調査は**採用**。ただし「1 step＋`extraFields`」だけでは不足 | 通常出力には更新後の保存量と更新前の勾配が混在する。これは[既存診断の注意書き](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/limiter_d.cu:256)にも明記されている。**同じ残差組立時点の原始量・勾配・リミッタ・面値**を対応付ける。保存後の`roUx/ro`から速度を作り直して混ぜない。 |
| **Major** | 「共有壁節点が冷点のLSQへ直接入る」は**要再検証** | m6_onの記録では、冷点517159→内部点517198→壁点517237で、壁まで2接続ある。[帳簿ブリーフ:9](/home/sano/work/forge-sern-design/notes/reviews/briefs/2026-09-27-sern3d-r5h-ledger-result.md:9)。現行LSQは内部面で直接接続した節点値だけを読む。[実装:870](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/calcGradient_d.cu:870)。壁ピンの**間接的な影響**と直接の寄与を分ける。 |
| **Major** | 「辺方向の比＞3、寄与＞50%なら勾配異常」は**却下** | 現行の正規行列は`M = Σ d̂d̂ᵀ`で、辺長そのものは消える。[重みの実装:709](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/calcGradient_d.cu:709)。非線形場では微分と辺の差分が一致する必要もない。比と寄与率は**追跡対象を選ぶ指標**に留め、係数の再現・線形場の再現・実際の面速度を判定に使う。 |
| **Major** | 既存`check_lsq_gradient.py`の無修正流用は**却下** | ツールは[`CELLS/centCoords`と境界疑似点を使用](/home/sano/work/forge-sern-design/solver_density_cuda/tools/check_lsq_gradient.py:88)し、[通常の逆行列](/home/sano/work/forge-sern-design/solver_density_cuda/tools/check_lsq_gradient.py:198)を取る。現行nodeは[節点座標へ置換](/home/sano/work/forge-sern-design/solver_density_cuda/mesh/mesh.cpp:357)し、境界疑似点を除外し、[スペクトル打切り](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/calcGradient_d.cu:650)を使う。現行作用素に合わせた局所診断が必要。 |
| **Major** | GG比較は**補助として採用**。GGへの切替・勾配修正との同時投入は**保留** | GGには非一様格子で線形場を再現しない既知の制約がある。[仕様:86](/home/sano/work/forge-sern-design/methods/gradient.md:86)。また速度勾配は[粘性流束も読む](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/viscousFlux_d.cu:130)。GGとの差や温度上昇だけで、LSQの誤りや修正の妥当性を認定しない。 |
| **Minor** | 「速度10→1333 m/s」の表現は**訂正** | 記録は`Ux: 10.3→1332.9`。3成分から計算した速度の大きさは約195.5→1372.8 m/s。[帳簿ブリーフ:17](/home/sano/work/forge-sern-design/notes/reviews/briefs/2026-09-27-sern3d-r5h-ledger-result.md:17)。比は同じ成分で計算する。相手の`|U|`を`Ux`の差に使わない。 |

結論: **次は、冷却中の共通状態を固定し、実際のGPU勾配と同じ近傍・重み・打切りを使う独立倍精度LSQを照合する、局所作用素のA/Bを一つ行う。**

第 1 仮説: **角部近傍の急変する速度分布に対するLSQの当てはめと成分別リミッタの組合せが、大きな速度外挿を許している。** 確度: **中**。現時点では、勾配の計算誤りよりこちらを優先する。ただし特定の短辺が主因という部分は未確認。

根拠: 現行LSQは`gᵢ = Σ cᵢⱼ(uⱼ−uᵢ)`という直接隣接値の線形演算である。[実装:875](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/calcGradient_d.cu:875)。リミッタは全近傍の成分別上下限を使い、各辺の再構成増分を検査するため、**辺の両端値だけによる制約ではない**。[実装:323](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/limiter_d.cu:323)。報告されたm6_onでは、密度が節点値のまま、`Ux`だけ約1323 m/s増えている。この組合せはコード上成立するが、正しい面状態である保証はない。

反証条件: 時点と入力を揃えた独立計算がその外挿を再現せず、係数生成・適用の誤差が外挿量の大部分を説明すること。その場合、当てはめの性質に帰属する説明を退ける。

第 2 仮説: **係数生成・float32適用・打切り方向の処理に局所的な問題がある。** 確度: **低・未確認**。固有値、保持モード、実係数を取得するまで除外しない。

判別 A/B: **変更するのは勾配の評価実装だけ。A＝実際のGPU評価、B＝独立した倍精度評価。ソルバの数値設定と状態は変えない。**

1. **対象状態を固定する。**  
   主対象は`case/46.sern_design/run_1055_r7b_m10_A_c/`から既存の単面介入で再生した、517199が冷却中かつ床到達前のcall 10相当。保存された集計では、この区間の実更新は−7.755e3 J/kg、EOS側変更は+4.188e−4 J/kgである。[集計:16](/home/sano/work/forge-sern-design/notes/investigations/2026-10-05-sern-r7a/R5H_M10_LEDGER100.txt:16)。元の`run_1055`の初回だけでは、移動先が冷却する状態を調べられない。**比較自体は同じ状態の残差組立1回分で、追加の時間発展は不要。**

2. **同時刻の入力と作用素を採る。**  
   517160・517199・後縁壁点・接続先について、実節点座標、内部面の接続、境界所属、速度3成分、勾配9成分、リミッタ、実際に流束へ渡した面速度、`cᵢⱼ`を記録する。Bは同じfloat32格納値をdoubleへ変換して計算し、`M`の固有値に対する**0.01λmaxの打切り**も揃える。

3. **寄与は目的の辺へ射影して測る。**  
   辺`i→k`について、各隣接点の寄与を  
   `tⱼ = (dᵢₖ·cᵢⱼ)(uⱼ−uᵢ)`  
   とする。符号付き和、絶対値和、相殺率を記録する。重み`1/|d|²`の大きさだけでは、下流辺への寄与を判定しない。成分差が丸め誤差程度の辺では比を出さない。

4. **判定を事前固定する。**  
   定数場と独立な3方向の線形場でも係数を確認する。打切りのない点では`Σ cᵢⱼdᵢⱼᵀ ≈ I`、打切り点では保持部分空間への射影を期待値にする。A/Bの辺方向増分差は、係数丸め・差分・積和から求めた誤差上限の10倍以内を整合条件とし、**その許容幅自体が問題の外挿量の1%未満**であることを要求する。

**→ 結果A:** 係数・線形場・実場の照合が通り、大きな外挿も再現するなら、「計算精度や逆行列実装の誤りが主因」を棄却する。残るのは局所場と再構成方式の適合性であり、短辺の大きな寄与だけではLSQ修正へ切り替えない。  
**→ 結果B:** 時点・入力・保持モードを揃えても許容を超える不一致があれば、先に係数生成・適用を調べる。LSQの設計上の性質という説明を保留する。  
入力欠落、時点不一致、誤差幅が大きすぎる場合は**判定不能**とする。

やらない方がよいこと: **GGを正解扱いすること、比＞3だけでリミッタを追加すること、共有勾配配列に`w`を掛けること、勾配修正と速度一次化を同時に有効化すること。** 計算誤りが実証された場合はその修正を先に行い、既存の速度抑制案は診断用として保持する。

呼び出し側の前提への異議: **「高ARだからLSQが悪条件」「大きな外挿だから勾配が誤り」は成立しない。** 実装の式による合成確認では、直交6近傍の辺長比1000でも`M`の固有値はすべて2、条件数は1だった。実格子では方向配置と方向別の増幅を測る必要がある。

`w`については、**端の抽出・領域指定・来歴管理は再利用できる**。しかし[現在のカーネル処理](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:284)は面速度の補間であり、LSQの近傍や重みを変更しない。勾配側へ用途を広げるなら別の作用として設計し、粘性への影響も検証する必要がある。

不足情報: 対象runの生HDF5・勾配・実係数・実効設定・VERDICT原本は、このcheckoutにはない。さらに[run台帳:38](/home/sano/work/forge-sern-design/case/46.sern_design/README.md:38)では旧m6_onの`run_0971`等の全場は削除済みである。保持された`run_1050_r7b_x_lsw08/`を補助測定に使えても、旧帳簿の「約6倍」の再現とは呼べない。今回はコードと持ち帰り記録のレビューであり、収束・準定常性を新たに認定していない。

ファイル変更・forge起動なし。**plan未反映**。呼び出し側の反映先は`plans/active/convection-zero-thickness-edge-reconstruction.md` §4.1・§5.1 #2b・測定の判定条件。
