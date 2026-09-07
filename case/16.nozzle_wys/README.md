# case/16.nozzle_wys — Wyslouzil 超音速ノズル

Wyslouzil et al., *J. Chem. Phys.* **113**, 7317 (2000)
("Binary condensation in a supersonic nozzle") の超音速 Laval ノズル。
凝縮実験用の矩形断面ノズルで、forge による 2D / 3D・非粘性 / 粘性 / 乱流の計算を行う。

## ジオメトリ (`mesh/nozzle_H.geo`, 単位 cm → ScalingFactor 0.01 で m)

論文との寸法照合 (おおむね一致):

| 量 | `.geo` 由来 | 論文 |
| --- | --- | --- |
| 入口全高 | 2×1.236 cm ≈ **24.7 mm** | flow straightener 25.4 mm |
| 断面幅 (z 押し出し) | 1.259 cm = **12.59 mm** | 12.7 mm |
| スロート半高さ | ≈ 0.225 cm = **2.25 mm** | — |
| 出口/スロート面積比 | **1.69** → 等エントロピー M≈2.0 | 設計 M≈1.95 (γ=1.4) |

矩形断面のため 2D 平面近似が妥当 (ユーザ確認済)。

> **⚠️ 重要 (2026-06): `mesh/nozzle_H.geo` は Fig.3 の実験ノズルと別形状**だった
> (`13.nozzle_H` とバイト同一の別 Wyslouzil 論文ノズル。スロート半 2.25mm・A/A\*=1.69 で過発散)。
> JCP 113,7317 (2000) Fig.3 の正しいノズルは **スロート全高 5mm・収束 38mm・発散 95mm 直線壁・
> 壁間 1.8°・A/A\*≈1.58** で `mesh/make_nozzle_fig3.py` が生成する (`nozzle_fig3_2d/3d`)。
> **凝縮検証は run_0047 以降 (修正形状) を使うこと**。run_0001〜0046 (nozzle_H) は旧形状で、
> dry が実験より過膨張・凝縮も相殺で「合って見えていた」だけ。詳細は下記「## 修正形状 (Fig.3) run」。
> **さらに 2026-08-19: `make_nozzle_fig3.py` の Fig.3 形状 (smoothstep 収縮 + throat から直線膨張) もユーザ提示の壁式と
> 不一致** (収縮部 |Δy| 最大 1.56 mm、throat 直後の 3 次ブレンド無し、A/A\* 1.597 vs 1.580)。**今後は `mesh/make_nozzle_user.py`
> (`nozzle_user_profile.py` の式) を正とし、run_0189 以降を使うこと** (下記「## ユーザ指定形状 (2026-08-19)」)。

## 条件

- 入口: `inlet_Pressure`, Pt=101325 Pa, Tt=293.15 K (亜音速圧力入口)。
- 出口: `outlet_statPress`, Ps=2000 Pa (超音速流出のため実質外挿)、逆流用 Pt/Tt 設定。
- 気体: γ=1.4, cp=1039 J/kgK (N2 相当)。時間積分は陰解法 (block-DPLUR, `timeIntegration:11`)。
- 初期場: `initial: nozzle_wys` (一様 M≈0.3, Pt/Tt そろえ。`setInitial.hpp` に追加)。
- 段階起動: 非粘性(slip) → 粘性層流(no-slip) → 乱流(SST) を引き継ぎ計算で段階的に。

### 2D / 3D メッシュの壁の扱い

- **2D** (`nozzle_2d*.geo`, 1層スラブ): Euler は全面 slip。粘性/乱流は輪郭壁=no-slip、
  front/back=slip (対称面) とするため `nozzle_2d_visc.geo` で front/back を別 physID(4) に分離。
- **3D** (`nozzle_3d.geo`, 12層): 矩形ダクトの4壁すべてが実壁。Euler は slip、粘性/乱流は全壁 no-slip。
  ※ **乱流(SST)では mesh を no-slip 壁付き bcond で `convertGmshToForge` する**こと
  (壁が slip だと `wall_dist=0` になり SST の ω が発散する)。

## 計算 run 一覧

すべて SLAU + 陰解法 (block-DPLUR)。中心線 exit Mach は等エントロピー値 2.00 を基準に比較。
詳細・後処理図は各 run の `centerline_compare.png` / `residual_history.png`。

| run | 物理 | 壁 | exit M (中心線) | exit Ps | 状態・備考 |
| --- | --- | --- | --- | --- | --- |
| `run_0001_slau_2d_imp` | 2D 非粘性 Euler | slip | **1.989** | 12.9 kPa | active. 等エントロピー 2.00 と一致 (0.6%) |
| `run_0002_slau_3d_imp` | 3D 非粘性 Euler | slip | **1.990** | 12.9 kPa | active. 2D と完全一致 (検証) |
| `run_0003_slau_2d_visc` | 2D 粘性 層流 | 輪郭no-slip / fb slip | **1.940** | 13.9 kPa | active. BL 変位で M 低下 |
| `run_0005_slau_3d_visc` | 3D 粘性 層流 | 全壁 no-slip | **1.937** | 14.0 kPa | active. 2D 層流とほぼ同等 |
| `run_0004_slau_2d_sst` | 2D 乱流 SST | 輪郭no-slip / fb slip | **1.883** | 15.3 kPa | active. 乱流 BL でさらに M 低下。**乱流熱伝導修正後 T≤Tt** |
| `run_0006_slau_3d_sst` | 3D 乱流 SST | 全壁 no-slip | **1.534** | 25.4 kPa | active. 4壁 no-slip の閉塞で M 大幅低下。T≤Tt |

exit Mach の単調序列 (非粘性 1.99 > 層流 1.94 > 2D乱流 1.88 > 3D乱流 1.53) は
境界層変位による有効面積減少として物理的に整合。中心線静温はいずれも全温 293 K 以下。

## 非平衡凝縮 (H2O) run 一覧 — Wyslouzil Fig.3 検証

Wyslouzil et al. JCP 113, 7317 (2000) **Fig.3 (pv0=1.0 kPa 水)** 条件で、N2 キャリア中の
希薄水蒸気 (Y_H2O=0.01095) の非平衡凝縮を 4 モーメント法で計算し、凝縮あり/なしの静圧比を比較する。
凝縮潜熱で中心線静圧が dry 等エントロピー線より上振れする現象 ([condensation plan](../../plans/active/condensation-nonequilibrium.md))。

- 入口: `inlet_Pressure` Pt=59070 Pa, Tt=286.65 K, Y(N2,H2O)=(0.98905, 0.01095)。
- **carrier=CPG**: 低温 (<200K) で NASA-9 TP が外挿不可・発散するため、N2 キャリアは CPG (γ=1.4) で解く
  (この温度域の N2 cp はほぼ一定で CPG が正確かつ頑健)。H2O は移流種 (vapor budget) として追跡し、
  凝縮潜熱・物性 (Murphy-Koop psat / Hertz-Knudsen 成長) は CPG エネルギー経路に組み込む。
- 検証図: `fig3_compare.png` (非粘性), `fig3_compare_visc.png` (非粘性+粘性SST 重ね合わせ)。

| run | 物理 | 凝縮 | 主要結果 | 状態 |
| --- | --- | --- | --- | --- |
| `run_0007_h2o_cond` | 2D 非粘性 (CPG+species) | off | dry baseline。収束 (rms_ro↓3.5dec) | active |
| `run_0008_h2o_cond_on` | 2D 非粘性 | **on** | 潜熱で T 143→185K, p/p0 +19〜28%。cond/dry 比が exp と ~5% 一致 | active |
| `run_0009_h2o_sst_dry` | 2D 粘性+SST 乱流 | off | dry SST baseline (Wyslouzil 条件)。T≤Tt | active |
| `run_0010_h2o_sst_cond` | 2D 粘性+SST 乱流 | **on** | 凝縮効果 (ユーザ要望)。等温CNT + Hertz-Knudsen (baseline) | active |

**核生成/成長モデル感度 (SST, 凝縮 on, `condKantrowitz`×`condGrowthModel`)** — `compare_models.png`:

| run | Kantrowitz | 成長則 | onset 位置 | 備考 |
| --- | --- | --- | --- | --- |
| `run_0010_h2o_sst_cond` | off | Hertz-Knudsen | x≈1.5cm | baseline。等温CNT は onset 早すぎ・bump 過大 |
| `run_0011_h2o_sst_kw_hk` | **on** | Hertz-Knudsen | x≈2.3cm | Kantrowitz が J 抑制 → onset 遅延 (実験寄り)、bump 控えめ |
| `run_0012_h2o_sst_nokw_gyar` | off | **Gyarmathy** | x≈1.3cm | Gyarmathy は成長急 → bump が最も鋭く過大 |
| `run_0013_h2o_sst_kw_gyar` | **on** | **Gyarmathy** | x≈2.0cm | **onset 位置・peak が実験に最も一致** (Kantrowitz 遅延 + Gyarmathy 急成長) |

知見: **Kantrowitz 非等温補正は核生成 J を抑え onset を下流へ遅らせる**(等温CNT は onset が早すぎ実験の
bump 位置を外す)。**Gyarmathy(熱伝導律速)は Hertz-Knudsen より成長が急**で p/p0 の立ち上がりが鋭い。
両者を併せた `run_0013` (Kantrowitz+Gyarmathy) が x≈2〜3.5cm の onset/peak を実験と最もよく再現
(x=3.1cm: forge 0.353 vs exp 0.354)。下流 x≳4cm は全モデルほぼ収束し、実験より ~0.02 低い (膨張側オフセット)。

**Gyarmathy Knudsen 係数 C 感度 (`condGyarmathyC`, base=Kantrowitz+Gyarmathy)** — `compare_gyarC.png`:

| run | C | onset | x=3.1cm p/p0 |
| --- | --- | --- | --- |
| `run_0014_h2o_sst_gyarC0p0` | 0 (連続体, Kn 補正なし) | 最も早い (x≈1.3) | 0.342 |
| `run_0015_h2o_sst_gyarC1p59` | 1.59 | 早い | 0.348 |
| `run_0013_h2o_sst_kw_gyar` | **3.18 (標準)** | x≈2.0 | **0.353** (exp 0.354) |
| `run_0016_h2o_sst_gyarC6p36` | 6.36 | 遅い | 0.328 |
| `run_0017_h2o_sst_gyarC12p72` | 12.72 | 最も遅い (x≈3.5) | 0.284 |

C は成長率の Knudsen 抑制を制御し、**小さいほど成長が速く onset が早まり、大きいほど遅く弱くなる**
(g プロファイルの立ち上がり位置が C で単調にシフト)。**標準値 3.18 が onset/peak を実験に最も良く再現**し、
ずらすと悪化 (C<3.18 で早すぎ・C>3.18 で遅すぎ)。下流 x≳5cm は最終的に全量凝縮するため C 依存は小。
→ **この系では Gyarmathy 標準係数 3.18 がほぼ最適**であることを確認。

**液滴温度 T_d 考慮 (一温度 vs 二温度, `condTwoTemp`, Hertz-Knudsen)** — `compare_2temp.png`:

| run | モデル | 液滴温度 |
| --- | --- | --- |
| `run_0010_h2o_sst_cond` | isoCNT+HK | 一温度 (T_d=T_g) |
| `run_0018_h2o_sst_nokw_hk_2T` | isoCNT+HK | **二温度 (Hill T_d)** |
| `run_0011_h2o_sst_kw_hk` | Kantrowitz+HK | 一温度 |
| `run_0019_h2o_sst_kw_hk_2T` | Kantrowitz+HK | **二温度** |

知見: **この希薄水/N2 系では液滴温度の影響は小さい**(p/p0 差 ≤0.008、概ね ~0.002 = <1%)。液滴超過温度
T_d−T_g は**成長活発な onset 前線で局所的に最大 ~15K**(過冷却 ~25–30K の一部)に達するが、中央値は ~0.2K
(下流は凝縮完了で蒸気枯渇→成長停止→T_d→T_g)。**N2 キャリアが潜熱を効率よく奪う**ため巨視的影響は
Kantrowitz/Gyarmathy のモデル選択より一桁小さい。方向は予想どおり成長を僅かに遅らせ onset を微小に下流へ
(g プロファイルで二温度=破線が僅かに遅延)。→ **純蒸気凝縮では T_d が必須だが、キャリア希薄凝縮では一温度近似で十分**、という結論。実装詳細は [methods/condensation.md §4](../../methods/condensation.md#理論)。

**複数 H2O 分圧スイープ (モデル汎化, Kantrowitz+Gyarmathy 固定)** — `compare_multicond.py` / `compare_multicond.png`:

入口 H2O 質量分率のみを変えて pv0 を 1.00 / 0.50 / 0.26 kPa とし、Wyslouzil Fig.3 multi-P (`wyslouzil_fig3_multiP.csv`, P0=59.07 kPa) と中心線 p/p0 を比較。**単一条件で較正したモデルを他分圧に汎化できるかの検証**。

| run | pv0 | 入口 Y_H2O | onset (exp) | x=3.1cm p/p0 (forge / exp) |
| --- | --- | --- | --- | --- |
| `run_0013_h2o_sst_kw_gyar` | 1.00 kPa | 0.01095 | x≈2cm | 0.353 / 0.354 |
| `run_0020_h2o_sst_kg_0p50kPa` | 0.50 kPa | 0.005461 | x≈3.5cm | 0.279 / 0.297 |
| `run_0021_h2o_sst_kg_0p26kPa` | 0.26 kPa | 0.002835 | x≈5cm | 0.274 / 0.293 |

知見: **分圧依存 (凝縮 onset の下流シフトと圧力上昇幅の序列) を 3 条件すべてで再現**。pv0 が下がるほど onset
が下流へ移り bump が弱まり dry に漸近する、という実験傾向を捉える。forge は絶対 p/p0 が一様に exp より
~5–8% 低い (既知の dry baseline オフセット = forge dry が非粘性等エントロピー、実験は粘性排除厚で上振れ) が、
**凝縮物理の相対挙動・分圧応答は単一条件較正なしで汎化**する。
(注: `run_0021` は step6000 で外部中断したが、p/p0・g とも step2000→6000 で 4 桁不変＝定常化済みのため発達場として採用。0.26 kPa は g max≈6e-4 と凝縮ごく僅か。)

**非平衡凝縮の検証結果**:

- **非粘性 (run_0008)**: 凝縮ありの中心線 p/p0 は dry 比で 1.19〜1.28 倍 (スロート下流 2.5〜8 cm)、
  実験の cond/isentrope 比 1.18〜1.23 と ~5% 以内で一致。**凝縮効果 (比) は定量一致**するが、絶対 p/p0 は
  非粘性 dry が実験等エントロピーよりやや低い (粘性排除厚さ未考慮)。
- **粘性+SST 乱流 (run_0010 vs run_0009, ユーザ要望)**: 粘性排除厚さで dry baseline が押し上がり、
  **絶対 p/p0 も実験と一致**する (下表)。凝縮 onset (x≈1.7cm の潜熱スパイク) も再現。x≳6cm の下振れは
  壁解像不足ではなく (y+ は良好、下記)、2D 近似/ノズル面積比など膨張側のオフセット (dry/cond 共通)。

  | x [cm] | SST dry | SST cond | exp isentrope | exp cond |
  | --- | --- | --- | --- | --- |
  | 2.5 | 0.300 | 0.381 | 0.305 | 0.375 |
  | 3.0 | 0.279 | 0.349 | 0.290 | 0.355 |
  | 5.0 | 0.223 | 0.267 | 0.240 | 0.290 |

  → **「粘性+乱流で凝縮あり/なしを計算すると Fig.3 の差が出る」というユーザ仮説を定量的に確認**。
  比較図は `fig3_compare_visc.png` (非粘性破線 + 粘性SST実線 + 実験点)。

## SU2 クロスチェック (dry, 3 物理モデル) — `compare_su2_forge.py` / `compare_su2_forge.png`

同一ノズルの中心線 p/p0 を forge (SLAU, セル中心) と SU2 v8.5 (ROE, 節点) で比較。

| 物理 | forge run | SU2 run | forge↔SU2 一致 |
| --- | --- | --- | --- |
| Euler (非粘性) | `run_0001_slau_2d_imp` | `run_0020_su2_euler` | 全域 ≤0.6% |
| Laminar (粘性) | `run_0003_slau_2d_visc` | `run_0023_su2_lam_conv` | 全域 ≤0.6% |
| Turbulent (SST) | `run_0004_slau_2d_sst` | `run_0024_su2_sst_conv` | 全域 ≤0.8% |

→ **forge は 3 モデルとも SU2 と全域 <1% 一致**。粘性/SST は下流で等エントロピー線の上 (境界層排除厚) に乗り、
forge が SU2 と同等の BL 変位を出していることを確認。

⚠️ **教訓**: SU2 の初回 run (`run_0021_su2_lam` / `run_0022_su2_sst`) は前セッション中断で 975/752 iter で
SIGTERM 終了 (`su2.log` に `Exit Success` は出るが `rms[RhoE]≈-0.2` の未収束) しており、その未収束場と比べると
下流で偽の 20-25% 差が出ていた。途中解から継続収束 (`run_0023`/`run_0024`, `rms[RhoE]≈-1.4`, 出口積分量
ドリフト <0.2%) させると上記のとおり一致。**SU2 は `Exit Success` でなく iter 数・`rms[RhoE]` で収束判定すること**
([procedures/su2-cross-check.md](../../procedures/su2-cross-check.md) の収束確認節)。

## TP gas (thermalMethod 2) + 凝縮の検証 — 真因 = sub-200K NASA-9 外挿

ユーザ要望で TP 気相 + 凝縮を wys で計算しようとしたが dry TP すら発散。系統的切り分け
(`run_0025`〜`run_0042`、大半は破棄可の診断用) で**入口/壁/出口 BC・凝縮・種数・粘性・乱流・
時間積分・convMethod・float 精度・初期値をすべて棄却**。**真因は wys が出口 T≈159K (元の凝縮ゾーン ~27K) で
NASA-9 有効下限 200K を割り、低温外挿が TP を不安定化**すること。全温 Tt スイープで確定:

| run | 全温 Tt | 出口 T | 結果 |
| --- | --- | --- | --- |
| `run_0041_tp_n2_hot500` | 500 K | 278 K (>200K) | **400 step 完走・残差安定** |
| `run_0042_tp_n2_Tt360` | 360 K | ~194 K | step24 で発散 |
| (Tt=286.65, 通常) | 286.65 K | 159 K | step12 で発散 |

→ TP 実装は壊れていない (T>200K では wys でも動く)。**wys の極低温凝縮で CPG を強制してきた方針は正しい**。
詳細・全棄却リストは [`.github/plans/condensation-nonequilibrium.md`](../../plans/active/condensation-nonequilibrium.md) の 2026-06-16 ログ。

## 壁面 y+ (SST メッシュ `nozzle_2d_sst.h5`)

`wall_yplus.py` で run_0009/run_0010 の壁第一層 y+ を収束場から実測 (y+=√(u_t·y1/ν), ν=分子粘性):

- 第一セル距離 y1 = **1.6〜9.0 µm** (スロート半高 2.25mm の ~0.1%)。
- **y+ ≤ 1 が壁面の 98.2%**、mean 0.61、max **1.88**、全面 **y+ ≤ 2**。
- 1 を超えるのは入口リーディングエッジ (x≈−6.3cm) のごく一部のみ。スロート下流 (超音速・凝縮領域) は
  y+ ≈ 0.5 で完全に解像 → **低Re SST の壁まで積分が成立**。`wall_yplus.png` に分布。
- **乱流粘性も正しく作動**: μ_turb/μ_lam は壁第一セルで ≈0 (粘性サブレイヤ、k→0/ω→∞ ゆえ)、対数層で
  median 4.5・max 29、コアで ~3。場全体で 79% のセルが μ_turb>μ_lam。第一セルが粘性サブレイヤ内ゆえ
  τ_w を有効粘性 (μ_lam+μ_turb) で評価しても y+ は不変 (= 壁解像が本物である証拠)。

(旧記述「近壁メッシュは Euler デモ用で y+~1 ではない」は、非粘性デモ用 `nozzle_2d.h5` に関するもので、
SST 専用メッシュ `nozzle_2d_sst.h5` には当てはまらない。上記実測のとおり壁解像済み。)

## 既知の課題

- 本ケース検証中に **RANS エネルギー方程式の乱流熱伝導欠落バグ**を発見・修正
  (`.github/plans/diffusion-turbulent-thermal-conductivity.md`)。修正前は近壁静温が
  449 K (全温 293 K 超過)、修正後 293 K に収束。

## 後処理スクリプト

- `postproc_centerline.py <run> <label>` — 中心線 Mach / 静圧を面積比等エントロピーと比較。
- `plot_residuals.py <run> <label>` — 全 rms 残差の片対数プロット。
- `wall_yplus.py <run> [--wall-physid 3]` — SST 壁第一層 y+ を収束場から算出・可視化。
- `compare_fig3.py` / `compare_fig3_visc.py` — Wyslouzil Fig.3 (凝縮あり/なし) と中心線 p/p0 を比較。
- `build_restart.py <src_res> <src_mesh> <dst_mesh> [--rok K --roomega W]` — 引き継ぎ初期場生成。

## 修正形状 (Fig.3) run 一覧 — run_0047 以降

正しい Fig.3 ノズル (`mesh/make_nozzle_fig3.py` → `nozzle_fig3_2d.h5`, 2D pseudo, front/back=slip) で再構築。
SLAU 陰解法 (timeIntegration:11, blockDPLUR, **nStepInner:5, cfl_pseudo:4**)。dry は実験 isentrope、
凝縮は実験 1kPa と比較。比較図: `dry_vs_exp.png` (laminar/SST vs isentrope)、`cond_models_compare.png` (4モデル vs 1kPa)。

| run | 物理 | 凝縮モデル | 主要結果 | 状態 |
| --- | --- | --- | --- | --- |
| `run_0047_fig3_2d_lam_dry` | 2D laminar | off | dry。中心線 p/p0 が実験 isentrope と −2〜5% | active |
| `run_0048_fig3_2d_sst_dry` | 2D SST | off | dry。**実験 isentrope と ±1.5% (laminar より良)** → 凝縮は SST 上で実施 | active |
| `run_0049_fig3_2d_sst_hk` | 2D SST | Hertz–Knudsen | 全凝縮 (g~90%)。onset 過早・overshoot (0.44 vs exp 0.36 @2cm) | active |
| `run_0050_fig3_2d_sst_kwhk` | 2D SST | **Kantrowitz+HK** | **onset 遅延で実験最良** (Fluent UDF 知見と一致) | active |
| `run_0051_fig3_2d_sst_gyar` | 2D SST | Gyarmathy | 全凝縮。onset 過早・overshoot | active |
| `run_0052_fig3_2d_sst_kwgyar` | 2D SST | Kantrowitz+Gyarmathy | onset 遅延で実験良好 | active |

知見: **凝縮成長停止バグ修正後** (核生成を $r_{\rm nuc}=1.01r_*$ で生成し液滴を $r_*$ から離脱させる;
[methods/condensation.md](../../methods/condensation.md#実装) §検証 case/16) に全モデルが全凝縮を再現。
下流 (x≳4cm) は全モデル実験 ±数%。**onset 域は Kantrowitz 有り (Kw+HK/Kw+Gyar) が overshoot を抑え実験に最良**。
残課題: x≈3cm のピークが実験よりやや高い (onset レート微調整)、3D 壁解像 (`make_nozzle_3d_wallres.py`) は未実施。

## TP split (MIXDRY 擬似種 + H₂O 独立種) と node 検証 — run_0170 以降 (2026-08-17)

計画 [tooling-nozzle-tp-split-h2o-condensation.md](../../plans/accepted/tooling-nozzle-tp-split-h2o-condensation.md)。
イソブタン燃焼ガス (case/42) の H₂O 凝縮に向け、「H₂O 以外を NASA-9 擬似種 `MIXDRY` に畳み、H₂O だけ独立種」
とする 2 種 TP を、既存 Fig.3 (2D SST, Kw+HK, `run_0050`) で検証した。組成は既存どおり **H₂O 1.095 %**。
投入: `run_tp_split_wys.py` (cell, 同メッシュ) / `run_tp_split_wys_node.py` (node, 平面メッシュ)、比較: `compare_tpsplit_wys.py` → `compare_tpsplit_wys.png`。

| run | 離散化 / メッシュ | 熱力学 | 物理 | 主要結果 | 状態 |
| --- | --- | --- | --- | --- | --- |
| `run_0170_fig3_2d_sst_kwhk_tpsplit` | cell / 既存 `nozzle_fig3_2d` (押し出し擬似2D) | **TP `[MIXDRY(=N2),H2O]`** thermoHrefTemp 298.15 (IC は run_0048 dry SST 場を `cpg_field_to_tp` で TP 化) | SST 粘性, Kw+HK | **run_0050 (CPG) を再現**: g_max 0.0106 vs 0.0107, onset 2.45 vs 2.35 cm, 中心線 p/p0 差 ≤3.1 % (TP の H₂O cp 差)。NaN 0、rms_ro 1.2e-9 (2.2 桁, still falling) | active (**TP split 検証の正本**) |
| `run_0184_wysinv_cell_cpg_cond` / `run_0185_wysinv_node_cpg_cond` / `run_0186_wysinv_node_tpsplit_cond` | cell / **node** / node、平面 2D 一様メッシュ `mesh/nozzle_fig3_2d_planar_inv` (301×120, 壁クラスタなし) | CPG / CPG / TP split | **非粘性 slip**, Kw+HK, 2 次 (段階起動 soft→mid→本段 12000) | **node = cell**: g 0.0109 = 0.0109, onset 1.75 cm 同一, p/p0 差 2 %; **node TP split** g 0.0108, onset 1.85 cm, p/p0 差 1.7 %。node rms_ro 3.4e-9 (3.3 桁)。非粘性なので onset は SST 基準 (2.35 cm) より上流 | active (**node の検証**) |
| `run_0187_wysinv_node_tpsplit_cond_regress_noevap` / `run_0188_wysinv_node_tpsplit_cond_evap` | run_0186 と同一入力 (node, 非粘性, TP split, Kw+HK) | TP split | 蒸発実装後バイナリの回帰 (0187 = `condEvaporation: 0` [当時の既定]; 現在の既定は 1) と蒸発 ON (0188 = 1) | **回帰同一**: 0187 vs 0186 max|ΔT| 3.8e-3 K・|ΔP/P| 2.4e-5・|Δg| 1.2e-7 (run 間ノイズ)。**蒸発 ON も同一** (0188 vs 0187 |ΔT| 7.8e-3 K): Wyslouzil は液相域が全域 S≥2.03 で蒸発分岐が発火しない。残差 3.1 桁 (0186 と同じプラトー) | active ([蒸発 plan](../../plans/accepted/condensation-evaporation.md) 回帰) |
| `run_0172_fig3_2d_sst_kwhk_tpsplit_restart_ctrl` / `run_0173_fig3_2d_sst_kwhk_tpsplit_evap` | run_0170 (cell, SST no-slip 断熱壁, TP split, Kw+HK) の res_15000 から restart 9000 step | TP split | 蒸発 OFF (対照) / ON | 壁近傍 (wall_dist<80 µm) に S 0.4–0.99 の液相 (g 中央値 1.7e-5, max 1.6e-3) が両者に残る = HK 有限速度 (dr/dt −1e-4 m/s) と上流からの移流供給の釣り合い (Tt=287 K の壁は熱くなく蒸発が遅い、λ 律速は非発火)。evap vs ctrl: max|ΔT| 0.22 K, |Δg| 2.6e-4。**restart 後に両 run とも壁 T が 290 K>Tt へドリフト** ([[forge-sst-restart-nonfidelity]] 既知)。高温壁デモは case/42 NS で実施 | 削除済み(2026-08-31) (蒸発 plan §6-3 の低温壁 対照) |
| `run_0171_fig3_2d_sst_kwhk_tpsplit_node` / `run_0179`, `run_0181`〜`0183_wysnode_*` | node / 平面 `nozzle_fig3_2d_planar` (Bump 0.004, y₁≈0.5–5 µm, AR 728) | CPG or TP | SST or laminar no-slip (0182 のみ slip) | **node × 平面壁クラスタ no-slip は不成立**: 壁ノード T が Tt を超え (出口コーナーから 326→400+ K)、衝撃列が遡上して 6000–9000 step で unstart (2 次・1 次・SST・laminar・nodeWallDirichlet 0/1・出口 Pt=Ps いずれも)。slip (0182) だけ健全 (T ≤ Tt)。cross-mesh IC の壁ジグザグは 1D 等エントロピー IC で除いたが本質でない → **node 粘性壁 (平面, 極薄壁セル) の申し送り** | 削除済み(2026-08-31) (**申し送りの根拠**) |

診断 run (0172–0178, 0180, 0187) は削除済 (上記と同じ結論の切り分け)。

## ユーザ指定形状 (2026-08-19) — 形状の確認と Euler/NS × 凝縮 off/on (TP split) — run_0189 以降

ユーザ提示の上壁式 (半高 y [mm], throat x=0; `mesh/nozzle_user_profile.py`):
12.7 (−60≤x≤−38) / 0.16−0.33x (−38<x<−21.27) / 2.5+0.015512821x²+0.000243078x³ (−21.27<x<0) /
2.5+0.00194105x²−7.99459e−5x³ (0<x<8.093) / 2.457620744+0.015709255x (8.093<x<95)。
**現行 `make_nozzle_fig3.py` 形状 (smoothstep 収縮 + 直線膨張) とは一致していなかった** (`check_geometry_user.py` → `geometry_user_vs_current.png`):

| 量 | ユーザ式 | 現行 mesh (`nozzle_fig3_2d*`) |
| --- | --- | --- |
| 入口直管 (x −60〜−38 mm, 半高 12.7) | あり | **無し** (x=−38 から開始) |
| 収縮部 (−38〜0) | 直線 (勾配 −0.33) + 3 次 (曲率半径 32 mm@throat) | smoothstep (曲率半径 24 mm@throat)。**最大 |Δy| 1.56 mm @x=−27** |
| 膨張部 (0〜95) | 3 次ブレンド (0〜8.09 mm, R=258 mm) → 直線 (切片 2.4576, 勾配 0.015709=0.9°) | throat から直線 (切片 2.5, 同勾配)。**|Δy| ≤0.042 mm** (出口で最大) |
| 出口半高 / A/A\* | **3.950 mm / 1.580** | 3.992 mm / 1.597 |
| 1D 等エントロピー p/p0 差 (γ=1.4) | — | 現行が x=5 mm で −8 %、x=20 mm −4 %、出口 −2 % (ユーザ形状の方が膨張が緩い) |

→ 収縮部と throat 直後のブレンドが違い、膨張側の p/p0 に 2〜8 % 効く。以下はユーザ式でメッシュを作り直して再計算した。

**メッシュ** `mesh/make_nozzle_user.py` (5 区間を別カーブ + Transfinite 複合辺、単位 mm): `nozzle_user_2d.geo` (cell 用 1 層押し出し, front/back slip)、
`nozzle_user_2d_planar.geo` (node 粘性用, 壁クラスタ)、`nozzle_user_2d_planar_inv.geo` (node 非粘性用, 一様)。
streamwise 347 節点 (Δx≈0.42 mm; 直管 0.63 mm) × 壁法線 120 (Bump 0.004, y₁≈0.5–5 µm)、41174 cells。品質 **PASS** (AR max 724, skew 0.47)。
**投入** `run_user_profile.py RUN --disc cell|node --phys euler|sst [--cond]`: TP split `[MIXDRY(=N2), H2O]` (Y_H2O 0.01095, `thermoHrefTemp` 298.15)、
Pt 59070 Pa / Tt 286.65 K、出口 Ps=Pt=2000 Pa、Euler = slip + `visc 0` `thermCond 0`、NS = SST 低 Re (`wallTreatmentSST 0`) no-slip 断熱壁 + Sutherland + const-Pr 0.72 (Prt 0.9)、
凝縮 = Kw+HK (`condModel 1, condKantrowitz 1, condGrowthModel 0`)、IC = 1D 等エントロピー (TP semi-perfect)、
起動 soft (1 次 cfl0.5, 3000) → mid (1 次 cfl1, 3000) → 本段 **全域 2 次 (convMethod 1/limiter 2) cfl_pseudo 2**, 12000 step (NS は `continue_run.py` で +24000 継続)。
比較 `compare_user_profile.py` → `compare_user_profile.png` / `compare_user_profile_centerline.csv` (中心線 x, M, P, T, g を 0.5 mm 刻みで保存)。

| run | 離散化 | 物理 | 凝縮 | 主要結果 (中心線, x≈94 mm 出口) | 収束/定常 | 状態 |
| --- | --- | --- | --- | --- | --- | --- |
| `run_0189_user_cell_euler_dry` | cell | Euler | off | M 1.912, p/p0 0.1465, T 165.8 K (1D 等エントロピー M 1.918 / p/p0 0.1448) | rms_ro 5.4e-10 (1.6 桁, still falling), `check_quasisteady` machmax/pmax **STEADY** | active |
| `run_0190_user_cell_euler_cond` | cell | Euler | **on** | onset (g>1e-3) **x=20.2 mm**, g_exit 0.0109 (全量凝縮), M 1.704, p/p0 0.1771, T 199.3 K | rms_ro 1.6e-9 plateau, 全 snapshot で onset/p/p0 4 桁不変, STEADY | active |
| `run_0191_user_cell_sst_dry` → **`run_0196_user_cell_sst_dry_cont`** (+24000) | cell | NS (SST) | off | M 1.792, p/p0 0.1762, T 174.8 K。**exp isentrope と x≥11 mm で +0.9〜+1.9 %** (throat 直後 x=0.9 mm は +4.6 %) | rms_ro 4e-8 plateau (2 次リミッタ), 継続 24000 step で p/p0 変化 ≤1e-4, **STEADY** | active (0191 は中継) |
| `run_0192_user_cell_sst_cond` → **`run_0197_user_cell_sst_cond_cont`** (+24000) | cell | NS (SST) | **on** | onset **x=23.2 mm**, g_exit 0.0108, M 1.600, p/p0 0.2090, T 208.0 K。exp cond 1 kPa と −5〜−3 % (x 21–32 mm, onset 帯) / +4〜+5 % (x 42–72 mm) | rms_ro 4e-8 plateau, 継続で不変, **STEADY** | active (0192 は中継) |
| `run_0193_user_node_euler_dry` / `run_0194_user_node_euler_cond` | **node** (平面一様) | Euler | off / on | **cell と一致**: p/p0 0.1460/0.1767 (cell 0.1465/0.1771), onset 20.2 mm 同一, g 0.0109 同一 | rms_ro 4.3e-9 (2.7 桁), STEADY | active (node 対照) |
| `run_0195_user_node_sst_dry` | node (平面壁クラスタ) | NS (SST) | off | **不成立 (2026-09-07 に解決: run_0212/0213 参照)** (既知の申し送りと同じ指紋: 壁ノード T 504 K > Tt, P 126 kPa > Pt が x=−40 mm から、本段 6000 step で出口側 unstart)。途中で停止 | — | 削除済み(2026-08-31) |

| `run_0198_user_node_sst_dry` / `run_0199_..._outfix` | node (平面壁クラスタ) | NS (SST) | off | **unstart** (出口列 P 45–57 kPa > 指定 2 kPa、全域 M≈0.3)。`nodeInletCornerWall: 1` で run_0195 の入口角 P>Pt は消えたが出口側が残る。0199 は TP 出口 γ 修正の効果確認 (bit 同一 = 無関係) | — | 破棄予定 (診断) |
| `run_0200`–`run_0211` (診断) | node | NS | off | 切り分け: TP/CPG・dilatation・壁関数・nodeWallDirichlet・出口 k/ω・壁第一層 0.6/2.4/8 µm は全て同じ失敗、層流 (0203) は健全。100 step 再現 (0201) で出口列の壁向き Uy・+5 % P バンプ→SST の k 膨張→剥離前線が上流へ、と特定。**`outflow` (全量外挿) で解消** (0210 CPG / 0211 TP) | — | 破棄予定 (診断) |
| `run_0214_ab_ps_matched` | node | NS (SST, CPG) | off | `outlet_statPress` のまま Ps を実出口圧 10.4 kPa に合わせても解消 (case/45 と同じ構成)。真因 = **Ps 指定 (2 kPa) ≪ 実出口圧 (10.4 kPa) × node の亜音速壁列** | mid 段 rms_ro 3.9e-7 (outflow と同一) | 破棄予定 (診断) |
| `run_0212_user_node_sst_dry_outflow` → **`run_0213_user_node_sst_dry_outflow_cont`** (+24000) | **node** (平面壁クラスタ, `nodeInletCornerWall: 1`, 出口 `outflow`, katoLaunder) | NS (SST) | off | **cell run_0196 と中心線 0.1 % 以内で一致** (出口 M 1.793 / p/p0 0.1759 / T 174.7 K)。`compare_node_vs_cell_sst.png` | `check_convergence`: rms_ro 4.3e-9 プラトー (0212 から計 2.7 桁↓), roe/roK/roOmega still falling; `check_quasisteady` machmax/pmax **STEADY** | active (node NS 生産) |

| `run_0219_user_node3d_sst_dry_half` (AWS g5) | **node 3D** 半幅+対称面 (`mesh/nozzle_user_3d.msh`: 12.7 mm 幅の半分 6.35 mm, z=0 側壁 no-slip / z=6.35 対称面 slip, 半幅 50 節点 z1 2 µm, **2.08 M 節点**, 品質 SOFT-PASS AR 1196@0.02 %) | NS (SST) | off | soft 段は収束したが mid 段 step 2 で k/ω 爆発 → 段間引き継ぎの `interp_field` が (x,y) 最近傍で z 列を混同していた (index コピーに修正) | — | 破棄予定 (h5 は 0220 が流用) |
| **`run_0220_user_node3d_sst_dry_half`** (AWS g5, res_*.h5 は AWS 側) | node 3D 半幅+対称面 (同上), 出口 `outflow`, IC = run_0213 の 2D 場を (x,y) 移植, soft/mid 2000+2000 → 本段 **cfl 6 + relax 0.7** 12000 | NS (SST) | off | **3D node SST 完走** (12000 step 14.6 min @A10G = 73 ms/step, GPU 3.2 GB)。対称面中心線: throat M 0.898, 出口 M **1.604** / p/p0 0.235 (2D は 0.951 / 1.793 / 0.176 → 側壁 BL の閉塞。旧形状 4 壁 3D run_0006 の 1.53 と同傾向)。`compare_2d_vs_3d_node_sst.png` / `centerline_sym.csv`。**注意**: 側壁∩輪郭壁の角線近傍 (x>15 mm, wd≤0.5 mm) で T が Tt+20 K まで上がる (壁温は node の既知アーティファクト [node-wall-entropy-checkerboard] の 3D 版、コアは無影響) | `check_convergence`: rms_ro 5.7e-11 (falling), roUx 7e-9, roK 4e-8, roOmega 0.028 plateau; `check_quasisteady` machmax/pmax **ALL STEADY**; NaN 0 | active (3D node SST 初完走。ただし run_0222 で出口角 unstart と判明 → 定常解ではない) |

| `run_0222_user_node3d_sst_dry_half_cont` (AWS) | run_0220 の継続 +24000 (index コピー restart) | NS (SST) | off | **出口角線から unstart が進行し step 25941 (通算) で NaN**: 壁 p/p0 @x=95 が 0.20 (2k) → 0.25 (12k) → 0.30 (18k) → 0.44 (24k)、角線 T−Tt +2 → +19 → +41 K。run_0220 の「+30 %」は過渡値で定常解ではない。2D では `outflow` で解消した出口列の問題が、3D では側壁∩輪郭壁∩出口の角線で残る (SST 固有、層流は無事) | `wall_series.py` (時系列表) | 破棄予定 (診断記録) |
| **`run_0221_user_node3d_lam_dry_half`** (AWS) | 3D 半幅 **層流** (turbulence none, 同メッシュ h5 流用, IC = 2D 層流 run_0203, soft/mid 2000+2000 → 本段 cfl 5 + relax 0.7 12000) | NS (層流) | off | **収束・健全**: T ≤ Tt (超過ノード 0)、角線 T−Tt −7 K。**壁 p/p0 は 2D SST / 実験 dry と ±1〜2 %** (x=16.4: 0.367 vs exp 0.364, x=45.6: 0.255 vs 0.254)。側壁 BL の閉塞は実験と整合する大きさ | `check_convergence`: 全列 2.0〜2.4 桁↓ still converging (rms_ro 4e-12); `check_quasisteady` pmax STEADY; 壁 p/p0 は 6000 step 以降 4 桁不変 | active (**3D node NS 実証・実験一致**) |

| **`run_0224_user_node3d_lam_cond_half`** (AWS) | 3D 半幅 層流 **凝縮 ON** (Kw+HK, 同メッシュ h5 流用, IC = run_0221 の index コピー, 段階起動なし, cfl 2 + relax 0.7, 12000) | NS (層流) | **on** | **完走・健全** (T ≤ Tt, NaN 0)。onset (中心 g>1e-3) **x=23.4 mm** (2D NS SST cond run_0197 は 23.2)、g_exit 0.0109 (全量凝縮)。壁 p/p0 の cond/dry 比 1.16/1.23/1.22/1.21/1.20 (x=32/42/52/62/72) が **実験の cond/isentrope 比 1.20/1.20/1.19/1.19/1.18 と一致**。絶対値は実験 cond と onset 帯 −2〜−5 %・下流 +1〜+3 %。`wall_pp0_3d_lam_cond_vs_exp.png` | `check_convergence`: 全列 1.6〜3.7 桁↓ still converging (rms_ro 1.5e-11); pmax STEADY; 壁 p/p0 は 6000→12000 で +1.4 % (x=42) とわずかに動く | active (**3D 凝縮の初結果**) |
| `run_0223_user_node3d_sst_dry_half_ext` (AWS) | 3D SST + **出口バッファ 20 mm** (`mesh/nozzle_user_3d_ext.msh` 2.37 M 節点, x>95 の輪郭壁・側壁 slip physID 6) | NS (SST) | off | soft 段 step 826 で roOmega NaN: 接合 x=95 で wall_dist が 0→数 mm に跳び、BL を運ぶ slip 壁近傍 (x=100, y=±3.93, z 0.1) で ω 1e21。→ `mesh.wallDistExtraPhysIDs: [6]` (slip 壁も壁距離に含める) を実装して再試行 (run_0225) | — | 破棄予定 (診断) |
| `run_0300_corner_node3d_sst_coarse` (ローカル RTX 3060) | **角線ノード仮説の判別**: 同じ生成器で x 半分・NY 60・z 25 (`mesh/nozzle_user_3d_coarse.msh` 257k 節点, 側壁 z₁ 2 µm 同一・輪郭 y₁ 1.8 µm, 角線構造は 2.08M と同一, QC SOFT-PASS AR 1214@0.04 %), node SST, run_0220 と同手順 (IC=run_0213 interp, soft/mid 2000+2000 → cfl 6 + relax 0.7 8000, outflow) | NS (SST) | off | **T>Tt のノード 0** (全 res で), T0 最大 +2.3 K (2D と同等), 出口角壁 T−Tt −5.7 K・P 11.68 kPa で頭打ち。角線ノード自体は 2.08M でも隣接と連続 (受動) → 角線の双対構造は加熱の原因でない (plan turbulence-sst-node-corner-heating §2.1) | `check_convergence`: rms_ro 3.0 桁↓ falling, roe 2.3 桁↓, roOmega 2.8 桁↓ (roUz 1.2e-9 plateau); `check_quasisteady` pmax/machmax STEADY | active (判別実験) |
| `run_0301_corner_cell3d_sst_coarse` (ローカル) | 同上メッシュの **cell** 離散化 (IC=run_0196 interp) | NS (SST) | off | **T>Tt のセル 0**, T0 最大 +2.4 K, node run_0300 と角部分布一致 | `check_convergence`: 全列 plateau (cell 既知の床, rms_ro 7.6e-9); pmax/machmax STEADY | active (判別実験) |
| `run_0302_corner_node3d_sst_finexy_z25` (ローカル) | (x,y) は 2.08M と同一 (y₁ 0.6 µm)・z 25 節点 (`mesh/nozzle_user_3d_finexy_z25.msh` 1.04M, **QC FAIL: AR>1000 が 0.66 %, max 2601** の診断 run), node SST, 同手順 | NS (SST) | off | **角部加熱が再現**: T>Tt 151k ノード (14.5 %), 出口角壁 T−Tt +2.1→+5.2→+6.7→+7.1 K (2000→8000 step, 上昇中), 角壁 P 12.2→13.5 kPa 上昇 (unstart ドリフト), x=60 角壁 +3.3 K。z を粗くしても再現するので引き金は輪郭壁 y₁ 0.6 µm × 側壁 2 µm の異方角セル (AR 700) 側 (角線構造・z 解像度は無関係) | `check_convergence`: rms_ro 1.3 桁↓, roe 0.6 桁↓, roOmega plateau 1.3e-1; `check_quasisteady` pmax/machmax STEADY だが角壁 T/P は単調上昇 | active (判別実験, H1 検証用の小型メッシュ) |
| `run_0303_sstdef_ref` / `run_0304_sstdef_keys0` / `run_0305_sstdef_new` (ローカル) | **SST 整合オプション既定 ON の回帰 (node 2D)**: run_0213/res_24000 から +6000 (index コピー)。ref=旧バイナリ (b12d74a4)・旧既定 / keys0=新バイナリで `sstOmegaProdFromPk:0, sstSigmaBlend:0` 明記 (F1 ラグ除去・等方項リミッタ前の効果のみ) / new=新既定 | NS (SST) | off | keys0−ref ≤4e-6 (ビット同等)。new−ref: 壁 p/p0 @16.4/45.6/85 = 0.3657/0.2574/0.1874 → 0.3658/0.2575/0.1876 (+0.03〜+0.1 %), k 最大 2.8 %・μt 6 % (出口側)。plan [turbulence-sst-consistency-options](../../plans/active/turbulence-sst-consistency-options.md) §3.6 | `check_convergence`: 参照 plateau 継続 (rms_ro 4.3e-9), NaN 0 | active (A/B 記録) |
| `run_0306_sstdef_cell_ref` / `run_0307_sstdef_cell_keys0` / `run_0308_sstdef_cell_new` (ローカル) | 同上の **cell** 版 (run_0196/res_24000 から +6000) | NS (SST) | off | keys0−ref ≤8e-5。new−ref: 壁 p/p0 0.3643/0.2578/0.1870 → 0.3644/0.2579/0.1872 (+0.03〜+0.1 %), k 3.5 %・μt 6 % | plateau 継続 (rms_ro 4.2e-8), NaN 0 | active (A/B 記録) |
| `run_0309_ek_node_off` / `run_0310_ek_node_on` / `run_0311_ek_cell_on` (ローカル) | **全エネルギーに ρk を含める (`sstEnergyIncludesK`, plan [turbulence-sst-energy-includes-k](../../plans/active/turbulence-sst-energy-includes-k.md))**: run_0213 (node) / run_0196 (cell) の res_24000 から +6000。off=キー 0 (前バイナリ run_0305 とのビット同等確認) / on=キー 1 | NS (SST) | off | off−run_0305 ≤7e-6 (同等)。**on−off: 壁 p/p0 @16.4/45.6/85 = node 0.3658/0.2575/0.1876 → 0.3659/0.2577/0.1877, cell 0.3644/0.2579/0.1872 → 0.3645/0.2581/0.1874 (+0.03〜+0.08 %)**, max(T0−Tt) 1.96→1.56 K (node) / 2.00→1.58 K (cell), k/μt ≤ 数 %。全温は T0 + k/c_p で比較 | `check_convergence`: node rms_ro 3.0 桁↓, roK 3.6 桁↓, roe 1.9 桁↓ still converging; cell plateau (既知)。NaN 0。切替直後 rms_roK 7e-3 のスパイク→6000 step で床 | active (A/B 記録) |

| `run_0225_user_node3d_sst_dry_half_ext_ps` (AWS) | 3D SST + 出口バッファ + `wallDistExtraPhysIDs: [6]` + **出口静圧固定 Ps 10.3 kPa** | NS (SST) | off | soft 段 step 11 で NaN: 接合 x=95 の壁ノードで `wall_y_eff` が wall_dist=0 の slip 隣接を拾い ω ピン 1e23 → `ransBoundary_d.cu` で wall_dist>0 の非壁隣接だけを使うよう修正 | — | 破棄予定 (診断) |
| **`run_0226_user_node3d_sst_dry_half_ext_ps`** → **`run_0227_..._cont`** (+24000, AWS) | 3D SST + バッファ 20 mm (slip, 壁距離に算入) + 静圧固定 10.3 kPa + wall_y_eff 修正, cfl 6 + relax 0.7 | NS (SST) | off | **完走・NaN 0。角線 T−Tt が通算 12000→36000 で +9.4→+0.7 (x=46) / +10.9→+1.8 (x=85) と単調減少、x=16 は −1.1 (Tt 未満)。壁 p/p0 は 0.3199→0.3282 (x=46) で頭打ち**。上流 (x≤85) の壁圧は run_0220 (バッファ無し) と 4 桁同一 = 出口処理は上流に無影響、出口角の暴走だけを止めた。実験 dry より +29 % (角線の SST k 未減衰、plan turbulence-sst-node-corner-heating) | `wall_series.txt`; rms_ro 5e-11 台 | active (0228 で更新) |

| **`run_0228_user_node3d_sst_dry_half_ext_ps_relguard`** (AWS) | run_0227 の最終場から index コピー、**k/ω 拡散の相対ガード** (`scalarTransport_d.cu`, 旧 1e-12 [m³] 絶対値 → 1e-6·|d||S|)、他は run_0226 と同一 (バッファ + 静圧固定 10.3 kPa + cfl 6) | NS (SST) | off | **角部加熱が消滅**: T > Tt+1 のノード 0 個 (旧 42896)、角線 T−Tt −9〜−17 K。角対角線の ω が壁漸近解 6ν/(β₁d²) に一致 (0.6 µm: 1.8e10 vs 1.7e10; 2.5 µm: 1.4e9 vs 1.0e9)、k は粘性底層で 1e-5 (旧 1.2e3)。壁 p/p0 @x=46 0.328 → **0.273** (2D SST 0.257, 層流 0.255, exp 0.254)、@x=85 0.251 → 0.206 (2D 0.187)。`wall_pp0_guard_fix.png` | `check_convergence`: roUx/roUy 3 桁↓, roOmega 5 桁↓, rms_ro 7e-11 falling; pmax STEADY; 壁 p/p0 は末尾 2000 step で −0.0006 (まだ僅かに下降) | active (**3D node SST の正本**) |

| **`run_0229_user_node3d_sst_cond_half_ext_ps`** (AWS) | 3D SST **凝縮 ON** (run_0228 の場から index コピー, バッファ + 静圧固定, cfl 2, 12000) | NS (SST) | **on** | 完走・NaN 0。onset **24.7 mm** (2D SST cond 23.0 / 3D 層流 cond 23.4)、g_exit 0.0109。壁 p/p0 は実験 cond より +9〜+11 % (dry の +7 % 上乗せ)、cond/dry 比 1.13〜1.24 は実験 1.18〜1.20 と整合。T max 297 K (潜熱, 2D と同じく物理)。`wall_pp0_cond_all.png` | rms_ro 2.6e-11 falling, roe/roK 2.7 桁↓; pmax STEADY | active (3D SST 凝縮) |
| `run_0230_user_node_sst_cond_outflow` (ローカル) | **2D node SST 凝縮** (run_0213 の場から index コピー, outflow, cfl 2, 12000) | NS (SST) | **on** | onset 23.0 mm、壁 p/p0 は実験 cond と onset 帯 −3〜−5 % / 下流 +4〜+5 % (cell run_0197 と同傾向) | — | active (2D node 凝縮対照) |
| `run_0231_ab2d_{base_kpin0,base,omegaPk,sigmaBlend,isoStress,energyK,all}` | SST 整合オプション A/B (run_0213 の場から 3000 step, plan turbulence-sst-consistency-options §3.1) | NS (SST) | off | 壁 p/p0 差は base 比 ≤0.25 % (isoStress +0.2, energyK −0.1, omegaPk +0.1, sigmaBlend 0)。壁 k ピンは 2D で無差 | — | 破棄予定 (A/B 記録) |
| `run_0232_ab3d_{base,all}` (AWS) | 同 A/B の 3D 版 (run_0228 の場から 6000 step) | NS (SST) | off | (plan §3.3) | — | 破棄予定 (A/B 記録) |
| `run_0233_ab2d_nodewd` | node 壁距離をノード座標基準に修正した変換 h5 で run_0213 の場から 3000 step (旧 = 双対重心間距離) | NS (SST) | off | 第一層 wall_dist 最大 +30 % (1.34→1.74 µm)、壁 p/p0 −0.1〜−0.3 %、x=46 の壁 ω ピン 2.0e11→1.2e11、Tmax ≤ Tt。plan turbulence-sst-consistency-options §3.5 | — | 破棄予定 (A/B 記録) |

**結果 (`compare_user_profile.png`)**:
- **NS dry が実験 isentrope に乗る** (x≥11 mm で +1〜2 %; 旧形状 run_0048 は 1 次精度で ±1.5 %)。Euler dry は排除厚がないので −4〜−13 % 下 (下流ほど乖離)。
- **凝縮の効果**: onset は Euler 20.2 mm → NS 23.2 mm (境界層で膨張が緩み過冷却到達が遅れる)。潜熱で中心線 M は出口で 1.91→1.70 (Euler) / 1.79→1.60 (NS)、静圧は +21 % (Euler) / +19 % (NS)、静温 +33 K。液滴質量分率は両者とも x≈50 mm で 0.010 に達しほぼ全量 (Y_H2O 0.01095) が凝縮。
- **実験 cond 1 kPa との比較 (NS)**: onset 帯 (21–32 mm) で −3〜−5 % (forge の bump 立ち上がりが遅く低い)、下流 42–72 mm で +4〜+5 % (旧形状 run_0170 は +5.6〜+7.4 %。形状修正で下流の過大が 2 pt 縮んだ)。onset 域の差は凝縮モデル (Kantrowitz 抑制量) 側の残課題。
- NS 凝縮 run の壁温は凝縮域下流で 291–297 K > Tt=286.65 K (壁 g≈0、コアの潜熱が乱流熱伝導で壁へ届く: 上限 g·L/cp≈26 K の一部)。dry は全域 T≤Tt。
- 軸 h0=(ρe+P)/ρ (Euler cond) は入口/出口で一致 (差 2 J/kg)、onset 帯で最大 −263 J/kg (cp·Tt の 0.09 %)。

## 多成分 TP 発散の再検証 run 一覧 — run_0069 以降 (2026-06-18)

「多成分 thermally perfect gas (thermalMethod:2) が収束しない」件の原因切り分け。**結論 (0acca05 の
「sub-200K のみ・種数無関係」を訂正): 真因は 2 つ**。
1. **IC 生成器の運動エネルギーバグ**: `gen_ic_2sp_from_n2.py:69` が `ek=0.5*(roUx²+..)/ro` (solver 規約は
   `/ro²`)。高速域 (出口 |u|~500, ro~0.1) で roe が ~10⁵ J/kg 過小 → 読み戻しで多数セル T<200K に床張り付き。
2. **nSpecies≥2 コードパス固有の不安定**: ek を修正した整合 IC・全>200K・species 完全保存でも 2 成分は
   step6-7 で局所 ro<0 → NaN。真の単成分は同条件で安定。Y_H2O≈0 でも発散 (step15)、実 H2O で加速 (step7)。
3. sub-200K NASA-9 外挿は別の第 3 要因 (cold で全部落ちる。Tt スイープ `run_0041_hot500`[単成分]安定/
   `run_0042_Tt360`/`run_0043_Tt286` は単成分の話)。

棄却済 (再調査不要): float32 精度 (double 同一)・CFL/剛性 (explicit cfl0.05 でも発散)・implicit 分離更新。
全 run native ビルド (`.build-native/`)。**全て診断用・破棄予定**。

| run | 変更点 | 結果 (初 NaN) | 判定 |
| --- | --- | --- | --- |
| `run_0069_ver_min2sp_float` | 2sp, inviscid+explicit+laminar (最小化) | step4。ro<0・sonic→0 が最冷出口で発生 | 発散機構を特定 |
| `run_0070_ver_min2sp_double` | =0069 を **double** ビルド | step4、残差軌跡が float と6桁一致 | **①float精度 棄却** |
| `run_0071_ver_min_n2o2` | =0069 で H2O→O2 (IC不整合, 参考) | step3 | 参考 (IC不整合) |
| `run_0072_ver_min2sp_implicit` | =0069 を implicit | step4 | **④分離更新 棄却** |
| `run_0073_ver_2sp_zeroH2O` | 整合IC・**Y_H2O≈0** (実質純N2を2sp機構で) | step5 | **②H2O/種数 棄却** |
| `run_0074_ver_2sp_realH2O_consistIC` | 整合IC・Y_H2O=0.01095 | step4 | 2sp は発散 |
| `run_0076_ver_1sp_N2_stablecfg` | **単成分N2**・run_0056 と同一 config | step4 | cold で単成分も発散(交絡: ek+sub200K) |
| `run_0077_ver_2sp_stablecfg` | =0076 を 2sp に | step4 | cold は交絡で切り分け不可 |

注: cold (Tt=286.65) では sub-200K と ek バグの両方で 1sp も落ちるため切り分け不能。**hot (Tt=500K, 全>200K) が
クリーンな対照** (下表)。hot で 1sp は安定・2sp は発散し、種数効果が分離できる。

| run (hot Tt=500K) | 変更点 | 結果 | 判定 |
| --- | --- | --- | --- |
| `run_0081_ver_1sp_hot500_control` | **単成分N2**, res_400 から restart | **exit0・400step 安定・残差フラット** | 単成分 hot は安定 |
| `run_0080_ver_2sp_hot500` | N2+H2O, ek バグ IC | step7。res_0 で 9466 セル T<200 (IC破壊) | 要因1顕在 |
| `run_0083_ver_2sp_n2o2_hot500` | N2+O2 (生成エンタルピー≈0), ek バグ IC | step18 (序盤は残差低下) | O2 でも発散・H2Oより遅い |
| `run_0084_ver_2sp_hot500_fixedIC` | N2+H2O, **ek 修正 IC** | res_0 整合(全>200K, 残差=1sp並)・なお step7 | **要因2: 修正後も発散** |
| `run_0085_ver_2sp_hot500_zeroH2O_fixedIC` | **Y_H2O≈0**, ek 修正 IC | step15 (局所 ro<0・species完全保存) | **要因2はコードパス側** |

### 要因2 の特定 (run_0089〜0094, hot Tt=500K, ek修正IC, condensation無)

O2 と H2O の対照で **H2O 固有**と判明、cfl_pseudo 依存で **implicit 緩和ミスマッチ**と確定:

| run | 2nd種/Y | cfl_pseudo | 結果 |
| --- | --- | --- | --- |
| `run_0089_cmp_1sp` | なし | 2 | 40step 安定 |
| `run_0089_cmp_o2_1e6` | O2 / 1e-6 | 2 | **40step 安定** (h_O2−h_N2≈0) |
| `run_0089_cmp_h2o_1e6` | H2O / 1e-6 | 2 | step13 発散。H2O が入ったセルの T が偽上昇 |
| `run_0089_cmp_h2o_real` | H2O / 0.01095 | 2 | step5 発散 |
| `run_0091_fix_*` | (energy補正カーネル配線) | 2 | **悪化** (O2もstep12)=**二重計上**。デッドカーネルは修正でない |
| `run_0092b_h2o_cflp10` | H2O / 1e-6 | 10 | step2 |
| `run_0092b_h2o_cflp1` | H2O / 1e-6 | **1** | **60step 安定** |
| `run_0092b_h2o_cflp0p1` | H2O / 1e-6 | 0.1 | 安定 |
| `run_0093_h2o_real_cflp1_long` | H2O / 0.01095 | **1** | **400step 完走・T[288,504]・Y1保持** (SST残差はプラトー=単成分と同) |
| `run_0094_cold_real_fixedIC_cflp1` | H2O / 0.01095 **cold Tt=286** | 1 | step13 発散=**要因3(sub-200K)** が cold で残存 |

**結論 (3 要因)**:
1. **IC ek バグ** (`gen_ic_2sp_from_n2.py:72` を `/ro²` に修正済)。
2. **多成分 implicit 結合不安定**: roe は block-DPLUR、roY は別 point-implicit で更新→擬似時間緩和がミスマッチ→roY 変化に roe が対応せず Newton で T ジャンプ。H2O の |h|≈1.3e7 が増幅 (O2 は無害)。**回避策: 多成分 TP は `cfl_pseudo ≤ 1`**。恒久修正は species を block 同梱 or 緩和整合 (未実装)。対流流束自体は整合 (energy 補正配線で二重計上→O2 発散が証拠)。
3. **sub-200K NASA-9 外挿** (cold 固有): 出口<200K の極低温では CPG を使うか thermo 低温拡張。

`run_0070` 用 double バイナリは `.build-native/double/` (要 `FORGE_CUDA_BLOCKSIZE=128`)。デッド energy 補正カーネルは
[`.github/plans/thermophysics-multicomponent-tpgas.md`](../../plans/accepted/thermophysics-multicomponent-tpgas.md) 参照。

### cfl_pseudo パラスタ + 等エントロピー IC + ramp (run_0095〜0100, hot Tt=500K, N2+H2O 1%)

「cfl を上げても収束するか」「初期値を等エントロピーに」の検証。**結論: cfl_pseudo≤1 は硬い上限**。

- **20step 刻み発達** (`run_0095`, cfl_pseudo=1, 600step): 場は定常・物理的 (T[288,504]・sub-200K=0・Y1=0.01095 保持・ro<0=0)。ただし残差はプラトー (rms_ro~1e-7, rms_roe~8e-3、機械ゼロには落ちない=単成分 run_0056 と同じ既存挙動)。cfl=1 vs 0.5 (`run_0096`) で最終場 max|ΔP|≈3.4kPa(~5%)・max|ΔT|≈13K=厳密収束前。
- **cfl_pseudo パラスタ** (warm-start IC, `run_0099_ws_cflp*`):

  | cfl_pseudo | 1 | 2 | 4 | 5 | 8 | 10 |
  |---|---|---|---|---|---|---|
  | 結果 | **安定(399)** | step6 | step3 | step2 | step2 | step2 |

- **等エントロピー IC** (`gen_ic_isentropic.py`: 混合 R/Cp/s0 から (Pt,Tt) 等エントロピー再構成。T[283,509]・sub-200K=0): 熱力学的には妥当だが **連続の式 (ρuA=const) 非整合** (軸速度を sqrt(2Δh) から出したため)。no-slip 粘性壁では壁せん断ショック、inviscid+slip でも質量不整合で **cfl=1 でも発散** (`run_0097`/`run_0098`)。warm-start (mass-consistent) の方が安定。
- **cfl ramp** (`run_0100_ramp_cflp*`: cfl=1 で落ち着いた `run_0099_ws_cflp1/res_400` から高 cfl 再起動): **cfl≥2 はやはり発散** (cfl=2→step7)。過渡限定でなく**構造的上限**。

**まとめ**: 多成分 TP は (1) CEA 有効域 (>200K) の高温条件 + (2) `cfl_pseudo≤1` で**安定・物理的な定常場**が得られる (打ち手は有効)。ただし残差はプラトー品質 (厳密収束・高 cfl には要因2の恒久修正=species を block 同梱/緩和整合が必要)。等エントロピー IC を使うなら ρuA 連続を満たす quasi-1D 構成が要る。

### エンタルピー基準オフセット (sensible-enthalpy datum) による安定化 (run_0101〜0109, 2026-06-19)

**要因2 (多成分 implicit 結合不安定) への打ち手**として、各化学種のエンタルピー基準を
`h_s(Tref=298.15K)=0` へ平行移動する **sensible-enthalpy datum** を実装 (config `physProp.thermoHrefTemp`、
solver は NASA-9 係数 `a7` に焼き込み・全経路整合、Python IC 生成器は `gen_ic_2sp_from_n2.py` 第6引数で同一基準)。
**狙い**: H2O の生成エンタルピー `h_H2O(298K)≈−13.4 MJ/kg` (N2/O2 は ≈0、`thermo_href_compare.png` 参照) が
roe(block-DPLUR)/roY(point-implicit) の擬似時間緩和ミスマッチを増幅して Newton 温度ジャンプを起こすため、
**桁違いの生成エンタルピーを除いて増幅を抑える**。非反応流では物理不変 (e_mix と Σh_s J_s* の基準移動が連続式で相殺)。

同一 hot N2+H2O(1%) 場 (`run_0081/res_400` 由来、全>200K)・同一メッシュ・同一 BC で、絶対基準 (abs) と
オフセット基準 (href) を `cfl_pseudo` スイープで比較 (全 native ビルド、要 `solverConfig.hpp` 変更後の full rebuild)。

| run | 基準 | cfl_pseudo | 結果 (初 NaN step) |
| --- | --- | --- | --- |
| `run_0101_abs_cflp1`   | 絶対       | 1  | **400step 完走** (control = run_0099 再現) |
| `run_0102_abs_cflp2`   | 絶対       | 2  | step6 発散 |
| `run_0103_abs_cflp4`   | 絶対       | 4  | step3 発散 |
| `run_0104_href_cflp1`  | オフセット | 1  | **400step 完走** |
| `run_0105_href_cflp2`  | オフセット | 2  | **400step 完走** ← abs は step6 で落ちる cfl で安定化 |
| `run_0106_href_cflp4`  | オフセット | 4  | step7 発散 (abs step3 より遅延) |
| `run_0107_href_cflp5`  | オフセット | 5  | step4 発散 |
| `run_0108_href_cflp8`  | オフセット | 8  | step3 発散 |
| `run_0109_href_cflp10` | オフセット | 10 | step3 発散 |

**結論 (打ち手は有効だが部分的)**:
1. **安定 `cfl_pseudo` 上限が 1→2 へ向上** (abs は cfl2 で step6 発散、href は cfl2 で 400step 完走)。高 cfl (4,5,8,10) でも発散ステップが一貫して後退 (cfl4: step3→7 等)。**生成エンタルピーの増幅が安定性に効いている直接証拠**。
2. **物理不変を確認**: 同一入力場の step0 再構成 (T,P,ρ,Ux) が abs と href で**機械精度一致** (T rel 1.2e-7, P 6.1e-8, ρ/Ux=0)。roe は基準分だけ平行移動 (`roe_abs−roe_href ≈ ρ·Y_H2O·h_ref(H2O)`, 整合 5%)。`roe/ρ` は abs の **−86〜−19 kJ/kg (負・桁落ち)** から href の **+61〜+127 kJ/kg (sensible・良条件)** へ。
3. **限界**: 上限は 2× 止まりで `cfl_pseudo≥4` は依然発散、残差はプラトー品質 (`check_convergence.py`=NOT CONVERGED、単成分 run_0056 と同等)。**要因2 の構造的ミスマッチ自体は未解消** — オフセットは増幅振幅を下げるだけ。厳密収束・高 cfl には species を 5×5 block 同梱 (full coupling) が引き続き必要。

図: `thermo_href_compare.png` (cp/h vs T・オフセット効果), `href_cfl_ceiling.png` (cfl 上限比較), 各 run の `residual_history.png`。
**run_0101〜0109 は診断用** (恒久保持は不要、結論は本表)。

### 化学種陰解法の緩和整合 (案B, matched-relaxation scalar-DPLUR) 検証 — run_0110 以降 (2026-06-19)

**要因2 の恒久修正候補「案B (緩和整合)」の切り分け**。化学種 `ρY_s` を流れ block-DPLUR と**同一緩和**
(同一凍結残差・`dt_local`・`implicitRelax`・`nStepInner` sweep) のスカラ DPLUR で前進させ、roe(block)/roY(点陰的)
の擬似時間緩和ミスマッチを解消できるかを検証 (config `time.deltaT.speciesImplicitCoupling`, 既定 0=従来 segregated・
ビット不変, 1=緩和整合)。実装は plan [`thermophysics-species-implicit-coupling.md`](../../plans/accepted/thermophysics-species-implicit-coupling.md)、
docs `thermophysics/{theory §3.1, implementation §4b}`。全 native full rebuild。hot N2+H2O(1%, href ON)・同一場/BC。

| run | 化学種更新 | cfl_pseudo | 結果 (初 NaN/最終 step) | 判定 |
| --- | --- | --- | --- | --- |
| `run_0110_match_cflp1`  | 緩和整合 (=1) | 1 | **400step 完走** (NaN無・T[288,504]・ΣY=1) | active (keeper) |
| `run_0111_match_cflp2`  | 緩和整合 (=1) | 2 | **400step 完走** | active (keeper) |
| `run_0112_match_cflp4`  | 緩和整合 (=1) | 4 | step7 発散 | 削除済み(2026-08-31) |
| `run_0113_match_cflp5`  | 緩和整合 (=1) | 5 | step4 発散 | 削除済み(2026-08-31) |
| `run_0114_match_cflp8`  | 緩和整合 (=1) | 8 | step3 発散 | 削除済み(2026-08-31) |
| `run_0115_match_cflp10` | 緩和整合 (=1) | 10 | step3 発散 | 削除済み(2026-08-31) |
| `run_0116_off_cflp1_regr` | 従来 (=0) | 1 | 400step 完走 (default-path 回帰) | 削除済み(2026-08-31) |
| `run_0119_match_cflp3`  | 緩和整合 (=1) | 3 | step20 発散 | 削除済み(2026-08-31) |
| `run_0120_off_cflp3`    | 従来 (=0) | 3 | step20 発散 | 削除済み(2026-08-31) |

**単成分 N2 対照 (H2O 無・案B は構造的に no-op)**:

| run | cfl_pseudo | 結果 | 用途 |
| --- | --- | --- | --- |
| `run_0117_1sp_n2_cflp2` | 2 | 400step 完走 | 単成分回帰 (新旧バイナリ ~atomicAdd 一致) |
| `run_0118_1sp_n2_cflp3` | 3 | step55 発散 | 単成分 cfl3 上限 |
| `run_0117_1sp_n2_cflp4` | 4 | step7 発散 | **単成分でも cfl4=step7** |

**結論 (案B は要因2 の打ち手として無効 — 切り分け完了)**:
1. **緩和整合は cfl 上限を一切動かさない**。matched (=1) と off (=0) の発散 step が**全 cfl で完全一致**
   (cfl≤2 完走 / cfl3 step20=step20 / cfl4 step7=step7 / cfl5 step4 / cfl8,10 step3)。cfl1/2 の残差プラトーも
   matched≡off (`check_convergence.py`=NOT CONVERGED, rms_roK/roOmega 上昇が律速で単成分と同型)。
2. **cfl4 の壁は化学種起因ではない**。**単成分 N2 (H2O 完全に無し) も cfl4 で step7** に発散 (run_0117_cflp4)
   = 2成分と同一。⚠️ **【2026-06-20 訂正】当初ここで「流れ+SST block-DPLUR 律速」と書いたが誤り** — その単成分 N2 も
   `thermalMethod=2` (TP/NASA-9) だった。後続の CPG/TP 対照 (下節「ノズル CFL 上限」run_0121〜0155) で
   **真の律速は `thermalMethod=2` (TP) であり、ノズル形状・乱流モデルではない**ことが判明
   (CPG-SST は同一ノズルで cfl100 まで rms 1e-12 収束)。
3. **案B が無効な機構**: 本ケースの組成は**完全に一様** (Y_N2=0.98905, Y_H2O=0.010950 が全セル同値, min=max)
   → `res_roY≈0` → 化学種補正は matched/segregated どちらでも ~0。要因2 の H2O ペナルティは化学種**移流**の緩和ではなく、
   **EOS エネルギー再構成** (Newton 反転の組成依存 h_mix) にある。移流緩和をそろえても触れない。
4. cfl3 では 2成分 (step20) が単成分 (step55) より早く落ちる=残存する組成-エネルギー結合ペナルティはあるが、
   **緩和率の統一では消えない** (cross-Jacobian ∂roe/∂roY が必要=案A 領域)。ただし cfl4 で単成分も壁に当たるため、
   案A の上積み余地は流れ/SST 上限に頭打ち。**次の高優先は流れ/SST 陰解法の cfl 上限引き上げ** (cf.
   `time_integration-implicit-stable-cfl.md` が flat_plate で cfl120 達成。hot ノズルで効かない理由の調査)。

`speciesImplicitCoupling` は**既定 0 で温存** (本ケースでは無益だが、組成勾配の大きいケースで移流緩和が効く可能性は残すため
削除せず flag-gated)。matched 経路の commit/再正規化/物理量は健全 (run_0110/0111 で NaN無・ΣY=1・T>200K 実証)。
**run_0112〜0120 は診断用** (恒久保持不要、結論は本表)。keeper は run_0110/0111 (緩和整合の安定動作例)。
図: 各 keeper の `residual_history.png`。

## ノズル block-DPLUR CFL 上限 — 真因は thermalMethod (TP) であってノズル/乱流ではない (run_0121〜0155, 2026-06-20)

「多成分/単成分 TP が cfl_pseudo 2〜3 で頭打ち」を、**ノズル形状そのものか・乱流モデルか・thermalMethod(TP)か**で
切り分けた。**初期値は quasi-1D 等エントロピーで生成** (`gen_ic_isentropic_cpg.py`: 面積-マッハ関係から M(x) を解き
P,T,ρ,u を等エントロピーで与える。ρuA=const を機械精度で満たし startup ショックが出ない)。1次風上・blockDPLUR・
nStepInner=20・3000step・detectNaN。binary は `.build-native/relwithdebinfo/forge`。

### Test 1/2: CPG (定数 cp) は乱流の有無に関係なく cfl 100〜200 まで深く収束

| cfl_pseudo | 2 | 5 | 10 | 20 | 50 | 100 | 200 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| **inviscid Euler** (run_0121〜0127) rms_ro fin | 5.3e-12 | 2.0e-12 | 1.9e-12 | 2.0e-12 | 2.2e-12 | 2.4e-12 | **2.7e-12** |
| **SST 粘性乱流** (run_0131〜0136) rms_ro fin | 1.0e-10 | 2.5e-11 | 5.9e-12 | 2.3e-12 | 2.1e-12 | **2.2e-12** | — |

- inviscid (`nozzle_2d.h5` 全slip, `ic_isen_cpg.h5`, Tt=293) は全 cfl で **rms_ro 5 桁低下→1e-12 (float32 floor)・NaN無**。
  exit M=1.97 (等エントロピー 2.0 と一致)。
- SST (`nozzle_2d_sst.h5` no-slip, `ic_isen_cpg_sst.h5`) も同様に **cfl100 まで深く収束** (rms_roOmega も 6.3 桁低下)。
  → **ノズル超音速流も SST 乱流も block-DPLUR の cfl 上限を律速しない**。cfl2-3 の壁は別要因。
- 注: `check_convergence.py` は floor 到達 (5 桁低下後フラット) を "NOT CONVERGED(plateau)" と誤判定する。
  実体は機械精度収束 (TP の 0.5 桁プラトーとは別物)。

### 対照 (decisive): 同一ノズル・同一 SST・同一 hot 等エントロピー IC で thermalMethod だけ変える

`nozzle_2d_sst.h5` + SST + hot 等エントロピー IC (Pt=59070, Tt=500 → 全 T>200K で sub-200K 外挿を排除)。
**唯一の差は `thermalMethod`** (0=CPG / 2=TP NASA-9・単成分 N2)。IC は `ic_isen_cpg_hot.h5` / `ic_isen_tpN2_hot.h5`
(後者は前者の ρ,u を保ち roe のみ NASA-9 N2 の e(T) に置換、`cpg_to_tpN2_roe.py`)。

| cfl_pseudo | 2 | 5 | 20 | 50 | 100 |
| --- | --- | --- | --- | --- | --- |
| **CPG (=0)** (run_0141〜0145) rms_ro fin | 5.4e-11 | 1.3e-11 | 2.5e-12 | 1.5e-12 | **1.4e-12** ✓ |
| **TP/NASA-9 (=2)** (run_0151〜0155) | 9.6e-8 *(plateau)* | **発散 step14** | 発散 step5 | 発散 step3 | 発散 step3 |

**結論 (真因 = thermalMethod=2 の TP 熱力学が block-DPLUR を律速)**:
1. **cfl 上限もプラトーも TP 固有**。同一ノズル・SST・IC で CPG は cfl100 まで rms 1e-12 収束、TP は **cfl2 でも
   1e-7 プラトー (収束しない)・cfl≥5 で発散**。ノズル形状・乱流モデル・初期値 (等エントロピー fresh IC で restart 過渡も排除) は
   すべて無罪。
2. **TP cfl2 のプラトー 9.6e-8 は元の「行き詰まり」runs (run_0093/0101 等) と同レベル** → あの残差プラトーも
   SST ではなく **TP 律速**だった (CPG-SST は機械ゼロまで落ちる)。
3. **機構 (推定)**: block-DPLUR の 5×5 Jacobian は CPG 閉形式 (`chi=(γ-1)/sonic`) を **per-cell 凍結 γ[ic]=γ_mix(T)** で
   流用している。CPG は γ 厳密一定ゆえ Jacobian が厳密で機械収束。TP は γ=γ(T) が変動し、かつ T が毎ステップ NASA-9
   Newton 反転で決まるため、凍結 γ の線形化が実 EOS 微分 (∂P/∂(ρe), ∂P/∂ρ) と不整合 → 線形化誤差が残差プラトーを作り、
   高 cfl で増幅して発散。**恒久修正の方向 = block-DPLUR Jacobian に TP EOS 微分を整合的に入れる** (frozen-γ 流用をやめ、
   ∂P/∂(ρe)=γ_eff−1 等を TP の実微分で評価する)。これは要因2 (組成-エネルギー結合) よりも上位の律速。

run keeper: run_0121/0126 (inviscid), run_0131/0136 (SST), run_0141/0145 (CPG hot), run_0151 (TP hot プラトー実証)。
他は cfl 掃引の診断用 (破棄可)。図: 各 keeper の `residual_history.png`。IC 生成: `gen_ic_isentropic_cpg.py`,
`cpg_to_tpN2_roe.py`。

## 一般EOSヤコビアンで TP 陰解法律速を根治 (run_0161〜0166, 2026-06-20)

上節で「TP の cfl≈2 頭打ち・残差プラトーの真因は block-DPLUR の固有系が CPG 専用 (`h=c²/(γ−1)`,
`∂P/∂ρ|_e=0` を仮定)」と確定したのを受け、**一般EOS固有系**を実装 (plan
[`time_integration-general-eos-jacobian.md`](../../plans/accepted/time_integration-general-eos-jacobian.md))。
音響右固有ベクトルのエネルギーに実全エンタルピー `Ht`、左密度成分に `χ=c²−κh` を入れる**閉形式**3項改変
(`accumulate_split_jacobian_cf`、`generalEOS=thermalMethod==2` で分岐、CPG ビット不変)。閉形式は数値 LU・FD と
機械精度一致を host 単体テスト (`tools/test_eos_jacobian.cpp`) で確認済。

同一ノズル (`nozzle_2d_sst.h5`)・SST・hot 等エントロピー IC (`ic_isen_tpN2_hot.h5`, Tt=500) で:

| cfl_pseudo | 2 | 5 | 20 | 50 | 100 |
| --- | --- | --- | --- | --- | --- |
| **TP 旧** (CPG 形 Jacobian) | 9.6e-8 *(plateau)* | 発散 s14 | 発散 s5 | 発散 s3 | 発散 s3 |
| **TP 新** (一般EOS, run_0161〜0165) rms_ro | 5.6e-11 | **3.3e-11** | 3.7e-11 | 3.9e-11 | **3.8e-11** |

**結論 (TP 律速を根治)**:
1. **cfl 上限 2→≥100**: 新 TP は cfl 2/5/20/50/100 すべて 3000step 完走・NaN 無 (旧は cfl5 で step14 発散)。
2. **残差プラトー突破**: 9.6e-8 → ~4e-11 (**3.6 桁改善**)。rms_roe 3.7 桁・rms_roOmega 5.6 桁低下。物理健全
   (T∈[296,505]K, 超音速, exitM≈1.15=SST 閉塞)。
3. **CPG 回帰** (`run_0166` vs 旧 `run_0143`): rms_ro 2.4816e-12 vs 2.4821e-12 = ~atomicAdd 一致
   (`generalEOS=0` は旧コードと同一・ビット不変)。
4. **残差床は別トラック**: 新 TP は ~4e-11 で下げ止まり (CPG は 1.4e-12)。これは Jacobian でなく NASA Newton
   反転・単精度 EOS 評価の床 (レビュー予測どおり)。後続課題。

keeper: run_0161/0163/0165 (TP 新, cfl2/20/100)・run_0166 (CPG 回帰)。図: 各 `residual_history.png`。
IC: `gen_ic_isentropic_cpg.py` + `cpg_to_tpN2_roe.py`。実装: `cuda_forge/eos_jacobian_d.cuh`,
`timeIntegration_d.cu`, 検証 `tools/test_eos_jacobian.cpp`。
