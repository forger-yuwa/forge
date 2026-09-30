# 諮問: 共役平板 B1 予備確認の「採取間隔の検定」FAIL をどう扱うか

AGENTS.md 条件 3 (事前に書いた判定が FAIL)・4 (plan に無い手順変更の前)。plan `plans/active/boundary-cht-conjugate-flat-plate.md` §4.1.1・§4.1.2 (commit d57e8309)。

## 観測事実 (AWS FP64 `86115cb1`、md5 159220e4…、run は case/65.conjugate_flat_plate/)
- 予備確認 `run_0013_pre_c1_n64_coupled` (warmup 5000) / `run_0014_pre_c1_n64_fixed_tw` (warmup 900000)。どちらも `run_0007_c1_n64` の res_600000.h5 (restart_field --keep-src-dtype、ビット一致 OK) + conjugate_state_5.h5 + conjugate_Tw_5.csv→wall_profile_5.csv で 200 step、毎 step 出力。ログ: 「conjugate_state_5.h5 から固体温度と更新位相を復元 (step 600000)」「applyWallProfiles … applied: Ts」。
- 判定ツール `case/65.conjugate_flat_plate/ab_series.py pre` (本 brief と同時に追加、未 commit)。出力 `run_0013_pre_c1_n64_coupled/AB_PRE_CHECK.txt`。
- 受入条件: 全 PASS — 初期壁温の A/B 差 0、B の初期壁温 − wall_profile 0、B の壁温の時間変化 0、A は更新 4 回 (step 0,50,100,150)・壁温変化 max 2.1e-3 K。
- 採取間隔の検定 (毎 step 200 点 vs 11 step 間引き 18〜19 点の、節点群ごとの標準偏差の節点最大、10 % 以内): **FAIL**。差は 0.1〜40 %。FAIL は le/te/up の P (14〜41 %)、up の T・Uy (14〜19 %)、B の iface_le q (10.1 %) 等。down 群や iface_te は 0〜7 %。
  - 18 点の標本標準偏差の相対標準誤差は約 1/√(2·17) ≈ 17 % なので、差の多くは標本数で説明できる規模 (エイリアシングの証拠とは限らない)。
- **目を引く事実**: 界面熱流束の前縁帯 (x/L≤0.02) の時間変動 (標準偏差/窓内平均 q) は A 0.1245、**B (壁温固定) 0.1206** でほぼ同じ。200 step だけの観測。
- ディスク: AWS 残り約 25 GB。流体 res は 1 枚 3 MB、壁ダンプ 75 KB。60000 step を 11 step 間隔で流体も出すと 1 run 16 GB で 2 本は載らない。壁ダンプは流体出力と同じ間隔でしか出ない。
- forge の point probe (`probe.yaml` の points、KD-tree で最寄り節点) は毎 step `point_probe_<i>.out` に追記できる (`solver_density_cuda/probe/point_probes.cu`)。節点座標を与えれば流体節点群 (164 点) の P・T・Uy は毎 step 採れる。界面熱流束 (iface_q_eff) は probe では採れない。

## 私の案
1. 流体節点群 (164 点) は point probe で**毎 step** 採る (採取間隔の問題を無くす)。
2. 界面熱流束は壁ダンプでしか採れない。流体+壁の出力間隔を 101 step (素数、50 と互いに素、594 枚 ≈ 1.8 GB/run) にする。検定は、別に 2020 step の毎 step 予備確認 (約 6 GB、判定後に削除) で、毎 step と 101 step 間引き (20 点) では同じ問題が出るので、**標準偏差の比較ではなく、毎 step 系列の自己相関時間 / スペクトルから 101 step 間引きの標準偏差推定の偏りを出す**形に検定を変える…が、事前登録を事後に変えることになる。
3. もしくは検定の基準を「差 ≤ 10 %」から「標本数から決まる 95 % 区間内」に改める。

## 質問
1. 予備確認の FAIL をどう扱うか (登録どおり本段に進まない。検定の再設計は事後改訂として記録)。
2. 本段の採取方法として上の 1・2 は妥当か。代案は。
3. 前縁帯の界面熱流束の変動が A/B でほぼ同じという 200 step の観測を、本段の判別に活かす形はあるか (これだけで結論しない)。
