# V-ax0 追加の負例と判別 A/B (2026-09-27、ローカル float ビルド build-ypls = commit 後の output.cpp を含む、FORGE_CUDA_BLOCKSIZE=128、各 1 step)

## 軸上の辺 (case/61 run_0005 の入力を y −10 mm 平行移動し、|y|<1e-9 を厳密に 0 に丸め = 界面 5 節点すべてが r=0)
```
118:[conjugateWall] ERROR: physID 4: 軸対称の界面節点 5 個が軸上 (r ≤ r_floor = 1e-20 m)。最初: 界面節点 0 (0.002, 0)。
rc=1
```
(丸め前は float32 の y=−2.2e-10 で「r<0 の節点が 5 個」の拒否に先に掛かった: 軸上判定の分岐を通すため丸めた)

## method 属性 A/B (case/62 run_0003_dry_r32u の入力、CHT なし・interfaceDiag 1、axisymMethod 0 / 1 だけ変更)
```
axisymMethod 0 res_wall_cj_4_1.h5 iface_q_eff: all finite=True all NaN=False status=None
axisymMethod 0 res_wall_cj_4_1.h5 iface_q_eff_raw: all finite=True all NaN=False status=None
axisymMethod 0 res_wall_cj_4_1.h5 iface_Qf_eff: all finite=True all NaN=False status=None
axisymMethod 0 res_wall_hot_3_1.h5 iface_q_eff: all finite=True all NaN=False status=None
axisymMethod 0 res_wall_hot_3_1.h5 iface_q_eff_raw: all finite=True all NaN=False status=None
axisymMethod 0 res_wall_hot_3_1.h5 iface_Qf_eff: all finite=True all NaN=False status=None
axisymMethod 1 res_wall_cj_4_1.h5 iface_q_eff: all finite=False all NaN=True status='unverified: axisymMethod != 0 (planar geometry with 1/y source); value is NaN'
axisymMethod 1 res_wall_cj_4_1.h5 iface_q_eff_raw: all finite=False all NaN=True status='unverified: axisymMethod != 0 (planar geometry with 1/y source); value is NaN'
axisymMethod 1 res_wall_cj_4_1.h5 iface_Qf_eff: all finite=False all NaN=True status='unverified: axisymMethod != 0 (planar geometry with 1/y source); value is NaN'
axisymMethod 1 res_wall_hot_3_1.h5 iface_q_eff: all finite=False all NaN=True status='unverified: axisymMethod != 0 (planar geometry with 1/y source); value is NaN'
axisymMethod 1 res_wall_hot_3_1.h5 iface_q_eff_raw: all finite=False all NaN=True status='unverified: axisymMethod != 0 (planar geometry with 1/y source); value is NaN'
axisymMethod 1 res_wall_hot_3_1.h5 iface_Qf_eff: all finite=False all NaN=True status='unverified: axisymMethod != 0 (planar geometry with 1/y source); value is NaN'
```
