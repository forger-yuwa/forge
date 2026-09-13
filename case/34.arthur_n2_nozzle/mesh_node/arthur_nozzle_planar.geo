// =====================================================================
// Arthur (1952) 2D source-flow hypersonic nozzle — diverging section only
//   再現対象: P.D. Arthur PhD thesis (Caltech/GALCIT, 1952), II.A Apparatus
//   - 2 次元 source-flow ノズル: スロート 0.010 in (=0.254 mm 全高),
//     出口 ~1 in (=25.4 mm 全高), 全開き角 11 deg (half-angle 5.5 deg)。
//   - 本メッシュは上半分のみ (中心線 y=0 を対称面 slip でモデル化)。
//   - 発散部のみを切り出し、入口をスロートに置く (超音速インレットで M=1.05
//     を固定 → 収束部・choking 不要)。
//   - 壁形状: 双曲線  y(x) = sqrt(yt^2 + (x*tan a)^2)
//       * x=0 (スロート) で dy/dx=0 (壁が水平接線) → 鋭い凸コーナーが無く
//         入口の膨張特異点を回避できる。
//       * x が増えると y -> x*tan a に漸近 = Arthur の source-flow 直線壁。
//       * スロート曲率半径 R_t = yt/tan^2 a ≈ 13.7 mm (原典に記載が無いため
//         この滑らかな双曲線で代表させる。下流の膨張は直線 5.5 deg と一致)。
//   - 単位: メートル (SI)。
// =====================================================================

yt   = 0.000127;     // スロート半高 [m]
ye   = 0.0127;       // 出口半高   [m]
tana = 0.0962992;    // tan(5.5 deg)
xe   = Sqrt(ye*ye - yt*yt)/tana;   // 出口軸位置 [m] (双曲線で y=ye となる x, ~0.13187 m)

nx = 400;            // 軸方向セル数
ny = 60;             // 横方向セル数
rp = 1.008;          // throat 側へ寄せる progression 比
nw = 400;            // 壁スプライン点数

Point(1) = {0.0, 0.0, 0.0};   // スロート・中心線
Point(2) = {xe , 0.0, 0.0};   // 出口・中心線
Point(3) = {xe , ye , 0.0};   // 出口・壁
Point(4) = {0.0, yt , 0.0};   // スロート・壁

// 壁スプライン (throat_wall P4 -> exit_wall P3), 双曲線
pW[0] = 4;
For i In {1 : nw-1}
  xq = xe * i/(nw-1);
  yq = Sqrt(yt*yt + (xq*tana)*(xq*tana));
  pW[i] = newp;
  Point(pW[i]) = {xq, yq, 0.0};
EndFor
pW[nw-1] = 3;   // 終点は厳密に exit_wall

Line(1)   = {1, 2};      // 中心線 (+x)
Line(2)   = {2, 3};      // 出口   (+y)
wall = newl; Spline(wall) = pW[];   // 壁 P4->P3 (+x)
Line(4)   = {4, 1};      // 入口   (-y)

Curve Loop(1) = {1, 2, -wall, 4};
Plane Surface(1) = {1};

Transfinite Curve {1}    = nx Using Progression rp;   // 中心線: throat で小
Transfinite Curve {wall} = nx Using Progression rp;   // 壁:     throat で小
Transfinite Curve {2, 4} = ny;
Transfinite Surface {1} = {1, 2, 3, 4};
Recombine Surface {1};

// ---- 平面 2D (node/median-dual 用: 押し出し無し) ----
//   node は 2 ノード spanwise 押し出しで MUSCL 散逸が消えるため平面メッシュ必須 (memory: node-2node-spanwise)
Physical Curve("inlet" , 1) = { 4 };
Physical Curve("outlet", 2) = { 2 };
Physical Curve("wall"  , 3) = { 1, wall };   // 中心線 (対称 slip) + 壁 (slip)
Physical Surface("fluid", 4) = { 1 };
