// =====================================================================
// 超音速冷却壁平板 (case/48) — 平面 2D 構造メッシュ (node 用, 押し出しなし)
//   領域: x in [-0.1, 1.0], y in [0, H]。x<0 は slip 助走, x>=0 が no-slip 平板
//   入口 x=-0.1 (超音速一様), 出口 x=1.0, 上面 y=H slip
//   壁垂直: 第一セル高さ Y1 (テンプレート置換), 等比 R_Y, セル数 NY (H に届く値を gen_mesh.py が決める)
//   streamwise: 平板 NX_PLATE セル (前縁側細かく r=R_PLATE), 助走 NX_UP セル
//   派生元: case/26 mesh/flat_plate_planar.geo
// =====================================================================
Geometry.PointNumbers = 0;
Mesh.ScalingFactor = 1.0;
lc = 0.05;
x0 = -0.1;  x1 = 0.0;  x2 = 1.0;  H = __H__;
nx_up = __NX_UP__;  r_up = 1.08;
nx_plate = __NX_PLATE__;  r_plate = __R_PLATE__;
ny = __NY__;  r_y = __R_Y__;
Point(1) = {x0, 0.0, 0.0, lc};
Point(2) = {x1, 0.0, 0.0, lc};
Point(3) = {x2, 0.0, 0.0, lc};
Point(4) = {x2, H,   0.0, lc};
Point(5) = {x1, H,   0.0, lc};
Point(6) = {x0, H,   0.0, lc};
Line(1) = {1, 2};   // bottom slip (助走)
Line(2) = {2, 3};   // bottom wall (平板)
Line(3) = {3, 4};   // outlet
Line(4) = {4, 5};   // top (平板上)
Line(5) = {5, 2};   // 前縁の垂直線 (上->下)
Line(6) = {5, 6};   // top (助走上)
Line(7) = {6, 1};   // inlet (上->下)
Transfinite Line {7} = ny Using Progression 1.0/r_y;
Transfinite Line {5} = ny Using Progression 1.0/r_y;
Transfinite Line {3} = ny Using Progression r_y;
Transfinite Line {1} = nx_up Using Progression 1.0/r_up;
Transfinite Line {6} = nx_up Using Progression r_up;
Transfinite Line {2} = nx_plate Using Progression r_plate;
Transfinite Line {4} = nx_plate Using Progression 1.0/r_plate;
Curve Loop(1) = {7, 1, -5, 6};
Plane Surface(1) = {1};
Transfinite Surface {1};
Recombine Surface(1);
Curve Loop(2) = {5, 2, 3, 4};
Plane Surface(2) = {2};
Transfinite Surface {2};
Recombine Surface(2);
Physical Curve("inlet",  1) = {7};
Physical Curve("outlet", 2) = {3};
Physical Curve("top",    3) = {4, 6};
Physical Curve("wall",   4) = {2};
Physical Curve("sym",    5) = {1};
Physical Surface("fluid", 8) = {1, 2};
