// V1 共役スラブ: 一様構造 quad (physID 1 左 / 2 右 / 3 下=共役壁 / 4 上=等温壁)
W = 0.02; H = 0.01; nx = 20; ny = 16;
Point(1) = {0, 0, 0, 1};  Point(2) = {W, 0, 0, 1};
Point(3) = {W, H, 0, 1};  Point(4) = {0, H, 0, 1};
Line(1) = {1, 2}; Line(2) = {2, 3}; Line(3) = {3, 4}; Line(4) = {4, 1};
Line Loop(1) = {1, 2, 3, 4}; Plane Surface(1) = {1};
Transfinite Line{1, 3} = nx + 1;
Transfinite Line{2, 4} = ny + 1;
Transfinite Surface{1}; Recombine Surface{1};
Physical Curve("side_left",  1) = {4};
Physical Curve("side_right", 2) = {2};
Physical Curve("wall_bot",   3) = {1};
Physical Curve("wall_top",   4) = {3};
Physical Surface("fluid",    5) = {1};
Mesh.MshFileVersion = 4.1;
