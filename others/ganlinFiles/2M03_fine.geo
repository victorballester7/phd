//==========================================
//+ 主要网格点
Point(1) = {-2, 0, 0, 1.0};
//+
Point(2) = {0, 0, 0, 1.0};
//+
Point(3) = {0, -1, 0, 1.0};
//+
Point(4) = {2, -1, 0, 1.0};
//+
Point(5) = {2, 0, 0, 1.0};
//+
Point(6) = {6, 0, 0, 1.0};
//+
Point(7) = {6, 6, 0, 1.0};
//+
Point(8) = {-2, 6, 0, 1.0};
//==========================================
//+辅助网格点
Point(9) = {4, 0, 0, 1.0};
//+
Point(10) = {-2, 0.6, 0, 1.0};
//+
Point(11) = {0, 0.6, 0, 1.0};
//+
Point(12) = {2, 0.6, 0, 1.0};
//+
Point(13) = {4, 0.6, 0, 1.0};
//+
Point(14) = {6, 0.6, 0, 1.0};
//+
Point(15) = {0, -0.5, 0, 1.0};
//+
Point(16) = {2, -0.5, 0, 1.0};
//+
Point(17) = {1, -1, 0, 1.0};
//+
Point(18) = {1, 0, 0, 1.0};
//+
Point(19) = {1, 0.6, 0, 1.0};
//==========================================
Line(1) = {2, 1};
//+
Line(2) = {5, 9};
//+
Line(3) = {2, 18};
//+
Line(4) = {5, 18};
//+
Line(5) = {2, 15};
//+
Line(6) = {5, 16};
//+
Line(7) = {3, 15};
//+
Line(8) = {4, 16};
//+
Line(9) = {3, 17};
//+
Line(10) = {4, 17};
//+
Line(11) = {9, 6};
//+
Line(12) = {1, 10};
//+
Line(13) = {6, 14};
//+
Line(14) = {10, 8};
//+
Line(15) = {14, 7};
//+
Line(16) = {7, 8};
//+
Line(17) = {11, 10};
//+
Line(18) = {12, 13};
//+
Line(19) = {11, 19};
//+
Line(20) = {12, 19};
//+
Line(21) = {13, 14};
//==========================================
//+
Curve Loop(1) = {9, -10, 8, -6, 4, -3, 5, -7};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {12, -17, 19, -20, 18, 21, -13, -11, -2, 4, -3, 1};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {14, -16, -15, -21, -18, 20, -19, 17};
//+
Plane Surface(3) = {3};
//==========================================
//+
Physical Curve("inflow") = {12, 14};
//+
Physical Curve("outflow") = {13, 15};
//+
Physical Curve("top") = {16};
//+
Physical Curve("wall_main") = {1, 2, 11};
//+
Physical Curve("gap_bottom") = {9, 10};
//+
Physical Curve("gap_sides") = {5, 7, 8, 6};
//==========================================
//+ Gap结构网格
Transfinite Surface {1} = {3, 4, 5, 2};
//+
Transfinite Surface {2} = {1, 6, 14, 10};
//+
Transfinite Curve {3, 4, 9, 10, 19, 20} = 46 Using Progression 1.05;
//+
Transfinite Curve {6, 8, 7, 5} = 35 Using Progression 1.12;
//+
Recombine Surface {1};
//==========================================
//+ 主壁面x方向网格分布
Transfinite Curve {1, 17} = 24 Using Progression 1.2;
//+
Transfinite Curve {2, 18} = 40 Using Progression 1.1;
//+
Transfinite Curve {21, 11} = 12 Using Progression 1.0;
//+ 主壁面y方向网格分布
Transfinite Curve {12, 13} = 25 Using Progression 1.2;
//+
Recombine Surface {2};
//==========================================
//+ 外部网格分布
//+ 顶部
Transfinite Curve {16} = 8 Using Progression 1;
//+ 侧边
Transfinite Curve {14, 15} = 10 Using Progression 1.2;
//==========================================
//+ Nektar网格定义
//+ Define surface for NekMesh
Physical Surface(100) = {1};
Physical Surface(101) = {2};
Physical Surface(102) = {3};
