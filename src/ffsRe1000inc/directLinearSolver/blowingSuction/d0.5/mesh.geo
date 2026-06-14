// Macro parameters
deltaStar = 1;
D = 0.5 * deltaStar; // upstream channel depth below y = 0
r = 16 * deltaStar; // refinement extent around the step corner (both sides)
BL = 0.25 * deltaStar; // first quad layer thickness
BL_upper = 3.5 * deltaStar; // second quad layer thickness
x0 = 100 * deltaStar; // inflow location is x = -x0
xStep = 100 * deltaStar; // forward-facing step location
lengthOutflow = 1000 * deltaStar; // length after the step
x3 = xStep + lengthOutflow; // outlet location
xA = xStep - r;
xB = xStep + r;
height = 150 * deltaStar; // top boundary y-coordinate
triagHeightRegion = height - BL - BL_upper;

// Micro parameters
dx = 0.2;
dy = 1.25 * dx;
size_triag_v_top = 64;
size_quads_h = 5;

// Progression parameters
p_corner_h = 0.8; // streamwise grading towards the corner from both sides
p_step_v = 0.9;
p_bl = 0.5;
p_bl_upper = 0.8;
p_upstream_sparse = 0.99;
p_downstream_sparse = 1;
p_triag_h_inflow = 0.8;
p_triag_h_outflow = p_triag_h_inflow;
enlarge_triag_inflow = 1.1;
enlarge_triag_outflow = 1.1;

// Automated parameters
L_upstream_sparse = xA + x0;
L_downstream_sparse = x3 - xB;

a_corner_h = dx / r;
N_corner_h = Max(2, Ceil(Log(a_corner_h / (1 - p_corner_h + a_corner_h * p_corner_h)) / Log(p_corner_h)));

a_step_v = dy / D;
N_step_v = Max(2, Ceil(Log(a_step_v / (1 - p_step_v + a_step_v * p_step_v)) / Log(p_step_v)));

a_bl = dx / BL;
N_bl = Max(2, Ceil(Log(a_bl / (1 - p_bl + a_bl * p_bl)) / Log(p_bl)) + 1);

a_bl_upper = (BL / BL_upper) * (1 - p_bl) / (1 - p_bl^N_bl);
N_bl_upper = Max(2, Ceil(Log(a_bl_upper / (1 - p_bl_upper + a_bl_upper * p_bl_upper)) / Log(p_bl_upper)));

N_upstream_sparse = Max(2, Ceil(L_upstream_sparse / size_quads_h));
N_downstream_sparse = Max(2, Ceil(L_downstream_sparse / size_quads_h));

a_triag_h = (BL_upper / triagHeightRegion) * (1 - p_bl_upper) / (1 - p_bl_upper^N_bl_upper);
N_triag_v_top = Max(2, Ceil((x3 + x0) / size_triag_v_top));
N_triag_h_inflow = Max(2, Ceil(Log(a_triag_h / (1 - p_triag_h_inflow + a_triag_h * p_triag_h_inflow)) / Log(p_triag_h_inflow)));
N_triag_h_outflow = Max(2, Ceil(Log(a_triag_h / (1 - p_triag_h_outflow + a_triag_h * p_triag_h_outflow)) / Log(p_triag_h_outflow)));
N_triag_h_inflow = Ceil(N_triag_h_inflow / enlarge_triag_inflow);
N_triag_h_outflow = Ceil(N_triag_h_outflow / enlarge_triag_outflow);

// Geometry
Point(1) = {-x0, 0, 0};
Point(2) = {xA, 0, 0};
Point(3) = {xStep, 0, 0};
Point(4) = {-x0, D, 0};
Point(5) = {xA, D, 0};
Point(6) = {xStep, D, 0};
Point(7) = {xB, D, 0};
Point(8) = {x3, D, 0};
Point(9) = {-x0, BL+D, 0};
Point(10) = {xA, BL+D, 0};
Point(11) = {xStep, BL+D, 0};
Point(12) = {xB, BL+D, 0};
Point(13) = {x3, BL+D, 0};
Point(14) = {-x0, BL+D + BL_upper, 0};
Point(15) = {xA, BL+D + BL_upper, 0};
Point(16) = {xStep, BL+D + BL_upper, 0};
Point(17) = {xB, BL+D + BL_upper, 0};
Point(18) = {x3, BL+D + BL_upper, 0};
Point(19) = {-x0, height+D, 0};
Point(20) = {x3, height+D, 0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 6};
Line(4) = {6, 7};
Line(5) = {7, 8};
Line(6) = {8, 13};
Line(7) = {13, 18};
Line(8) = {18, 20};
Line(9) = {20, 19};
Line(10) = {19, 14};
Line(11) = {14, 9};
Line(12) = {9, 4};
Line(13) = {4, 1};
Line(14) = {4, 5};
Line(15) = {5, 6};
Line(16) = {2, 5};
Line(17) = {5, 10};
Line(18) = {10, 15};
Line(19) = {6, 11};
Line(20) = {11, 16};
Line(21) = {7, 12};
Line(22) = {12, 17};
Line(23) = {9, 10};
Line(24) = {10, 11};
Line(25) = {11, 12};
Line(26) = {12, 13};
Line(27) = {14, 15};
Line(28) = {15, 16};
Line(29) = {16, 17};
Line(30) = {17, 18};

// lower quads (upstream only)
Curve Loop(1) = {1, 16, -14, 13};
Plane Surface(1) = {1};

Curve Loop(2) = {2, 3, -15, -16};
Plane Surface(2) = {2};

// first quad layer
Curve Loop(3) = {14, 17, -23, 12};
Plane Surface(3) = {3};

Curve Loop(4) = {15, 19, -24, -17};
Plane Surface(4) = {4};

Curve Loop(5) = {4, 21, -25, -19};
Plane Surface(5) = {5};

Curve Loop(6) = {5, 6, -26, -21};
Plane Surface(6) = {6};

// second quad layer
Curve Loop(7) = {23, 18, -27, 11};
Plane Surface(7) = {7};

Curve Loop(8) = {24, 20, -28, -18};
Plane Surface(8) = {8};

Curve Loop(9) = {25, 22, -29, -20};
Plane Surface(9) = {9};

Curve Loop(10) = {26, 7, -30, -22};
Plane Surface(10) = {10};

// upper triangular region
Curve Loop(11) = {30, 8, 9, 10, 27, 28, 29};
Plane Surface(11) = {11};

// quad surfaces
Transfinite Surface {1} = {1, 2, 5, 4};
Transfinite Surface {2} = {2, 3, 6, 5};
Transfinite Surface {3} = {4, 5, 10, 9};
Transfinite Surface {4} = {5, 6, 11, 10};
Transfinite Surface {5} = {6, 7, 12, 11};
Transfinite Surface {6} = {7, 8, 13, 12};
Transfinite Surface {7} = {9, 10, 15, 14};
Transfinite Surface {8} = {10, 11, 16, 15};
Transfinite Surface {9} = {11, 12, 17, 16};
Transfinite Surface {10} = {12, 13, 18, 17};

// streamwise refinement: corner from both sides
Transfinite Curve {1, 14, 23, 27} = N_upstream_sparse Using Progression p_upstream_sparse;
Transfinite Curve {2, 15, 24, 28} = N_corner_h Using Progression p_corner_h;
Transfinite Curve {-4, -25, -29} = N_corner_h Using Progression p_corner_h;
Transfinite Curve {5, 26, 30} = N_downstream_sparse Using Progression p_downstream_sparse;

// wall-normal and step-normal refinement
Transfinite Curve {3, 16, -13} = N_step_v Using Progression p_step_v;
Transfinite Curve {12, -17, -19, -21, -6} = N_bl Using Progression p_bl;
Transfinite Curve {11, -18, -20, -22, -7} = N_bl_upper Using Progression p_bl_upper;

Recombine Surface {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

// triangular cap
Transfinite Curve {9} = N_triag_v_top Using Progression 1;
Transfinite Curve {10} = N_triag_h_inflow Using Progression p_triag_h_inflow;
Transfinite Curve {-8} = N_triag_h_outflow Using Progression p_triag_h_outflow;

// boundaries
Physical Curve(1) = {10, 11, 12, 13}; // inlet
Physical Curve(2) = {6, 7, 8}; // outlet
Physical Curve(3) = {9}; // top
Physical Curve(4) = {1, 2, 3, 4, 5}; // wall

Physical Surface(100) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
Physical Surface(101) = {11};
