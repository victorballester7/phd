// Macro parameters
deltaStar = 1;
D = 3.5 * deltaStar; // depth of the gap 
W = 10 * deltaStar; // width of the gap
r = 16; // length of non-constant quads upstream and downstream of the gap
s = 8; // length of non-constant quads left and right inside the gap
BL = 0.25 * deltaStar; // height of the second layer of quad elements
BL_upper = 3.5 * deltaStar; // height of the third layer of quad elements
x0 = 100 * deltaStar; // x distance from inflow to gap
lengthOutflow = 1000 * deltaStar; // length after the gap
x3 = W + lengthOutflow; // last point of the domain
height = 150 * deltaStar; // height of the triangular region
triagHeightRegion_inflow = height - BL - BL_upper; // height of the triangular region
triad_eps = 2 * deltaStar;

// Micro parameters
dx = 0.2; // dx of square elements near the leading and trailing edge of the gap (most resolution-demanding points in the domain)
dy = 1.25 * dx; 
size_triag_v_top = 64; // size of the elements in the upper triangular region
size_quads_v_downstream = 5; 

// All the value below are between 0 and 1. =1  means the elements are all equal in size. Avoid using values very close to 1, because we dividing by log(p) in the formula.
dy_in_v_right = 1.5; // heighest quad in the upper layer of quads, in the right part of the gap
dy_in_v_left = 1.75; // heighest quad in the upper layer of quads, in the left part of the gap

// fit (check the file scripts/fit_gmsh.py for more details)
p_in_v_left = 1.0335663036230878 - 1.0234149444433507/(W - 0.8098798097181781)^0.5991729377115732; // densitiy concentration of elements inside the gap, vertically, near the left wall.
p_in_v_right = 1.008073572191136 - 1.1471875292840257/(W - 1.2941255994826328)^0.7890392764815325; // densitiy concentration of elements inside the gap, vertically, near the right wall.
p_in_h_bottom = 0.9; // densitiy concentration of elements inside the gap, horizontally. 
p_in_h_top = 0.9; // densitiy concentration of elements inside the gap, horizontally. 

p_out_h = 0.5; // densitiy concentration of elements outside the gap, horizontally.
p_out_h_upper = 0.8; // densitiy concentration of elements outside the gap, horizontally on the upper layer of quads
p_out_v_inflow_dense = 0.8; // densitiy concentration of elements outside the gap, vertically, near the inflow and the gap
p_out_v_inflow_sparse = 0.99; // densitiy concentration of elements outside the gap, vertically, near the inflow and far from the gap
p_out_v_outflow = 1; // densitiy concentration of elements outside the gap, vertically, near the outflow and the gap
p_triag_h_inflow = 0.8; // densitiy concentration of elements in the triangular region of the domain, horizontally.
p_triag_h_outflow = p_triag_h_inflow; // densitiy concentration of elements in the triangular region of the domain, horizontally.

enlarge_triag_inflow = 1.1; // heigher the value, larger the elements in the triangular region near the inflow
enlarge_triag_outflow = 1.1; // heigher the value, larger the elements in the triangular region near the outflow

// Automated parameters 
a_in_v_right = dx / (W / 2);
a_in_v_left = dx / (W / 2);
a_in_h_bottom = dy / (D / 2);
a_in_h_top = dx / (D / 2);
a_out_h = dx / BL;
a_out_v_inflow_dense = dx / r;
a_out_v_outflow_dense = dx / r;

// for the derivation of the formula, see file 20250123_gmshFormula.md in the docs folder
N_in_v_right = Ceil( Log (a_in_v_right / (1 - p_in_v_right + a_in_v_right * p_in_v_right)) / Log(p_in_v_right) );
N_in_v_left = Ceil( Log (a_in_v_left / (1 - p_in_v_left + a_in_v_left * p_in_v_left)) / Log(p_in_v_left) );
N_in_h_bottom = Ceil( Log (a_in_h_bottom / (1 - p_in_h_bottom + a_in_h_bottom * p_in_h_bottom)) / Log(p_in_h_bottom) );
N_in_h_top = Ceil( Log (a_in_h_top / (1 - p_in_h_top + a_in_h_top * p_in_h_top)) / Log(p_in_h_top) );
N_out_h = Ceil( Log (a_out_h / (1 - p_out_h + a_out_h * p_out_h)) / Log(p_out_h) ) + 1; // + 1 is to correct the fact that the distance is too small to properly approximate the dx length
N_out_v_inflow_dense = Ceil( Log (a_out_v_inflow_dense / (1 - p_out_v_inflow_dense + a_out_v_inflow_dense * p_out_v_inflow_dense)) / Log(p_out_v_inflow_dense) );

p_for_a_out_v_inflow_sparse_computation = p_out_v_inflow_dense;
a_out_v_inflow_sparse = r / (x0 - r) * (1 - p_for_a_out_v_inflow_sparse_computation) / (1 - p_for_a_out_v_inflow_sparse_computation^(N_out_v_inflow_dense + 1));
N_out_v_inflow_sparse = Ceil( Log(a_out_v_inflow_sparse / (1 - p_out_v_inflow_sparse + a_out_v_inflow_sparse * p_out_v_inflow_sparse)) / Log(p_out_v_inflow_sparse) );

N_out_v_outflow = Ceil( (x3 -W/2 - triad_eps) / size_quads_v_downstream );

// second layer of quads
a_out_h_upper = (BL) / (BL_upper) * (1 - p_out_h) / (1 - p_out_h^(N_out_h + 1));
N_out_h_upper = Ceil( Log(a_out_h_upper / (1 - p_out_h_upper + a_out_h_upper * p_out_h_upper)) / Log(p_out_h_upper) );

// ---------------------------
// for triangular region
// ---------------------------
a_triag_h_inflow = (BL_upper) / (triagHeightRegion_inflow) * (1 - p_out_h_upper) / (1 - p_out_h_upper^(N_out_h_upper + 1));
a_triag_h_outflow = a_triag_h_inflow;
N_triag_v_top = Ceil( (x3 + x0) / size_triag_v_top ); // we assume p_triag_v = 1, otherwise the formula is different
N_triag_h_inflow = Ceil( Log(a_triag_h_inflow / (1 - p_triag_h_inflow + a_triag_h_inflow * p_triag_h_inflow)) / Log(p_triag_h_inflow) );
N_triag_h_outflow = Ceil( Log(a_triag_h_outflow / (1 - p_triag_h_outflow + a_triag_h_outflow * p_triag_h_outflow)) / Log(p_triag_h_outflow) );
N_triag_h_inflow = Ceil( N_triag_h_inflow / enlarge_triag_inflow );
N_triag_h_outflow = Ceil( N_triag_h_outflow / enlarge_triag_outflow );

Point(1) = {-x0, 0, 0};
Point(2) = {0,0,0};
Point(3) = {0,-D,0};
Point(4) = {W/2+triad_eps,-D,0};
Point(5) = {W/2+triad_eps,BL+BL_upper,0};
Point(6) = {x3,0,0};
Point(8) = {x3,height,0};
Point(9) = {-x0,height,0};
Point(10) = {-x0,BL,0};
Point(11) = {0,BL,0};
Point(13) = {W/2,0,0};
Point(14) = {W/2,BL,0};
Point(15) = {0,-D/2,0};
Point(16) = {W/2,-D,0};
Point(18) = {-r,0,0};
Point(19) = {-r,BL,0};
Point(22) = {-x0,BL+BL_upper,0};
Point(23) = {-r,BL+BL_upper,0};
Point(25) = {x3,BL+BL_upper,0};
Point(26) = {0,BL+BL_upper,0};
Point(27) = {W/2,BL+BL_upper,0};
Point(28) = {W/2,-D/2,0};

Line(1) = {1,18};
Line(2) = {2,15};
Line(3) = {3,16};
Line(6) = {6,25};
Line(8) = {8,9};
Line(9) = {9,22};
Line(10) = {10,19};
Line(11) = {10,1};
Line(12) = {11,14};
Line(14) = {2,13};
Line(15) = {15,3};
Line(16) = {16,4};
Line(18) = {4,5};
Line(20) = {18,2};
Line(21) = {14,13};
Line(22) = {13,28};
Line(23) = {19,11};
Line(24) = {22,10};
Line(25) = {25,8};
Line(27) = {22,23};
Line(26) = {4,6};
Line(28) = {5,25};
Line(30) = {23,26};
Line(31) = {26,27};
Line(32) = {27,14};
Line(33) = {28,16};

// outer quads 1
Curve Loop(1) = {-21, -12, -23, -10, 11, 1, 20, 14};
Plane Surface(1) = {1};

// inside gap
Curve Loop(2) = {2, 15, 3, -33, -22, -14};
Plane Surface(2) = {2};

// triangular region
Curve Loop(3) = {25, 8, 9, 27, 30, 31,32,21,22,33,16,18,28};
Plane Surface(3) = {3};

// 3rd layer of quads
Curve Loop(4) = {-24, 27, 30, 31,32,-12,-23,-10};
Plane Surface(4) = {4};

// outer quads downstream
Curve Loop(5) = {18, 28,-6, -26};
Plane Surface(5) = {5};



// This defines the surfaces that will be meshed with quad elements
Transfinite Surface {1} = {14,13, 1, 10};
Transfinite Surface {2} = {3, 13, 16, 2}; // number '2' must match the number '2' of plane surface
Transfinite Surface {4} = {27, 14, 22, 10};
Transfinite Surface {5} = {25, 6, 4, 5};

Transfinite Curve {-14, -3, -12, -31} = N_in_v_left Using Progression p_in_v_left;

Transfinite Curve {15, 33} = N_in_h_bottom Using Progression p_in_h_bottom;
Transfinite Curve {-2, -22} = N_in_h_top Using Progression p_in_h_top;

Transfinite Curve {11, -21} = N_out_h Using Progression p_out_h;
Transfinite Curve {24, 32} = N_out_h_upper Using Progression p_out_h_upper;
Transfinite Curve {10, 1, 27} = N_out_v_inflow_sparse Using Progression p_out_v_inflow_sparse;
Transfinite Curve {20, 23, 30} = N_out_v_inflow_dense Using Progression p_out_v_inflow_dense;
Transfinite Curve {28, 26} = N_out_v_outflow Using Progression p_out_v_outflow;
Transfinite Curve {-18, -6} = N_out_h_upper Using Progression p_out_h_upper;

Recombine Surface {1, 2, 4, 5};

// for triangluar region
Transfinite Curve {8} = N_triag_v_top Using Progression 1;
Transfinite Curve {-25} = N_triag_h_outflow Using Progression p_triag_h_outflow;
Transfinite Curve {9} = N_triag_h_inflow Using Progression p_triag_h_inflow;


// defining boundary
Physical Curve(1) = {9,24,11}; // inlet
Physical Curve(2) = {25,6}; // outlet
Physical Curve(3) = {8}; // top
Physical Curve(4) = {1, 20, 2, 15, 3, 16, 26}; // wall


Physical Surface(100) = {1, 2, 4, 5};
Physical Surface(101) = {3};
