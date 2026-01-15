// Macro parameters
deltaStar = 1;
D = 4 * deltaStar; // depth of the gap 
W = 15 * deltaStar; // width of the gap
BL = 0.5 * deltaStar; // height of the second layer of quad elements
x0 = 100 * deltaStar; // x distance from inflow to gap
lengthOutflow = 150 * deltaStar; // length after the gap
x1 = 2 * BL; // first point of quads in the domain
lengthQuadsOutflow = 3 * BL;
x2 = W + lengthQuadsOutflow; // last point of quads in the domain
x3 = W + lengthOutflow; // last point of the domain
triagHeightRegion = 120 * deltaStar; // height of the triangular region
h = triagHeightRegion; // height of the domain


// Micro parameters
dx = 0.05; // dx of square elements near the leading and trailing edge of the gap (most resolution-demanding points in the domain)
size_triag_v_top = 8; // size of the elements in the upper triangular region
size_triag_v_inflow = 1; // size of the elements in the lower triangular region (inflow part)
size_triag_v_outflow = 4; // size of the elements in the lower triangular region (outflow part)
size_triag_h_inflow = 4; // size of the elements in the horizontal direction in the triangular region (inflow part)
size_triag_h_outflow = 6; // size of the elements in the horizontal direction in the triangular region (outflow part)

// All the value below are between 0 and 1. =1  means the elements are all equal in size. Avoid using values very close to 1, because we dividing by log(p) in the formula.
p_in_v_right = 0.87; // densitiy concentration of elements inside the gap, vertically, near the right wall.
p_in_v_left = 0.8; // densitiy concentration of elements inside the gap, vertically, near the left wall.
p_in_h = 0.9; // densitiy concentration of elements inside the gap, horizontally. 

p_out_h = 0.75; // densitiy concentration of elements outside the gap, horizontally.
p_out_v_inflow = 0.9; // densitiy concentration of elements outside the gap, vertically, near the inflow.
p_out_v_outflow = 0.9; // densitiy concentration of elements outside the gap, vertically, near the outflow.
p_triag_h_inflow = 0.9; // densitiy concentration of elements in the triangular region of the domain, horizontally.
p_triag_h_outflow = p_triag_h_inflow; // densitiy concentration of elements in the triangular region of the domain, horizontally.
p_triag_v_inflow = 0.96; // densitiy concentration of elements in the triangular region of the domain, vertically, near the inflow.
p_triag_v_outflow = 0.93; // densitiy concentration of elements in the triangular region of the domain, vertically, near the outflow.


// Automated parameters 
a_in_v_right = dx / (W / 2);
a_in_v_left = dx / (W / 2);
a_in_h = dx / (D / 2);
a_out_h = dx / BL;
a_out_v_inflow = dx / x1;
a_out_v_outflow = dx / lengthQuadsOutflow;

// for the derivation of the formula, see file 20250123_gmshFormula.md in the docs folder
N_in_v_right = Ceil( Log (a_in_v_right / (1 - p_in_v_right + a_in_v_right * p_in_v_right)) / Log(p_in_v_right) );
N_in_v_left = Ceil( Log (a_in_v_left / (1 - p_in_v_left + a_in_v_left * p_in_v_left)) / Log(p_in_v_left) );
N_in_h = Ceil( Log (a_in_h / (1 - p_in_h + a_in_h * p_in_h)) / Log(p_in_h) );
N_out_h = Ceil( Log (a_out_h / (1 - p_out_h + a_out_h * p_out_h)) / Log(p_out_h) );
N_out_v_inflow = Ceil( Log (a_out_v_inflow / (1 - p_out_v_inflow + a_out_v_inflow * p_out_v_inflow)) / Log(p_out_v_inflow) );
N_out_v_outflow = Ceil( Log (a_out_v_outflow / (1 - p_out_v_outflow + a_out_v_outflow * p_out_v_outflow)) / Log(p_out_v_outflow) );

// for triangular region
N_triag_v_top = Ceil( (x3 - x0) / size_triag_v_top ); // we assume p_triag_v = 1, otherwise the formula is different
N_triag_v_inflow = Ceil( Abs(x1 - x0) / size_triag_v_inflow );
N_triag_v_outflow = Ceil( Abs(x3 - x2) / size_triag_v_outflow );
N_triag_h_inflow = Ceil( h / size_triag_h_inflow );
N_triag_h_outflow = Ceil( h / size_triag_h_outflow );


Point(1) = {-x0, 0, 0};
Point(2) = {0,0,0};
Point(3) = {0,-D,0};
Point(4) = {W,-D,0};
Point(5) = {W,0,0};
Point(6) = {x3,0,0};
Point(8) = {x3,h,0};
Point(9) = {-x0,h,0};
Point(11) = {0,BL,0};
Point(12) = {W,BL,0};
Point(13) = {W/2,0,0};
Point(14) = {W/2,BL,0};
Point(15) = {0,-D/2,0};
Point(16) = {W/2,-D,0};
Point(17) = {W,-D/2,0};
Point(18) = {-x1,0,0};
Point(19) = {-x1,BL,0};
Point(20) = {x2,0,0};
Point(21) = {x2,BL,0};

Line(1) = {1,18};
Line(2) = {2,15};
Line(3) = {3,16};
Line(4) = {4,17};
Line(5) = {5,20};
Line(6) = {6,8};
Line(8) = {8,9};
Line(9) = {9,1};
Line(11) = {19,11};
Line(12) = {11,14};
Line(13) = {12,21};
Line(14) = {2,13};
Line(15) = {15,3};
Line(16) = {16,4};
Line(17) = {17,5};
Line(18) = {13,5};
Line(19) = {14,12};
Line(20) = {20,6};
Line(21) = {20,21};
Line(22) = {18,2};
Line(23) = {18,19};

Curve Loop(1) = {5, 21, -13, -19, -12, -11, -23, 22, 14, 18};
Plane Surface(1) = {1};

Curve Loop(2) = {2, 15, 3, 16, 4, 17, -18, -14};
Plane Surface(2) = {2};

Curve Loop(3) = {1, 23, 11, 12, 19, 13, -21, 20, 6, 8, 9};
Plane Surface(3) = {3};


// This defines the surfaces that will be meshed with quad elements
Transfinite Surface {1} = {18, 20, 21, 19}; // number '2' must match the number '2' of plane surface
Transfinite Surface {2} = {3, 4, 5, 2}; // number '1' must match the number '1' of plane surface

Transfinite Curve {18, 16, 19} = N_in_v_right Using Progression p_in_v_right;
Transfinite Curve {-14, -3, -12} = N_in_v_left Using Progression p_in_v_left;
Transfinite Curve {-2, 15, 17, -4} = N_in_h Using Progression p_in_h;


Transfinite Curve {-21,-23} = N_out_h Using Progression p_out_h;
Transfinite Curve {22, 11} = N_out_v_inflow Using Progression p_out_v_inflow;
Transfinite Curve {-13, -5} = N_out_v_outflow Using Progression p_out_v_outflow;

Recombine Surface {1, 2};

Transfinite Curve {8} = N_triag_v_top Using Progression 1;
Transfinite Curve {1} = N_triag_v_inflow Using Progression p_triag_v_inflow;
Transfinite Curve {-20} = N_triag_v_outflow Using Progression p_triag_v_outflow;
Transfinite Curve {9} = N_triag_h_inflow Using Progression p_triag_h_inflow;
Transfinite Curve {-6} = N_triag_h_outflow Using Progression p_triag_h_outflow;

// defining boundary
Physical Curve(1) = {9}; // inlet
Physical Curve(2) = {6}; // outlet
Physical Curve(3) = {8}; // top
Physical Curve(4) = {1, 22, 2, 15, 3, 16, 4, 17, 5, 20}; // wall


Physical Surface(100) = {1};
Physical Surface(101) = {2};
Physical Surface(102) = {3};
