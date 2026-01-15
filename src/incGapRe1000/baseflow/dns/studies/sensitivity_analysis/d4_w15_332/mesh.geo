// Macro parameters
deltaStar = 1;
D = 4 * deltaStar; // depth of the gap 
W = 15 * deltaStar; // width of the gap
BLwidth = 2.85 * deltaStar; // boundary layer thickness, fixed value
BL = 3 * D; // height of the second layer of quad elements
x0 = 5 * W; // x distance from inflow to gap
lengthOutflow = 20 * W; // length after the gap
x1 = W + lengthOutflow; // last point of the domain
triagHeightRegion = 4.5 * BL; // height of the triangular region
h = BL + triagHeightRegion; // height of the domain


// Micro parameters
dx = 0.05; // dx of square elements near the leading and trailing edge of the gap (most resolution-demanding points in the domain)

// All the value below are between 0 and 1. =1  means the elements are all equal in size. Avoid using values very close to 1, because we dividing by log(p) in the formula.
p_in_v = 0.95; // densitiy concentration of elements inside the gap, vertically.
p_in_h = 0.95; // densitiy concentration of elements inside the gap, horizontally. 

p_out_h = 0.9; // densitiy concentration of elements outside the gap, horizontally.
p_out_v_inflow = 0.94; // densitiy concentration of elements outside the gap, vertically, near the inflow.
p_out_v_outflow = 0.98; // densitiy concentration of elements outside the gap, vertically, near the outflow.
p_triag_h = 1.2; // densitiy concentration of elements in the triangular region of the domain, horizontally.


// Automated parameters 
a_in_v = dx / (W / 2);
a_in_h = dx / (D / 2);
a_out_h = dx / BL;
a_out_v_inflow = dx / x0;
a_out_v_outflow = dx / lengthOutflow;

// for the derivation of the formula, see file 20250123_gmshFormula.md in the docs folder
N_in_v = Ceil( Log (a_in_v / (1 - p_in_v + a_in_v * p_in_v)) / Log(p_in_v) );
N_in_h = Ceil( Log (a_in_h / (1 - p_in_h + a_in_h * p_in_h)) / Log(p_in_h) );
N_out_h = Ceil( Log (a_out_h / (1 - p_out_h + a_out_h * p_out_h)) / Log(p_out_h) );
N_out_v_inflow = Ceil( Log (a_out_v_inflow / (1 - p_out_v_inflow + a_out_v_inflow * p_out_v_inflow)) / Log(p_out_v_inflow) );
N_out_v_outflow = Ceil( Log (a_out_v_outflow / (1 - p_out_v_outflow + a_out_v_outflow * p_out_v_outflow)) / Log(p_out_v_outflow) );

// for triangular region
dx_triag_h = BL * (1 - p_out_h) / (1 - p_out_h ^ (N_out_h + 1));
a_triag_h = dx_triag_h / triagHeightRegion;
N_triag_h = Ceil( Log (1 - (1 - p_triag_h) / a_triag_h) / Log(p_triag_h) );
dx_triag_v = triagHeightRegion * p_triag_h ^ N_triag_h * (1 - p_triag_h) / (1 - p_triag_h ^ (N_triag_h + 1));
a_triag_v = dx_triag_v / (x0 + W + lengthOutflow);
N_triag_v = Ceil( 1/(a_triag_v) - 1 ); // we assume p_triag_v = 1, otherwise the formula is different


Point(1) = {-x0, 0, 0};
Point(2) = {0,0,0};
Point(3) = {0,-D,0};
Point(4) = {W,-D,0};
Point(5) = {W,0,0};
Point(6) = {x1,0,0};
Point(7) = {x1,BL,0};
Point(8) = {x1,h,0};
Point(9) = {-x0,h,0};
Point(10) = {-x0,BL,0};
Point(11) = {0,BL,0};
Point(12) = {W,BL,0};
Point(13) = {W/2,0,0};
Point(14) = {W/2,BL,0};
Point(15) = {0,-D/2,0};
Point(16) = {W/2,-D,0};
Point(17) = {W,-D/2,0};

Line(1) = {1,2};
Line(2) = {2,15};
Line(3) = {3,16};
Line(4) = {4,17};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,9};
Line(9) = {9,10};
Line(10) = {10,1};
Line(11) = {10,11};
Line(12) = {11,14};
Line(13) = {12,7};
Line(14) = {2,13};
Line(15) = {15,3};
Line(16) = {16,4};
Line(17) = {17,5};
Line(18) = {13,5};
Line(19) = {14,12};

Curve Loop(1) = {2, 15, 3, 16, 4, 17, -18, -14};
Plane Surface(1) = {1};

Curve Loop(2) = {1,14,18,5,6,-13,-19,-12,-11,10};
Plane Surface(2) = {2};

Curve Loop(3) = {11, 12, 19, 13, 7, 8, 9};
Plane Surface(3) = {3};


Transfinite Surface {1} = {3, 4, 5, 2}; // number '1' must match the number '1' of plane surface
Transfinite Surface {2} = {1, 6, 7, 10}; // number '2' must match the number '2' of plane surface

Transfinite Curve {-14, 18, -3, 16, -12, 19} = N_in_v Using Progression p_in_v;
Transfinite Curve {-2, 15, 17, -4} = N_in_h Using Progression p_in_h;


Transfinite Curve {1, 11} = N_out_v_inflow Using Progression p_out_v_inflow;
Transfinite Curve {10, -6} = N_out_h Using Progression p_out_h;
Transfinite Curve {-13, -5} = N_out_v_outflow Using Progression p_out_v_outflow;

Recombine Surface {1, 2};

Transfinite Curve {8} = N_triag_v Using Progression 1;
Transfinite Curve {7, -9} = N_triag_h Using Progression p_triag_h;

// defining boundary
Physical Curve(1) = {9, 10}; // inlet
Physical Curve(2) = {6, 7}; // outlet
Physical Curve(3) = {8}; // top
Physical Curve(4) = {1, 3, 16, 5, 2, 15, 4, 17}; // wall


Physical Surface(100) = {1};
Physical Surface(101) = {2};
Physical Surface(102) = {3};
