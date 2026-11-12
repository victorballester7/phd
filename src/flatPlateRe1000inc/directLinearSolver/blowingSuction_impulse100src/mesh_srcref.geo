// Macro parameters
deltaStar = 1;
D = 4 * deltaStar; // depth of the gap 
BL = 3.75 * deltaStar; // height of the second layer of quad elements
x0 = 100 * deltaStar; // x distance from inflow to gap
lengthOutflow = 1000 * deltaStar; // length after the gap
x3 = lengthOutflow; // last point of the domain
triagHeightRegion = 150 * deltaStar; // height of the triangular region
h = triagHeightRegion; // height of the domain


// Micro parameters
size_triag_v_top = 32; // size of the elements in the upper triangular region
//size_triag_h_inflow = 4; // size of the elements in the horizontal direction in the triangular region (inflow part)
//size_triag_h_outflow = 8;
size_quad_h = BL / 6;
size_quad_v = 2.25;

// All the value below are between 0 and 1. =1  means the elements are all equal in size. Avoid using values very close to 1, because we dividing by log(p) in the formula.
p_quad_v = 1; // densitiy concentration of quad elements, vertically, near the inflow.
p_quad_h = 0.45;  
p_triag_h_outflow = 0.8; // densitiy concentration of elements in the triangular region of the domain, horizontally.
p_triag_h_inflow = p_triag_h_outflow; // density concentration of elements in the triangular region of the domain, horizontally.


// for the derivation of the formula, see file 20250123_gmshFormula.md in the docs folder
N_quad_h = Ceil( BL / size_quad_h);
N_quad_v = Ceil ((x0 + x3) / size_quad_v);

// ---------------------------
// for triangular region
// ---------------------------
a_triag_h_inflow = (BL) / (h - BL) * (1 - p_quad_h) / (1 - p_quad_h^(N_quad_h + 1));
a_triag_h_outflow = a_triag_h_inflow;
N_triag_v_top = Ceil( (x3 + x0) / size_triag_v_top ); // we assume p_triag_v = 1, otherwise the formula is different
//N_triag_h_inflow = Ceil( h / size_triag_h_inflow );
//N_triag_h_outflow = Ceil( h / size_triag_h_outflow );
N_triag_h_inflow = Ceil( Log(a_triag_h_inflow / (1 - p_triag_h_inflow + a_triag_h_inflow * p_triag_h_inflow)) / Log(p_triag_h_inflow) );
N_triag_h_outflow = Ceil( Log(a_triag_h_outflow / (1 - p_triag_h_outflow + a_triag_h_outflow * p_triag_h_outflow)) / Log(p_triag_h_outflow) );


j_patch = Round(12 * (N_quad_v - 1) / (x0 + x3));
xa = -x0 + j_patch * (x0 + x3) / (N_quad_v - 1);
Printf("source patch: j = %g, xa = %g", j_patch, xa);
Point(1) = {-x0, 0, 0};
Point(101) = {xa, 0, 0};
Point(102) = {xa, BL, 0};
Point(6) = {x3,0,0};
Point(7) = {x3,BL,0};
Point(8) = {x3,h,0};
Point(9) = {-x0,h,0};
Point(10) = {-x0,BL,0};

Line(1) = {101,6};
Line(101) = {1,101};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,9};
Line(9) = {9,10};
Line(11) = {10,1};
Line(12) = {7, 102};
Line(103) = {102, 10};
Line(104) = {102, 101};


// outer quads
Curve Loop(1) = {6, 12, 104, 1};
Plane Surface(1) = {1};

// triangular region
Curve Loop(3) = {9, -103, -12, 7, 8};
Plane Surface(3) = {3};

// source patch
Curve Loop(5) = {101, -104, 103, 11};
Plane Surface(5) = {5};


// This defines the surfaces that will be meshed with quad elements
Transfinite Surface {1} = {101, 6, 7, 102};

Transfinite Curve {11, 104, -6} = N_quad_h Using Progression p_quad_h;
Transfinite Curve {1, -12} = N_quad_v - j_patch Using Progression p_quad_v;
Transfinite Curve {103} = j_patch + 1 Using Progression p_quad_v; // top of the patch keeps the original nodes (triangle region unchanged)

Recombine Surface {1, 2};

// for triangluar region
Transfinite Curve {8} = N_triag_v_top Using Progression 1;
Transfinite Curve {9} = N_triag_h_inflow Using Progression p_triag_h_inflow;
Transfinite Curve {-7} = N_triag_h_outflow Using Progression p_triag_h_outflow;


// defining boundary
Physical Curve(1) = {9,11}; // inlet
Physical Curve(2) = {7,6}; // outlet
Physical Curve(3) = {8}; // top
Physical Curve(4) = {101, 1}; // wall


Physical Surface(100) = {1};
Physical Surface(102) = {3};

// ---- source refinement patch (scripts/refineSourceMesh.py) ----
Point(200) = {-95.0, 0.5, 0};
Point{200} In Surface{5};
Field[1] = Distance; Field[1].PointsList = {200};
Field[2] = Threshold; Field[2].InField = 1; Field[2].SizeMin = 0.05; Field[2].SizeMax = 1.0;
Field[2].DistMin = 0.35; Field[2].DistMax = 4.1499999999999995;
Field[3] = Distance; Field[3].CurvesList = {101}; Field[3].Sampling = 400;
Field[4] = Threshold; Field[4].InField = 3; Field[4].SizeMin = 0.12; Field[4].SizeMax = 1.0;
Field[4].DistMin = 0; Field[4].DistMax = 2.5142857142857147;
Field[5] = Min; Field[5].FieldsList = {2, 4};
Field[6] = Restrict; Field[6].InField = 5; Field[6].SurfacesList = {5};
Background Field = 6;
MeshAlgorithm Surface {5} = 6;
Physical Surface(103) = {5};
