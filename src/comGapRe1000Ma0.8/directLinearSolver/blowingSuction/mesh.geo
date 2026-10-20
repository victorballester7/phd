// Keep coincident points/curves of the non-conforming interfaces as distinct entities
Geometry.AutoCoherence = 0;

// Macro parameters
deltaStar = 1;
D = 0.5 * deltaStar; // depth of the gap 
W = 30 * deltaStar; // width of the gap
r = 16; // length of non-constant quads upstream and downstream of the gap
s = 8; // length of non-constant quads left and right inside the gap
BL = 0.25 * deltaStar; // height of the second layer of quad elements
x0 = 100 * deltaStar; // x distance from inflow to gap
lengthOutflow = 1000 * deltaStar; // length after the gap
x3 = W + lengthOutflow; // last point of the domain
height = 150 * deltaStar; // height of the domain

// Micro parameters
dx = 0.15; // dx of square elements near the leading and trailing edge of the gap (most resolution-demanding points in the domain)
dy = 1.25 * dx; 

// Non-conforming quad strips above the first layer of quads (y > BL).
// Each interface is a horizontal line where the mesh below and above do not share nodes.
// In strip k, the horizontal size is h_min = dx * coarsen_h^(k+1) at the two edges of the gap (x = 0 and x = W),
// grows geometrically (ratio q_h) away from them and is capped at h_far * coarsen_h^(k+1), and everything is capped at h_top.
// Hence, once dx * coarsen_h^(k+1) >= h_top, the strip is uniform in x.
y_if[] = {BL, 0.75 * deltaStar, 1.75 * deltaStar, 3.75 * deltaStar, 8 * deltaStar, 16 * deltaStar, 35 * deltaStar}; // heights of the non-conforming interfaces (first one must be BL)
coarsen_h = 2; // horizontal coarsening factor across each interface
q_h = 1.15; // horizontal growth ratio of the elements away from the edges of the gap
h_far = 4; // horizontal size far from the gap in the first strip (multiplied by coarsen_h in each strip)
h_top = 30; // horizontal size of the elements in the uniform upper strips (maximum horizontal size)
q_v = 1.3; // vertical growth ratio of the elements in the strips (continuous across interfaces)

// All the value below are between 0 and 1. =1  means the elements are all equal in size. Avoid using values very close to 1, because we dividing by log(p) in the formula.
dy_in_v_right = 1.5; // heighest quad in the upper layer of quads, in the right part of the gap
dy_in_v_left = 1.75; // heighest quad in the upper layer of quads, in the left part of the gap

// fit (check the file scripts/fit_gmsh.py for more details)
p_in_v_left = 1.0335663036230878 - 1.0234149444433507/(W - 0.8098798097181781)^0.5991729377115732; // densitiy concentration of elements inside the gap, vertically, near the left wall.
p_in_v_right = 1.008073572191136 - 1.1471875292840257/(W - 1.2941255994826328)^0.7890392764815325; // densitiy concentration of elements inside the gap, vertically, near the right wall.
p_in_h_bottom = 0.8; // densitiy concentration of elements inside the gap, horizontally. 
p_in_h_top = 0.8; // densitiy concentration of elements inside the gap, horizontally. 

p_out_h = 0.5; // densitiy concentration of elements outside the gap, horizontally.
p_out_v_inflow_dense = 0.8; // densitiy concentration of elements outside the gap, vertically, near the inflow and the gap
p_out_v_inflow_sparse = 0.99; // densitiy concentration of elements outside the gap, vertically, near the inflow and far from the gap
p_out_v_outflow_dense = 0.9; // densitiy concentration of elements outside the gap, vertically, near the outflow and the gap
p_out_v_outflow_sparse = 0.995; // densitiy concentration of elements outside the gap, vertically, near the outflow and far from the gap

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
N_out_v_outflow_dense = Ceil( Log (a_out_v_outflow_dense / (1 - p_out_v_outflow_dense + a_out_v_outflow_dense * p_out_v_outflow_dense)) / Log(p_out_v_outflow_dense) );

p_for_a_out_v_inflow_sparse_computation = p_out_v_inflow_dense;
p_for_a_out_v_outflow_sparse_computation = (p_out_v_outflow_dense + p_out_v_outflow_dense) / 2;
a_out_v_inflow_sparse = r / (x0 - r) * (1 - p_for_a_out_v_inflow_sparse_computation) / (1 - p_for_a_out_v_inflow_sparse_computation^(N_out_v_inflow_dense + 1));
a_out_v_outflow_sparse = r / (lengthOutflow - r) * (1 - p_for_a_out_v_outflow_sparse_computation) / (1 - p_for_a_out_v_outflow_sparse_computation^(N_out_v_outflow_dense + 1));

N_out_v_inflow_sparse = Ceil( Log(a_out_v_inflow_sparse / (1 - p_out_v_inflow_sparse + a_out_v_inflow_sparse * p_out_v_inflow_sparse)) / Log(p_out_v_inflow_sparse) );
N_out_v_outflow_sparse = Ceil( Log(a_out_v_outflow_sparse / (1 - p_out_v_outflow_sparse + a_out_v_outflow_sparse * p_out_v_outflow_sparse)) / Log(p_out_v_outflow_sparse) );


Point(1) = {-x0, 0, 0};
Point(2) = {0,0,0};
Point(3) = {0,-D,0};
Point(4) = {W,-D,0};
Point(5) = {W,0,0};
Point(6) = {x3,0,0};
Point(7) = {x3,BL,0};
Point(10) = {-x0,BL,0};
Point(11) = {0,BL,0};
Point(12) = {W,BL,0};
Point(13) = {W/2,0,0};
Point(14) = {W/2,BL,0};
Point(15) = {0,-D/2,0};
Point(16) = {W/2,-D,0};
Point(17) = {W,-D/2,0};
Point(18) = {-r,0,0};
Point(19) = {-r,BL,0};
Point(20) = {W + r,0,0};
Point(21) = {W + r,BL,0};

Line(1) = {1,18};
Line(2) = {2,15};
Line(3) = {3,16};
Line(4) = {4,17};
Line(5) = {5,20};
Line(6) = {6,7};
Line(10) = {10,19};
Line(11) = {10,1};
Line(12) = {11,14};
Line(13) = {12,21};
Line(14) = {2,13};
Line(15) = {15,3};
Line(16) = {16,4};
Line(17) = {17,5};
Line(18) = {13,5};
Line(19) = {14,12};
Line(20) = {18,2};
Line(21) = {20,6};
Line(22) = {21,7};
Line(23) = {19,11};

// outer quads
Curve Loop(1) = {5, 21, 6, -22, -13, -19, -12, -23, -10, 11, 1, 20, 14, 18};
Plane Surface(1) = {1};

// inside gap
Curve Loop(2) = {2, 15, 3, 16, 4, 17, -18, -14};
Plane Surface(2) = {2};


// This defines the surfaces that will be meshed with quad elements
Transfinite Surface {1} = {1, 6, 7, 10};
Transfinite Surface {2} = {3, 4, 5, 2}; // number '2' must match the number '2' of plane surface

Transfinite Curve {18, 16, 19} = N_in_v_right Using Progression p_in_v_right;
Transfinite Curve {-14, -3, -12} = N_in_v_left Using Progression p_in_v_left;

Transfinite Curve {15, -4} = N_in_h_bottom Using Progression p_in_h_bottom;
Transfinite Curve {-2, 17} = N_in_h_top Using Progression p_in_h_top;


Transfinite Curve {11, -6} = N_out_h Using Progression p_out_h;
Transfinite Curve {10, 1} = N_out_v_inflow_sparse Using Progression p_out_v_inflow_sparse;
Transfinite Curve {20, 23} = N_out_v_inflow_dense Using Progression p_out_v_inflow_dense;
Transfinite Curve {-21, -22} = N_out_v_outflow_sparse Using Progression p_out_v_outflow_sparse;
Transfinite Curve {-5, -13} = N_out_v_outflow_dense Using Progression p_out_v_outflow_dense;

Recombine Surface {1, 2};

// ---------------------------
// non-conforming quad strips
// ---------------------------
// Horizontal intervals, from left to right, and whether the fine end (an edge of the gap) is at their left (1) or right (0)
x_int[] = {-x0, 0, W/2, W, x3};
fine_left[] = {0, 1, 0, 1};

// Macro creating a horizontal chain of lines at height y_chain, from -x0 to x3, with size h_min_k at the edges
// of the gap, growing with ratio q_h up to h_max_k. Output: chain[] (lines), p_first, p_last (end points)
Macro HorizontalChain
  chain[] = {};
  pa = newp; Point(pa) = {x_int[0], y_chain, 0};
  p_first = pa;
  For i In {0:#fine_left[]-1}
    L = x_int[i + 1] - x_int[i];
    pb = newp; Point(pb) = {x_int[i + 1], y_chain, 0};
    // geometric part: n_g elements from h_min_k up to (just below) h_max_k
    n_g = 0;
    If (h_max_k > h_min_k)
      n_g = Ceil( Log(h_max_k / h_min_k) / Log(q_h) );
    EndIf
    L_g = h_min_k * (q_h^n_g - 1) / (q_h - 1);
    n_u = Round((L - L_g) / h_max_k); // elements of the uniform part
    If (n_g == 0)
      // uniform only
      l = newl; Line(l) = {pa, pb};
      Transfinite Curve {l} = Max(Round(L / h_max_k), 1) + 1 Using Progression 1;
      chain[] += {l};
    ElseIf (n_u < 1)
      // geometric only (the interval is too short to fit a uniform part)
      n = Max( Round( Log(1 + L * (q_h - 1) / h_min_k) / Log(q_h) ), 1 );
      l = newl; Line(l) = {pa, pb};
      Transfinite Curve {l} = n + 1 Using Progression q_h^(2 * fine_left[i] - 1);
      chain[] += {l};
    Else
      // geometric part next to the gap edge, uniform part away from it
      If (fine_left[i])
        pm = newp; Point(pm) = {x_int[i] + L_g, y_chain, 0};
        lg = newl; Line(lg) = {pa, pm};
        lu = newl; Line(lu) = {pm, pb};
        Transfinite Curve {lg} = n_g + 1 Using Progression q_h;
        Transfinite Curve {lu} = n_u + 1 Using Progression 1;
        chain[] += {lg, lu};
      Else
        pm = newp; Point(pm) = {x_int[i + 1] - L_g, y_chain, 0};
        lu = newl; Line(lu) = {pa, pm};
        lg = newl; Line(lg) = {pm, pb};
        Transfinite Curve {lu} = n_u + 1 Using Progression 1;
        Transfinite Curve {lg} = n_g + 1 Using Progression 1/q_h;
        chain[] += {lu, lg};
      EndIf
    EndIf
    pa = pb;
  EndFor
  p_last = pa;
Return

nStrips = #y_if[];
y_if[] += {height};
lower_if[] = {10, 23, 12, 19, 13, 22}; // lower side of the first interface (top of the first layer of quads)
h_v = BL; // height of the element below the first strip
inlet_strips[] = {};
outlet_strips[] = {};
For k In {0:nStrips-1}
  m = coarsen_h^(k + 1);
  h_max_k = Min(h_far * m, h_top);
  h_min_k = Min(dx * m, h_max_k);

  // horizontal lines (new points also at the bottom, so that the interface is non-conforming)
  y_chain = y_if[k];
  Call HorizontalChain;
  bottom[] = {chain[]}; pb0 = p_first; pb1 = p_last;
  y_chain = y_if[k + 1];
  Call HorizontalChain;
  top[] = {chain[]}; pt0 = p_first; pt1 = p_last;

  // vertical lines, geometric growth q_v starting from the size of the last element below
  H = y_if[k + 1] - y_if[k];
  E_v = Max( Round( Log(1 + H * (q_v - 1) / (q_v * h_v)) / Log(q_v) ), 1 );
  h_v = H * (q_v - 1) / (q_v^E_v - 1) * q_v^(E_v - 1); // actual size of the last element of the strip
  ll = newl; Line(ll) = {pb0, pt0};
  lr = newl; Line(lr) = {pb1, pt1};
  Transfinite Curve {ll, lr} = E_v + 1 Using Progression q_v;
  inlet_strips[] += {ll};
  outlet_strips[] += {lr};

  // surface
  cl = newcl; Curve Loop(cl) = {bottom[], lr, -top[{#top[]-1:0:-1}], -ll};
  s = news; Plane Surface(s) = {cl};
  Transfinite Surface {s} = {pb0, pb1, pt1, pt0};
  Recombine Surface {s};
  Physical Surface(102 + k) = {s};

  // interface k: lower side 200 + 2k, upper side 201 + 2k
  Physical Curve(200 + 2*k) = {lower_if[]};
  Physical Curve(201 + 2*k) = {bottom[]};
  lower_if[] = {top[]};
EndFor
top_strips[] = {lower_if[]};


// defining boundary
Physical Curve(1) = {11, inlet_strips[]}; // inlet
Physical Curve(2) = {6, outlet_strips[]}; // outlet
Physical Curve(3) = {top_strips[]}; // top
Physical Curve(4) = {1, 20, 2, 15, 3, 16, 4, 17, 5, 21}; // wall


Physical Surface(100) = {1};
Physical Surface(101) = {2};
// Physical Surface(102), (103), ... are the non-conforming strips (defined in the loop above)
