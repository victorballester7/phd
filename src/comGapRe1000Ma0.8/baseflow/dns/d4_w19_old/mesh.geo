// Keep coincident points/curves of the non-conforming interfaces as distinct entities
Geometry.AutoCoherence = 0;

// Macro parameters
deltaStar = 1;
D = 4 * deltaStar; // depth of the gap
W = 19 * deltaStar; // width of the gap
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

// Every interface is split in x into segments of length ~seg_len (plus the fixed breakpoints -r, 0, W/2, W, W + r).
// Each segment is a separate Nektar++ interface. Nektar++ pairs every rank holding one side of an interface with
// every rank holding the other side, so short segments keep the interface communication local (big speed-up in parallel).
seg_len = 50;

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
Line(23) = {19,11};

// ---------------------------
// top of the first layer of quads, split into segments
// ---------------------------
// Macro splitting the line between the existing points pl (left) and pr (right) at height y_split into ~seg_len long
// sub-lines. The nodes are exactly those of a single line with N_split nodes and progression p_split starting at x_start
// (one of the two ends), so the first layer of quads is the same as without splitting.
// Output: split[] (sub-lines, from left to right, all oriented to the right), split_x[] (x of the inner split nodes)
Macro SplitLine
  x_end = x_start == x_left ? x_right : x_left;
  L_split = x_end - x_start;
  n_seg = Max( Round(Fabs(L_split) / seg_len), 1 );
  // node indices (counted from x_start) at the split points
  j_split[] = {0};
  For m In {1:n_seg-1}
    t = m / n_seg;
    j = Round( Log(1 + t * (p_split^(N_split - 1) - 1)) / Log(p_split) );
    If (j > j_split[#j_split[]-1] && j < N_split - 1)
      j_split[] += {j};
    EndIf
  EndFor
  j_split[] += {N_split - 1};
  // from left to right
  If (x_start != x_left)
    j_split[] = {j_split[{#j_split[]-1:0:-1}]};
  EndIf
  split[] = {};
  split_x[] = {};
  pa = pl;
  For m In {0:#j_split[]-2}
    If (m == #j_split[]-2)
      pb = pr;
    Else
      xs = x_start + L_split * (p_split^j_split[m + 1] - 1) / (p_split^(N_split - 1) - 1);
      pb = newp; Point(pb) = {xs, y_split, 0};
      split_x[] += {xs};
    EndIf
    l = newl; Line(l) = {pa, pb};
    Transfinite Curve {l} = Fabs(j_split[m + 1] - j_split[m]) + 1 Using Progression (x_start == x_left ? p_split : 1/p_split);
    split[] += {l};
    pa = pb;
  EndFor
Return

y_split = BL;
// upstream of the gap (old line 10): progression starts at x = -x0
pl = 10; pr = 19; x_left = -x0; x_right = -r; x_start = -x0; N_split = N_out_v_inflow_sparse; p_split = p_out_v_inflow_sparse;
Call SplitLine;
split_in[] = {split[]}; split_in_x[] = {split_x[]};
// downstream of the gap (old line 22): progression starts at x = x3
pl = 21; pr = 7; x_left = W + r; x_right = x3; x_start = x3; N_split = N_out_v_outflow_sparse; p_split = p_out_v_outflow_sparse;
Call SplitLine;
split_out[] = {split[]}; split_out_x[] = {split_x[]};

// breakpoints of the segments of all interfaces, and the curves on top of the first layer of quads in each segment (one per segment)
x_int[] = {-x0, split_in_x[], -r, 0, W/2, W, W + r, split_out_x[], x3};
lower_if[] = {split_in[], 23, 12, 19, 13, split_out[]};
lower_cnt[] = {};
For i In {0:#lower_if[]-1}
  lower_cnt[] += {1};
EndFor
nInt = #x_int[] - 1;

// outer quads
Curve Loop(1) = {5, 21, 6, -split_out[{#split_out[]-1:0:-1}], -13, -19, -12, -23, -split_in[{#split_in[]-1:0:-1}], 11, 1, 20, 14, 18};
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
Transfinite Curve {1} = N_out_v_inflow_sparse Using Progression p_out_v_inflow_sparse;
Transfinite Curve {20, 23} = N_out_v_inflow_dense Using Progression p_out_v_inflow_dense;
Transfinite Curve {-21} = N_out_v_outflow_sparse Using Progression p_out_v_outflow_sparse;
Transfinite Curve {-5, -13} = N_out_v_outflow_dense Using Progression p_out_v_outflow_dense;

Recombine Surface {1, 2};

// ---------------------------
// non-conforming quad strips
// ---------------------------
// Macro creating a horizontal chain of lines at height y_chain, from -x0 to x3, with points at all the breakpoints x_int[]
// and size h_min_k at the edges of the gap, growing with ratio q_h up to h_max_k.
// Output: chain[] (lines), chain_cnt[] (number of lines in each segment), p_first, p_last (end points)
Macro HorizontalChain
  chain[] = {};
  chain_cnt[] = {};
  pa = newp; Point(pa) = {x_int[0], y_chain, 0};
  p_first = pa;
  d_c = (h_max_k - h_min_k) / (q_h - 1); // distance to the gap edge where the size reaches h_max_k
  For i In {0:nInt-1}
    pb = newp; Point(pb) = {x_int[i + 1], y_chain, 0};
    // distance of both ends to the closest edge of the gap (the breakpoints include 0, W/2 and W)
    xm = (x_int[i] + x_int[i + 1]) / 2;
    If (xm < 0)
      da = -x_int[i]; db = -x_int[i + 1];
    ElseIf (xm < W/2)
      da = x_int[i]; db = x_int[i + 1];
    ElseIf (xm < W)
      da = W - x_int[i]; db = W - x_int[i + 1];
    Else
      da = x_int[i] - W; db = x_int[i + 1] - W;
    EndIf
    d1 = Min(da, db); d2 = Max(da, db);
    h1 = Min(h_min_k + (q_h - 1) * d1, h_max_k);
    If (d2 <= d_c + 0.5 * h_max_k || d1 >= d_c - 0.5 * h1)
      // a single line: geometric (fine end at the gap side) or uniform
      If (d1 >= d_c - 0.5 * h1)
        n = Max( Round((d2 - d1) / h_max_k), 1 ); prog = 1;
      Else
        h2 = Min(h_min_k + (q_h - 1) * d2, h_max_k);
        n = Max( Round(Log(h2 / h1) / Log(q_h) + (d2 - Min(d2, d_c)) / h_max_k), 1 ); prog = q_h;
      EndIf
      l = newl; Line(l) = {pa, pb};
      Transfinite Curve {l} = n + 1 Using Progression (da < db ? prog : 1/prog);
      chain[] += {l};
      chain_cnt[] += {1};
    Else
      // geometric part next to the gap side, uniform part away from it
      n_g = Max( Round(Log(h_max_k / h1) / Log(q_h)), 1 );
      n_u = Max( Round((d2 - d_c) / h_max_k), 1 );
      pm = newp;
      If (da < db)
        Point(pm) = {x_int[i] + (d_c - d1), y_chain, 0};
        lg = newl; Line(lg) = {pa, pm};
        lu = newl; Line(lu) = {pm, pb};
        Transfinite Curve {lg} = n_g + 1 Using Progression q_h;
        Transfinite Curve {lu} = n_u + 1 Using Progression 1;
        chain[] += {lg, lu};
      Else
        Point(pm) = {x_int[i + 1] - (d_c - d1), y_chain, 0};
        lu = newl; Line(lu) = {pa, pm};
        lg = newl; Line(lg) = {pm, pb};
        Transfinite Curve {lu} = n_u + 1 Using Progression 1;
        Transfinite Curve {lg} = n_g + 1 Using Progression 1/q_h;
        chain[] += {lu, lg};
      EndIf
      chain_cnt[] += {2};
    EndIf
    pa = pb;
  EndFor
  p_last = pa;
Return

nStrips = #y_if[];
y_if[] += {height};
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
  bottom[] = {chain[]}; bottom_cnt[] = {chain_cnt[]}; pb0 = p_first; pb1 = p_last;
  y_chain = y_if[k + 1];
  Call HorizontalChain;
  top[] = {chain[]}; top_cnt[] = {chain_cnt[]}; pt0 = p_first; pt1 = p_last;

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

  // interface k, segment i: lower side 200 + 2*(k*nInt + i), upper side 201 + 2*(k*nInt + i)
  cl_lo = 0; cl_up = 0;
  For i In {0:nInt-1}
    Physical Curve(200 + 2*(k*nInt + i)) = {lower_if[{cl_lo:cl_lo + lower_cnt[i] - 1}]};
    Physical Curve(201 + 2*(k*nInt + i)) = {bottom[{cl_up:cl_up + bottom_cnt[i] - 1}]};
    cl_lo += lower_cnt[i];
    cl_up += bottom_cnt[i];
  EndFor
  lower_if[] = {top[]};
  lower_cnt[] = {top_cnt[]};
EndFor
top_strips[] = {lower_if[]};

Printf("nStrips = %g, nSegments = %g, nInterfaces = %g", nStrips, nInt, nStrips * nInt);


// defining boundary
Physical Curve(1) = {11, inlet_strips[]}; // inlet
Physical Curve(2) = {6, outlet_strips[]}; // outlet
Physical Curve(3) = {top_strips[]}; // top
Physical Curve(4) = {1, 20, 2, 15, 3, 16, 4, 17, 5, 21}; // wall


Physical Surface(100) = {1};
Physical Surface(101) = {2};
// Physical Surface(102), (103), ... are the non-conforming strips (defined in the loop above)
// Physical Curve(200), (201), ... are the two sides of the segments of the non-conforming interfaces (defined in the loop above)
