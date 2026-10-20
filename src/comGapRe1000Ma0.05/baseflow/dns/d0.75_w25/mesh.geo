// Keep coincident points/curves of the non-conforming interfaces as distinct entities
Geometry.AutoCoherence = 0;

// Macro parameters
deltaStar = 1;
D = 0.75 * deltaStar; // depth of the gap
W = 25 * deltaStar; // width of the gap
BL = 0.25 * deltaStar; // height of the first layer of quad elements (conforming with the gap)
x0 = 100 * deltaStar; // x distance from inflow to gap
lengthOutflow = 1000 * deltaStar; // length after the gap
x3 = W + lengthOutflow; // last point of the domain
height = 150 * deltaStar; // height of the domain

// Micro parameters
dx = 0.15; // horizontal size of the elements at the leading edge of the gap (x = 0)
dx_te = dx; // horizontal size of the elements at the trailing edge of the gap (x = W, most resolution-demanding point in the domain)
dy = 1.25 * dx;

// Horizontal sizes. The first layer of quads (layer L = 0) has size dx (dx_te) at the leading (trailing) edge of the gap, growing with
// ratio q_h[0] away from them up to h_max[0]. Layer L and L + 1 (L = k + 1 is the strip k) are separated by a non-conforming interface
// (except inside the gap, where the first layer is conforming with it).
// NESTED interfaces: the nodes of strip k are a subset of the nodes of the layer below. Two neighbouring faces of the layer below are
// merged into one face (2:1) only if the merged face is not larger than the target size of strip k,
// h_target_k = Min(dx * coarsen_h^(k + 1) + (q_h[k + 1] - 1) * d, h_max[k + 1]) (d = distance to the closest edge of the gap);
// otherwise the face is kept (1:1). Hence every vertex of the upper side of an interface is also a vertex of the lower side, and no
// face of the lower side straddles a vertex of the upper side (a staggered, non-nested interface makes P4 unstable).
y_if[] = {BL, 0.75 * deltaStar, 1.75 * deltaStar, 3.75 * deltaStar, 8 * deltaStar, 16 * deltaStar, 35 * deltaStar}; // heights of the non-conforming interfaces (first one must be BL)
coarsen_h = 2; // horizontal coarsening factor of h_min across each interface
q_h[] = {1.1, 1.15, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4}; // horizontal growth ratio away from the edges of the gap, for each layer
h_max[] = {8, 10, 20, 50, 100, 200, 300, 500}; // maximum horizontal size, for each layer
q_v[] = {1.3, 1.3, 1.3, 1.4, 1.6, 2.0, 2.5}; // vertical growth ratio of the elements in each strip
n_v_min[] = {1, 1, 1, 1, 1, 1, 1, 1};
k_fine_end = 4; // from strip k_fine_end on, the 2:1 merging in each interval starts at the end closest to the gap (otherwise the last
                // face before the leading edge is never merged and becomes a spike of very large aspect ratio in the upper strips)

// Every interface is split in x into segments of length ~seg_len (plus the breakpoints 0, W/2, W). Each segment is a separate
// Nektar++ interface: Nektar++ pairs every rank holding one side of an interface with every rank holding the other side,
// so short segments keep the interface communication local.
seg_len = 50;

// fit (check the file scripts/fit_gmsh.py for more details)
p_in_v_left = 1.0335663036230878 - 1.0234149444433507/(W - 0.8098798097181781)^0.5991729377115732; // densitiy concentration of elements inside the gap, vertically, near the left wall.
p_in_v_right = 1.008073572191136 - 1.1471875292840257/(W - 1.2941255994826328)^0.7890392764815325; // densitiy concentration of elements inside the gap, vertically, near the right wall.
p_in_h_bottom = 0.8; // densitiy concentration of elements inside the gap, horizontally.
p_in_h_top = 0.8; // densitiy concentration of elements inside the gap, horizontally.
p_out_h = 0.5; // densitiy concentration of elements of the first layer, vertically (towards the wall).

// Automated parameters
a_in_v_right = dx_te / (W / 2);
a_in_v_left = dx / (W / 2);
a_in_h_bottom = dy / (D / 2);
a_in_h_top = dx / (D / 2);
a_out_h = dx / BL;

// for the derivation of the formula, see file 20250123_gmshFormula.md in the docs folder
N_in_v_right = Ceil( Log (a_in_v_right / (1 - p_in_v_right + a_in_v_right * p_in_v_right)) / Log(p_in_v_right) );
N_in_v_left = Ceil( Log (a_in_v_left / (1 - p_in_v_left + a_in_v_left * p_in_v_left)) / Log(p_in_v_left) );
N_in_h_bottom = Ceil( Log (a_in_h_bottom / (1 - p_in_h_bottom + a_in_h_bottom * p_in_h_bottom)) / Log(p_in_h_bottom) );
N_in_h_top = Ceil( Log (a_in_h_top / (1 - p_in_h_top + a_in_h_top * p_in_h_top)) / Log(p_in_h_top) );
N_out_h = Ceil( Log (a_out_h / (1 - p_out_h + a_out_h * p_out_h)) / Log(p_out_h) ) + 1; // + 1 is to correct the fact that the distance is too small to properly approximate the dx length
N_out_h = Max(N_out_h, n_v_min[0] + 1);

// ---------------------------
// gap
// ---------------------------
Point(2) = {0,0,0};
Point(3) = {0,-D,0};
Point(4) = {W,-D,0};
Point(5) = {W,0,0};
Point(13) = {W/2,0,0};
Point(15) = {0,-D/2,0};
Point(16) = {W/2,-D,0};
Point(17) = {W,-D/2,0};
Line(2) = {2,15};
Line(3) = {3,16};
Line(4) = {4,17};
Line(14) = {2,13};
Line(15) = {15,3};
Line(16) = {16,4};
Line(17) = {17,5};
Line(18) = {13,5};
Curve Loop(2) = {2, 15, 3, 16, 4, 17, -18, -14};
Plane Surface(2) = {2};
Transfinite Surface {2} = {3, 4, 5, 2};
Transfinite Curve {18, 16} = N_in_v_right Using Progression p_in_v_right;
Transfinite Curve {-14, -3} = N_in_v_left Using Progression p_in_v_left;
Transfinite Curve {15, -4} = N_in_h_bottom Using Progression p_in_h_bottom;
Transfinite Curve {-2, 17} = N_in_h_top Using Progression p_in_h_top;
Recombine Surface {2};

// ---------------------------
// breakpoints of the segments of the interfaces (common to all the horizontal lines)
// ---------------------------
x_int[] = {-x0};
n_seg = Max(Round((x0) / seg_len), 1);
For m In {1:n_seg-1}
  x_int[] += {-x0 + m * x0 / n_seg};
EndFor
x_int[] += {0, W/2, W};
n_seg = Max(Round((x3 - W) / seg_len), 1);
For m In {1:n_seg-1}
  x_int[] += {W + m * (x3 - W) / n_seg};
EndFor
x_int[] += {x3};
nInt = #x_int[] - 1;

// Macro creating a horizontal chain of lines at height y_chain, from -x0 to x3, with points at all the breakpoints x_int[]
// and size h_min_k at the edges of the gap, growing with ratio q_k up to h_max_k.
// If gap_lines == 1, the intervals [0, W/2] and [W/2, W] use the lines gap_l and gap_r instead (conforming with the gap).
// Output: chain[] (lines), chain_cnt[] (number of lines in each interval), p_first, p_last (end points)
// Macro appending to Xc[] the n_ nodes (excluding the first one) of a transfinite line [xa_, xa_ + L_] with n_ elements and progression r_
Macro AddNodes
  For j In {1:n_}
    If (Fabs(r_ - 1) < 1e-12)
      Xc[] += {xa_ + L_ * j / n_};
    Else
      Xc[] += {xa_ + L_ * (r_^j - 1) / (r_^n_ - 1)};
    EndIf
  EndFor
Return

// Macro setting isb_ = 1 if xq_ is one of the breakpoints x_int[]
Macro IsBreakpoint
  isb_ = 0;
  For ib_ In {0:#x_int[]-1}
    If (Fabs(xq_ - x_int[ib_]) < 1e-9)
      isb_ = 1;
    EndIf
  EndFor
Return

// Macro creating the nested chain of strip k from the nodes X_lo[] of the layer below, interval by interval (between breakpoints):
// two neighbouring faces are merged (2:1) if the merged face is not larger than the target size of the strip; otherwise the face
// is kept (1:1). The faces are paired from the left end of each interval, except from strip k_fine_end on in the intervals whose
// fine end (closest edge of the gap) is on the right (x < 0 and W/2 < x < W), which are paired from the right. Output: X_up[]
Macro NestedNodes
  X_up[] = {X_lo[0]};
  j0_ = 0;
  For ii_ In {0:nInt-1}
    // nodes j0_..j1_ of X_lo[] in the interval [x_int[ii_], x_int[ii_ + 1]]
    found_ = 0;
    For jj_ In {j0_ + 1:#X_lo[]-1}
      If (found_ == 0 && Fabs(X_lo[jj_] - x_int[ii_ + 1]) < 1e-9)
        j1_ = jj_; found_ = 1;
      EndIf
    EndFor
    xmi_ = (x_int[ii_] + x_int[ii_ + 1]) / 2;
    rev_ = (k >= k_fine_end && (xmi_ < 0 || (xmi_ > W/2 && xmi_ < W)));
    seg_[] = {};
    If (rev_)
      For jj_ In {j1_:j0_:-1}
        seg_[] += {X_lo[jj_]};
      EndFor
    Else
      For jj_ In {j0_:j1_}
        seg_[] += {X_lo[jj_]};
      EndFor
    EndIf
    // greedy pairing along seg_[] (no breakpoint inside it)
    out_[] = {seg_[0]};
    skip = 0;
    For j In {0:#seg_[]-2}
      If (skip)
        skip = 0;
      Else
        merge = 0;
        If (j + 2 <= #seg_[] - 1)
          xm_ = (seg_[j] + seg_[j + 2]) / 2;
          d_ = (xm_ < W/2 ? Fabs(xm_) : Fabs(xm_ - W));
          hm_ = (xm_ > W/2 ? dx_te : dx) * coarsen_h^(k + 1);
          ht_ = Min(hm_ + (q_h[k + 1] - 1) * d_, h_max[k + 1]);
          If (Fabs(seg_[j + 2] - seg_[j]) <= ht_)
            merge = 1;
          EndIf
        EndIf
        If (merge)
          out_[] += {seg_[j + 2]}; skip = 1;
        Else
          out_[] += {seg_[j + 1]};
        EndIf
      EndIf
    EndFor
    // append the new nodes of the interval (without its first node, already in X_up[]) in increasing x
    If (rev_)
      For jj_ In {#out_[]-2:0:-1}
        X_up[] += {out_[jj_]};
      EndFor
    Else
      For jj_ In {1:#out_[]-1}
        X_up[] += {out_[jj_]};
      EndFor
    EndIf
    j0_ = j1_;
  EndFor
Return

// Macro creating a chain of straight lines (one element each) through the nodes X_up[] at height y_chain.
// Output: chain[] (lines), chain_cnt[] (number of lines in each interval between breakpoints), p_first, p_last
Macro ExplicitChain
  chain[] = {}; chain_cnt[] = {};
  pa = newp; Point(pa) = {X_up[0], y_chain, 0}; p_first = pa;
  cnt = 0;
  For j In {1:#X_up[]-1}
    pb = newp; Point(pb) = {X_up[j], y_chain, 0};
    l = newl; Line(l) = {pa, pb};
    Transfinite Curve {l} = 2;
    chain[] += {l}; cnt += 1;
    xq_ = X_up[j]; Call IsBreakpoint;
    If (isb_)
      chain_cnt[] += {cnt}; cnt = 0;
    EndIf
    pa = pb;
  EndFor
  p_last = pa;
Return

Macro HorizontalChain
  chain[] = {};
  chain_cnt[] = {};
  Xc[] = {x_int[0]}; // x of the nodes of the chain
  pa = newp; Point(pa) = {x_int[0], y_chain, 0};
  p_first = pa;
  For i In {0:nInt-1}
    xm = (x_int[i] + x_int[i + 1]) / 2;
    h_min_i = (xm > W/2 ? Min(h_min_k * dx_te / dx, h_max_k) : h_min_k); // size at the closest edge of the gap
    d_c = (h_max_k - h_min_i) / (q_k - 1); // distance to the gap edge where the size reaches h_max_k
    If (gap_lines && xm > 0 && xm < W)
      If (xm < W/2)
        chain[] += {gap_l}; pb = gap_pm;
        For j In {1:N_in_v_left-1}
          m = N_in_v_left - 1 - j;
          Xc[] += {W/2 - W/2 * (p_in_v_left^m - 1) / (p_in_v_left^(N_in_v_left - 1) - 1)};
        EndFor
      Else
        chain[] += {gap_r}; pb = gap_pr;
        For j In {1:N_in_v_right-1}
          Xc[] += {W/2 + W/2 * (p_in_v_right^j - 1) / (p_in_v_right^(N_in_v_right - 1) - 1)};
        EndFor
      EndIf
      chain_cnt[] += {1};
    Else
      If (gap_lines && x_int[i + 1] == 0)
        pb = gap_pl;
      Else
        pb = newp; Point(pb) = {x_int[i + 1], y_chain, 0};
      EndIf
      // distance of both ends to the closest edge of the gap (the breakpoints include 0, W/2 and W)
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
      h1 = Min(h_min_i + (q_k - 1) * d1, h_max_k);
      If (d2 <= d_c + 0.5 * h_max_k || d1 >= d_c - 0.5 * h1)
        // a single line: geometric (fine end at the gap side) or uniform
        If (d1 >= d_c - 0.5 * h1)
          n = Max( Round((d2 - d1) / h_max_k), 1 ); prog = 1;
        Else
          h2 = Min(h_min_i + (q_k - 1) * d2, h_max_k);
          n = Max( Round(Log(h2 / h1) / Log(q_k) + (d2 - Min(d2, d_c)) / h_max_k), 1 ); prog = q_k;
        EndIf
        l = newl; Line(l) = {pa, pb};
        Transfinite Curve {l} = n + 1 Using Progression (da < db ? prog : 1/prog);
        xa_ = x_int[i]; L_ = x_int[i + 1] - x_int[i]; n_ = n; r_ = (da < db ? prog : 1/prog); Call AddNodes;
        chain[] += {l};
        chain_cnt[] += {1};
      Else
        // geometric part next to the gap side, uniform part away from it
        n_g = Max( Round(Log(h_max_k / h1) / Log(q_k)), 1 );
        n_u = Max( Round((d2 - d_c) / h_max_k), 1 );
        pm = newp;
        If (da < db)
          Point(pm) = {x_int[i] + (d_c - d1), y_chain, 0};
          lg = newl; Line(lg) = {pa, pm};
          lu = newl; Line(lu) = {pm, pb};
          Transfinite Curve {lg} = n_g + 1 Using Progression q_k;
          Transfinite Curve {lu} = n_u + 1 Using Progression 1;
          chain[] += {lg, lu};
          xa_ = x_int[i]; L_ = d_c - d1; n_ = n_g; r_ = q_k; Call AddNodes;
          xa_ = x_int[i] + (d_c - d1); L_ = x_int[i + 1] - xa_; n_ = n_u; r_ = 1; Call AddNodes;
        Else
          Point(pm) = {x_int[i + 1] - (d_c - d1), y_chain, 0};
          lu = newl; Line(lu) = {pa, pm};
          lg = newl; Line(lg) = {pm, pb};
          Transfinite Curve {lu} = n_u + 1 Using Progression 1;
          Transfinite Curve {lg} = n_g + 1 Using Progression 1/q_k;
          chain[] += {lu, lg};
          xa_ = x_int[i]; L_ = x_int[i + 1] - (d_c - d1) - xa_; n_ = n_u; r_ = 1; Call AddNodes;
          xa_ = x_int[i + 1] - (d_c - d1); L_ = d_c - d1; n_ = n_g; r_ = 1/q_k; Call AddNodes;
        EndIf
        chain_cnt[] += {2};
      EndIf
    EndIf
    pa = pb;
  EndFor
  p_last = pa;
Return

// ---------------------------
// first layer of quads (y in [0, BL]), conforming with the gap
// ---------------------------
Point(11) = {0,BL,0};
Point(14) = {W/2,BL,0};
Point(12) = {W,BL,0};
Line(12) = {11,14};
Line(19) = {14,12};
Transfinite Curve {19} = N_in_v_right Using Progression p_in_v_right;
Transfinite Curve {-12} = N_in_v_left Using Progression p_in_v_left;

h_min_k = dx; h_max_k = h_max[0]; q_k = q_h[0];
gap_lines = 1;
y_chain = 0; gap_l = 14; gap_r = 18; gap_pl = 2; gap_pm = 13; gap_pr = 5;
Call HorizontalChain;
wall0[] = {chain[]}; pw0 = p_first; pw1 = p_last;
y_chain = BL; gap_l = 12; gap_r = 19; gap_pl = 11; gap_pm = 14; gap_pr = 12;
Call HorizontalChain;
top0[] = {chain[]}; top0_cnt[] = {chain_cnt[]}; pt0 = p_first; pt1 = p_last;
X_lo[] = {Xc[]};
gap_lines = 0;

l_in0 = newl; Line(l_in0) = {pt0, pw0};
l_out0 = newl; Line(l_out0) = {pw1, pt1};
Transfinite Curve {l_in0} = N_out_h Using Progression p_out_h;
Transfinite Curve {-l_out0} = N_out_h Using Progression p_out_h;
cl = newcl; Curve Loop(cl) = {wall0[], l_out0, -top0[{#top0[]-1:0:-1}], l_in0};
s0 = news; Plane Surface(s0) = {cl};
Transfinite Surface {s0} = {pw0, pw1, pt1, pt0};
Recombine Surface {s0};
wall[] = {};
For i In {0:#wall0[]-1}
  If (wall0[i] != 14 && wall0[i] != 18)
    wall[] += {wall0[i]};
  EndIf
EndFor

// ---------------------------
// non-conforming quad strips
// ---------------------------
nStrips = #y_if[];
y_if[] += {height};
h_v = BL * (1 - p_out_h) / (1 - p_out_h^(N_out_h - 1)) * p_out_h^0; // largest element of the first layer (top)
inlet_strips[] = {};
outlet_strips[] = {};
lower_if[] = {top0[]};
lower_cnt[] = {top0_cnt[]};
For k In {0:nStrips-1}
  // nodes of the strip: nested in the nodes of the layer below (1:1 or 2:1)
  Call NestedNodes;
  // horizontal lines (new points also at the bottom, so that the interface is non-conforming, but nested)
  y_chain = y_if[k];
  Call ExplicitChain;
  bottom[] = {chain[]}; bottom_cnt[] = {chain_cnt[]}; pb0 = p_first; pb1 = p_last;
  y_chain = y_if[k + 1];
  Call ExplicitChain;
  top[] = {chain[]}; top_cnt[] = {chain_cnt[]}; pt0 = p_first; pt1 = p_last;
  X_lo[] = {X_up[]};

  // vertical lines, geometric growth q_v starting from the size of the last element below
  H = y_if[k + 1] - y_if[k];
  E_v = Max( Round( Log(1 + H * (q_v[k] - 1) / (q_v[k] * h_v)) / Log(q_v[k]) ), n_v_min[k + 1] );
  h_v = H * (q_v[k] - 1) / (q_v[k]^E_v - 1) * q_v[k]^(E_v - 1); // actual size of the last element of the strip
  ll = newl; Line(ll) = {pb0, pt0};
  lr = newl; Line(lr) = {pb1, pt1};
  Transfinite Curve {ll, lr} = E_v + 1 Using Progression q_v[k];
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

Printf("nStrips = %g, nSegments = %g, nInterfaces = %g (nested)", nStrips, nInt, nStrips * nInt);


// defining boundary
Physical Curve(1) = {l_in0, inlet_strips[]}; // inlet
Physical Curve(2) = {l_out0, outlet_strips[]}; // outlet
Physical Curve(3) = {top_strips[]}; // top
Physical Curve(4) = {wall[], 2, 15, 3, 16, 4, 17}; // wall


Physical Surface(100) = {s0};
Physical Surface(101) = {2};
// Physical Surface(102), (103), ... are the non-conforming strips (defined in the loop above)
// Physical Curve(200), (201), ... are the two sides of the segments of the non-conforming interfaces (defined in the loop above)
