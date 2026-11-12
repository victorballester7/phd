"""Local refinement around the blowing/suction source for the incompressible gap and flat-plate meshes (2026-11-10).

The production source is a Gaussian vortex exp(-w r^2) at (xs, ys) = (-95, 0.5). With w = 100 (the noise runs'
envelope) its e-folding radius is 0.1, while the structured mesh there has a 4.4 x 0.5 element (gaps, P5) or a
2.25 x 0.43 one (flat plate, P6): the vortex is sampled at a single quadrature column. This script rewrites a mesh.geo
so that the upstream end of the boundary-layer quad block, x in [-x0, xa] (xa ~ -x0 + LPATCH) and y in [0, top of the
quad layers], becomes an unstructured triangle patch (Physical Surface 103) graded around the source:

    h = min(HMAX, HMIN + G * max(0, r - R0),  HW0 + GW * y)      r = distance to the source

xa is put on an existing node of the structured streamwise distribution, and the remaining structured curves keep
the same node count and progression, so every node of the structured mesh downstream of xa is unchanged. The patch
centre is a mesh vertex. Sessions must add composite 103 to their EXPANSIONS.

Because gmsh re-triangulates the free-stream triangle region whenever its boundary curves are split (flat plate),
the final mesh is SPLICED (splice()): the original mesh (regenerated from the original geo, msh2) minus its quads
and wall edges in the patch region, plus the patch triangles and patch wall edges of the refined mesh, merged on
the coinciding interface nodes. Everything outside the patch is then the original mesh, node for node.

usage: python3 refineSourceMesh.py gap|flat|bfs <mesh.geo> <out.geo> [--hmin 0.05]
       python3 refineSourceMesh.py splice <orig.msh2> <refined.msh2> <xa> <ytop> <out.msh2>
"""
import argparse
import re

XS, YS = -95.0, 0.5
LPATCH = 12.0  # patch length from the inlet
R0, G, HMAX = 0.35, 0.25, 1.0
HW0, GW = 0.12, 0.35


def sub1(s, old, new):
    assert s.count(old) == 1, old
    return s.replace(old, new)


def subre(s, pat, new):
    s2, n = re.subn(pat, new, s, flags=re.M)
    assert n == 1, pat
    return s2


def fields(hmin, wall_curves, surf=5):
    dmax_src = R0 + (HMAX - hmin) / G
    dmax_wall = (HMAX - HW0) / GW
    return f"""
// ---- source refinement patch (scripts/refineSourceMesh.py) ----
Point(200) = {{{XS}, {YS}, 0}};
Point{{200}} In Surface{{{surf}}};
Field[1] = Distance; Field[1].PointsList = {{200}};
Field[2] = Threshold; Field[2].InField = 1; Field[2].SizeMin = {hmin}; Field[2].SizeMax = {HMAX};
Field[2].DistMin = {R0}; Field[2].DistMax = {dmax_src};
Field[3] = Distance; Field[3].CurvesList = {{{wall_curves}}}; Field[3].Sampling = 400;
Field[4] = Threshold; Field[4].InField = 3; Field[4].SizeMin = {HW0}; Field[4].SizeMax = {HMAX};
Field[4].DistMin = 0; Field[4].DistMax = {dmax_wall};
Field[5] = Min; Field[5].FieldsList = {{2, 4}};
Field[6] = Restrict; Field[6].InField = 5; Field[6].SurfacesList = {{{surf}}};
Background Field = 6;
MeshAlgorithm Surface {{{surf}}} = 6;
Physical Surface(103) = {{{surf}}};
"""


def gap(s, hmin):
    # split position on a node of the inflow-sparse progression (lines 1, 10, 27: x = -x0 .. -r, increasing x)
    s = sub1(s, "Point(1) = {-x0, 0, 0};", """L_sp = x0 - r;
l0_sp = L_sp * (1 - p_out_v_inflow_sparse) / (1 - p_out_v_inflow_sparse^(N_out_v_inflow_sparse - 1));
j_patch = Round(Log(1 - %g * (1 - p_out_v_inflow_sparse) / l0_sp) / Log(p_out_v_inflow_sparse));
xa = -x0 + l0_sp * (1 - p_out_v_inflow_sparse^j_patch) / (1 - p_out_v_inflow_sparse);
Printf("source patch: j = %%g, xa = %%g", j_patch, xa);
Point(1) = {-x0, 0, 0};
Point(101) = {xa, 0, 0};
Point(102) = {xa, BL, 0};
Point(103) = {xa, BL + BL_upper, 0};""" % LPATCH)
    s = sub1(s, "Line(1) = {1,18};", "Line(1) = {101,18};\nLine(101) = {1,101};")
    s = sub1(s, "Line(10) = {10,19};", "Line(10) = {102,19};")
    s = sub1(s, "Line(27) = {22,23};", "Line(27) = {103,23};\nLine(103) = {22,103};\nLine(104) = {102,101};\nLine(105) = {103,102};")
    s = sub1(s, "Curve Loop(1) = {5, 21, 6, -22, -13, -19, -12, -23, -10, 11, 1, 20, 14, 18};",
             "Curve Loop(1) = {5, 21, 6, -22, -13, -19, -12, -23, -10, 104, 1, 20, 14, 18};")
    s = sub1(s, "Curve Loop(3) = {-26, 25, 8, 9, 27, 30, 31,32,33};", "Curve Loop(3) = {-26, 25, 8, 9, 103, 27, 30, 31,32,33};")
    s = sub1(s, "Curve Loop(4) = {-24, 27, 30, 31,32,33, -26,-7,-22,-13,-19,-12,-23,-10};",
             "Curve Loop(4) = {-105, 27, 30, 31,32,33, -26,-7,-22,-13,-19,-12,-23,-10};")
    s = sub1(s, "Plane Surface(4) = {4};", "Plane Surface(4) = {4};\n\n// source patch\nCurve Loop(5) = {101, -104, -105, -103, 24, 11};\nPlane Surface(5) = {5};")
    s = sub1(s, "Transfinite Surface {1} = {1, 6, 7, 10};", "Transfinite Surface {1} = {101, 6, 7, 102};")
    s = sub1(s, "Transfinite Surface {4} = {7, 25, 22, 10};", "Transfinite Surface {4} = {7, 25, 103, 102};")
    s = sub1(s, "Transfinite Curve {11, -6} = N_out_h Using Progression p_out_h;",
             "Transfinite Curve {11, 104, -6} = N_out_h Using Progression p_out_h;")
    s = sub1(s, "Transfinite Curve {24, 28, 29, -7} = N_out_h_upper Using Progression p_out_h_upper;",
             "Transfinite Curve {24, 105, -7} = N_out_h_upper Using Progression p_out_h_upper;")
    s = sub1(s, "Transfinite Curve {10, 1, 27} = N_out_v_inflow_sparse Using Progression p_out_v_inflow_sparse;",
             "Transfinite Curve {10, 1, 27} = N_out_v_inflow_sparse - j_patch Using Progression p_out_v_inflow_sparse;\n"
             "Transfinite Curve {103} = j_patch + 1 Using Progression p_out_v_inflow_sparse; // top of the patch keeps the original nodes (triangle region unchanged)")
    s = sub1(s, "Physical Curve(4) = {1, 20,", "Physical Curve(4) = {101, 1, 20,")
    return s + fields(hmin, "101")


def flat(s, hmin):
    s = sub1(s, "Point(1) = {-x0, 0, 0};", """j_patch = Round(%g * (N_quad_v - 1) / (x0 + x3));
xa = -x0 + j_patch * (x0 + x3) / (N_quad_v - 1);
Printf("source patch: j = %%g, xa = %%g", j_patch, xa);
Point(1) = {-x0, 0, 0};
Point(101) = {xa, 0, 0};
Point(102) = {xa, BL, 0};""" % LPATCH)
    s = sub1(s, "Line(1) = {1,6};", "Line(1) = {101,6};\nLine(101) = {1,101};")
    s = sub1(s, "Line(12) = {7, 10};", "Line(12) = {7, 102};\nLine(103) = {102, 10};\nLine(104) = {102, 101};")
    s = sub1(s, "Curve Loop(1) = {6, 12, 11, 1};", "Curve Loop(1) = {6, 12, 104, 1};")
    s = sub1(s, "Curve Loop(3) = {9, -12, 7, 8};", "Curve Loop(3) = {9, -103, -12, 7, 8};")
    s = sub1(s, "Plane Surface(3) = {3};", "Plane Surface(3) = {3};\n\n// source patch\nCurve Loop(5) = {101, -104, 103, 11};\nPlane Surface(5) = {5};")
    s = sub1(s, "Transfinite Surface {1} = {1, 6, 7, 10};", "Transfinite Surface {1} = {101, 6, 7, 102};")
    s = sub1(s, "Transfinite Curve {11, -6} = N_quad_h Using Progression p_quad_h;",
             "Transfinite Curve {11, 104, -6} = N_quad_h Using Progression p_quad_h;")
    s = sub1(s, "Transfinite Curve {1, -12} = N_quad_v Using Progression p_quad_v;",
             "Transfinite Curve {1, -12} = N_quad_v - j_patch Using Progression p_quad_v;\n"
             "Transfinite Curve {103} = j_patch + 1 Using Progression p_quad_v; // top of the patch keeps the original nodes (triangle region unchanged)")
    s = sub1(s, "Physical Curve(4) = {1}; // wall", "Physical Curve(4) = {101, 1}; // wall")
    return s + fields(hmin, "101")


def bfs(s, hmin):
    """Backward-facing step geo (bfsRe1000inc): same upstream block as the gap geo (lines 1, 10, 11, 24, 27, surfaces 1
    and 4, same orientations), but other curve loops, and Surface(5) is taken, so the patch is Surface(6) (2026-11-12)."""
    s = sub1(s, "Point(1) = {-x0, 0, 0};", """L_sp = x0 - r;
l0_sp = L_sp * (1 - p_out_v_inflow_sparse) / (1 - p_out_v_inflow_sparse^(N_out_v_inflow_sparse - 1));
j_patch = Round(Log(1 - %g * (1 - p_out_v_inflow_sparse) / l0_sp) / Log(p_out_v_inflow_sparse));
xa = -x0 + l0_sp * (1 - p_out_v_inflow_sparse^j_patch) / (1 - p_out_v_inflow_sparse);
Printf("source patch: j = %%g, xa = %%g", j_patch, xa);
Point(1) = {-x0, 0, 0};
Point(101) = {xa, 0, 0};
Point(102) = {xa, BL, 0};
Point(103) = {xa, BL + BL_upper, 0};""" % LPATCH)
    s = sub1(s, "Line(1) = {1,18};", "Line(1) = {101,18};\nLine(101) = {1,101};")
    s = sub1(s, "Line(10) = {10,19};", "Line(10) = {102,19};")
    s = sub1(s, "Line(27) = {22,23};", "Line(27) = {103,23};\nLine(103) = {22,103};\nLine(104) = {102,101};\nLine(105) = {103,102};")
    s = sub1(s, "Curve Loop(1) = {-21, -12, -23, -10, 11, 1, 20, 14};", "Curve Loop(1) = {-21, -12, -23, -10, 104, 1, 20, 14};")
    s = sub1(s, "Curve Loop(3) = {25, 8, 9, 27, 30, 31,32,21,22,33,16,18,28};",
             "Curve Loop(3) = {25, 8, 9, 103, 27, 30, 31,32,21,22,33,16,18,28};")
    s = sub1(s, "Curve Loop(4) = {-24, 27, 30, 31,32,-12,-23,-10};", "Curve Loop(4) = {-105, 27, 30, 31,32,-12,-23,-10};")
    s = sub1(s, "Plane Surface(5) = {5};", "Plane Surface(5) = {5};\n\n// source patch\nCurve Loop(6) = {101, -104, -105, -103, 24, 11};\nPlane Surface(6) = {6};")
    s = sub1(s, "Transfinite Surface {1} = {14,13, 1, 10};", "Transfinite Surface {1} = {14,13, 101, 102};")
    s = sub1(s, "Transfinite Surface {4} = {27, 14, 22, 10};", "Transfinite Surface {4} = {27, 14, 103, 102};")
    s = sub1(s, "Transfinite Curve {11, -21} = N_out_h Using Progression p_out_h;",
             "Transfinite Curve {11, 104, -21} = N_out_h Using Progression p_out_h;")
    s = sub1(s, "Transfinite Curve {24, 32} = N_out_h_upper Using Progression p_out_h_upper;",
             "Transfinite Curve {24, 105, 32} = N_out_h_upper Using Progression p_out_h_upper;")
    s = sub1(s, "Transfinite Curve {10, 1, 27} = N_out_v_inflow_sparse Using Progression p_out_v_inflow_sparse;",
             "Transfinite Curve {10, 1, 27} = N_out_v_inflow_sparse - j_patch Using Progression p_out_v_inflow_sparse;\n"
             "Transfinite Curve {103} = j_patch + 1 Using Progression p_out_v_inflow_sparse; // top of the patch keeps the original nodes (triangle region unchanged)")
    s = sub1(s, "Physical Curve(4) = {1, 20,", "Physical Curve(4) = {101, 1, 20,")
    return s + fields(hmin, "101", surf=6)


def read_msh2(f):
    L = open(f).read().split("\n")
    i = L.index("$Nodes"); n = int(L[i + 1])
    nodes = {int(l.split()[0]): tuple(map(float, l.split()[1:4])) for l in L[i + 2:i + 2 + n]}
    i = L.index("$Elements"); n = int(L[i + 1])
    elems = []
    for l in L[i + 2:i + 2 + n]:
        v = list(map(int, l.split()))
        nt = v[2]
        elems.append((v[1], v[3:3 + nt], v[3 + nt:]))  # type, tags (physical, elementary), node ids
    return nodes, elems


def splice(orig, refined, xa, ytop, out):
    import numpy as np
    from scipy.spatial import cKDTree
    no, eo = read_msh2(orig)
    nr, er = read_msh2(refined)
    ids = sorted(no)
    tree = cKDTree(np.array([no[i][:2] for i in ids]))
    cen = lambda nodes, e: np.mean([nodes[k][:2] for k in e[2]], axis=0)
    inpatch = lambda c: c[0] < xa - 1e-9 and c[1] < ytop - 1e-9
    keep = [e for e in eo if not ((e[0] in (2, 3) or (e[0] == 1 and e[1][0] == 4)) and inpatch(cen(no, e)))]
    removed = len(eo) - len(keep)
    add = [e for e in er if (e[0] in (2, 3) and e[1][0] == 103) or (e[0] == 1 and e[1][0] == 4 and inpatch(cen(nr, e)))]
    nmap, nxt = {}, max(ids) + 1
    nodes = dict(no)
    for e in add:
        for k in e[2]:
            if k in nmap:
                continue
            d, j = tree.query(nr[k][:2])
            if d < 1e-6:
                nmap[k] = ids[j]
            else:
                nmap[k] = nxt; nodes[nxt] = nr[k]; nxt += 1
    add = [(t, tags, [nmap[k] for k in nn]) for t, tags, nn in add]
    els = keep + add
    used = {k for e in els for k in e[2]}
    with open(out, "w") as f:
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n%d\n" % len(used))
        for k in sorted(used):
            f.write("%d %.16g %.16g %.16g\n" % ((k,) + nodes[k]))
        f.write("$EndNodes\n$Elements\n%d\n" % len(els))
        for i, (t, tags, nn) in enumerate(els, 1):
            f.write(" ".join(map(str, [i, t, len(tags)] + tags + nn)) + "\n")
        f.write("$EndElements\n")
    print("splice: removed %d original elements, added %d (patch), %d new nodes, %d unused original nodes dropped"
          % (removed, len(add), nxt - max(ids) - 1, len(no) - len(used & set(no))))


if __name__ == "__main__":
    import sys
    if sys.argv[1] == "splice":
        splice(sys.argv[2], sys.argv[3], float(sys.argv[4]), float(sys.argv[5]), sys.argv[6])
        sys.exit()
    ap = argparse.ArgumentParser()
    ap.add_argument("kind", choices=["gap", "flat", "bfs"])
    ap.add_argument("geo")
    ap.add_argument("out")
    ap.add_argument("--hmin", type=float, default=0.05)
    a = ap.parse_args()
    s = open(a.geo).read()
    s = {"gap": gap, "flat": flat, "bfs": bfs}[a.kind](s, a.hmin)
    open(a.out, "w").write(s)
