"""Base flow on a spliced source-refined mesh (scripts/refineSourceMesh.py), exact outside the patch (2026-11-10).

FieldConvert interpfield is not exact even between identical meshes (~3e-3 in u for the incompressible gap base
flows, which are elementwise P3 polynomials that are slightly discontinuous across elements; the solver uses each
element's own polynomial). So: every element of the new mesh whose ordered edge/vertex sequence is identical to an
element of the old mesh gets the old element's coefficients copied verbatim; only the remaining (patch) elements take
their coefficients from an interpfield result on the new mesh.

usage: python3 spliceBaseflow.py <old mesh.xml> <old baseflow.fld> <new mesh.xml> <interpfield.fld on new mesh> <out.fld>
"""
import base64
import re
import sys
import zlib

import numpy as np

I8 = "<i8"


def dec(s, dt):
    return np.frombuffer(zlib.decompress(base64.b64decode(s.strip())), dtype=dt)


def enc(a):
    return base64.b64encode(zlib.compress(np.ascontiguousarray(a, dtype="<f8").tobytes())).decode()


def mesh(path):
    m = open(path).read()
    v = dec(re.search(r"<VERTEX[^>]*>(.*?)</VERTEX>", m, re.S).group(1), np.dtype([("id", I8), ("x", "<f8"), ("y", "<f8"), ("z", "<f8")]))
    xy = {int(r["id"]): (round(float(r["x"]), 7), round(float(r["y"]), 7)) for r in v}
    e = dec(re.search(r"<EDGE[^>]*>(.*?)</EDGE>", m, re.S).group(1), np.dtype([("id", I8), ("v0", I8), ("v1", I8)]))
    ed = {int(r["id"]): (xy[int(r["v0"])], xy[int(r["v1"])]) for r in e}
    sig = {}
    for tag, n in (("T", 3), ("Q", 4)):
        g = re.search(r"<%s [^>]*>(.*?)</%s>" % (tag, tag), m, re.S)
        if not g:
            continue
        a = dec(g.group(1), np.dtype([("id", I8)] + [("e%d" % k, I8) for k in range(n)]))
        for r in a:
            es = [set(ed[int(r["e%d" % k])]) for k in range(n)]
            # element vertex traversal (Nektar Tri/QuadGeom): v_k = e_{k-1} & e_k ... rotated so that v0 = e0 & e_{n-1};
            # the local coordinates, hence the coefficient layout, depend on this order only (not on edge directions)
            vs = tuple(next(iter(es[k - 1] & es[k])) for k in range(n))
            sig[int(r["id"])] = (tag, vs)
    return sig


def ids(s):
    out = []
    for t in s.split(","):
        if "-" in t:
            a, b = map(int, t.split("-"))
            out += range(a, b + 1)
        else:
            out.append(int(t))
    return out


def fld(path):
    """{element id: (shape, nummodes, fields, coeff array [nfields, ncoeff])}, header text"""
    f = open(path).read()
    out = {}
    for g in re.finditer(r"<ELEMENTS ([^>]*)>(.*?)</ELEMENTS>", f, re.S):
        at = dict(re.findall(r'(\w+)="([^"]*)"', g.group(1)))
        el = ids(at["ID"])
        fields = at["FIELDS"].split(",")
        d = dec(g.group(2), "<f8").reshape(len(fields), len(el), -1)
        for k, i in enumerate(el):
            out[i] = (at["SHAPE"], at["NUMMODESPERDIR"], at["FIELDS"], at["BASIS"], d[:, k, :])
    meta = re.search(r"<Metadata>.*?</Metadata>", f, re.S).group(0)
    return out, meta


def main(old_mesh, old_fld, new_mesh, interp_fld, out):
    so, sn = mesh(old_mesh), mesh(new_mesh)
    inv = {}
    for i, s in so.items():
        inv.setdefault(s, i)
    fo, meta = fld(old_fld)
    fi, _ = fld(interp_fld)
    groups = {}
    ncopy = ninterp = 0
    for i in sorted(sn):
        j = inv.get(sn[i])
        if j is not None and j in fo:
            rec = fo[j]; ncopy += 1
        else:
            rec = fi[i]; ninterp += 1
        key = rec[:4]
        groups.setdefault(key, []).append((i, rec[4]))
    with open(out, "w") as f:
        f.write('<?xml version="1.0" encoding="utf-8" ?>\n<NEKTAR>\n    %s\n' % meta)
        for (shape, nm, fields, basis), lst in groups.items():
            data = np.stack([c for _, c in lst], axis=1)  # [nfields, nel, ncoeff] -> field-major
            f.write('    <ELEMENTS FIELDS="%s" SHAPE="%s" BASIS="%s" NUMMODESPERDIR="%s" ID="%s" COMPRESSED="B64Z-LittleEndian" BITSIZE="64">%s</ELEMENTS>\n'
                    % (fields, shape, basis, nm, ",".join(str(i) for i, _ in lst), enc(data.ravel())))
        f.write("</NEKTAR>\n")
    print("%s: %d elements copied from the old base flow, %d interpolated (patch)" % (out, ncopy, ninterp))


if __name__ == "__main__":
    main(*sys.argv[1:6])
