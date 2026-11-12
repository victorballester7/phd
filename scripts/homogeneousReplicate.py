"""Replicate a quasi-3D (3DH1D) Nektar++ field n times in the spanwise direction.

The field on Lz is mapped exactly onto n Lz: Fourier mode k of the input becomes
mode n k of the output, and the modes in between are set to zero. The number of
planes grows by the same factor, so Delta z is unchanged. (FieldConvert's
homstretch instead stretches the field and keeps the number of planes.)

Nektar++ stores the homogeneous coefficients plane by plane: plane 0 is mode 0,
plane 1 the Nyquist cosine, and planes 2k, 2k + 1 the cosine and sine of mode k.
Within an ELEMENTS block the data are ordered [field][element][plane][2D coeffs].

Works on a single .fld file or a parallel .chk directory (Info.xml plus
P*.fld). Only the homogeneous attributes and the data change.

Usage: python3 homogeneousReplicate.py in.chk out.chk [n=2]
"""

import base64
import os
import re
import shutil
import sys
import zlib

import numpy as np


def _expand_ids(s):
    n = 0
    for part in s.split(","):
        a, _, b = part.partition("-")
        n += int(b) - int(a) + 1 if b else 1
    return n


def _block(m, n):
    tag, data = m.group(1), m.group(2)
    nf = len(re.search(r'FIELDS="([^"]*)"', tag)[1].split(","))
    ne = _expand_ids(re.search(r' ID="([^"]*)"', tag)[1])
    nm = re.search(r'NUMMODESPERDIR="UNIORDER:([^"]*)"', tag)[1].split(",")
    npl = int(nm[-1])
    assert "B64Z-LittleEndian" in tag and 'BITSIZE="64"' in tag
    a = np.frombuffer(zlib.decompress(base64.b64decode(data)), "<f8")
    nc, r = divmod(a.size, nf * ne * npl)
    assert r == 0, "data length is not fields x elements x planes x coeffs"
    a = a.reshape(nf, ne, npl, nc)
    b = np.zeros((nf, ne, n * npl, nc))
    b[:, :, 0] = a[:, :, 0]                 # mode 0
    b[:, :, 1] = a[:, :, 1]                 # Nyquist: mode npl/2 -> n npl/2
    for k in range(1, npl // 2):
        b[:, :, 2 * n * k:2 * n * k + 2] = a[:, :, 2 * k:2 * k + 2]
    lz = float(re.search(r'HOMOGENEOUSLENGTHS="([^"]*)"', tag)[1])
    tag = re.sub(r'HOMOGENEOUSLENGTHS="[^"]*"', f'HOMOGENEOUSLENGTHS="{n * lz:g}"', tag)
    tag = re.sub(r'HOMOGENEOUSZIDS="[^"]*"',
                 'HOMOGENEOUSZIDS="' + ",".join(map(str, range(n * npl))) + '"', tag)
    tag = re.sub(r'NUMMODESPERDIR="UNIORDER:[^"]*"',
                 'NUMMODESPERDIR="UNIORDER:' + ",".join(nm[:-1] + [str(n * npl)]) + '"', tag)
    enc = base64.b64encode(zlib.compress(b.astype("<f8").tobytes())).decode()
    return f"{tag}{enc}</ELEMENTS>"


def convert_file(src, dst, n):
    s = open(src).read()
    s = re.sub(r"(<ELEMENTS[^>]*>)([^<]*)</ELEMENTS>", lambda m: _block(m, n), s)
    with open(dst, "w") as f:
        f.write(s)


def main():
    src, dst = sys.argv[1], sys.argv[2]
    n = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    if os.path.isdir(src):
        os.makedirs(dst, exist_ok=True)
        for f in sorted(os.listdir(src)):
            if f.endswith(".fld"):
                convert_file(os.path.join(src, f), os.path.join(dst, f), n)
            else:
                shutil.copy(os.path.join(src, f), dst)
    else:
        convert_file(src, dst, n)


if __name__ == "__main__":
    main()
