"""Write the TikZ wrapper of a render.py image: projected wall faces, gap outline
and fades of the plates where the frame cuts them.

usage: python3 gen_tex.py <png> <out.tex> [--fade-strength s] [--fade-start a]

--fade-strength  how much the downstream plate fades: 1 vanishes completely
                 (default), 0 does not fade
--fade-start     where the downstream fade starts, as a fraction of the distance
                 from the gap to the frame edge (default FADE_START)
--front-faces    hide what lies behind the front (z = Lz) cut face of the wall below
                 the upstream plate (painted white, no outline): the floor edge seen
                 through the open front no longer makes the upstream step look inclined
"""
import os, sys
import numpy as np
from PIL import Image
from camera import W, D, LZ, X0, X1, RES_X, RES_Y, TIKZ_W, ndc, depth

args = sys.argv[1:]
png, out = args[:2]
strength = float(args[args.index('--fade-strength') + 1]) if '--fade-strength' in args else 1.0
name = os.path.splitext(os.path.basename(png))[0]

im = Image.open(png)
pw, ph = im.size
LEFT = float(im.text.get('left', 0))  # frame extension to the left (render.py --left)
S = TIKZ_W / ((1 + LEFT) * RES_X)  # cm per pixel of the reference frame
H0 = RES_Y * S  # height of the reference frame (cm)
H = TIKZ_W * ph / pw  # height of this (possibly expanded) frame (cm)

# fades of the plates along the streamwise direction, away from the gap: from
# FADE_START of the way between the gap edge and the frame edge until they vanish
# at FADE_END of that distance (downstream; upstream uses FADE_START_UP)
FADE_START, FADE_END = 0.2, 0.97
FADE_START_UP = 0.15
if '--fade-start' in args:
    FADE_START = float(args[args.index('--fade-start') + 1])
NLEV = 40  # grey levels of the fades


def tikz(p):
    """3D points -> TikZ coordinates (cm); the reference frame sits at the bottom."""
    q = ndc(np.asarray(p, float))
    return np.c_[((q[:, 0] + 1) / 2 + LEFT) * RES_X * S, (q[:, 1] + 1) / 2 * H0]


def clip3d(poly, near=1.0):
    """Clip a 3D polygon to the half-space in front of the camera."""
    out = []
    for a, b in zip(poly, poly[1:] + poly[:1]):
        da, db = depth(a)[0] - near, depth(b)[0] - near
        if da >= 0:
            out.append(a)
        if da * db < 0:
            out.append(a + da / (da - db) * (b - a))
    return out


def clip2d(poly, lo, hi):
    """Sutherland-Hodgman clipping of a 2D polygon to a rectangle."""
    for ax in (0, 1):
        for sgn, lim in ((1, lo[ax]), (-1, hi[ax])):
            out = []
            for a, b in zip(poly, poly[1:] + poly[:1]):
                da, db = sgn * (a[ax] - lim), sgn * (b[ax] - lim)
                if da >= 0:
                    out.append(a)
                if da * db < 0:
                    out.append(a + da / (da - db) * (b - a))
            poly = out
            if not poly:
                return []
    return poly


def face(poly, hi=(TIKZ_W, H)):
    poly = clip3d([np.array(p, float) for p in poly])
    return clip2d(list(tikz(poly)), (-0.5, -0.5), (hi[0] + 0.5, hi[1] + 0.5)) if len(poly) >= 3 else []


def path(pts, cycle=True):
    s = ' -- '.join(f'({x:.4f},{y:.4f})' for x, y in pts)
    return s + (' -- cycle' if cycle else '')


def smoothstep(x):
    x = np.clip(x, 0, 1)
    return 3 * x**2 - 2 * x**3


def smoothstep_inv(y):
    return 0.5 - np.sin(np.arcsin(1 - 2 * y) / 3)


def clip_half(poly, n, c):
    """Clip a 2D polygon to the half-plane n.p >= c."""
    out = []
    for a, b in zip(poly, poly[1:] + poly[:1]):
        da, db = n @ a - c, n @ b - c
        if da >= 0:
            out.append(a)
        if da * db < 0:
            out.append(a + da / (da - db) * (b - a))
    return out


def faded(poly, edges, halfplanes, opacity, fill_pct, edge_pct):
    """TikZ code of a polygon that fades to white without transparency.

    The opacity field is quantised in NLEV levels; level k covers the region
    where opacity >= k/NLEV, i.e. the intersection of halfplanes(k/NLEV), and is
    painted over the previous one in a lighter grey, so there are no seams.
    """
    lines = []
    for k in range(1, NLEV + 1):
        o = k / NLEV
        reg = poly
        for n, c in halfplanes(o):
            reg = clip_half(reg, np.asarray(n, float), c)
            if len(reg) < 3:
                break
        if len(reg) >= 3:
            lines.append(f'\\fill[black!{fill_pct * o:.3f}] {path(reg)};')
    # edges: pieces of constant (quantised) opacity
    for a, b in edges:
        ts = np.linspace(0, 1, 401)
        pts = a + ts[:, None] * (b - a)
        mid = (pts[1:] + pts[:-1]) / 2
        lev = np.floor(np.array([opacity(m) for m in mid]) * NLEV + 1e-9).astype(int)
        i = 0
        while i < len(lev):
            j = i
            while j + 1 < len(lev) and lev[j + 1] == lev[i]:
                j += 1
            if lev[i] > 0:
                lines.append(f'\\draw[wall edge, draw=black!{edge_pct * lev[i] / NLEV:.3f}] '
                             f'{path([pts[i], pts[j + 1]], cycle=False)};')
            i = j + 1
    return '\n\t\t'.join(lines)


def edge(a, b):
    """Projected 3D segment, clipped to the near plane and to the frame."""
    seg = clip3d([np.array(a, float), np.array(b, float)])
    if len(seg) < 2:
        return None
    a, b = tikz(seg[:2])
    # Liang-Barsky clipping to the frame (with the same margin as the faces)
    lo, hi = np.array([-0.5, -0.5]), np.array([TIKZ_W + 0.5, H + 0.5])
    t0, t1, d = 0.0, 1.0, b - a
    for ax in (0, 1):
        for p, q in ((-d[ax], a[ax] - lo[ax]), (d[ax], hi[ax] - a[ax])):
            if p == 0:
                if q < 0:
                    return None
                continue
            r = q / p
            if p < 0:
                t0 = max(t0, r)
            else:
                t1 = min(t1, r)
    return (a + t0 * d, a + t1 * d) if t0 < t1 else None


floor = face([(0, -D, 0), (W, -D, 0), (W, -D, LZ), (0, -D, LZ)])
wall = face([(W, -D, 0), (W, 0, 0), (W, 0, LZ), (W, -D, LZ)])
plate_up = face([(X0, 0, 0), (0, 0, 0), (0, 0, LZ), (X0, 0, LZ)])
plate_down = face([(W, 0, 0), (X1, 0, 0), (X1, 0, LZ), (W, 0, LZ)])
front_edge = tikz([(0, 0, LZ), (0, -D, LZ)])  # front edge of the upstream gap wall



def streamwise_fade(plate, xg, xf, edges, start, strength=1.0, fill_pct=6, edge_pct=80):
    """Fade a plate along its projected streamwise direction, from the gap edge
    x = xg towards x = xf. The distances are measured on the reference frame, so
    every panel fades at the same place, and end before the plate meets the frame."""
    pg, pf = tikz([(xg, 0, LZ / 2), (xf, 0, LZ / 2)])
    u = (pf - pg) / np.linalg.norm(pf - pg)
    xe = X1 if xf > xg else X0  # far end of the plate
    cut = face([(xg, 0, 0), (xe, 0, 0), (xe, 0, LZ), (xg, 0, LZ)], hi=(TIKZ_W, H0))
    cut = clip2d(cut, (0, 0), (TIKZ_W, H0))
    on_edge = [p for p in cut if min(abs(p[0]), abs(p[0] - TIKZ_W), abs(p[1]), abs(p[1] - H0)) < 1e-6]
    t_edge = min((p - pg) @ u for p in on_edge if (p - pg) @ u > 1e-6)
    t0, t1 = start * t_edge, FADE_END * t_edge
    return faded(
        plate,
        [e for e in edges if e],
        # opacity 1 - strength*smoothstep(x): it is >= o where x <= x(o)
        lambda o: [] if (1 - o) / strength >= 1 else
        [(-u, -(pg @ u + t0 + (t1 - t0) * smoothstep_inv((1 - o) / strength)))],
        lambda p: 1 - strength * smoothstep(((p - pg) @ u - t0) / (t1 - t0)),
        fill_pct, edge_pct)


down = streamwise_fade(
    plate_down, W, X1,
    (edge((W, 0, 0), (X1, 0, 0)), edge((W, 0, LZ), (X1, 0, LZ))),
    FADE_START, strength)
# upstream the far end of the plate lies behind the camera: a point close to
# the gap gives the projected streamwise direction
up = streamwise_fade(
    plate_up, 0.0, -5.0,
    (edge((X0, 0, 0), (0, 0, 0)), edge((X0, 0, LZ), (0, 0, LZ)), edge((0, 0, 0), (0, 0, LZ))),
    FADE_START_UP)

# front (z = Lz) cut face of the wall below the upstream plate, painted white so that
# it only hides the part of the floor that lies behind it
fronts = ''
if '--front-faces' in args:
    fronts = ('% front (z = Lz) cut face of the wall below the upstream plate\n\t\t'
              '\\fill[white] ' + path(face([(X0, -D, LZ), (0, -D, LZ), (0, 0, LZ), (X0, 0, LZ)])) + ';')
    # lower edge of that face, fading upstream with the plate, so the step reads as a block
    fronts += '\n\t\t' + streamwise_fade(
        [], 0.0, -5.0, (edge((X0, -D, LZ), (0, -D, LZ)),), FADE_START_UP)

tex = r"""\documentclass{standalone}

\usepackage{tikz}

\def\image{../pngs/NAME.png}
\def\figurewidth{FIGW}

\begin{document}
\begin{tikzpicture}[
		line join=round,
		line cap=round,
		wall edge/.style={draw=black!80, line width=0.6pt},
		floor/.style={wall edge, fill=black!15},
		wall/.style={wall edge, fill=black!25},
	]
	\def\xmin{0}
	\def\xmax{\figurewidth}
	\def\ymin{0}
	\def\ymax{YMAX}
	\path[use as bounding box]
	(\xmin,\ymin)
	rectangle
	(\xmax,\ymax);

	\begin{scope}
		\clip (\xmin,\ymin) rectangle (\xmax,\ymax);
		\filldraw[floor] FLOOR;
		\filldraw[wall] WALL;
		% downstream plate: fades out along the streamwise direction
		PLATEDOWN
		FRONTFACES
		% upstream plate: fades out along the streamwise direction
		PLATEUP
		\draw[wall edge] FRONTEDGE;
	\end{scope}
	\node[
		inner sep=0pt,
		anchor=south west
	] at (0,0)
	{
		\includegraphics[
			width=\figurewidth cm
		]{\image}
	};

\end{tikzpicture}
\end{document}
"""
for k, v in {
    'NAME': name,
    'FIGW': f'{TIKZ_W:g}',
    'YMAX': f'{H:.4f}',
    'FLOOR': path(floor),
    'WALL': path(wall),
    'PLATEDOWN': down,
    'PLATEUP': up,
    'FRONTFACES': fronts,
    'FRONTEDGE': path(front_edge, cycle=False),
}.items():
    tex = tex.replace(k, v)
open(out, 'w').write(tex)
