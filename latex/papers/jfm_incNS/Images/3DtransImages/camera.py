"""Camera and geometry shared by render.py (pvbatch) and gen_tex.py (python3).

The case is selected with the GAP_CASE environment variable (default d4w19).
"""
import math, os
import numpy as np

CASE = os.environ.get('GAP_CASE', 'd4w19')

if CASE == 'd4w19':
    # geometry (mesh.geo / session.xml)
    W, D, LZ = 19.0, 4.0, 20.0
    X0, X1 = -100.0, 1019.0  # streamwise extent of the domain

    # reference (paper) resolution: ndc_y in [-1, 1] spans this frame
    RES_X, RES_Y = 7632, 3748
    TIKZ_W = 6.62  # figure width in cm

    # camera from Images/3dTransParaview/camaraviewTransd4w19.pvcc
    CAM_POS = np.array([-18.82844914780593, 37.554888112401606, 52.46179957477808])
    CAM_FOC = np.array([82.62030503429251, -39.90063129684653, -29.137813977026564])
    CAM_UP = np.array([0.3027752863017271, 0.8502292903816899, -0.4306242907476947])
    CAM_ANGLE = 30.0  # vertical view angle (deg)
    # the .pvcc was set for Lz=16: shift along z so the front (z=Lz) edge keeps its place
    CAM_POS[2] += 4.0
    CAM_FOC[2] += 4.0

    L2_RANGE = (5e-7, 1e-3)  # colour scale of -lambda_2
    L2_THR = 5e-7  # isovolume -lambda_2 >= L2_THR
elif CASE == 'd1.5w80':
    W, D, LZ = 80.0, 1.5, 5.0
    X0, X1 = -100.0, 1080.0

    RES_X, RES_Y = 7632, 1668
    TIKZ_W = 13.5

    # view fitted to the gap outline and the roller positions of the earlier version
    # of the figure (pdfs/d1.5w80Lz5_3d_L2crit_31.tex, git history), which was an
    # off-centre crop: the principal point lies NDC_OFFSET below the frame centre
    CAM_POS = np.array([-76.83060219231919, 51.38042258190629, 157.34463293792976])
    CAM_FOC = np.array([38.983260227456185, -79.78749931434987, 60.48326015790474])
    CAM_UP = np.array([0.43396177154842874, 0.7508440309139675, -0.4979060373960275])
    CAM_ANGLE = 11.59619976896764
    NDC_OFFSET = (0.0, -5.41253590332756)

    L2_RANGE = (5e-5, 1e-1)
    L2_THR = 3e-4
elif CASE == 'd1.5w70':
    # forced region-III case (seed amplitude 1e-3), same view as d1.5w80; Lz = 6 instead
    # of 5, so the camera is shifted along z to keep the front (z = Lz) edge in place
    W, D, LZ = 70.0, 1.5, 6.0
    X0, X1 = -100.0, 1070.0

    RES_X, RES_Y = 7632, 1668
    TIKZ_W = 13.5

    CAM_POS = np.array([-76.83060219231919, 51.38042258190629, 158.34463293792976])
    CAM_FOC = np.array([38.983260227456185, -79.78749931434987, 61.48326015790474])
    CAM_UP = np.array([0.43396177154842874, 0.7508440309139675, -0.4979060373960275])
    CAM_ANGLE = 11.59619976896764
    NDC_OFFSET = (0.0, -5.41253590332756)

    L2_RANGE = (5e-5, 1e-1)
    L2_THR = 3e-4
else:
    raise ValueError(f'unknown GAP_CASE {CASE!r}')
if 'NDC_OFFSET' not in globals():
    NDC_OFFSET = (0.0, 0.0)  # NDC of the principal point (vtkCamera window centre)

LOOK = (CAM_FOC - CAM_POS) / np.linalg.norm(CAM_FOC - CAM_POS)
RIGHT = np.cross(LOOK, CAM_UP)
RIGHT /= np.linalg.norm(RIGHT)
UP = np.cross(RIGHT, LOOK)  # view-up orthogonalised as in VTK


def ndc(p):
    """Normalised device coordinates (x, y) of the reference frame for points p (n, 3)."""
    r = np.atleast_2d(p) - CAM_POS
    depth = r @ LOOK
    t = math.tan(math.radians(CAM_ANGLE / 2))
    return np.c_[(r @ RIGHT) / (depth * t * RES_X / RES_Y) + NDC_OFFSET[0],
                 (r @ UP) / (depth * t) + NDC_OFFSET[1]]


def depth(p):
    return (np.atleast_2d(p) - CAM_POS) @ LOOK
