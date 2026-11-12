"""Render lambda_2 IsoVolume images of a gap DNS (or, with --mode, Re(u_x) of an
eigenmode). The case (geometry, camera, lambda_2 range) is set in camera.py and
the walls are drawn afterwards in TikZ by gen_tex.py.

usage: pvbatch render.py <in.vtu> <out.png> [--expand] [--left p] [--scale s] [--mode [--thr t]]

--expand  extend the frame upwards (same camera, same pixel scale) until the
          whole isovolume fits vertically
--left    extend the frame to the left by a fraction p of its width (same camera);
          p is stored in the PNG ('left' text chunk) for gen_tex.py
--scale   resolution factor w.r.t. the paper resolution (e.g. 0.2 for tests)
--mode    eigenmode: IsoVolumes of |Re(u_x)|/max|Re(u_x)| >= 10^t (default t=-1.5),
          coloured on the two-sided log scale of Images/2D_3Dmode/pdfs/colorbar.pdf
          (10^-2.5 ... 1); a 2D mode (z-independent) is extruded over [0, LZ]
"""
import os, sys, math
import numpy as np
from PIL import Image
from PIL.PngImagePlugin import PngInfo
from paraview.simple import *
from vtkmodules.util.numpy_support import vtk_to_numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from camera import *

args = sys.argv[1:]
vtu, png = args[:2]
expand = '--expand' in args
scale = float(args[args.index('--scale') + 1]) if '--scale' in args else 1.0
left = float(args[args.index('--left') + 1]) if '--left' in args else 0.0
mode = '--mode' in args
logthr = float(args[args.index('--thr') + 1]) if '--thr' in args else -1.5
LOGMIN = -2.5  # lower end of the two-sided log colour scale
resX, resY0 = int(RES_X * scale), int(RES_Y * scale)
thr = L2_THR

# colorbrewer RdBu (11 classes), as pgfplots' colormap/RdBu
RDBU = ['67001f', 'b2182b', 'd6604d', 'f4a582', 'fddbc7', 'f7f7f7',
        'd1e5f0', '92c5de', '4393c3', '2166ac', '053061']



def extrude(src, nz=2):
    """Extrude a 2D unstructured grid (z = 0) over z in [0, LZ] into a 3D one."""
    from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid, VTK_WEDGE, VTK_HEXAHEDRON
    from vtkmodules.vtkCommonCore import vtkPoints
    from vtkmodules.util.numpy_support import numpy_to_vtk
    g = servermanager.Fetch(src)
    p2 = vtk_to_numpy(g.GetPoints().GetData())[:, :2]
    n = len(p2)
    zs = np.linspace(0, LZ, nz + 1)
    pts = vtkPoints()
    pts.SetData(numpy_to_vtk(np.vstack([np.c_[p2, np.full(n, z)] for z in zs]), deep=True))
    out = vtkUnstructuredGrid()
    out.SetPoints(pts)
    for c in range(g.GetNumberOfCells()):
        ids = [g.GetCell(c).GetPointId(i) for i in range(g.GetCell(c).GetNumberOfPoints())]
        typ = VTK_WEDGE if len(ids) == 3 else VTK_HEXAHEDRON
        for k in range(nz):
            out.InsertNextCell(typ, len(ids) * 2, [i + k * n for i in ids] + [i + (k + 1) * n for i in ids])
    u = vtk_to_numpy(g.GetPointData().GetArray('u'))
    arr = numpy_to_vtk(np.tile(u, nz + 1), deep=True)
    arr.SetName('u')
    out.GetPointData().AddArray(arr)
    tp = TrivialProducer()
    tp.GetClientSideObject().SetOutput(out)
    tp.UpdatePipeline()
    return tp


reader = XMLUnstructuredGridReader(FileName=[vtu])
reader.PointArrayStatus = ['u'] if mode else ['L2']
reader.UpdatePipeline()
x0, x1, y0, y1, z0, z1 = reader.GetDataInformation().GetBounds()

if mode:
    src = reader
    if z1 - z0 < 1e-9:
        src = extrude(reader)
    umax = max(abs(r) for r in reader.PointData['u'].GetRange())
    # r = u/umax is thresholded (linear, so the surfaces are smooth); the colour is
    # s = sign(u) (log10|r| - LOGMIN)/(-LOGMIN), in [-1, 1] for |r| >= 10^LOGMIN
    calc = Calculator(Input=src)
    calc.ResultArrayName = 'r'
    calc.Function = f'u/{umax!r}'
    calc = Calculator(Input=calc)
    calc.ResultArrayName = 's'
    calc.Function = f'sign(r)*max(log10(abs(r)+1e-30)-({LOGMIN}),0)/({-LOGMIN})'
    isos = []
    for rng in ([10**logthr, 2.0], [-2.0, -10**logthr]):
        iv = IsoVolume(Input=calc)
        iv.InputScalars = ['POINTS', 'r']
        iv.ThresholdRange = rng
        isos.append(iv)
    iso = AppendDatasets(Input=isos)
    color_array = 's'
else:
    calc = Calculator(Input=reader)
    calc.ResultArrayName = 'mL2'
    calc.Function = '-L2'
    iso = IsoVolume(Input=calc)
    iso.InputScalars = ['POINTS', 'L2']
    iso.ThresholdRange = [-1e30, -thr]
    color_array = 'mL2'
iso.UpdatePipeline()

# top of the frame, in units of the original half-height above the image centre
top = 1.0
if expand:
    pts = vtk_to_numpy(iso.GetClientSideObject().GetOutputDataObject(0).GetPoints().GetData())
    pts = pts[depth(pts) > 1.0]
    q = ndc(pts)
    q = q[(q[:, 0] >= -1.0 - 2 * left) & (q[:, 0] <= 1.0)]
    top = max(1.0, q[:, 1].max() + 0.03)
f = (top + 1) / 2  # height factor w.r.t. the original image
print(f'height factor {f:.3f}')

view = CreateView('RenderView')
view.OrientationAxesVisibility = 0
view.UseColorPaletteForBackground = 0
view.Background = [1, 1, 1]
# render symmetrically about the camera axis with a larger view angle (same pixel
# scale) and keep the upper part, so the original frame stays at the bottom
# (the view angle is vertical, so a wider view adds the same margin on both sides:
# the right one is cropped afterwards)
resYs = int(round(resY0 * top))
resXs = int(round(resX * (1 + 2 * left)))
view.ViewSize = [resXs, resYs]

disp = Show(iso, view)
ColorBy(disp, ('POINTS', color_array))
lut = GetColorTransferFunction(color_array)
if mode:
    # s = -1 ... 1 -> RdBu index 10 (dark blue) ... 0 (dark red), evenly spaced
    pts = []
    for k, c in enumerate(reversed(RDBU)):
        pts += [-1 + 2 * k / 10] + [int(c[i:i + 2], 16) / 255 for i in (0, 2, 4)]
    lut.RGBPoints = pts
    lut.ColorSpace = 'RGB'
    lut.AutomaticRescaleRangeMode = 'Never'
else:
    lut.ApplyPreset('Viridis', True)
    lut.InvertTransferFunction()
    lut.UseLogScale = 1
    lut.RescaleTransferFunction(*L2_RANGE)
    lut.AutomaticRescaleRangeMode = 'Never'
disp.SetScalarBarVisibility(view, False)

# walls, slightly sunk into the solid and rendered flat white: they hide the
# structures that lie behind them and are made transparent afterwards
eps = 0.02
walls = [
    ((x0, -eps, 0), (0, -eps, 0), (x0, -eps, LZ)),          # upstream plate
    ((0, -D - eps, 0), (W, -D - eps, 0), (0, -D - eps, LZ)),  # gap floor
    ((W + eps, -D, 0), (W + eps, 0, 0), (W + eps, -D, LZ)),   # downstream gap wall
    ((W, -eps, 0), (x1, -eps, 0), (W, -eps, LZ)),             # downstream plate
]
for o, p1, p2 in walls:
    pl = Plane(Origin=o, Point1=p1, Point2=p2)
    d = Show(pl, view)
    d.AmbientColor = d.DiffuseColor = [1, 1, 1]
    d.Ambient, d.Diffuse, d.Specular = 1.0, 0.0, 0.0

# Show()/Render() may reset the camera, so it is applied right before the screenshot
Render(view)
view.CameraPosition = list(CAM_POS)
view.CameraFocalPoint = list(CAM_FOC)
view.CameraViewUp = list(CAM_UP)
view.CameraViewAngle = 2 * math.degrees(math.atan(top * math.tan(math.radians(CAM_ANGLE / 2))))
view.CameraParallelProjection = 0
# off-centre principal point: vtkCamera shifts the image by minus the window
# centre, which is in NDC of the rendered window ([-top, top] vertically and
# 1 + 2 left times the width of the reference frame)
view.GetActiveCamera().SetWindowCenter(-NDC_OFFSET[0] / (1 + 2 * left), -NDC_OFFSET[1] / top)
SaveScreenshot(png, view, ImageResolution=[resXs, resYs], TransparentBackground=1)

im = np.array(Image.open(png).convert('RGBA'))
im = im[: int(round(resY0 * f)), : int(round(resX * (1 + left)))]  # rows run from the top
white = (im[..., :3] >= 250).all(axis=-1)
im[white, 3] = 0
info = PngInfo()
info.add_text('left', repr(left))
Image.fromarray(im).save(png, pnginfo=info)
print('done', png)
