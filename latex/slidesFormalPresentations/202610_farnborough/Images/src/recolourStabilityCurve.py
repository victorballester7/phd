"""Recolour the experimental bypass curve (tab:brown) of the slide copies of
3d2dStabilityCurveRe1000.pdf and its NoInset variant to black, without re-exporting the figure.
(scripts/pp/plots/3d2dStabilityCurveRe1000.py already uses black.)

Run from the Images/ folder:  python3 src/recolourStabilityCurve.py
"""
import pypdf
from pypdf.generic import ContentStream

BROWN = "0.5490196078 0.337254902 0.2941176471"

def recolour(src, dst, n):
    reader = pypdf.PdfReader(src)
    writer = pypdf.PdfWriter(clone_from=reader)
    page = writer.pages[0]
    data = page.get_contents().get_data().decode("latin-1")
    assert data.count(BROWN) == n, data.count(BROWN)
    stream = ContentStream(None, writer)
    stream.set_data(data.replace(BROWN, "0 0 0").encode("latin-1"))
    page.replace_contents(stream)
    page.compress_content_streams()
    writer.write(dst)


recolour("3d2dStabilityCurveRe1000.pdf", "3d2dStabilityCurveRe1000Black.pdf", 3)  # main axes, inset, legend
recolour("3d2dStabilityCurveRe1000NoInset.pdf", "3d2dStabilityCurveRe1000NoInsetBlack.pdf", 2)  # main axes, legend
