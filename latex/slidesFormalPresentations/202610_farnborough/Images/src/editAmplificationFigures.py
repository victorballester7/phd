"""Derive slide variants of the qualitative n-factor sketches.

The original TikZ sources of amplificationQualititiveNfactor{,OnlyExp}.pdf are
lost, so the PDFs are edited directly:
  * amplificationQualititiveNfactorCurves.pdf : OnlyExp without the purple
    (x^trans / Delta n_exp) annotations -> only the n_flat and n_gap curves.
  * amplificationQualititiveNfactorOnlyExpLabels.pdf : OnlyExp with the
    same n_gap / n_flat labels.
  * amplificationQualititiveNfactorLabels.pdf : full sketch with n_gap / n_flat
    labels next to the top of each curve.
Labels reuse the subset fonts embedded in the figures (CMMI12 'n', CMSS10
'gap' and the 'fl' ligature + 'at'), so they match the existing text.

Run from the Images/ folder:  python3 src/editAmplificationFigures.py
"""
import pypdf
from pypdf.generic import ContentStream

BLUE = "0 0.44707 0.74118"
ORANGE = "0.85098 0.32549 0.09804"
PURPLE_LIGHT = "0.84824 0.7553 0.86705"


def label(x, y, colour, sub):
    """n_<sub> with baseline-left at (x, y) pt. The stream ends inside the
    picture's '44.096 33.385 cm' shift, so these are relative to that."""
    return (f"q {colour} rg BT /F44 14.3462 Tf {x} {y} Td [(n)]TJ "
            f"/F41 9.9626 Tf 8.4 -2.151 Td [({sub})]TJ ET Q\n")


GAP = "gap"
FLAT = "\\015at"  # 'fl' ligature, the only way the subset font spells it


def rewrite(src, dst, edit):
    reader = pypdf.PdfReader(src)
    writer = pypdf.PdfWriter(clone_from=reader)
    page = writer.pages[0]
    data = page.get_contents().get_data().decode("latin-1")
    stream = ContentStream(None, writer)
    stream.set_data(edit(data).encode("latin-1"))
    page.replace_contents(stream)
    page.compress_content_streams()
    writer.write(dst)


def drop_purple(data):
    lines = data.split("\n")
    first = next(i for i, l in enumerate(lines) if l.startswith(PURPLE_LIGHT))
    start = first - 1  # the 'q' opening the first purple scope
    gap = next(i for i, l in enumerate(lines) if l.startswith("/pgf@CA0.5 gs"))
    end = gap - 3  # stop before 'q', 'q', '0 g 0 G' that open the gap sketch
    block = lines[start:end]
    assert lines[start].strip() == "q"
    # Inline 'q ... Q' lines balance themselves; whole-line q/Q must too.
    depth = sum(l.strip() == "q" for l in block) - sum(l.strip() == "Q" for l in block)
    assert depth == 0, depth
    return "\n".join(lines[:start] + lines[end:])


def add_labels(data):
    return (data + "\n" + label(90.3, 134.7, ORANGE, GAP)
            + label(162.3, 136.6, BLUE, FLAT))


rewrite("amplificationQualititiveNfactorOnlyExp.pdf",
        "amplificationQualititiveNfactorCurves.pdf", drop_purple)
rewrite("amplificationQualititiveNfactorOnlyExp.pdf",
        "amplificationQualititiveNfactorOnlyExpLabels.pdf", add_labels)
rewrite("amplificationQualititiveNfactor.pdf",
        "amplificationQualititiveNfactorLabels.pdf", add_labels)
