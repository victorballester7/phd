#!/usr/bin/env bash
# Build the JFM submission files:
#   submission/main_review.pdf     -> "Main Document" (must be < 6 MB)
#   submission/cover_letter.pdf    -> "Cover Letter"
# The manuscript is compiled in a temporary copy so that it does not clash with a
# running `latexmk -pvc` (VimTeX) in the project directory.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$ROOT/submission"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

cp "$ROOT"/{main.tex,jfm.cls,jfm.bst,jfm.bib,lineno-FLM.sty} "$TMP"/
ln -s "$ROOT/Images" "$TMP/Images"
(cd "$TMP" && latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex >/dev/null)

# Downsample raster content (the lambda_2 renderings, natively ~1400-3400 dpi) and re-encode
# as high-quality JPEG; vector content is untouched. Note: -dJPEGQ is ignored by pdfwrite, so
# the quality is set through the distiller ColorImageDict: QFactor 0.15 (~JPEG q95) and no
# chroma subsampling (gs default is QFactor 0.9 with 2x2 subsampling, which blurs the figures).
# Tries the highest resolution first and falls back if the file exceeds ~5.5 MB (limit: 6 MB).
compress() {
  gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -dCompatibilityLevel=1.6 \
     -dDetectDuplicateImages=true \
     -dDownsampleColorImages=true -dColorImageDownsampleType=/Bicubic -dColorImageResolution="$1" \
     -dColorImageDownsampleThreshold=1.0 -dAutoFilterColorImages=false -dColorImageFilter=/DCTEncode \
     -dDownsampleGrayImages=true -dGrayImageResolution="$1" \
     -sOutputFile="$OUT/main_review.pdf" \
     -c "<< /ColorImageDict << /QFactor 0.15 /Blend 1 /HSamples [1 1 1 1] /VSamples [1 1 1 1] >> >> setdistillerparams" \
     -f "$TMP/main.pdf"
}
for dpi in 750 600 450 300; do
  compress "$dpi"
  [ "$(stat -c %s "$OUT/main_review.pdf")" -lt 5500000 ] && break
done
echo "raster figures downsampled to $dpi dpi"

(cd "$OUT" && latexmk -pdf -interaction=nonstopmode -halt-on-error cover_letter.tex >/dev/null && latexmk -c cover_letter.tex >/dev/null)

size=$(stat -c %s "$OUT/main_review.pdf")
echo "main_review.pdf: $((size / 1024)) KB ($(pdfinfo "$OUT/main_review.pdf" | awk '/Pages/{print $2}') pages)"
[ "$size" -lt 6000000 ] || echo "WARNING: main_review.pdf exceeds the 6 MB JFM limit"
if pdftotext "$OUT/main_review.pdf" - | grep "TODO" >/dev/null; then
  echo "WARNING: main_review.pdf still contains [TODO: ...] placeholders"
fi
if pdftotext "$OUT/cover_letter.pdf" - | grep "TODO" >/dev/null; then
  echo "WARNING: cover_letter.pdf still contains [TODO: ...] placeholders"
fi
words=$(wc -w < "$OUT/statement_of_significance.txt")
echo "statement_of_significance.txt: $words words (JFM: 50-100)"
