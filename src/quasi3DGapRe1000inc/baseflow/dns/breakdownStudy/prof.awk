# per-x wall-normal integrals (trapezoid, y <= 150) of uvw=|u'|^2 and ww=w'^2 of a pointsPO_n600 file
NR <= 3 { next }
{
  x = $1; y = $2; a = $7; b = $8; if (a < 0) a = -a; if (b < 0) b = -b
  if (x != px) { if (NR > 4) printf "%.6g %.6e %.6e\n", px, Sa, Sb; Sa = 0; Sb = 0; py = ""; }
  else if (y <= 150) { Sa += 0.5*(a+pa)*(y-py); Sb += 0.5*(b+pb)*(y-py) }
  px = x; py = y; pa = a; pb = b
}
END { printf "%.6g %.6e %.6e\n", px, Sa, Sb }
