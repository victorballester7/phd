# zeta(x) = (int_0^150 <w'^2>_z dy)^(1/2) and eta(x) = (int <|u'|^2>_z dy)^(1/2) from a
# pointsPO_n600_chk<n>.dat snapshot (columns x y u v w p uvw ww, x-major, y increasing).
# Prints x_b - w (first station past x - w = 25 with zeta >= 1e-2), eta there, and
# zeta at x - w = 40, 60, 80.
FNR <= 3 { next }
{
  x = $1; y = $2
  if (x != xc) { if (nx) { Z[nx] = sqrt(sz); E[nx] = sqrt(se) } nx++; X[nx] = x; xc = x; sz = 0; se = 0; py = "" }
  if (y <= 150) {
    if (py != "") { sz += 0.5 * (y - py) * (($8 < 0 ? -$8 : $8) + pw); se += 0.5 * (y - py) * (($7 < 0 ? -$7 : $7) + pu) }
    py = y; pw = ($8 < 0 ? -$8 : $8); pu = ($7 < 0 ? -$7 : $7)
  }
}
END {
  Z[nx] = sqrt(sz); E[nx] = sqrt(se)
  xb = "nan"; eb = "nan"
  for (i = 1; i <= nx; i++) if (X[i] - W >= 25 && Z[i] >= 1e-2) { xb = X[i] - W; eb = E[i]; break }
  for (i = 1; i <= nx; i++) { d = X[i] - W; if (d >= 40 && !z40) z40 = Z[i]; if (d >= 60 && !z60) z60 = Z[i]; if (d >= 80 && !z80) z80 = Z[i] }
  printf "%s %.3f %.2e %.2e %.2e\n", xb, eb, z40, z60, z80
}
