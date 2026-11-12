Checks of the b=8d beta scan (2026-11-12): are the stationary (real) d4_w15 modes and the scan values resolved?
Production: NUMMODES 6 (p 5), dt 0.007, tau 70, base flow interpolated from the P=8 DNS onto the short mesh at NUMMODES 4 (P=3!).
  d4_k15_bfP6   production, base flow re-interpolated at NUMMODES 6 (same as perturbation)
  d4_k15_tau98  production (P4 base flow), tau = 98 instead of 70 (aliasing check of the omega = 0 mode)
  d4_k15_P8     NUMMODES 8 (p 7), dt 0.0035, base flow at NUMMODES 8 (resolution + time step)
  d3_k07_P8     same for the d3_w21 oscillatory maximum
Base flows: FieldConvert -m interpfield:fromxml=<dns>/mesh.xml:fromfld=<dns>/mesh_<chk>.chk mesh_P<N>.xml bf.fld
