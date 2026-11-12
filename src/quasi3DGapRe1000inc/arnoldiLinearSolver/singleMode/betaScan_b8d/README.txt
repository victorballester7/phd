Single-mode Arnoldi beta scan with b = 8d (beta_k = 2 pi k/b, k = 1..31), Re = 1000, three 2D-stable/3D-unstable gaps:
  d4_w15, d3_w21, d1.75_w40 (base flows: incGapRe1000/baseflow/dns/<case>/mesh_{110,55,100}.chk, steady to rel. change < 2e-5,
  interpolated with FieldConvert interpfield onto the mesh here).
Mesh: quasi-3D dx = 0.2 family (singleMode/d4_w15/mesh.geo) with the outflow shortened 1000 -> 200 (~2000-2450 elements, half of
  the full domain); perturbation BCs are homogeneous. d4_w15_L1000_k16 = same run on the full-length mesh, to validate the cut
  (old full-domain reference, tau = 1750: sigma = 5.02e-3 at lambda_z = 0.5d, i.e. k = 16).
Arnoldi: tau = NumSteps*dt = 1e4*0.007 = 70 (pi/tau = 0.045 > largest 3D frequency seen, 0.044), kdim 16, nvec 2, evtol 1e-5,
  nits 500, no checkpoints (IO_CheckSteps = 0).
Submit: cd <case> && qsub betaScan.job   (array 1-31);   cd d4_w15_L1000_k16 && qsub pbspro.job
Created 2026-11-11.
