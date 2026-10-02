---

Part 1 – Nonlinear incompressible Navier–Stokes (IncNavierStokesSolver)

1.0 Entry point and object creation

1. solvers/IncNavierStokesSolver/IncNavierStokesSolver.cpp builds a Driver (Standard, ModifiedArnoldi, Arpack, SteadyState, …).
2. Driver::v_InitObject (library/SolverUtils/Driver.cpp:83) reads EvolutionOperator (default Nonlinear). It sets the session tag AdvectiveType (Nonlinear→Convective, Direct→Linearised, Adjoint→Adjoint, TransientGrowth→ both, SkewSymmetric, AdaptiveSFD). It then creates the EquationSystem named by SolverType, usually VelocityCorrectionScheme.
3. IncNavierStokes::v_InitObject (EquationSystems/IncNavierStokes.cpp:89):
   - maps u,v,w to m_velocity;
   - reads EqType and Kinvis (ν);
   - chooses the default advection form (:129–146), which the AdvectiveType tag can override (:151);
   - creates m_advObject from the Advection factory (:166);
   - loads m_forcing, and sets up Robin/radiation and Womersley BCs.
4. VelocityCorrectionScheme::v_InitObject (EquationSystems/VelocityCorrectionScheme.cpp:73):
   - pressure is the last field, and only the velocities are time-integrated (m_intVariables);
   - SetUpExtrapolation() (:152) creates the high-order pressure BC object (StandardExtrapolate, SubStepping, …);
   - SetUpSVV() and GJP stabilisation;
   - registers the operators with the time integrator: DefineOdeRhs → EvaluateAdvection_SetPressureBCs (explicit) and DefineImplicitSolve → SolveUnsteadyStokesSystem (implicit).

1.1 Continuous problem and discretisation

∂u/∂t = N(u) − ∇p + ν∇²u + f, with ∇·u = 0 and N(u) = −(u·∇)u. The spatial discretisation is C0 spectral/hp. Time integration uses a stiffly-stable IMEX scheme (IMEXOrder1..4); StandardExtrapolate::v_SubSteppingTimeIntegration rejects anything that isn't IMEX. The splitting is the Karniadakis–Israeli–Orszag velocity-correction scheme.

1.2 One time step, in the order the code runs it

UnsteadySystem::v_DoSolve (library/SolverUtils/UnsteadySystem.cpp:206) calls m_intScheme->TimeIntegrate. For each step the IMEX GLM does the following.

Step 1 – Explicit part: VelocityCorrectionScheme::v_EvaluateAdvection_SetPressureBCs (:758)
- a. Advection: IncNavierStokes::EvaluateAdvectionTerms (IncNavierStokes.cpp:317) copies the velocity, lets forcings PreApply (for example, a moving frame), then calls m_advObject->Advect.
  - NavierStokesAdvection::v_Advect (AdvectionTerms/NavierStokesAdvection.cpp:80) computes PhysDeriv of each component, then outarray[n] = u_j ∂u_n/∂x_j, then negates it.
  - Optional 3/2-rule dealiasing (SPECTRALHPDEALIASING) uses PhysInterp1DScaled and PhysGalerkinProjection1DScaled. Homogeneous directions use HomogeneousBwdTrans.
  - Other forms: SkewSymmetricAdvection, AlternateSkewAdvection, NoAdvection.
- b. Optional SmoothAdvection.
- c. Forcing (explicit): x->Apply(...) for every m_forcing (body force, sponge, moving body, …).
- d. High-order pressure BCs: m_extrapolation->EvaluatePressureBCs → StandardExtrapolate::v_EvaluatePressureBCs (StandardExtrapolate.cpp:69):
  - Extrapolate::v_CalcNeumannPressureBCs (Extrapolate.cpp:114) computes Q = N(uⁿ) − ν∇×∇×uⁿ on boundary elements (CurlCurl at :164, MountHOPBCs at :169), then NormVectorIProductWRTBase gives n·Q.
  - ExtrapolateArray (:969) extrapolates to tⁿ⁺¹ with the stiffly-stable β coefficients.
  - AddDuDt (:72) subtracts n·∂u_bc/∂t using AccelerationBDF.
  - CopyPressureHBCsToPbndExp writes the result into the pressure Neumann BCs. CalcOutflowBCs handles the convective and outflow BCs.
  - m_IncNavierStokesBCs->Update(...) handles the special BC classes (StaticWall, TransMovingWall, MRF*).

Step 2 – The GLM forms the explicit combination û = Σα_q u^{n−q} + Δt Σβ_q N^{n−q} and calls the implicit solve with aii_Dt (= Δt/γ₀).

Step 3 – Implicit part: v_SolveUnsteadyStokesSystem (:799)
- a. SubStepSetPressureBCs, used only with substepping.
- b. Pressure forcing: v_SetUpPressureForcing (:884) computes F = (∇·û)/aii_Dt.
- c. Pressure Poisson: v_SolvePressure (:949) calls HelmSolve with λ=0 to solve ∇²p = F with the HOPBCs, then AddPressureToOutflowBCs.
- d. SolveSolid → UpdateVelocityBCs (time-dependent velocity BCs).
- e. Viscous forcing: v_SetUpViscousForcing (:907) computes F_i = (∂p/∂x_i − û_i/aii_Dt)/ν.
- f. Helmholtz for each component: v_SolveViscous (:966) solves ∇²u_i − λu_i = F_i with λ = 1/(ν·aii_Dt). SVV factors come from AppendSVVFactors and the GJP jump term from ComputeGJPNormalVelocity. BwdTrans gives uⁿ⁺¹. This is equivalent to u/Δt′ − ν∇²u = û/Δt′ − ∇p.
- g. Optional flow-rate correction (MeasureFlowrate, which adds α·u_Stokes).

Step 4 – Back in v_DoSolve: phys is copied into the fields (VCS returns v_RequireFwdTrans()=false because HelmSolve already wrote the coefficients), then filters run and checkpoints are written.

1.3 Where each term is computed

┌──────────┬────────────────────────────────────────────────────────────────────┐
│   Term   │                              Function                              │
├──────────┼────────────────────────────────────────────────────────────────────┤
│ ∂u/∂t    │ IMEX GLM (γ₀, α_q) + aii_Dt in Step 3                              │
├──────────┼────────────────────────────────────────────────────────────────────┤
│ −(u·∇)u  │ NavierStokesAdvection::v_Advect                                    │
├──────────┼────────────────────────────────────────────────────────────────────┤
│ f        │ Forcing::Apply inside v_EvaluateAdvection_SetPressureBCs           │
├──────────┼────────────────────────────────────────────────────────────────────┤
│ −∇p      │ PhysDeriv of p in v_SetUpViscousForcing                            │
├──────────┼────────────────────────────────────────────────────────────────────┤
│ ∇·u = 0  │ Pressure Poisson (v_SetUpPressureForcing, v_SolvePressure) + HOPBC │
├──────────┼────────────────────────────────────────────────────────────────────┤
│ ν∇²u     │ HelmSolve in v_SolveViscous                                        │
├──────────┼────────────────────────────────────────────────────────────────────┤
│ ∂p/∂n BC │ Extrapolate::v_CalcNeumannPressureBCs, ExtrapolateArray, AddDuDt   │
└──────────┴────────────────────────────────────────────────────────────────────┘

---

Part 2 – Linearised (and adjoint) incompressible solver

2.1 How it is switched on

Set EvolutionOperator=Direct. Driver::v_InitObject sets AdvectiveType=Linearised, which makes IncNavierStokes::v_InitObject create LinearisedAdvection. Alternatively use EqType=UnsteadyLinearisedNS (default "Linearised", IncNavierStokes.cpp:142). Nothing else in the VCS changes.

2.2 Base flow: LinearisedAdvection::v_InitObject (AdvectionTerms/LinearisedAdvection.cpp:58)

- BaseFlow is mandatory. It can be a file (ImportFldBase, :410, which handles HalfMode/SingleMode planes) or analytic.
- Time-periodic base flows (Floquet): set N_slices>1. DFT (:661) takes an FFT of the slices, and UpdateBase (:499) reconstructs Ū(t) by Fourier series, or by Lagrange interpolation if BaseFlow_interporder≥2.
- UpdateGradBase (:570) precomputes ∂U_i/∂x_j once for a steady base flow, or at every call for a periodic one.

2.3 The linearised term: LinearisedAdvection::v_Advect (:264)

Per component i:
1. PhysDeriv(u′_i) (brought back to physical Fourier space for MultipleModes).
2. U_j ∂u′_i/∂x_j (:344).
3. + u′_j ∂U_i/∂x_j from m_gradBase (:352). HalfMode/SingleMode reduce the sum because W=0 and ∂/∂z=0 for the base.
4. Negate (:369), so the result is N_L(u′) = −(U·∇u′ + u′·∇U).

Everything downstream is the same linear machinery: forcing, then HOPBCs with N ← N_L(u′) (so ∂p′/∂n is computed from the linearised advection), then pressure Poisson, then Helmholtz. Perturbation BCs are not generated automatically. You write homogeneous BCs (u′=0 on walls and inflow, p′=0 on the outflow, …) in the session file.

2.4 Adjoint: AdjointAdvection::v_Advect (AdjointAdvection.cpp:51)

It returns +U·∇u* − (∇U)ᵀu* (plus scalar-coupling terms), with the base flow evaluated at m_period − time. This is the adjoint LNS written in reversed time τ = T − t, so the same forward integrator can be used.

2.5 Eigenvalue drivers

- DriverModifiedArnoldi (library/SolverUtils/DriverModifiedArnoldi.cpp, time-stepper approach):
  - EV_evolve (:461): CopyArnoldiArrayToField → TransCoeffToPhys → SetTime(0) → DoSolve (T = NumSteps·TimeStep) → CopyFieldToArnoldiArray. This applies the monodromy operator Φ_T = e^{LT} (or the adjoint-of-direct product for TransientGrowth, via CopyFwdToAdj).
  - v_Execute (:159) builds the Krylov sequence. EV_small (:492) does modified Gram–Schmidt, builds the Hessenberg matrix and calls dgeev. EV_test checks convergence and EV_post writes the eigenvectors.
  - DriverArnoldi::WriteEvs converts each eigenvalue μ to σ = ln|μ|/T and ω = arg μ/T.
  - Optional masks come from MaskInit.
- DriverArpack does the same with ARPACK.
- Steady coupled solver: CoupledLinearNS with EqType=SteadyLinearisedNS assembles the full (u,p) matrix with static condensation. SetUpCoupledMatrix (CoupledLinearNS.cpp:402, IsLinearNSEquation blocks at :755/:961/:1160) includes the U·∇u′ + u′·∇U terms and the shift λ. It is used in shift-invert mode (v_NegatedOp()=true, m_timeSteppingAlgorithm=false). The same matrix is the Newton Jacobian in SolveSteadyNavierStokes.
- Base flows come from DriverSteadyState (selective frequency damping, SFD) or from Newton in CoupledLinearNS.

---

Part 3 – Nonlinear compressible solver (CompressibleFlowSolver)

3.0 Classes and state

- State: q = (ρ, ρu, ρv, ρw, E). There is no pressure field.
- Discretisation: DG only (InitAdvection asserts eDiscontinuous).
- Class hierarchy: CompressibleFlowSystem is the base, with EulerCFE and NavierStokesCFE for explicit runs, and CFSImplicit → EulerImplicitCFE / NavierStokesImplicitCFE for implicit runs.

CompressibleFlowSystem::v_InitObject (EquationSystems/CompressibleFlowSystem.cpp:53) sets up:
- VariableConverter (equation of state: ideal gas, van der Waals, Peng–Robinson, …);
- InitAdvection (:164), which creates AdvectionWeakDG, GetFluxVector and the Riemann solver from UpwindType;
- artificial diffusion and forcing (sponge, …);
- one CFSBndCond object per user-defined BC;
- DefineOdeRhs(DoOdeRhs) and DefineProjection(DoOdeProjection).

NavierStokesCFE::InitObject_Explicit (NavierStokesCFE.cpp:64) adds Cp, Cv, μ (constant or Sutherland), κ or Pr, Twall, and DiffusionLDGNS or interior penalty.

3.1 One right-hand-side evaluation: DoOdeRhs (:214)

Step 1 – Trace states. GetFwdBwdTracePhys gives Fwd/Bwd on every trace. On a physical boundary, Bwd is the ghost state stored in the boundary expansion by the last BC application. Then PeriodicBwdRot.

Step 2 – Advection: DoAdvection (:360) → AdvectionWeakDG::AdvectCoeffs (library/SolverUtils/Advection/AdvectionWeakDG.cpp:117)
- Volume flux GetFluxVector (:452): F_j = (ρu_j, ρu_iu_j + pδ_ij, (E+p)u_j), with p from m_varConv->GetPressure. GetFluxVectorDeAlias is the dealiased variant.
- IProductWRTDerivBase computes (∇φ, F).
- Trace flux AdvectTraceFlux (:195) → RiemannSolver::Solve (rotates to the normal frame, RiemannSolver.cpp:100–131) → CompressibleSolver::v_Solve → v_PointSolve. Available solvers: Roe (RoeSolver.cpp:84), HLLC, AUSM, LaxFriedrichs, Exact, …
- AddTraceIntegral adds ⟨φ, F̂·n⟩, then MultiplyByElmtInvMass and BwdTrans.

Step 3 – Negate, so out = −∇·F.

Step 4 – Diffusion: NavierStokesCFE::v_DoDiffusion (:147)
- Convert to primitive variables [u,v,w,T] in the volume and on Fwd/Bwd (GetVelocityVector, GetTemperature).
- DiffusionLDGNS::v_DiffuseCoeffs (Diffusion/DiffusionLDGNS.cpp:229):
  a. DiffuseCalcDerivative (:315): LDG gradient. NumericalFluxO1 (:392) takes the Fwd value; ApplyBCsO1 (:430) handles wall velocity, T_wall or adiabatic walls, symmetry, and T from the equation of state on Dirichlet boundaries. The result is q = ∇(u,T).
  b. DiffuseVolumeFlux → v_GetViscousFluxVector (NavierStokesCFE.cpp:276): τ = μ(∇u+∇uᵀ) − ⅔μ(∇·u)I, energy flux u·τ + κ∇T, with μ and κ from GetViscosityAndThermalCondFromTemp (:672).
  c. DiffuseTraceFlux → NumericalFluxO2 (:691): the Bwd (downwind) q·n plus the penalty (1/h)·v_GetFluxPenalty (:630) = {μ}[u] and {κ}[T]; then ApplyBCsO2.
  d. IProductWRTDerivBase, AddTraceIntegral, MultiplyByElmtInvMass, and the result is added to out.
- Optional ArtificialDiffusion or physical shock-capturing artificial viscosity.

Step 5 – Forcing x->Apply (sponge / Absorption, axisymmetric, quasi-1D, synthetic eddies).

Step 6 – Local time stepping scaling, if enabled.

3.2 Projection and BCs: DoOdeProjection (:313)

Copy, then optional ExponentialFilter, then SetBoundaryConditions (:396). This calls ExtractTracePhys and then CFSBndCond::Apply for each BC (WallBC, RiemannInvariantBC::v_Apply at :67, PressureOutflow*, StagnationInflow, …), which writes the ghost states into the boundary expansions used in Step 1. The projection is called at every RK stage.

3.3 Time integration

- Explicit (RK, …): UnsteadySystem::v_DoSolve. Δt comes from the CFL via v_GetTimeStep (:650) when CFL is set.
- Implicit (CFSImplicit, CompressibleFlowSystemImplicit.cpp, SDIRK/BDF): DoImplicitSolve (:377) → DoImplicitSolveCoeff → Newton–GMRES:
  - residual NonlinSysEvaluatorCoeff (:216): x − λR(x) − b, where DoOdeProjection is followed by DoOdeRhsCoeff;
  - Jacobian-free J·v MatrixMultiplyMatrixFreeCoeff (:524), with ε = JacobiFreeEps·√((√‖x‖²+1)/‖v‖²) and optional central differences (GMRESCentralDifference);
  - block-Jacobi preconditioner from analytic element Jacobians: GetFluxVectorJacPoint (:1764), PointFluxJacobianPoint (:1904, split A± flux Jacobian), the viscous GetdFlux_dU_2D/3D and GetdFlux_dQx/y/z in NavierStokesImplicitCFE.cpp, and a finite-difference trace Jacobian NumCalcRiemFluxJac (:1249) with per-variable scaling m_magnitdEstimat.

These Jacobians are exactly the building blocks a linearised compressible solver needs.

---

Part 4 – Linearising the compressible solver

4.1 Target equations (conservative variables, frozen base q̄)

∂q′/∂t + ∂_j[A_j(q̄)q′] = ∂_j[(∂F^v_j/∂q)(q̄,∇q̄) q′ + (∂F^v_j/∂∇q)(q̄) ∇q′] + (forcing)′

with A_j = ∂F_j/∂q. For an ideal gas, p′ = (γ−1)[E′ − ū·m′ + ½|ū|²ρ′]. The viscous part includes μ′ = (dμ/dT)T′ when viscosity is variable. The BCs are the linearised versions of the ghost-state maps in §3.2.

4.2 Three ways to do it

A. Jacobian-free linearisation of the time-stepper (your current diff).
- A·v ≈ [Φ_T(q̄+εv) − Φ_T(q̄)]/ε, where Φ_T is the nonlinear solver run for time T.
- Pros: no solver changes; every term (Riemann solver, LDG, BCs, sponge) is linearised consistently at the discrete level.
- Cons:
  - the ε error builds up over the whole of T through a nonlinear trajectory;
  - the result depends on the base flow being a discrete fixed point;
  - it can't give an adjoint, so no transient growth or resolvent.

B. Jacobian-free linearisation at the right-hand-side level (my recommendation as the next step).
- A new operator mode inside CompressibleFlowSystem, triggered by AdvectiveType=Linearised. The driver already sets that tag for Direct; today the compressible solver ignores it.
// init: load BaseFlow into m_baseFlow (phys), then
//   DoOdeProjection(m_baseFlow) and DoOdeRhs(m_baseFlow) -> m_baseRhs
void CompressibleFlowSystem::DoOdeRhsLinearised(in, out, time)
{
    eps = m_linEps * normW(m_baseFlow) / normW(in);    // weighted norm, scalar eps
    qp  = m_baseFlow + eps * in;
    DoOdeProjection(qp, qp, time);   // nonlinear BCs on the perturbed state
    DoOdeRhs(qp, out, time);         // R(q̄ + ε q')
    out = (out - m_baseRhs) / eps;   // or central: [R(q̄+εq') − R(q̄−εq')]/2ε
}
// the ODE projection for q' becomes a plain copy (BCs are inside the finite difference)
- Register it with m_ode.DefineOdeRhs / DefineProjection when linearised. Keep filtering and local time stepping off, and freeze any artificial viscosity at its base-flow value.
- Why it's more robust than A:
  - the finite-difference error is O(ε) for each right-hand-side call instead of building up over T;
  - R(q̄) ≠ 0 is removed exactly, so you get the frozen base-flow linearisation even if q̄ isn't perfectly converged;
  - periodic base flows can reuse the UpdateBase/DFT logic from LinearisedAdvection (worth moving into a shared SolverUtils helper);
  - the same wrapper in DoOdeRhsCoeff gives an implicit linear solver for free through CFSImplicit (Newton converges in about one iteration);
  - it plugs into the existing ModifiedArnoldi / Arpack drivers unchanged.
- Cons: still no adjoint; it costs one nonlinear right-hand-side per linear one (two with central differences).

C. Analytic linearised operator (the true compressible counterpart of LinearisedAdvection). Term by term, mapped onto existing code:
- Volume flux: A_j(q̄)q′, reusing GetFluxVectorJacPoint / GetFluxVectorJacDirElmt, then the same IProductWRTDerivBase path.
- Trace flux: a linear upwind flux F̂′ = A⁺ₙ(q̄)q′_L + A⁻ₙ(q̄)q′_R, from PointFluxJacobianPoint with fsw=±1. This is the Roe or Steger–Warming linearisation about a continuous base state.
- Viscous: convert primitive perturbations pointwise (u′ = (m′−ρ′ū)/ρ̄, T′ from E′). The LDG gradient step is already linear, so it can be reused as is. Only the flux function passed through SetFluxVectorNS needs a linearised version (τ′, u′·τ̄ + ū·τ′, κ̄∇T′ + κ′∇T̄), and GetdFlux_dU / GetdFlux_dQ already provide its Jacobians. The penalty term uses {μ(T̄)}[u′].
- BCs: linearised CFSBndCond variants: wall with reflected m′ and T′=0 or ∂T′/∂n=0; far field with zero incoming characteristic perturbations.
- Pros: exact, cheapest per step, and the only route to an adjoint (discrete transpose or continuous adjoint), and so to transient growth and resolvent analysis.
- Cons: the most code, and each BC needs its own linearisation.

Suggested path: make A work for validation, build B as the production linearised solver, and move to C when you need the adjoint.

4.3 Problems I found in your current DriverModifiedArnoldi diff for the compressible case

1. The assert fails for the compressible solver. m_timeSteppingAlgorithm is only true for SolverType=VelocityCorrectionScheme (DriverArnoldi.cpp:60), so JacobianFreeInit asserts for NavierStokesCFE. The same flag also sets m_period=1 and switches WriteEvs to shift-invert output. It needs a general "time-stepping" test.
2. m_nfields is size()−1 when time-stepping (to drop p). For the compressible solver that would drop E, so it has to depend on whether a pressure field exists.
3. Critical: the Arnoldi vector never reaches the compressible solver. CompressibleFlowSystem doesn't override v_TransCoeffToPhys, and the EquationSystem default is empty (EquationSystem.cpp:1185). EV_evolve writes the Arnoldi vector into the coefficients, but v_DoSolve integrates from UpdatePhys(), so it would start from the stale physical state. Add a BwdTrans override (and a FwdTransLocalElmt one for v_TransPhysToCoeff).
4. ε scaling. The global ‖q̄‖ is dominated by E in dimensional units (p∞=101325), and the Arnoldi dot product has the same bias. Non-dimensionalise (ρ∞=u∞=1, p∞=1/(γM²)), or use a weighted or Chu-energy norm. The implicit solver scales per variable (CalcRefValues) for the same reason.
5. The time step must be fixed (CFL = 0, and no IO_CheckTime truncation). Otherwise Δt, and even the number of steps, depends on the state, and Φ_T isn't smooth.
6. Things that aren't smooth break the finite difference: shock capturing, exponential filtering, the Roe entropy fix, and wave-speed switches in HLL/HLLC. Use smooth subsonic cases first.
7. Forward difference over a long T. Offer a central-difference option (two evolves per iteration); then m_baseEvolved isn't needed.
8. Incompressible only: the pressure-BC history (Extrapolate::m_pressureCalls, m_pressureHBCs) is never reset between DoSolve calls (it is only zeroed in GenerateHOPBCMap, Extrapolate.cpp:690), while the IMEX scheme restarts every time. The base-flow evolve in init therefore runs with a fresh history, and later evolves run with the previous vector's history. The finite difference then isn't a clean function of v. The existing Direct Arnoldi has the same leak, but dividing by ε makes it worse.

4.4 Plan and checks

1. Fix points 1–3 above, or put B in a new EquationSystem mode.
2. Base flow: SFD through DriverSteadyState already works with the compressible solver (see Tests/implicitSolverCallsSFD_session.xml), or converge it with the implicit solver. Check the "Base flow drift" line your code prints.
3. Non-dimensionalise the case.
4. Checks:
   - linearity: ‖A(2v) − 2Av‖ and ‖A(v+w) − Av − Aw‖ should be small;
   - ε sweep: the leading eigenvalue should stay flat over several decades of ε;
   - A against B on the same case.
5. Physics check: a low-Mach (M≈0.1) cylinder near the incompressible Hopf point (Re≈47) should reproduce the incompressible Direct Arnoldi result.

If you want, I can make this into a shareable page, or start implementing option B or the driver fixes.
