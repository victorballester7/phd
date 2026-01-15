Great — let’s dig in properly. I’ll give a mathematically careful explanation in three parts:

1. the linearized problem and an energy identity that isolates the role of boundary terms;
2. the pressure Poisson / compatibility issue for incompressible flow;
3. an operator/functional-analytic view (how boundary conditions change the domain of the linear operator and therefore the spectrum, with a simple 1-D example to make the mechanism explicit).

I’ll try to be concise but rigorous. Tell me if you want any step expanded or written more formally (e.g. as a theorem + proof).

---

### 1) Linearized incompressible Navier–Stokes and the energy identity

Let $\Omega\subset\mathbb{R}^d$ be the fluid domain with outward unit normal $ \mathbf n$ on $\partial\Omega$. Let $ \mathbf U(\mathbf x)$ be a steady baseflow satisfying the steady NSE, and consider infinitesimal perturbations $\mathbf u'(\mathbf x,t), p'(\mathbf x,t)$ that satisfy the linearized incompressible Navier–Stokes equations (in primitive variables):

$$
\begin{aligned}
\partial_t \mathbf u' &= -(\mathbf U\cdot\nabla)\mathbf u' - (\mathbf u'\cdot\nabla)\mathbf U - \nabla p' + \nu \Delta \mathbf u',\\
\nabla\cdot\mathbf u' &= 0,
\end{aligned}
$$

with boundary conditions prescribed on $\partial\Omega$. For clarity denote the viscous stress operator by $ \mathcal{V}(\mathbf u') = \nu \Delta \mathbf u'$.

Take the $L^2$-inner product of the momentum equation with $\mathbf u'$ and integrate over $\Omega$. Using $\langle a,b\rangle=\int_\Omega a\cdot b\,d\Omega$ we get

$$
\frac{1}{2}\frac{d}{dt}\|\mathbf u'\|_{L^2}^2
= -\int_\Omega \mathbf u'_i \mathbf u'_j \partial_j U_i \,d\Omega
- \underbrace{\int_\Omega (\mathbf U\cdot\nabla)\left(\tfrac12|\mathbf u'|^2\right)\,d\Omega}_{\text{convective flux term}}
- \int_\Omega \mathbf u'\cdot\nabla p' \,d\Omega + \int_\Omega \mathbf u'\cdot \mathcal{V}(\mathbf u') \,d\Omega.
$$

Now integrate the divergence term by parts and the pressure/viscous terms by parts (using $\nabla\cdot\mathbf u'=0$ where needed). One obtains the exact balance with boundary integrals:

$$
\boxed{\;
\frac{1}{2}\frac{d}{dt}\|\mathbf u'\|_{L^2}^2
= -\int_\Omega \mathbf u'_i \mathbf u'_j \partial_j U_i \,d\Omega
- \nu\int_\Omega |\nabla \mathbf u'|^2 \,d\Omega
- \int_{\partial\Omega} \Big( \tfrac12(\mathbf U\cdot\mathbf n)\,|\mathbf u'|^2 + p'\,(\mathbf u'\cdot\mathbf n) - \nu\,\mathbf u'\cdot\partial_n\mathbf u' \Big)\,dS.
\;}
$$

Key point: **the whole effect of the outlet boundary condition on the perturbation energy enters via the boundary integral**. The sign and size of that boundary integral control whether energy is lost through the boundary (outgoing wave absorption), gained (incoming energy) or reflected back.

Compare two cases at an outflow portion $\Gamma_\text{out}$ where the baseflow has $ \mathbf U\cdot\mathbf n>0$ (flow exits domain).

* **Neumann / stress-free (or do-nothing)** type BC: typically $\partial_n \mathbf u' = 0$ or stress·n = 0 and no constraint on $\mathbf u'\cdot\mathbf n$. Then the boundary term includes $\tfrac12(\mathbf U\cdot\mathbf n)|\mathbf u'|^2$ which for positive $U\cdot n$ is an *outgoing positive flux* and thus represents energy leaving the domain (damping from viewpoint of interior). Physically outgoing perturbations can carry energy out.

* **Dirichlet $\mathbf u'|_{\Gamma_\text{out}} = 0$**: both $\mathbf u'$ and hence $|\mathbf u'|^2$ and $\mathbf u'\cdot\mathbf n$ vanish on $\Gamma_\text{out}$, so the boundary integral over $\Gamma_\text{out}$ is *zero* there. That removes the natural energy-leaving mechanism: outgoing convective flux $\tfrac12(\mathbf U\cdot\mathbf n)|\mathbf u'|^2$ is forced to zero at the boundary. Physically this behaves like a rigid constraint: perturbations cannot carry their energy out through that face — instead they must either vanish before reaching the face or be reflected and redistributed inside the domain.

Thus mathematically: Dirichlet at outflow eliminates the positive-definite outgoing flux term; it replaces an energetic boundary outflow with a hard kinematic constraint. This tends to produce reflected energy and change the global energy budget, which in turn shifts growth rates in the eigenvalue problem.

---

### 2) Pressure Poisson / solvability and compatibility

For incompressible flows the pressure is not an independent evolution variable but determined (instantaneously) from the divergence-free constraint. Take divergence of linearized momentum:

$$
\nabla\cdot(\partial_t\mathbf u') = \partial_t(\nabla\cdot\mathbf u') = 0
$$

so divergence of the RHS must vanish, giving the pressure Poisson equation (dropping $\partial_t$ when seeking normal-mode solutions $e^{\lambda t}$):

$$
\Delta p' = -\nabla\cdot\big((\mathbf U\cdot\nabla)\mathbf u' + (\mathbf u'\cdot\nabla)\mathbf U\big).
$$

To solve this Poisson problem you need boundary conditions for $p'$. Those are obtained from the normal component of the momentum equation on $\partial\Omega$:

$$
\mathbf n\cdot\big(-\nabla p' + \nu\Delta \mathbf u' - (\mathbf U\cdot\nabla)\mathbf u' - (\mathbf u'\cdot\nabla)\mathbf U\big)
= \partial_t(\mathbf n\cdot\mathbf u') = 0 \quad\text{(for steady/normal-mode analysis)}.
$$

If you impose $\mathbf u' = 0$ (Dirichlet) on the outflow, that enforces $\mathbf n\cdot\mathbf u'=0$ there and simplifies the normal momentum condition to a Neumann condition for $p'$. But if you instead impose a stress-free (Neumann) velocity BC, the pressure BC is different. In short:

* **Different velocity BCs imply different boundary conditions on $p'$** for the Poisson equation, which change solvability and the regularity of $p'$.
* If the BCs are inconsistent (for example specifying all velocity components plus an incompatible pressure BC), the Poisson problem can be ill-posed or lead to pressure singularities near the outlet.

So mathematically, Dirichlet outflow not only kills the outgoing energy flux but also changes the elliptic subproblem for pressure — this modifies the coupling between pressure and velocity and the structure of the linear operator.

---

### 3) Functional-analytic / spectral viewpoint

Let $L$ be the (linear) operator defined by the right-hand side of the linearized equations acting on divergence-free velocity fields, with domain $D(L)$ determined by the regularity and *boundary conditions*.

Concretely, think of a linear eigenproblem $\partial_t \mathbf u' = L\mathbf u'$ and normal-modes $\mathbf u'(\mathbf x,t)=\hat{\mathbf u}(\mathbf x)e^{\lambda t}$ leading to

$$
\lambda \hat{\mathbf u} = L\hat{\mathbf u},\qquad \hat{\mathbf u}\in D(L).
$$

Important general facts:

* The operator $L$ depends on the BCs because $D(L)$ (the set of admissible eigenfunctions) is defined by those BCs. Changing BCs changes $D(L)$ and so generally changes the spectrum $\sigma(L)$. This is a standard result in linear operator theory: an unbounded operator is determined by both its action and its domain.
* For transport-dominated problems (high Reynolds number, large convective terms), the continuous spectrum / convective spectrum associated with translation and radiation conditions is sensitive to whether you allow *outgoing wave solutions* that do not decay at the outlet. Imposing Dirichlet at the outlet kills those nondecaying outgoing modes and typically converts part of what would be a continuous spectrum into discrete eigenvalues that are strongly influenced by the boundary (and often localized near it).

A simple illustrative 1-D example (to see the reflection mechanism explicitly):

Consider the linear advection–diffusion eigenproblem on $x\in[0,L]$

$$
\lambda u = -c u_x + \nu u_{xx},\qquad c>0,
$$

with say $u(0)=0$ at inflow and at outflow either (A) Dirichlet $u(L)=0$, or (B) Neumann $u_x(L)=0$. Seek separated-mode solutions $u(x)=e^{\mu x}$. The characteristic equation is

$$
\nu\mu^2 - c\mu - \lambda = 0,
$$

so $\mu$ depends on $\lambda$. Applying the two boundary conditions leads to different transcendental equations for $\lambda$. Concretely, you end up with conditions like

$$
\text{(A) } e^{\mu L}=0 \quad\text{(impossible for real }\mu),\quad\text{so only discrete complex roots exist determined by the matching at }x=0,L,
$$

versus

$$
\text{(B) } \mu e^{\mu L}=0\quad\text{allowing different families of roots.}
$$

The precise algebra depends on the two roots $\mu_1(\lambda),\mu_2(\lambda)$; the point is: **the allowed modal shapes and hence the spectrum are different**. For the advection-dominated limit ($\nu\to0$), a physically natural outflow condition is the radiation condition $u_x=0$ is better than forcing $u(L)=0$: the latter forces incoming characteristics (from $x=L$) to be zero and reflects outgoing waves.

Although the 1-D case is easier to solve explicitly, the qualitative lesson holds in multi-D incompressible flows: Dirichlet at the outlet disallows outgoing waveforms (modes that would have been nonzero at the boundary) and therefore forces the eigenfunctions to "fit" imposed zeros — this generates artificial standing-wave–like modes between reflecting boundaries and can substantially alter eigenvalues.

---

### 4) Adjoint / sensitivity point of view

If you compute adjoint eigenmodes (for sensitivity, wavemaker, optimal forcing), the adjoint operator $L^*$ has a domain that depends on the primal BCs. Integration-by-parts formulae show boundary terms that must vanish for $L^*$ to be the formal adjoint. Changing primal BCs (e.g. to Dirichlet at outflow) changes the adjoint BCs; consequently, **eigenvalue sensitivities and receptivity are different**. Practically, that means boundary-induced reflection is not only changing eigenvalues but also changing which regions in the flow are most receptive to forcing.

---

### 5) Consequences stated as rigorous claims (informal theorem-like)

1. *Domain dependence:* Let $L_{D}$ and $L_{N}$ denote the linearized Navier–Stokes operators with identical interior differential expressions but with Dirichlet and Neumann/outflow BCs on $\Gamma_{\text{out}}$, respectively. Then $D(L_D)\neq D(L_N)$ and in general $\sigma(L_D)\neq\sigma(L_N)$. (This follows because an unbounded operator is determined by its action and domain; change of domain changes spectrum.)

2. *Energy flux statement:* For solutions of the linear problem,

$$
\frac{1}{2}\frac{d}{dt}\|\mathbf u'\|_{L^2}^2 + \nu\|\nabla\mathbf u'\|_{L^2}^2
= -\int_\Omega \mathbf u'_i\mathbf u'_j\partial_j U_i \,d\Omega
- \int_{\partial\Omega} \Big( \tfrac12(\mathbf U\cdot\mathbf n)|\mathbf u'|^2 + \cdots\Big)\,dS.
$$

Hence any BC that enforces $\mathbf u'|_{\Gamma_{\text{out}}}=0$ removes the outgoing energy flux term $\tfrac12(\mathbf U\cdot\mathbf n)|\mathbf u'|^2$ on $\Gamma_{\text{out}}$. Thus Dirichlet outflow eliminates boundary energy dissipation via convection, and may convert outgoing-transport damping into internal reflection.

3. *Pressure compatibility:* The pressure Poisson equation and its Neumann boundary data depend on the prescribed velocity BCs. Thus different velocity BCs produce different elliptic problems for $p'$ and different solvability/regularity properties, which feed back into the spectral problem.

Those are mathematically precise consequences (they can be formalized into a full theorem with hypotheses about regularity and baseflow smoothness if you want).

---

### 6) Practical mathematical signatures you can check (rigorous diagnostics)

* If $\hat{\mathbf u}$ is an eigenfunction that does not vanish in a neighborhood of $\Gamma_\text{out}$ but the Dirichlet BC forces $\hat{\mathbf u}=0$ on $\Gamma_\text{out}$, then the eigenfunction must have strong gradients near the boundary; this produces high viscous dissipation $\nu\int|\nabla\hat{\mathbf u}|^2$ and tends to move eigenvalues left (more damping) as $\nu>0$. You can make this precise by bounds on the Rayleigh quotient for the operator.

* Compute the energy flux through a cross-section near the outlet: if the flux is non-negligible but the BC forces zero on the boundary, then reflection must occur. Quantitatively compare $\int_{S(x)} \tfrac12(\mathbf U\cdot\mathbf n)|\mathbf u'|^2$ as $x$ approaches the outlet.

* Domain-extension limit: prove (or numerically test) that as you increase the downstream length $L\to\infty$, the spectrum for the problem with Dirichlet at $x=L$ converges to that of the problem with an outgoing radiation condition on the infinite domain, provided the modes decay sufficiently fast. If convergence is slow or fails, boundary is influencing spectrum.

---

### 7) Summary (short, maths-focused)

* Dirichlet at the outflow changes the domain $D(L)$ of the linear operator; therefore spectrum and eigenfunctions change.
* Energy identity shows Dirichlet removes the outgoing convective energy flux term at the outlet and thus prevents natural energy leaving the domain — mathematically this manifests as a boundary integral that vanishes under Dirichlet, producing internal reflection or enhanced gradients.
* For incompressible flow the pressure Poisson equation and its boundary conditions depend on the chosen velocity BC; imposing Dirichlet at outflow modifies the elliptic subproblem and can create compatibility/regularity issues.
* In operator-theory language, physically transparent BCs (radiation/convective or sponge) yield operators whose spectrum is closer to the open-domain problem; Dirichlet at the outlet typically produces spurious boundary-localized eigenmodes and spectrum shifts.

---

If you want, I can:

* write a formal lemma+proof about the energy identity and the sign of the boundary flux term; or
* derive the pressure Poisson BC explicitly for a rectangular outflow face; or
* work through the 1-D advection–diffusion eigenvalue problem fully (explicit characteristic roots and boundary-matching condition) so you see algebraically how the eigenvalues shift between Dirichlet and Neumann.

Which of those would you like next?

