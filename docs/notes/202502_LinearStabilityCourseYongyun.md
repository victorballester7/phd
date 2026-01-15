# Notes on linear stability course

## Squire's theorem. 

Consider a 3D parallel viscous flow of the form $(U(y),0,0)$. Given $Re_L$ as the critical Reynolds number for the onset of linear instability for a given wavenumbers $\alpha$ (in $x$ direction) and $\beta$ (in $z$ direction), the Reynolds number $\Re_c$ below which no exponential instabilities exist for any wave number satisfies

$$
\Re_c = \min_{\alpha,\beta} \Re_L(\alpha,\beta) = \min_{\alpha} \Re_L(\alpha,0)
$$

_Remark_: The most unstable linear instability is always 2D.

_Remark_: Viscosity can destabilize the flow (TS waves). For example, in parallel inviscid flows, Rayleigh inflection theorem tells as that if there exists a perturbation that generates linear instability then, $d^2U/dy^2$ must change sign (i.e. must have a zero). However, in the presence of viscosity, the flow can be linearly unstable even if $d^2U/dy^2$ does not change sign (example in Poiseuille flow, with $U(y) = 1-y^2$).

_Remark_: Usually the modes $\alpha$ and $\beta$ that generate the instability are low, because as they tend to infinity, the effect in the viscous term becomes more and more important, so instabilities are wiped out (think in Fourier transform).


## Disturbance growth in linearly stable flows (non-modal analysis)

 - The problem with linear stability analysis is that it only examines the evolution of the most unstable mode. But it *fails* to capture the linear interactions between eigenmodes.

 - There can be a short term growth in energy due to non-orthogonality of two eigenvectors, even for linearly stable flows. To see it more clear, consider the matrix 
 [ -1/Re   0  ]
 [   1   -2/Re]
 which has eigenvalues $-1/Re$ and $-2/Re$ and eigenvectors [1 Re] and [0 1]. 
   * In the limit $Re \to 0$, the eigenvalues are orthogonal, leading to a monotonic decay of the initial condition. 
   * In the limit $Re \to \infty$, the eigenvalues tend to the same one, yielding to an algebraically growing solution.

- Eigenvalues alone thus fail to capture transient effects.
- If there is a transient amplification, it is due to the non-orthogonality of the eigenvectors. This non-orthogonal eigenfunctions are the typical nature of the non-normal linear operator (due to the advection term). 
- *Convective instabilities are associated with this linear transient growth*.

_Definition_: A normal operator $L$ satisfies $L L^* = L^* L$, where $L^*$ is the complex-conjugate transpose of $L$. Linearized NS with non-zero advection term is a non-normal linear operator.


## Relation between absolute and global stability

- Intuitively in parallel flows, fixed x, we take the U profile in y and we seek for the saddle point that gives us the maximum growth rate. Depending where this point is located (hemisphere north in complex plane or south) we have convective or absolute instability. Now we vary x. If the region we have absolute instability is large enough, then we have global instability.


## Transient growth and Local convective instability

- Since we are simulating a finite domain, the total growth of the energy will increase within a short interval of time until the vortex governing the instability is advected out of the domain. This is the transient growth.
