# N factor of a gap flow from local stability theory

Code: `scripts/pp/src/pp/orrSommerfeld.py`, `scripts/pp/src/pp/localNfactor.py`,
driver `scripts/pp/scripts/predictNfactorLocal.py`,
benchmarks `scripts/pp/tests/testOrrSommerfeld.py`.

## The model

At every streamwise station the base flow profile `U(y)` is frozen and treated
as parallel. The spatial Orr-Sommerfeld problem is solved for a band of **real**
frequencies, giving `alpha(x; omega)`, and the perturbation is modelled as a
superposition of those modes,

```
u'(x, y, t) = sum_j a_j exp(N_j(x)) phi_j(x, y) exp(-i omega_j t) + c.c.
N_j(x)      = - int_{x0}^{x} alpha_i(xi; omega_j) dxi
```

with every eigenfunction `phi_j(x, .)` **normalised to unit L2 norm at its own
station**. That convention is what closes the model: the shape of the mode then
carries no amplitude information, all of it sits in `exp(N_j)`, and everything
local theory cannot know -- receptivity, how efficiently the forcing excites
each frequency, the non-parallel and non-modal corrections near the gap -- is
pushed into the constant weights `a_j`.

### The modes add in energy, not in amplitude

The DNS amplitude the prediction is compared against is

```
A(x)^2 = int <u'^2 + v'^2> dy
```

a **time averaged** quantity (`pp.DeltaN_computation.computeAmplitude`, which
reads the Reynolds stresses written by Nektar's `ReynoldsStresses` filter).
Modes of distinct real frequency are uncorrelated in time, so every cross term
averages to zero, and with `||phi_j|| = 1`

```
A(x)^2 = sum_j a_j^2 exp(2 N_j(x))
N(x)   = 1/2 log( A(x)^2 / A(x0)^2 )
```

`x0` is picked exactly as `computeNx` picks it -- the station of smallest
amplitude upstream of the gap -- so the two curves are directly comparable.

**The eigenfunctions drop out of the answer entirely.** Only `alpha_i` survives.
That is a consequence of the unit norm convention, and it is why the result does
not depend on the quadrature used for the norm.

### Weights

The DNS forcing is a localised divergence free Gaussian blob at `x = -95`
modulated by `awgn(1.0)`, i.e. **white in time**, so a flat `a_j` is the
zeroth order model and is the default. `fit_weights` instead solves for
`a_j^2 >= 0` by non negative least squares against the DNS amplitude; because
`A^2` is *linear* in `a_j^2` that fit has no local minima.

## Conventions

Everything is in the global non dimensionalisation of the DNS (lengths in
`delta*_le`, velocity in `U_inf`, `re` = `Re_delta*_le`). No rescaling of the
profile by the local `delta*` is done, so `alpha` and `omega` need no conversion
between stations and `N_j` is a plain integral in `x`. Frequencies are labelled
by the reduced frequency

```
F = omega * 1e6 / re
```

which is the invariant of a mode travelling downstream: `omega` and `Re` both
scale with the local `delta*`, so their ratio does not.

## Numerics, and the three things that actually mattered

**1. The quartic must be split.** Collecting powers of `alpha` in the
Orr-Sommerfeld equation gives a quartic polynomial eigenvalue problem whose
companion pencil mixes `D4` with the identity -- entries spanning `n^8`. On the
Poiseuille benchmark that lost eleven digits: the eigenvalue came out with an
error of `6e-4` where the temporal branch gave `1e-9`. Splitting off
`q = (D2 - alpha^2) v` (Bridges & Morris 1984) drops the highest derivative to
`D2`, and the *same sized* pencil then reproduces the temporal eigenvalue to
`1.4e-8`.

**2. The profile must be fitted, not interpolated.** The section files carry
about `1e-4` of noise in `u` (spectral element interpolation onto the sample
points). That is invisible in `u`, is ~1% of `dU/dy`, and leaves `U''` -- which
enters the operator directly, next to the critical layer -- meaningless. With a
spline interpolant the eigenvalue did not converge with `n` at all: at one gap
station it moved between `0.25 - 0.13i`, `0.14 + 0.06i` and "no mode" as `n`
went 61, 71, 81. Fitting a truncated Chebyshev series in the mapped coordinate
(degree 40, where the residual has reached the noise floor) restored clean
spectral convergence, `1e-8` from `n = 61` upward. Because the degree stays
below `n`, the collocation `D2` differentiates the fit exactly, so `U''` is
consistent with the discretisation of `v''` and costs nothing extra.

**3. Every station and frequency is solved independently.** A Newton
continuation marching one eigenvalue is ~100x cheaper, and it is what a first
implementation naturally does, but it assumes the branch it started on stays the
relevant one -- which is exactly what fails over the gap, where the shear layer
mode appears and overtakes the Tollmien-Schlichting one. A full eigensolve at
each (station, frequency) makes no such assumption; the stations are
embarrassingly parallel, so a case costs a few minutes.

The mode is picked by `localNfactor.select_mode`. The phase speed is the
discriminator: the discrete instability runs at `c_r = 0.2..0.5` while the
discretised continuous spectrum piles up at `c_r -> 1`. The window is
`0.15 < c_r < 0.75`, `|alpha_i| < 0.9 alpha_r`, `|alpha| < 3`, and the least
damped survivor is taken.

### Resolution

`n = 81` collocation nodes, free stream truncated at `y = 45`. Against `n = 101`
over a full case the synthesised `N(x)` agrees to four decimals; where the mode
is amplified `alpha_i` agrees to `6e-5`. The pairs that do differ are all
strongly damped high `F` modes, which contribute nothing to the sum.

## Validation

`tests/testOrrSommerfeld.py`, all against numbers this code did not produce:

| check | reference | obtained |
|---|---|---|
| Poiseuille, Re = 10000, alpha = 1 (Orszag 1971) | `c = 0.23752649 + 0.00373967i` | to `3e-9` |
| spatial branch fed the temporal eigenvalue | `alpha = 1` | to `1.4e-8` |
| Blasius branch I at `Re_delta* = 884.5`, vs the C++ solver (temporal + Gaster 2) | `omega = 0.061221`, `alpha = 0.183792` | `0.061052`, `0.183410` |
| Blasius critical Reynolds number | `Re_delta* = 519.4` | `519.5` |

The branch I difference is 0.3%, which is the expected error of the Gaster
transformation the C++ result relies on, not of the spatial solve.

## What the local theory says about the gap

For `d = 1.5`, `w = 30` at `Re = 1000` the tracked mode behaves as:

| region | `alpha` | `-alpha_i` |
|---|---|---|
| upstream, `x < 0` | `0.236 - 0.005i` | `+0.005` (TS) |
| inside the gap | `0.20 - 0.126i` | `+0.126` (shear layer) |
| just downstream | `0.21 + 0.004i` | damped |
| `x > 60` | `0.238 - 0.006i` | `+0.005` (TS again) |

The gap replaces the Tollmien-Schlichting mode by a cavity shear layer mode
about 25 times more amplified, over the gap width only. That jump is the whole
of `Delta N`.

The driver also reports `Delta N` with the gap interval removed from the
integral (`synthesise(..., freeze=(0, w))`, dotted in the figures). A parallel
flow analysis of a station *inside* a short gap is the weakest link in the
model -- the shear layer there is a few gap widths long, nowhere near the
slowly varying limit -- so the two numbers separate "the gap amplifies" from
"the local theory says the gap amplifies by this much".

## Results

Twelve cases, two Reynolds numbers, **flat weights -- nothing fitted**:

| Re | d | w | `dN` local | local, gap frozen | `dN` DNS |
|---|---|---|---|---|---|
| 800 | 1.50 | 70 | 6.522 | -1.030 | 5.848 |
| 800 | 1.75 | 33 | 3.537 | -0.678 | 3.146 |
| 800 | 2.00 | 41 | 4.303 | -0.866 | 3.875 |
| 800 | 3.50 | 21 | 2.091 | -0.450 | 1.926 |
| 1000 | 0.25 | 40 | 0.359 | -0.619 | 0.330 |
| 1000 | 0.50 | 40 | 1.495 | -0.716 | 1.246 |
| 1000 | 0.75 | 40 | 2.717 | -0.729 | 2.401 |
| 1000 | 1.00 | 40 | 3.572 | -0.774 | 3.314 |
| 1000 | 1.25 | 40 | 3.935 | -0.830 | 3.599 |
| 1000 | 1.50 | 30 | 3.107 | -0.658 | 2.758 |
| 1000 | 1.50 | 40 | 4.024 | -0.889 | 3.717 |
| 1000 | 1.75 | 40 | 4.043 | -0.911 | 3.918 |

Fitting a single line through the origin over all twelve:

```
dN_local = 1.098 * dN_DNS,   residual scatter 0.10,   correlation 0.9977
```

Two things are worth drawing out.

**One scalar covers everything.** The same slope holds at both Reynolds
numbers, for gap depths from 0.25 to 3.5 and widths from 21 to 70, over a range
of `dN` from 0.3 to 5.8, and after it the residual is 0.10 in `dN`. The local
model therefore has the mechanism; what it is missing is a constant efficiency,
not a trend. The natural reading is that a disturbance of finite streamwise
extent collects only about 91% of the shear layer's local growth while it
crosses the gap -- precisely the correction the parallel flow assumption cannot
contain.

**The gap interval is the whole of `dN`.** Freezing it leaves `dN` between
-0.45 and -1.03, i.e. slightly *below* the flat plate: downstream of the gap
the boundary layer is marginally fuller and marginally more stable, and every
bit of the excess N factor is produced in the few gap widths where the cavity
shear layer replaces the Tollmien-Schlichting mode.

### Where the prediction and the DNS part company

Comparing the two panels of `Nfactor_local_vs_dns_Re1000.pdf`: the jump across
the gap, the overshoot at `x ~ 40`, the dip at `x ~ 90` and the flat plate
curve itself are all reproduced closely. Downstream of `x ~ 200` the predicted
curves for the deeper gaps keep climbing while the DNS ones flatten and bunch
together near `N ~ 4.5`. The local theory has no mechanism for that: it carries
each frequency forward independently, so the envelope keeps following whichever
mode is still amplified. It is the region to be careful about if the model is
used far downstream; the `Delta N` window at `w + 80` is upstream of it.

## Outputs

* `data/localNfactor/modes_Re*_*.npz` -- cached `alpha(x, F)` per case
* `data/localNfactor/deltaN_local_vs_dns_Re*.dat` -- `d w dN_local dN_local_nogap dN_dns`
* `images/localNfactor/Nfactor_local_vs_dns_Re*.pdf` -- prediction and DNS side by side
* `images/localNfactor/deltaN_local_vs_dns.pdf` -- the twelve cases against the 1:1 line
* `images/localNfactor/modes_Re*_*.pdf` -- growth rate map in `(x, F)` and the per frequency `N_j(x)`
