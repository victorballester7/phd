# 13/12/2024

## Integration and differentiation in Nektar

We start in 1D. The idea is always to use Lagrange interpolating polynomial, that is, writing $u(\xi) = \sum_{i=0}^{Q-1} u(\xi_i) L_i(\xi)$, with $L_i(\xi)$ being the Lagrange interpolating polynomial. Integrating and differentiation formulas follow from integrating and differentiating the Lagrange polynomial. Different nodes $\xi_i$ give different formulas. Here is where we find all the common names (Gauss-Legendre, Gauss-Radau-Legendre or Gauss-Lobatto-Legendre), depending on whether we want to include no endpoints in the interval, only one or both.
