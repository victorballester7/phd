# derivation of the formula in gmsh

## Formula for quad elements

Consider a 1D line from $x_0$ to $x_f$ of length $L$, that is $x_f = x_0 + L$. We want to split the segment in $n$ elements using a geometric progression of ratio $p$. We are given this $p$ and let $dx$ be the distance between the last two elements. We want to compute an approximate number of the number of elements $N$ that we need.

Let $x_1$ be the first element in the domain. Thus:

$$
x_0 + (x_1 - x_0) * (1 + p + p^2 + \ldots + p^n) = x_f
$$

Which implies:

$$
(x_1 - x_0) * (1 + p + p^2 + \ldots + p^n) = (x_1 - x_0) (1 - p^{n+1}) / (1 - p) = L
$$

Moreover, we know that $(x_1 - x_0)p^n = dx$. Thus:

$$
dx / p^n = L * (1 - p) / (1 - p^{n+1})
$$

Let $z = p^n$ and $a = dx/L$. Then:

$$
d * (1 - p * z) = (1 - p) * z \implies
z = a / (1 - p + a * p)
$$

Finally:

$$
n = \log (z) / \log (p) = \log (a / (1 - p + a * p)) / \log (p)
$$

## Formula for triangular elements

Here, in order to simplify the notation, '_q' will mean that the quantity is related to the quads, and '_t' to the triangles.

To determine the first $dx_t$ based on the one from the quads (the largest one on a quads sequence), we just need $x_1_q - x_0_q$, which we have already seen that is equal to $x_1_q - x_0_q = L_q * (1 - p_q) / (1 - p_q ^ {n_q + 1})$.
Now:
<!-- TO BE REVIEWED -->
 <!-- - If we want to compute $n$ based on the greatest distance, as opposed to the last one (where we were based on the smallest distance between elements). So, since in this case we have $x_1_q - x_0_q = dx_t = x_1_t - x_0_t$, we have:  -->

 <!--    $$ -->
 <!--    dx_t = L_t * (1 - p_t) / (1 - p_t ^ {n_t + 1}) \implies -->
 <!--    n_t = \log (1 - (1 - p_t) / a_t) / \log (p_t) - 1 -->
 <!--    $$ -->

 - If we want to compute $n$ based on the smallest distance, we have that $dx_t = (x_1_q - x_0_q)$ is known and so we can apply the previous section to conclude:
 $$
   n_t = \log (a_t / (1 - p_t + a_t * p_t)) / \log (p_t)
 $$
 with $a_t = dx_t / L_t = L_q / L_t * (1 - p_q) / (1 - p_q ^ {n_q + 1})$.
