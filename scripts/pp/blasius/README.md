# Blasius profile (2D)

This solves the boundary value problem to get a blasius profile for an 2d inflow profile.

The equations that define the problem are:

$$
\begin{align}
f''' + 0.5 * ff'' &= 0 \\
f(0) &= 0\\
f'(0) &= 0\\
f'(\infty) &= 1\\
\end{align}
$$

## Usage

The code is written in Python and uses the `scipy` library to solve the boundary value problem. Run it with:

```python
python solve_bvp.py
```

It will output the coefficients of a polynomial approximation of order 8 (default), that then you can use in nektar session files.

