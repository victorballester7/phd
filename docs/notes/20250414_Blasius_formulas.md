# Blasius profile formulas.

A deduction of all the formulas used throughout the codes is given below.

## Blasius profile
Let's define a self similar variable $\eta$ as:

$$
\eta = \frac{y}{\sqrt{\nu x / U_\infty}} = \frac{y}{x}\sqrt{Re_x}
$$

Where $y$ is the distance from the wall, $x$ is the distance from the leading edge, $\nu$ is the kinematic viscosity, $U_\infty$ is the free stream velocity and $Re_x$ is the Reynolds number based on the distance from the leading edge of the plate. Then one defines $f(\eta)$ as:
$$
f(\eta) = \frac{\psi}{\nu \sqrt{Re_x}}
$$
Where $\psi$ is the stream function. The Blasius profile is the solution of the following third order ordinary differential equation:
$$
2f''' + ff'' = 0
$$
with boundary conditions $f(0) = 0$, $f'(0) = 0$ and $f'(\infty) = 1$. From here we can express $u(x,y)$ and $v(x,y)$ as:
$$
u(x,y) = U_\infty f'(\eta)
$$
$$
v(x,y) = \frac{1}{2} U_\infty \frac{1}{\sqrt{Re_x}} (\eta f'(\eta) - f(\eta))
$$

## $\delta^*$ for the Blasius profile
One can check that:
$$
\delta^* = \int_0^\infty \left( 1 - f'(\eta) \right) Y'(\eta) d\eta
$$
Where $y=Y(\eta)=\eta \frac{x}{\sqrt{Re_x}}$. Expanding the integral gives:
$$
\delta^* = \frac{x}{\sqrt{Re_x}} \int_0^\infty \left( 1 - f'(\eta) \right) d\eta=\frac{x}{\sqrt{Re_x}} \left( \eta - f(\eta) \right)\bigg|_0^\infty=C \frac{x}{\sqrt{Re_x}}
$$
with $C\simeq 1.72078$.

## Relations between Reynolds and $\delta^*$'s
We have
$$
\frac{C}{\sqrt{Re_x}}=\frac{\delta^*(x)}{x} = \frac{\frac{U_\infty \delta^*(x)}{\nu}}{\frac{U_\infty x}{\nu}} = \frac{Re_{\delta^*(x)}}{Re_x}
$$
which implies
$$
\sqrt{Re_x} = \frac{1}{C} Re_{\delta^*(x)}
$$
Moreover we have the clear equality:
$$
  Re_x=\frac{x}{\delta^*(x)} Re_{\delta^*(x)}
$$
Using the previous two relation we get:
$$
  \frac{x}{\delta^*(x)} = \frac{1}{C^2} Re_{\delta^*(x)}
$$

Now, we would like to find and expression of $\delta^*(x)$ in terms of $\delta^*(x_0)$ (reference length). Assume $x=x_0*\ell \delta^*(x_0)$, then we have:

$$
    \frac{x}{x_0}=1+\frac{\ell}{x_0/\delta^*(x_0)} = 1+\frac{\ell C^2}{Re_{\delta^*(x_0)}}
$$

Finally:
$$
    \frac{\delta^*(x)}{\delta^*(x_0)} = \sqrt{\frac{x}{x_0}} = \sqrt{1+\frac{\ell C^2}{Re_{\delta^*(x_0)}}}
$$

## Factor term for $v(x,y)$

We have that $V_\infty(x) = \frac{1}{2} U_\infty \frac{1}{\sqrt{Re_x}} \lim_{\eta\to\infty} (\eta f'(\eta) - f(\eta)) =\frac{1}{2} U_\infty \frac{C}{\sqrt{Re_x}}$. Thus:
$$
V_\infty(x) = \frac{1}{2} U_\infty \frac{C^2}{Re_{\delta^*(x)}} = \frac{1}{2} U_\infty \frac{C^2}{Re_{\delta^*(x_0)} \sqrt{1+\frac{\ell C^2}{Re_{\delta^*(x_0)}}}}
$$

