# Boundary Layer and Blasius profile

For the boundary layer theory, just reference some book of fluid mechanics, for example the one from McGrawHill.


## Blasius profile, some notes

When integrating the equations the v component cannot be neglected, so the boundary condition on the top has to be changed appropiatly with U_inf and V_inf parameters. V_inf is small compared to U_inf = 1, around 0.00147 I get when the flow is 2D and incompressible. Also, there is an expression analystic for it which is:

Let $a = lim_{\eta \to \infty} \eta - f(\eta) =  lim_{\eta \to \infty} \eta f'(\eta) - f(\eta) \simeq 1.72$. Then, since $Re_{\delta^*}/Re_x = \delta^*/x = a/\sqrt{Re_x}$, we deduce that $Re_{\delta^*} = a \sqrt{Re_x}$ and so:

$$
    v = Uinf/2 * 1/\sqrt{Re_x} * (\eta * f'(\eta) - f(\eta)) --> Vinf = Uinf/2 * a^2/Re_{\delta^*}
$$
