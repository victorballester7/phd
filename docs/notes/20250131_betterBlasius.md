# Non-Polynomial Fitting for Blasius Boundary Layer Profile

You're right—polynomials aren't the best choice for fitting the **Blasius boundary layer profile** because it has an asymptotic behavior and an inflection point, making polynomial approximations prone to oscillations (Runge's phenomenon). Instead, better approaches include:

### **1. Self-Similar Function Approximation (Recommended)**
Since the Blasius solution is given by \( u/U = f'(\eta) \), where \( \eta = y \sqrt{U / (\nu x)} \), a natural way to fit it is using a function of the form:
\[
u/U = A_1 \tanh(A_2 \eta) + A_3 \eta e^{-A_4 \eta}
\]
or other sigmoid-based functions.

✅ **Pros**: 
- Captures asymptotics correctly.
- Smooth and stable.

### **2. Exponential or Logarithmic Fits**
An alternative is an exponential series, such as:
\[
u/U = 1 - A_1 e^{-A_2 \eta} - A_3 e^{-A_4 \eta}
\]
or a log-like function.

✅ **Pros**: 
- Good for capturing the rapid change near the wall.
- Avoids high-degree polynomial oscillations.

### **3. Spline Interpolation**
If you're just trying to get a smooth numerical representation, **B-splines** or **cubic splines** can work well.

✅ **Pros**:
- Keeps local control over curvature.
- No artificial oscillations like polynomials.

Would you be fitting this numerically from CFD data or using theoretical values?
