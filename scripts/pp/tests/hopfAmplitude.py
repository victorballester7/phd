import numpy as np
import matplotlib.pyplot as plt

# x = 25
# y = 0.5
w = np.array([25.55, 25.6, 25.8, 26, 27, 30])
meanUL2 = np.array([1.911525e-01,2.125251e-01,2.954721e-01, 3.648846e-01,5.669644e-01,1.001902e+00])
meanVL2 = np.array([1.610369e-01,1.766622e-01,2.443431e-01,2.943190e-01,4.439065e-01,7.250598e-01])
wc = 25.41


x = np.linspace(np.min(w), np.max(w), 100)
y = 0.6 * np.sqrt(x - wc)

print(np.log(y) / np.log(x - wc))

plt.plot(np.log(w - wc), np.log(meanUL2), "o", label="DNS u L2")
plt.plot(np.log(w - wc), np.log(meanVL2), "o", label="DNS v L2")
plt.plot(np.log(w - wc), np.log(np.sqrt(meanUL2**2 + meanVL2**2)), "o", label="DNS sqrt(u^2 + v^2) L2")
plt.plot(np.log(x - wc), np.log(y), "-", label="~sqrt(w - wc)")
plt.xlabel("log(w - wc)")
plt.ylabel("log(||sqrt(u^2 + v^2)||_L2)")
plt.legend()
plt.title("Hopf Bifurcation Amplitude")
plt.grid()
plt.show()
