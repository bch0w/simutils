"""
Visualize what each of these spectral wavenumber filters looks like in 1D to
get a more intuitve sense of how they are used to modify the wavenumber domain
"""
import numpy as np
import matplotlib.pyplot as plt

a = 40
k = np.linspace(0, 1, 1000)

gaussian = np.pi ** 0.5 * a * np.exp(-k**2 * a**2 / 4)
exponential = (2 * a) / (1 + k**2 * a**2)
vonkarman = a / (1 + k**2 * a**2) ** 0.5

plt.plot(k, gaussian/gaussian.max(), c="C1", label=f"Gaussian")
plt.plot(k, exponential/exponential.max(), c="C2", label="Exponential")
plt.plot(k, vonkarman/vonkarman.max(), c="C3", label="von Karman")
plt.xlabel("K (wavenumber)")
plt.ylabel("Amplitude")
plt.legend()
plt.grid()
plt.title(f"1D Spectral Filters (a={a})")
plt.show()
