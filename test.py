import numpy as np
import matplotlib.pyplot as plt


def fp(rr: np.ndarray, p: int, n: int) -> np.ndarray:

    # Summation index
    nn = np.arange(n)

    # Mesh
    mesh_nn, mesh_rr = np.meshgrid(nn, rr)

    # Compute
    ff = np.real(np.sum((mesh_nn + 0.5 + 1j * mesh_rr)**(-p), axis=1))

    return ff



# Input values
p = 5
n = 2000

rmin = 0
rmax = 0.5
rste = 100

rr = np.linspace(rmin, rmax, rste)

ff = fp(rr, p, n)


plt.plot(rr, ff)
plt.show()