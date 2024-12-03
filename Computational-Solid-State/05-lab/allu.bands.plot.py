### Install libraries

import numpy as np
import matplotlib.pyplot as plt
bands = np.genfromtxt( 'allu.bands.dat.gnu' )
## Shift all energies by Fermi energy (can be extracted from QE output files)
Ef = 7.6703   # Remember to update the value!
bands[:,1] = bands[:,1] - Ef
fig = plt.figure(1, [10, 6], dpi=100)

plt.title('Al (fcc) Band structure')
plt.xlim(0.0,  3.0731)
plt.ylim(-16,10)
## ## Special k points
plt.xticks([0.0000, 1.0000, 1.5000, 2.2071, 3.0731],[r'$\Gamma$', 'X', 'W', 'L', r'$\Gamma$']) 
## Coordinates can be extracted from the output of bands.x
plt.yticks(np.linspace(-16, 10.0, num=14, dtype=float, endpoint=True))
plt.grid(axis='x', color='k', linestyle='-')
plt.ylabel(r'Energy ($eV$)')

#Number of bands and total number of k-points + 1 along the chosen path
num_bands = 10
num_pts = 161

for i in range(num_bands):
    plt.plot(bands[i*num_pts:i*num_pts+num_pts,0], bands[i*num_pts:i*num_pts+num_pts,1], linestyle = '-',color = 'k', linewidth=2.0)
## Plot Fermi level
plt.axhline(y = 0.0, color = 'b', linestyle = '--')
plt.savefig('allu.bands.plot.png', transparent=False, dpi=300, bbox_inches='tight')
plt.show()
