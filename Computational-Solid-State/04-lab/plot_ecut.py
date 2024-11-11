import matplotlib.pyplot as plt

# Read the data from energy_ecut.dat
ecut, etot = [], []
with open('energy_ecut.dat', 'r') as file:
    next(file)  # Skip the header line
    for line in file:
        e, t = line.split()
        ecut.append(float(e))
        etot.append(float(t))

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(ecut, etot, marker='o', linestyle='-', color='b')
plt.title('Total Energy (Ryd) vs Ecut (Ryd)', fontsize=18)
plt.xlabel('Ecut (Ryd)', fontsize=16)
plt.ylabel('Etot (Ryd)', fontsize=16)
plt.grid(True)
plt.xlim(min(ecut) - 10, max(ecut) + 10)  # Adjust x-axis range
plt.ylim(min(etot) - 1, max(etot) + 1)  # Adjust y-axis range
plt.savefig('etot_ecut.png')
plt.show()
