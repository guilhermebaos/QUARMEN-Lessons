import numpy as np
import matplotlib.pyplot as plt


filename = "/home/guilhermebaos/QUARMEN-Lessons/01-labs/01-lab/output.txt"
column = 1

# Load the data
data = np.loadtxt(filename, delimiter=',', dtype=float)
# data = data[:, column]

N = len(data)
k1 = int(0.05 * N)
k2 = N-1

Neq = k2 - k1 + 1


# Calculate Stats
ave = np.mean(data)
std = np.std(data)
var = np.var(data, mean=ave, ddof=1)

autocorrtotal = np.array([np.sum((data[0: N-i]-ave) * (data[0 + i: N]-ave)) for i in range(0, N-1)]) / (var * (N - 1))
print(autocorrtotal)

autocorr = np.array([np.sum((data[k1: k2-i]-ave) * (data[k1 + i: k2]-ave)) for i in range(1, Neq)]) / (var * (Neq - 1))
timecorr = 1 + 2 * np.sum(autocorr)

Neff = Neq / timecorr

# Prepare the figure
rows = 1
cols = 2
plt.figure(figsize=(8 * cols, 5 * rows))


# Distribution of the Data
heights, bin_edges = np.histogram(data, bins="auto")
bin_width = np.diff(bin_edges)

# Create a bar plot
plt.subplot(rows, cols, 1)
plt.title("Distribution of the Data")
plt.bar(bin_edges[:-1], heights / N, width=bin_width, edgecolor='black', align='edge')
plt.xticks(bin_edges, rotation=45)


# Cumulative distribution
# cumdata = np.cumsum(heights / N)
# cumaxis = np.array([(bin_edges[i] + bin_edges[i+1]) / 2 for i in range(len(heights))])
# plt.subplot(rows, cols, 2)
# plt.title("Cumulative Distribution of the Data")
# plt.plot(cumaxis, cumdata)
# plt.xticks(rotation=45)

# Autocorrelation time
plt.subplot(rows, cols, 2)
plt.title("Autocorrelation of Data")
plt.plot(range(0, N-1), autocorrtotal)


# Show Stats
left_column_width = 14
right_column_width = 14

print(f"""
{'Stat':<{left_column_width}} {'Value':>{right_column_width}}
{'-' * (left_column_width + right_column_width + 1)}""")
print(f"""
{'Mean':<{left_column_width}} {ave:>{right_column_width}.5f}
{'Std':<{left_column_width}} {std:>{right_column_width}.5f}
{'Var':<{left_column_width}} {var:>{right_column_width}.5f}
{'Std of Mean':<{left_column_width}} {std / np.sqrt(Neff):>{right_column_width}.5f}
{'Corr Time':<{left_column_width}} {timecorr:>{right_column_width}.5f}
\n
{'N':<{left_column_width}} {N:>{right_column_width}.5f}
{'Neq':<{left_column_width}} {Neq:>{right_column_width}.5f}
{'Neff':<{left_column_width}} {Neff:>{right_column_width}.5f}
""")

# Save the figure
plt.tight_layout()
plt.savefig('output.png')