import numpy as np
import matplotlib.pyplot as plt


filename = "/home/guilhermebaos/QUARMEN-Lessons/01-labs/01-lab/output.txt"
column = 1

# Load the data
data = np.loadtxt(filename, delimiter=',', dtype=float)
data = data[:, column]
samples = len(data)

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
plt.bar(bin_edges[:-1], heights / samples, width=bin_width, edgecolor='black', align='edge')
plt.xticks(bin_edges, rotation=45)


# Cumulative distribution
cumdata = np.cumsum(heights / samples)
cumaxis = np.array([(bin_edges[i] + bin_edges[i+1]) / 2 for i in range(len(heights))])
plt.subplot(rows, cols, 2)
plt.title("Cumulative Distribution of the Data")
plt.plot(cumaxis, cumdata)
plt.xticks(rotation=45)


# Stats
left_column_width = 14
right_column_width = 14

print(f"""
{'Stat':<{left_column_width}} {'Value':>{right_column_width}}
{'-' * (left_column_width + right_column_width + 1)}""")
print(f"""
{'Mean':<{left_column_width}} {np.mean(data):>{right_column_width}.5f}
{'Std':<{left_column_width}} {np.std(data):>{right_column_width}.5f}
{'Var':<{left_column_width}} {np.var(data):>{right_column_width}.5f}
{'Std of Mean':<{left_column_width}} {np.std(data) / np.sqrt(samples):>{right_column_width}.5f}
""")

# Save the figure
plt.tight_layout()
plt.savefig('output.png')