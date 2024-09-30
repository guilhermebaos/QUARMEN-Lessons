import numpy as np

# Sample data
data = np.array([1, 2, 1, 3, 2, 3, 4, 5, 4, 4, 5])

# Compute histogram
heights, bin_edges = np.histogram(data, bins='auto')  # Use 'auto' for optimal bin size

# Print the results
print("Heights:", heights)
print("Bin Edges:", bin_edges)
