# randomness - Script for exploring random number generation.

# Set up configuration options and special features
import numpy as np
import matplotlib.pyplot as plt

seed_in = 80634  # Dr. Hurst's childhood zipcode
num_samples = 10000
np.random.seed(seed_in)
x = np.random.rand(num_samples)

fig, ax = plt.subplots()
ax.hist(x, alpha=0.5, ec='C0', lw=2)
ax.set_xlabel('x')
ax.set_ylabel('N(x)')

plt.show()