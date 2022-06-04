import numpy as np
import matplotlib.pyplot as plt


mass = [1, 0.8181818182, 0.66666666667, 0.5384615385, 0.4285714286, 0.33333333333, 0.25, 0.1764705882, 0.1111111111, 0.05263157895]
life = [477, 472, 248, 199, 55, 18, 450, 22, 45, 220]



xData = mass
yData = life
plt.figure()
plt.xlabel('Mass Ratio', size=16)
plt.ylabel('lifetime (YR)', size=16)
plt.scatter(xData, yData, color='r', marker='o', label='Lifetime')
plt.legend()

plt.show()
