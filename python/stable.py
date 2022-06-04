import numpy as np
import matplotlib.pyplot as plt


mass = [1, 0.8181818182, 0.66666666667, 0.5384615385, 0.4285714286, 0.33333333333, 0.25, 0.1764705882, 0.1111111111, 0.05263157895]
sep1 = [0.5216504, 0.52013, 0.519, 0.5203, 0.51972, 0.523, 0.52309, 0.525, 0.552, 0.593]
sep2 = [0.5216, 0.52, 0.515, 0.5, 0.515, 0.52, 0.52, 0.52, 0.53, 0.58]
sep3 = [0.52, 0.515, 0.51, 0.495, 0.5, 0.5, 0.5, 0.5, 0.52, 0.56]
sep4 = [0.515, 0.51, 0.5, 0.49, 0.48, 0.48, 0.48, 0.48, 0.5, 0.54]
sep5 = [0.51, 0.5, 0.48, 0.485, 0.47, 0.46, 0.46, 0.46, 0.48, 0.52]
sep6 = [0.505, 0.49, 0.46, 0.484, 0.46, 0.44, 0.44, 0.44, 0.46, 0.5]
sep7 = [0.5, 0.485, 0.44, 0.482, 0.45, 0.42, 0.42, 0.42, 0.44, 0.48]
sep8 = [0.5216505, 0.52014, 0.52, 0.5204, 0.51973, 0.524, 0.52311, 0.526, 0.551, 0.594]
sep9 = [0.5217, 0.53, 0.53, 0.52, 0.52, 0.54, 0.54, 0.53, 0.6, 0.6]
sep10 = [0.522, 0.54, 0.56, 0.53, 0.54, 0.55, 0.55, 0.54, 0.62, 0.62]
sep11 = [0.53, 0.55, 0.58, 0.54, 0.56, 0.56, 0.56, 0.56, 0.64, 0.64]
sep12 = [0.54, 0.56, 0.6, 0.56, 0.58, 0.58, 0.58, 0.58, 0.66, 0.66]
sep13 = [0.55, 0.57, 0.61, 0.58, 0.6, 0.6, 0.6, 0.6, 0.68, 0.68]
sep14 = [0.56, 0.58, 0.62, 0.6, 0.62, 0.62, 0.62, 0.62, 0.7, 0.7]

stab = [0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36]


xData = mass
yData = sep1
plt.figure()
plt.xlabel('Mass Ratio', size=16)
plt.ylabel('Binary Separation (AU)', size=16)
plt.scatter(xData, yData, color='r', marker='o', label='Stable')
plt.scatter(xData, sep2, color='r', marker='o')
plt.scatter(xData, sep3, color='r', marker='o')
plt.scatter(xData, sep4, color='r', marker='o')
plt.scatter(xData, sep5, color='r', marker='o')
plt.scatter(xData, sep6, color='r', marker='o')
plt.scatter(xData, sep7, color='r', marker='o')
plt.scatter(xData, sep8, color='b', marker='o', label='Unstable')
plt.scatter(xData, sep9, color='b', marker='o')
plt.scatter(xData, sep10, color='b', marker='o')
plt.scatter(xData, sep11, color='b', marker='o')
plt.scatter(xData, sep12, color='b', marker='o')
plt.scatter(xData, sep13, color='b', marker='o')
plt.scatter(xData, sep14, color='b', marker='o')
plt.plot(xData, stab, color='m', label='Stable limit')
plt.legend()

plt.show()
