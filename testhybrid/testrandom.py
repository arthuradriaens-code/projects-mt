import numpy as np
import matplotlib.pyplot as plt

xpoints = np.random.uniform(low=-4000, high=4000.0, size=10000)
zpoints = np.random.uniform(low=-3000, high=0.0, size=10000)
points = [xpoints,zpoints]

for i in range(10000-1):
    plt.plot(points[0][i],points[1][i+1],marker='.')
plt.xlabel('$x_i$')
plt.ylabel('$z_{i+1}$')
plt.show()
