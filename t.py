mean = [52, 36]
cov = [[.1, 0],
       [0, .1]]  # diagonal covariance

import matplotlib.pyplot as plt
import numpy as np
np.random.seed(1)
x, y = np.random.multivariate_normal(mean, cov, 10).T
plt.plot(x, y, 'x')
plt.axis('equal')
plt.show()
