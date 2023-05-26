import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('data/epoch.dat')
phi, theta, loss, phi_grad, theta_grad = data.T

plt.plot(theta)
plt.show()
