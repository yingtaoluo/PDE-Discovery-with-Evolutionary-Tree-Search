from PDE_find import FiniteDiff
import numpy as np
import scipy.io as scio
import itertools
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pylab
pylab.rcParams['figure.figsize'] = (12, 8)
import pdb

data = scio.loadmat('./data/burgers.mat')
u = np.real(data['usol'])
x = np.real(data['x'][0])
t = np.real(data['t'])
dt = t[1]-t[0]
dx = x[2]-x[1]
n, m = u.shape

ut = np.zeros((n, m), dtype=np.complex64)
ux = np.zeros((n, m), dtype=np.complex64)
uxx = np.zeros((n, m), dtype=np.complex64)

for i in range(n):
    ut[i, :] = FiniteDiff(u[i, :], dt, 1)
for i in range(m):
    ux[:, i] = FiniteDiff(u[:, i], dx, 1)
    uxx[:, i] = FiniteDiff(u[:, i], dx, 2)

left = np.zeros((n, m), dtype=np.complex64)
numerator = u * ux
for i in range(m):
    left[:, i] = FiniteDiff(numerator[:, i], dx, 1)
right = ux**2 + u*uxx

print(left)
print('-----------------')
print(right)
print('-----------------')
error_relative = abs((left-right)/(right+1e-40))
print(error_relative)
print(error_relative.mean())

nx,nt = len(x), len(t)
X, T = np.meshgrid(x, t)
fig1 = plt.figure()
ax = fig1.gca(projection='3d')
surf = ax.plot_surface(X, T, left.real.T, rstride=1, cstride=1, cmap=plt.cm.coolwarm,
    linewidth=0, antialiased=False)
plt.title('Burgers Equation', fontsize = 20)
plt.xlabel('x', fontsize = 16)
plt.ylabel('t', fontsize = 16)
plt.show()

pdb.set_trace()
