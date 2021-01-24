import numpy as np
from PDE_find import Diff, FiniteDiff
import scipy.io as scio
from requests import get
from inspect import isfunction
import math
import pdb
import random

cheat = True
see_tree = None


def get_random_int(max_int):
    random_result = get('https://www.random.org/integers/?num=1&min=0&max={0}&col=1&base=10&format=plain&rnd=new'.format(max_int)).content
    try:
        int(random_result)
    except:
        print(random_result)
    return int(random_result)


rand = get_random_int(1e6)
print(rand)
# 237204
np.random.seed(rand)
random.seed(rand)


data = scio.loadmat('./data/burgers.mat')
u = np.real(data['usol'])
x = np.real(data['x'][0])
t = np.real(data['t'])[:,0]
n, m = u.shape
dt = t[1]-t[0]
dx = x[2]-x[1]
ut = np.zeros((n, m), dtype=np.complex64)
for idx in range(n):
    ut[idx, :] = FiniteDiff(u[idx, :], dt)

# 扩充维度使得与u的size相同
x = np.tile(x, (101, 1)).transpose((1, 0))
t = np.tile(t, (256, 1))

# calculate the error of correct cofs & correct terms
ux = np.zeros((n, m), dtype=np.complex64)
uxx = np.zeros((n, m), dtype=np.complex64)
for idx in range(m):
    ux[:, idx] = FiniteDiff(u[:, idx], dx)

for idx in range(m):
    uxx[:, idx] = FiniteDiff(ux[:, idx], dx)

right = np.reshape(-u*ux+0.1*uxx, (n*m, 1))
left = np.reshape(ut, (n*m, 1))
diff = np.linalg.norm(left-right, 2)
print(diff)

if cheat:
    # ALL = np.array([['sin', 1, np.sin], ['cos', 1, np.cos], ['+', 2, np.add], ['-', 2, np.subtract],
    #                 ['*', 2, np.multiply], ['d', 2, Diff], ['u', 0, u], ['t', 0, t], ['x', 0, x]])
    # OPS = np.array([['sin', 1, np.sin], ['cos', 1, np.cos],
    #                 ['+', 2, np.add], ['-', 2, np.subtract], ['*', 2, np.multiply], ['d', 2, Diff]])
    # OP1 = np.array([['sin', 1, np.sin], ['cos', 1, np.cos]])

    ALL = np.array([['+', 2, np.add], ['-', 2, np.subtract],
                    ['*', 2, np.multiply], ['d', 2, Diff], ['u', 0, u], ['t', 0, t], ['x', 0, x]])
    OPS = np.array([['+', 2, np.add], ['-', 2, np.subtract], ['*', 2, np.multiply], ['d', 2, Diff]])
    OP1 = np.array([])

    OP2 = np.array([['+', 2, np.add], ['-', 2, np.subtract], ['*', 2, np.multiply], ['d', 2, Diff]])
    VARS = np.array([['u', 0, u], ['t', 0, t], ['x', 0, x]])
    den = np.array([['x', 0, x]])

else:
    ALL = np.array([['sin', 1, np.sin], ['cos', 1, np.cos], ['log', 1, np.log], ['+', 2, np.add], ['-', 2, np.subtract],
                    ['*', 2, np.multiply], ['/', 2, np.divide], ['d', 2, Diff], ['u', 0, u], ['t', 0, t], ['x', 0, x]])
    OPS = np.array([['sin', 1, np.sin], ['cos', 1, np.cos], ['log', 1, np.log],
                    ['+', 2, np.add], ['-', 2, np.subtract], ['*', 2, np.multiply], ['/', 2, np.divide],
                    ['d', 2, Diff]])
    OP1 = np.array([['sin', 1, np.sin], ['cos', 1, np.cos], ['log', 1, np.log]])
    OP2 = np.array(
        [['+', 2, np.add], ['-', 2, np.subtract], ['*', 2, np.multiply], ['/', 2, np.divide], ['d', 2, Diff]])
    VARS = np.array([['u', 0, u], ['t', 0, t], ['x', 0, x]])
    den = np.array([['t', 0, t], ['x', 0, x]])


