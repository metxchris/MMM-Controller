import numpy as np
import numba
from numba import jit
from timeit import timeit, repeat

print(numba.__version__)

@jit(nopython=True)
def go_fast(a): # Function is compiled to machine code when called the first time
    trace = 0.0
    # assuming square input matrix
    for i in range(a.shape[0]):   # Numba likes loops
        trace += np.tanh(a[i, i]) # Numba likes NumPy functions
    return a + trace              # Numba likes NumPy broadcasting

def go_numpy(a):
    return a + np.tanh(np.diagonal(a)).sum()

x = np.arange(1000).reshape(100, 10)
t1 = np.array(repeat(lambda: go_fast(x), number=10000))
t2 = np.array(repeat(lambda: go_numpy(x), number=10000))
print(t1)
print(t2)

print("speed ratio", t1[1:] / t2[1:])
# go_fast(2*x)
