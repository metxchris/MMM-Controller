import sys; sys.path.insert(0, '../../')
import time

# 3rd Party Packages
import numpy as np
from scipy.linalg import eig

timing_loops = 1

# CSV info
file_dir = 'data/'
file_num = '00001'
f = f'{file_dir}{file_num}'

# Number formats for printing values
nfmt = '.2e'
sfmt = '>13'

# Load matrix A
AR = data_array = np.genfromtxt(f'{f}_AR.csv', delimiter=',', dtype=complex)
AI = data_array = np.genfromtxt(f'{f}_AI.csv', delimiter=',', dtype=complex) * 1j
A = AR + AI

# Load matrix B
BR = data_array = np.genfromtxt(f'{f}_BR.csv', delimiter=',', dtype=complex)
BI = data_array = np.genfromtxt(f'{f}_BI.csv', delimiter=',', dtype=complex) * 1j
B = BR + BI

# Print matrix A
print('\nMatrix AR:')
for i in range(A.shape[0]):
    for j in range(A.shape[1]):
        s = f'  {A[i, j].real:{nfmt}}'
        print(f'{s:{sfmt}}', end='')
    print('')

print('\nMatrix AI:')
for i in range(A.shape[0]):
    for j in range(A.shape[1]):
        s = f'  {A[i, j].imag:{nfmt}}'
        print(f'{s:{sfmt}}', end='')
    print('')

# Print matrix B if not identity matrix
if not (B == np.identity(B.shape[0], dtype=complex)).all():
    print('\nMatrix BR:')
    for i in range(B.shape[0]):
        for j in range(B.shape[1]):
            s = f'  {B[i, j].real:{nfmt}}'
            print(f'{s:{sfmt}}', end='')
        print('')

    print('\nMatrix BI:')
    for i in range(B.shape[0]):
        for j in range(B.shape[1]):
            s = f'  {B[i, j].imag:{nfmt}}'
            print(f'{s:{sfmt}}', end='')
        print('')

# Solve eigensystem
t0 = time.perf_counter_ns()
for i in range(timing_loops):
    w, v = eig(A, B)
t1 = time.perf_counter_ns()

# Sort eigenvalues and eigenvectors by largest to smallest w.imag
isorted = np.flip(np.argsort(w.imag))
w = w[isorted]
v = v[:, isorted]

print('\nEigenvector Ratios (Sorted):')
for i in range(v.shape[0] - 1):
    for j in range(v.shape[1]):
        s = f'  {np.abs(v[i, j]/v[v.shape[0] - 1, j]):{nfmt}}'
        print(f'{s:{sfmt}}', end='')
    print('')

print('\nEigenvalues (Sorted):')
for i in range(w.shape[0]):
    sr = f'{w[i].real:{nfmt}}'
    si = f'{w[i].imag:{nfmt}}'
    print(f'{sr:{sfmt}}{si:{sfmt}}')
print('')

if timing_loops > 1:
    print('\nTime (s):')
    print(f'{(t1-t0) * 1e-9 / timing_loops:{nfmt}}')
