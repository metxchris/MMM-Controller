import numpy as np
import fib4, fib5
import time
# import mmm

# print(fib1.fib.__doc__)
t1 = time.perf_counter()
a=np.zeros(100000000)
t = 0
# fib4.fib(a, len(a))
fib5.fib5(a)
t2 = time.perf_counter()
print(t2-t1)
print(t)
