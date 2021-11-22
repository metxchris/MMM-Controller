import numpy as np # 3rd Party Packages
import matplotlib.pyplot as plt

def test():
    data = {'a': np.arange(50),
            'c': np.random.randint(0, 50, 50),
            'd': np.random.randn(50)}
    data['b'] = data['a'] + 10 * np.random.randn(50)
    data['d'] = np.abs(data['d']) * 100

    plt.scatter('a', 'b', c='c', s='d', data=data)
    plt.xlabel('entry a')
    plt.ylabel('entry b')
    plt.show()

def plot2d(x1, y1, x2=None, y2=None):
    plt.figure()
    plt.plot(x1,y1,lw=2)
    if x2 is not None and y2 is not None:
        plt.plot(x2, y2, 'o', mew=1.5, fillstyle='none', markersize=4)
    plt.show()

if __name__ == '__main__':
    test()
