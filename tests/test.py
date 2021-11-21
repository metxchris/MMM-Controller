import numpy as np # 3rd Party Packages
import matplotlib.pyplot as plt

def basic_plot():
    data = {'a': np.arange(50),
            'c': np.random.randint(0, 50, 50),
            'd': np.random.randn(50)}
    data['b'] = data['a'] + 10 * np.random.randn(50)
    data['d'] = np.abs(data['d']) * 100

    plt.scatter('a', 'b', c='c', s='d', data=data)
    plt.xlabel('entry a')
    plt.ylabel('entry b')
    plt.show()

if __name__ == '__main__':
    basic_plot()

    # plt.figure()
    # plt.plot(x,y1,lw=2)
    # plt.plot(x0,y0, 'o', mew=1.5, fillstyle='none', markersize=4)
    # plt.show()
