from cycler import cycler
import matplotlib.pyplot as plt

line_cycle = (cycler(color=[
                '#008fd5',
                '#fc4f30',
                '#e5ae38',
                '#6d904f',
                '#8b8b8b',
                '#810f7c'])
            + cycler(dashes=[
                (1, 0), 
                (6, 2, 1, 2), 
                (5, 2), 
                (1.5, 1), 
                (2, 1), 
                (4, 1)])
            + cycler(lw=[
                1.75, 
                1.75, 
                1.75, 
                2.0, 
                2.0, 
                2.0]))

plt.rcParams.update({'axes.prop_cycle': line_cycle})