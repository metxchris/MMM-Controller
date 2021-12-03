from cycler import cycler
import matplotlib.pyplot as plt

line_cycle = (cycler(color=[
                (0.094, 0.353, 0.663), # Blue
                (0.933, 0.180, 0.184), # Red
                (0.957, 0.490, 0.137), # Orange
                (0, 0.549, 0.282), # Green
                (0.4, 0.173, 0.569), # Purple
                (0.984, 0.722, 0.153)]) # Yellow
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