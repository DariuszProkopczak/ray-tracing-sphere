import csv
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
# code is taken from UChicago HPC Summer 2019 Slack Thread

with open('animation.ascii') as file:

    data = np.loadtxt('animation.ascii')

m, n = data.shape

x = range(0, m, n)

fig, ax = plt.subplots(1, figsize=(8,8))

fig.subplots_adjust(left=0, right=1, bottom=0, top=1)

g = ax.imshow(data[0:n], cmap='inferno');

ax.axis('off')

def animate(i, x, data):

    g.set_array(data[x[i]:x[i+1]])

    return g,

a = FuncAnimation(fig, animate, frames=len(x)-1, fargs=(x, data))

plt.show(block=True)