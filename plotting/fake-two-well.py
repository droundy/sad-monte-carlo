#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt

# Ideas for Henry:

# Run and test Zeno's algorithm on LJ clusters.

# Develop new challenging double-well or multiple-well potential we can solve
# analytically.

# 




# I want a fake potential that models a phase transition
#
# a) is in many dimensions
#
# b) has at least two minima
#
# c) has a global minimum that has a configurably smaller funnel than a local
#    minimum
#
# d) we can solve for the density of states of!


x = np.arange(-1,1,0.01)
y = np.arange(-1,1,0.01)

X,Y = np.meshgrid(x,y)

x0 = 0
y0 = x0
R0 = 1

V0 = -1*(1 - ((X-x0)**2+(Y-y0)**2)/R0**2)

# plt.figure('V0')
# plt.pcolormesh(X, Y, V0, shading='auto')
# # plt.contour(X, Y, V0, levels = [0])
# plt.colorbar()

x1 = 0.82
y1 = x1
R1 = 0.18

V1 = -1.1*(1 - ((X-x1)**2+(Y-y1)**2)/R1**2)

# plt.figure('V1')
# plt.pcolormesh(X, Y, V1, shading='auto')
# plt.colorbar()

V2 = np.zeros_like(V0)

plt.figure('Vmin')
plt.pcolormesh(X, Y, np.minimum(V2, np.minimum(V1, V0)), shading='auto')
plt.colorbar()
plt.gca().set_aspect('equal')

plt.show()