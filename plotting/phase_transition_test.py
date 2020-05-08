import numpy as np
import matplotlib.pyplot as plt


row = 10
col = 10
N_sites = 10

pressure = np.zeros((row, col))
gibbs_free = np.zeros((row, col))
"""
T = np.zeros((row, col))
s_excess = np.zeros((row, col))
chem_pot = np.zeros((row, col))
N = np.zeros((row, col))
energy_data = np.zeros((row, col))
E = np.zeros((row, col))
for i in np.arange(0,col-1,1):
    for j in np.arange(0,row-1,1):
        T[i][j] = 2
        s_excess[i][j] = 3
        chem_pot[i][j] = 4
        N[i][j] = 5
        energy_data[i][j] = 6
        E = 7
"""


x = np.arange(1, 10)
y = x.reshape(-1, 1)
h = x * y

cs = plt.contourf(h, levels=[10, 30, 50],
    colors=['#808080', '#A0A0A0', '#C0C0C0'], extend='both')
cs.cmap.set_over('red')
cs.cmap.set_under('blue')
cs.changed()


T = np.random.rand(row, col)
s_excess = np.random.rand(row, col)
chem_pot = np.random.rand(row, col)
N = np.random.rand(row, col)
energy_data = np.random.rand(row, col)
E = np.random.rand(row, col)

#print((2 * 3 + 4*5 - 6)/10**2)

for i in np.arange(0,col,1):
    for j in np.arange(0,row,1):
        pressure[i][j] = (T[i][j] * s_excess[i][j] + chem_pot[i][j] * N[i][j] - energy_data[i][j])/N_sites**2
        #pressure[i][j] = 1
        gibbs_free[i][j] = energy_data[i][j] - 2* T[i][j] * s_excess[i][j]
    #print(pressure)
"""
    plt.figure('pressure')
    plt.clf()
    #pressure = T(ds/dv)U,N
    plt.pcolor(N,E,pressure,)
    plt.xlabel('$N$')
    plt.ylabel('$E$')
    plt.colorbar()

    plt.figure('gibbs')
    plt.clf()
    plt.pcolor(N,E,gibbs_free)
    plt.xlabel('$N$')
    plt.ylabel('$E$')
    plt.colorbar()

    plt.pause(1)

plt.ioff()
plt.show()
"""
"""
def f(x, y):
    return np.sin(x) ** 10 + np.cos(10 + y * x) * np.cos(x)

x = np.linspace(0, 5, 50)
y = np.linspace(0, 5, 40)

X, Y = np.meshgrid(x, y)
Z = f(X, Y)

plt.contour(X, Y, Z, colors='black');

plt.contour(X, Y, Z, 20, cmap='RdGy');

plt.colorbar();

plt.imshow(Z, extent=[0, 5, 0, 5], origin='lower',
           cmap='RdGy')
plt.colorbar()
plt.axis(aspect='image');

contours = plt.contour(X, Y, Z, 3, colors='black')
plt.clabel(contours, inline=True, fontsize=8)

plt.imshow(Z, extent=[0, 5, 0, 5], origin='lower',
           cmap='RdGy', alpha=0.5)
plt.colorbar();
"""
plt.contour(N, E, gibbs_free, 20, cmap='RdGy');
plt.imshow(gibbs_free, extent=[0, 5, 0, 5], origin='lower',
           cmap='RdGy', alpha=0.5)
plt.colorbar();
plt.show()
