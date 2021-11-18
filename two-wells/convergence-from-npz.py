import os
import numpy as np
import glob
import matplotlib.pyplot as plt
import system
import styles

lowest_interesting_E = -1.07
highest_interesting_E = -0.5

exact = np.load(system.name()+'.npz')
correct_S=exact['correct_S']
paths = glob.glob('*.npz')


for fname in paths:
    if ['sad', 'z', 'wl', 'itwl'] in fname:
        base = fname[:-4]
        method = base[:base.find('-')]
        data = np.load(fname)

        E = data['E']

        mean_e = data['mean_e']
        mean_which=data['mean_which']
        hist=data['hist']

        plt.figure('fraction-well')
        mean_which = np.loadtxt(f'{base}-which.dat')
        plt.plot(mean_e, mean_which, label=base)
    

        if hist is not None:
            plt.figure('histogram')
            plt.plot(mean_e, hist, label=base)


        S=data['S']
        moves=data['moves']
        errors_S=data['errors_S']

        plt.figure('latest-entropy')

        if method in {'wl','itwl','sad'}:
            plt.plot(E, S, label=base, marker = styles.marker(base),
                    color = styles.color(base), linestyle= styles.linestyle(base), markevery=100)
        elif method == 'z':
            plt.plot(E, S, label=base, color = styles.color(base), linestyle= styles.linestyle(base))

        plt.figure('convergence')
        if method in {'wl','itwl','sad'}:
            plt.loglog(moves, errors_S, label=base, marker = styles.marker(base), color = styles.color(base), linestyle= styles.linestyle(base), markevery=2)
        elif method == 'z':
            plt.loglog(moves, errors_S, label=base, color = styles.color(base), linestyle= styles.linestyle(base), linewidth = 3)

        


plt.figure('latest-entropy')
plt.xlabel(r'$E$')
plt.ylabel(r'$S(E)$')
plt.legend()
plt.savefig(system.system+'.svg')

plt.figure('fraction-well')
plt.xlabel(r'E')
plt.ylabel(r'Proportion in Small Well')
plt.legend()
plt.savefig(system.system+'-which.svg')

if hist is not None:
    plt.figure('histogram')
    plt.xlabel(r'$E$')
    plt.ylabel(r'# of Visitors')
    plt.legend()
    plt.savefig(system.system+'-histogram.svg')

plt.figure('convergence')
plt.xlabel(r'# of Moves')
plt.ylabel(rf'max error in entropy between {lowest_interesting_E} and {highest_interesting_E}')
plt.ylim(1e-2, 1e2)
plt.legend()
#make diagonal lines for convergence
x = np.linspace(1e-10,1e20,2)
y = 1/np.sqrt(x)
for i in range(20):
    plt.loglog(x,y*10**(4*i/5-2), color = 'y',alpha=0.5)
plt.savefig(system.system+'-convergence.svg')
plt.savefig(system.system+'-convergence.pdf')

plt.figure('convergence_heat_capacity')
plt.xlabel(r'# of Moves')
plt.ylabel(rf'max error in entropy between {lowest_interesting_E} and {highest_interesting_E}')
plt.ylim(1e-2, 1e2)
plt.legend()
#make diagonal lines for convergence
x = np.linspace(1e-10,1e20,2)
y = 1/np.sqrt(x)
for i in range(20):
    plt.loglog(x,y*10**(4*i/5-2), color = 'y',alpha=0.5)
plt.savefig(system.system+'-heat-capacity-convergence.svg')
plt.savefig(system.system+'-heat-capacity-convergence.pdf')

# TODO: Implement convergence saving in npz
# plt.figure('latest heat capacity')
# ax.legend()
# plt.savefig(system.system+'-heat-capacity.svg')
# plt.savefig(system.system+'-heat-capacity.pdf')

plt.show()