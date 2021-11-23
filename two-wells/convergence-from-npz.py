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
E=exact['E']
plt.plot
paths = glob.glob('*.npz')

plt.figure('latest-entropy')
plt.plot(E, correct_S, ':', label='exact', linewidth=2)

#Combines two strings, truncating the longer one
#to the length of the shorter one and adding them
def combine_data(a,b):
    if type(b) is int:
        return a
    elif type(a) is int:
        return b
    elif len(a) < len(b):
        return combine_data(b,a)
    elif len(a) == len(b):
        return a+b
    else:
        return a[:len(b)] + b

for fname in paths:
    if any( [method in fname[:-3] for method in ['sad', 'z', 'wl', 'itwl']] ) and 'seed-1+' in fname:
        tail = fname[:fname.find('seed')]
        b = fname[fname.find('seed'):]
        i= fname.find('seed') + b.find('+')
        front = fname[i:]

        mean_e=0
        mean_which=0
        hist=0
        S=0
        moves=0
        errors_S=0


        i=0
        for seed in ['1','12','123','1234','12345','123456','1234567','12345678']:
            print(fname)
            base = fname[:-4]
            i+=1
            fname = tail + 'seed-' + seed + front
            
            de_ind = base.rfind('de')
            precision = base[de_ind:]
            method = base[:base.find('+')]
            data = np.load(fname)

            mean_e=combine_data(mean_e,data['mean_e'])
            mean_which=combine_data(mean_which,data['mean_which'])
            try:
                hist=combine_data(hist,data['hist'])
            except:
                hist=None
            S=combine_data(S,data['S'])
            moves=combine_data(moves,data['moves'])
            errors_S=combine_data(errors_S,data['errors_S'])

        #finish average
        mean_e=mean_e/i
        mean_which=mean_which/i
        hist=hist/i
        S=S/i
        moves=moves/i
        errors_S=errors_S/i

        plt.figure('fraction-well')
        plt.plot(mean_e, mean_which, label=base)
    

        if hist is not None:
            plt.figure('histogram')
            plt.plot(mean_e, hist, label=base)

        plt.figure('latest-entropy')

        if method in {'wl','itwl','sad'}:
            plt.plot(E, S, label=precision, marker = styles.marker(base),
                    color = styles.color(base), linestyle= styles.linestyle(base), markevery=50)
        elif method == 'z':
            plt.plot(E, S, label=precision, color = styles.color(base), linestyle= styles.linestyle(base))

        plt.figure('convergence')
        if method in {'wl','itwl','sad'}:
            plt.loglog(moves, errors_S, label=precision, marker = styles.marker(base), color = styles.color(base), linestyle= styles.linestyle(base), markevery=2)
        elif method == 'z':
            plt.loglog(moves, errors_S, label=precision, color = styles.color(base), linestyle= styles.linestyle(base), linewidth = 3)

        


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