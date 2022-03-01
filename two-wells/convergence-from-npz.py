import os
import numpy as np
import glob
import matplotlib.pyplot as plt
import system
import styles
import heat_capacity
from results_plotting import Results

def arithmetic_mean(list_of_arr):
    mu=list_of_arr[0]
    i=1
    for arr in list_of_arr[1:]:
        try:
            print(len(mu))
            mu+=arr
            i+=1 
        except:
            pass
    return mu/i

def geometric_mean(list_of_arr):
    mu=list_of_arr[0]
    i=1
    for arr in list_of_arr[1:]:
        try:
            mu*=arr
            i+=1 
        except:
            pass
    return mu**(1/i)

lowest_interesting_E = -1.07
highest_interesting_E = -0.5

lowest_interesting_T=0.008

exact = np.load(os.path.join('.','thesis-data',system.name()+'.npz'))
correct_S=exact['correct_S']
E=exact['E']
dE = E[1] - E[0]
paths = glob.glob(os.path.join('thesis-data','*.npz'))

plt.figure('latest-entropy')
plt.plot(E, correct_S, ':', label='exact', linewidth=2)

fig, ax = plt.subplots(figsize=[5, 4], num='latest-heat-capacity')
axins = ax.inset_axes( 0.5 * np.array([1, 1, 0.47/0.5, 0.47/0.5]))#[0.005, 0.012, 25, 140])

heat_capacity.plot_from_data(exact['T'],exact['correct_C'], ax=ax, axins=axins)
ax.indicate_inset_zoom(axins, edgecolor="black")

#Combines two strings, truncating the longer one
#to the length of the shorter one and adding them
def combine_data(a,b, replace = False):
    if replace:#return only b--to replace a
        return b
    if type(b) is int:
        return a
    elif type(a) is int:
        return b
    elif len(a) < len(b):
        return np.concatenate(b[:len(a)] + a, b[len(a):])
    elif len(a) == len(b):
        return a+b
    else:
        return np.concatenate(a[:len(b)] + b, a[len(b):])

for fname in paths:
    if any( [method in fname[:-3] for method in [ 'sad', 'itwl']] ) and ('seed-1+' in fname and 'de-1e-05+step-0.0001'  in fname):
        tail = fname[:fname.find('seed')]
        b = fname[fname.find('seed'):]
        i= fname.find('seed') + b.find('+')
        front = fname[i:]

        mean_e=[]
        mean_which=[]
        hist=[]
        E=[]
        S=[]
        T=[]
        C=[]
        moves=[]
        errors_S=[]
        errors_C=[]

        i=0
        results = Results()
        for seed in [1, 12, 123, 1234, 12345, 123456, 1234567, 12345678]:
            seed = str(seed)
            try:
                base = fname[:-4]
                i+=1
                fname = tail + 'seed-' + seed + front
                #print(fname)
                
                de_ind = base.rfind('de')
                precision = base[de_ind:]
                method = base[12:base.find('+')]
                results.add_npz(fname)
            except:
                print(f'skipping file {fname}')
                pass
        results.median_method(ax, axins)

        


plt.figure('latest-entropy')
plt.xlabel(r'$E$')
plt.ylabel(r'$S(E)$')
plt.legend()
plt.ylim(-40,0)
plt.xlim(-1.15,-0.85)
plt.savefig(system.system+'latest-entropy.svg')
plt.savefig(system.system+'latest-entropy.pdf')

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
plt.xlim(1e8, 1e12)
plt.legend()
#make diagonal lines for convergence
x = np.linspace(1e-30,1e40,2)
y = 1/np.sqrt(x)
for i in range(50):
    plt.loglog(x,y*10**(4*i/5-2), color = 'y',alpha=0.5)
plt.savefig(system.system+'-convergence.svg')
plt.savefig(system.system+'-convergence.pdf')

plt.figure('convergence-heat-capacity')
plt.xlabel(r'# of Moves')
plt.ylabel(rf'max error in heat capacity above $T=${lowest_interesting_T}')
plt.ylim(1e-2, 1e2)
plt.xlim(1e8, 1e12)
plt.legend()
#make diagonal lines for convergence
x = np.linspace(1e-10,1e20,2)
y = 1/np.sqrt(x)
for i in range(20):
    plt.loglog(x,y*10**(4*i/5-2), color = 'y',alpha=0.5)
plt.savefig(system.system+'-heat-capacity-convergence.svg')
plt.savefig(system.system+'-heat-capacity-convergence.pdf')


plt.figure('latest-heat-capacity')
ax.legend(loc='lower right')
plt.savefig(system.system+'-heat-capacity.svg')
plt.savefig(system.system+'-heat-capacity.pdf')

plt.show()