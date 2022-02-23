import os
import numpy as np
import glob
import matplotlib.pyplot as plt
import system
import styles
import heat_capacity

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
    if any( [method in fname[:-3] for method in [ 'sad', 'itwl']] ) and ('seed-1+' in fname and 'de-1e-05+step-0.01'  in fname):
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
        for seed in ['1']:
            try:
                base = fname[:-4]
                i+=1
                fname = tail + 'seed-' + seed + front
                print(fname)
                
                de_ind = base.rfind('de')
                precision = base[de_ind:]
                method = base[12:base.find('+')]
                data = np.load(fname)
                

                mean_e.append(data['mean_e'])
                mean_which.append(data['mean_which'])
                try:
                    hist.append(data['hist'])
                except:
                    hist=None
                if len(E) < len(data['E']): E = data['E']
                S.append(data['S'])
                if len(T) < len(data['T']): T = data['T']
                
                C.append(data['C'])
                moves.append(data['moves'])
                
                errors_S.append(data['errors_S'])
                
                errors_C.append(data['errors_C'])
            except:
                print(f'skipping file {fname}')
                pass

            #finish average
            mean_e=geometric_mean(mean_e)
            mean_which=geometric_mean(mean_which)
            hist=geometric_mean(hist)
            errors_S=geometric_mean(errors_S)
            
            errors_C=geometric_mean(errors_C)

            S=S[0]
            C=C[0]
            if method == 'itwl':
                label = r'$1/t$-WL' + r'-$E_{barr}$=0.'+styles.get_barrier(base)[8]
            if method == 'sad':
                label = r'SAD' + r'-$E_{barr}$=0.'+styles.get_barrier(base)[8]
            print(label)


            plt.figure('fraction-well')
            plt.plot(mean_e[:len(mean_which)], mean_which[:len(mean_e)], label=label)
        
            if hist is not None:
                plt.figure('histogram')
                plt.plot(mean_e[:len(hist)], hist[:len(mean_e)], label=label)

            plt.figure('latest-entropy')

            if method in {'wl','itwl','sad'}:
                plt.plot(E[:len(S)], S[:len(E)], label=label, marker = styles.marker(base),
                        color = styles.color(base), linestyle= styles.linestyle(base), markevery=250)
            elif method == 'z':
                plt.plot(E[:len(S)], S[:len(E)], label=label, color = styles.color(base), linestyle= styles.linestyle(base))
            
            heat_capacity.plot_from_data(T[:len(C)], C[:len(T)], fname=fname,ax=ax, axins=axins)

            plt.figure('convergence')
            if method in {'wl','itwl','sad'}:
                plt.loglog(moves[0][:len(errors_S)], errors_S[:len(moves[0])], label=label, marker = styles.marker(base), color = styles.color(base), linestyle= styles.linestyle(base), markevery=2)
            elif method == 'z':
                plt.loglog(moves[0][:len(errors_S)], errors_S[:len(moves[0])], label=label, color = styles.color(base), linestyle= styles.linestyle(base), linewidth = 3)

            plt.figure('convergence-heat-capacity')
            if method in {'wl','itwl','sad'}:
                plt.loglog(moves[0][:len(errors_C)], errors_C[:len(moves[0])], label=label, marker = styles.marker(base), color = styles.color(base), linestyle= styles.linestyle(base), markevery=2)
            elif method == 'z':
                plt.loglog(moves[0][:len(errors_C)], errors_C[:len(moves[0])], label=label, color = styles.color(base), linestyle= styles.linestyle(base), linewidth = 3)

        


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
x = np.linspace(1e-30,1e40,2)
y = 1/np.sqrt(x)
for i in range(50):
    plt.loglog(x,y*10**(4*i/5-2), color = 'y',alpha=0.5)
plt.savefig(system.system+'-convergence.svg')
plt.savefig(system.system+'-convergence.pdf')

plt.figure('convergence-heat-capacity')
plt.xlabel(r'# of Moves')
plt.ylabel(rf'max error in heat capacity between {lowest_interesting_E} and {highest_interesting_E}')
plt.ylim(1e-2, 1e2)
plt.legend()
#make diagonal lines for convergence
x = np.linspace(1e-10,1e20,2)
y = 1/np.sqrt(x)
for i in range(20):
    plt.loglog(x,y*10**(4*i/5-2), color = 'y',alpha=0.5)
plt.savefig(system.system+'-heat-capacity-convergence.svg')
plt.savefig(system.system+'-heat-capacity-convergence.pdf')


plt.figure('latest-heat-capacity')
ax.legend()
plt.savefig(system.system+'-heat-capacity.svg')
plt.savefig(system.system+'-heat-capacity.pdf')

plt.show()