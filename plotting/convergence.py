#!/usr/bin/python3

import numpy as np
import yaml, cbor, argparse, sys, os, glob
import scipy.constants as const
import scipy.optimize as optimize
import matplotlib.pyplot as plt
import colorcet as cc

import compute

parser = argparse.ArgumentParser(description="fake energies analysis")
parser.add_argument('base', nargs='*', help = 'the yaml or cbor files')
parser.add_argument('--intensive', action='store_true')

args = parser.parse_args()

prop_cycle = plt.rcParams['axes.prop_cycle']

colors = cc.glasbey_dark

def beautiful_name(base):
    name = ''
    if 'r-' == base[:2]:
        name += "Zeno's "
        base = base[2:]
    elif 'wl-' == base[:3]:
        name += 'WL '
        base = base[3:]
    elif 'sad-' == base[:4]:
        name += 'SAD '
        base = base[4:]
    elif 'itwl-' == base[:5]:
        name += r'$t^{-1}$-WL '
        base = base[5:]
    if base == 'erfinv' or base == 'quadratic':
        name += ''
        base = ''
    elif base[:7] == 'erfinv-':
        base = base[7:]
        name += rf' $\Delta E = {base}$'
        base = ''
    elif base[:10] == 'quadratic-':
        base = base[10:]
        name += rf' $\Delta E = {base}$'
        base = ''
    return name + base

print('''
Changes to do:

''')

def linear_density_of_states(E):
    return np.heaviside(E, 0.5)*np.heaviside(1-E, 0.5)
def quadratic_density_of_states(E):
    return 1.5*np.sqrt(E)*np.heaviside(E, 0)*np.heaviside(1-E, 0)
def gaussian_density_of_states(E):
    return (np.pi*(sigma**3)*np.sqrt(32*np.log(-1/E))) / -E / (4*np.pi/3)
def other_density_of_states(E):
    return np.zeros_like(E)
def erfinv_density_of_states(E):
    return np.sqrt(1/(np.pi*N))*np.exp(-(E - mean_erfinv_energy)**2/N)
def erfinv_E_from_T(T):
    return -(N/2)/T
def piecewise_density_of_states(E):
    # D = np.zeros_like(E)
    # D[E>-e2] = (a**3 * np.sqrt(E[E>-e2]+e1) / 2) + ( (b - (b - a)*np.sqrt(E[E>-e2]/e2 + 1))**2 / (2*e2*np.sqrt(E[E>-e2]/e2+1)) ) + ( (b + (b - a)*np.sqrt(E[E>-e2]/e2 + 1))**2 / (2*e2*np.sqrt(E[E>-e2]/e2+1)) )
    # D[E>0] = (b + (b - a)*np.sqrt(E[E>0]/e2 + 1))**2 / (2*e2*np.sqrt(E[E>0]/e2+1))
    # return D
    if E > 0:
        return (b + (b - a)*np.sqrt(E/e2 + 1))**2 / (2*e2*np.sqrt(E/e2+1))
    elif E > -e2:
        return (a**3 * np.sqrt(E+e1) / 2) + ( (b - (b - a)*np.sqrt(E/e2 + 1))**2 / (2*e2*np.sqrt(E/e2+1)) ) + ( (b + (b - a)*np.sqrt(E/e2 + 1))**2 / (2*e2*np.sqrt(E/e2+1)) )
    return 0
piecewise_density_of_states = np.vectorize(piecewise_density_of_states)

#Read Data
moves = {}
error = {} #store the max error in each move
bases = []
print('base', args.base)

energy_boundaries = {}
entropy_boundaries = {}
mean_energy = {}
lnw = {}
systems = {}
#each file has different path (including extension) so concatenating is easy
for base in args.base:
    #change base to have the cbor files. currently has the directory
    if '.cbor' in base or '.yaml' in base:
        base = base[:-5]
    bases.append(base)

for base in bases:
    print('reading', base)
    with open(base+'-system.dat') as f:
        systems[base] = yaml.safe_load(f)


print('done reading bases')
sigma = 1

E = np.linspace(0.01, 0.4, 1000) # FIXME make this depend on which system we have
if 'pieces' in bases[0]:
    E = np.linspace(-systems[bases[0]]['e2']+0.02, 31, 1000)
elif 'erfinv' in bases[0]:
    N = 3
    E = np.linspace(erfinv_E_from_T(0.11), 0, 1000)

#Analysis
exact_density_of_states = linear_density_of_states
if 'linear' in bases[0]:
    print('\n\n\nusing the linear_density_of_states\n\n\n')
    exact_density_of_states = linear_density_of_states
elif 'quadratic' in bases[0]:
    print('\n\n\nusing the quadratic_density_of_states\n\n\n')
    exact_density_of_states = quadratic_density_of_states
elif systems[bases[0]]['kind'] == 'Gaussian':
    print('\n\n\nusing the gaussian_density_of_states\n\n\n')
    sigma = systems[bases[0]]['sigma']
    exact_density_of_states = gaussian_density_of_states
elif systems[bases[0]]['kind'] == 'Erfinv':
    print('\n\n\nusing the erfinv\n\n\n')
    mean_erfinv_energy = systems[bases[0]]['mean_energy']
    N = systems[bases[0]]['N']
    exact_density_of_states = erfinv_density_of_states
elif systems[bases[0]]['kind'] == 'Pieces':
    print('\n\n\nusing the pieces\n\n\n')
    a = systems[bases[0]]['a']
    b = systems[bases[0]]['b']
    e1 = systems[bases[0]]['e1']
    e2 = systems[bases[0]]['e2']
    exact_density_of_states = piecewise_density_of_states
else:
    print('\n\n\nusing the most bogus density of states\n\n\n', systems[bases[0]]['kind'])
    exact_density_of_states = other_density_of_states

# We can compute the exact entropy now, at our energies E
print('computing exact density of states')
exact_entropy = np.log(exact_density_of_states(E))
if 'pieces' in base:
    # FIXME we should properly normalize piecewise_density_of_states
    exact_entropy -= np.log(sum(exact_density_of_states(E))*(E[1]-E[0]))
print('done computing exact density of states')

for frame in range(10000):
    plt.clf()
    which_color = 0
    plotted_something = False
    if 'erfinv' in base:
        plt.plot(E, exact_entropy, '--', label='exact')
        plt.ylabel('$S(E)$')
        plt.ylim(exact_entropy.min()*1.1 - exact_entropy.max()*0.1, -exact_entropy.min()*0.1 + exact_entropy.max()*1.1)
    else:
        plt.plot(E, np.exp(exact_entropy), '--', label='exact')
        plt.ylabel('density of states')

    for base in bases:
        color = colors[which_color]
        which_color += 1

        if base not in moves:
            moves[base] = []
            error[base] = []
        frames = sorted(glob.glob(base+'/*.cbor'))
        if frame >= len(frames):
            continue
        f = os.path.splitext(frames[frame])[0]
        mymove = float(os.path.basename(f))
        moves[base].append(mymove)
        print(f'working on {base} with moves {mymove} which is {f}')

        energy_b = np.loadtxt(f+'-energy-boundaries.dat')
        mean_e = np.loadtxt(f+'-mean-energy.dat')
        my_lnw = np.loadtxt(f+'-lnw.dat')
        
        if energy_b.ndim == 0: #in case of a single value
            energy_b = np.array([energy_b.item()])

        if energy_b[0] < energy_b[-1]:
            energy_b = np.flip(energy_b)
            mean_e = np.flip(mean_e)
            my_lnw = np.flip(my_lnw)

        # Create a function for the entropy based on this number of moves:
        l_function, _, _ = compute.linear_entropy(energy_b, mean_e, my_lnw)
        # l_function, _, _ = compute.step_entropy(energy_b, mean_e, my_lnw)
        entropy_here = l_function(E)
        if 'erfinv' in base:
            plt.plot(E, entropy_here, label=beautiful_name(f))
        else:
            plt.plot(E, np.exp(entropy_here), label=beautiful_name(f))
        max_error = np.max(np.abs(entropy_here - exact_entropy))
        error[base].append(max_error)
        plotted_something = True
    if not plotted_something:
        break
    plt.xlabel('E')
    plt.legend()
    plt.draw_if_interactive()
    plt.pause(0.1)

#Plotting
mins = 1e10
maxs = -1e10
mint = 10
maxt = 0
plt.figure('convergence')
for base in bases:
    mins = min(error[base]+ [mins])
    maxs = max([e for e in error[base] if e < 100]+ [maxs])
    maxt = max([maxt]+ moves[base])
    plt.loglog(moves[base], error[base], label=beautiful_name(base))
t = np.linspace(mint, maxt, 3)
for t0 in 10.0**np.arange(-3, 14, 2.0):
    plt.loglog(t, np.sqrt(t0/t), color='xkcd:gray', alpha=0.2)

plt.ylim(mins, maxs)
plt.xlim(mint, maxt)
plt.xlabel('Moves')
plt.ylabel('Error (S - S$_{exact}$)')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('convergence.pdf')
plt.show()
