#!/usr/bin/python3

import numpy as np
import yaml, cbor, argparse, sys, os, glob
import scipy.constants as const
import scipy.optimize as optimize
import matplotlib.pyplot as plt
import colorcet as cc

parser = argparse.ArgumentParser(description="fake energies analysis")
parser.add_argument('base', nargs='*', help = 'the yaml or cbor files')
parser.add_argument('--intensive', action='store_true')

args = parser.parse_args()


prop_cycle = plt.rcParams['axes.prop_cycle']

colors = cc.glasbey_dark

print('''
Changes to do:

- Loop over all the movie files (so we can see the convergence process)

- Compare the entropy with the exact entropy (ln of exact DOS) and find the
  biggest error by comparing at a regular set of energies.  Store that biggest
  error for each number of moves.

- Plot error versus number of moves, on a log-log plot.

- (optionally) Create a movie showing how the entropy evolves over time.

        exact_entropy = np.log(linear_density_of_states) #for convergence analysis
        #select energies evenly; exclusively
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
    return np.sqrt(np.pi/N)*np.exp(-(E - mean_erfinv_energy)**2/N)
def piecewise_density_of_states(E):
    if E > 0:
        return (b + (b - a)*np.sqrt(E/e2 + 1))**2 / (2*e2*np.sqrt(E/e2+1))
    elif E > -e2:
        return (a**3 * np.sqrt(E+e1) / 2) + ( (b - (b - a)*np.sqrt(E/e2 + 1))**2 / (2*e2*np.sqrt(E/e2+1)) ) + ( (b + (b - a)*np.sqrt(E/e2 + 1))**2 / (2*e2*np.sqrt(E/e2+1)) )
piecewise_density_of_states = np.vectorize(piecewise_density_of_states)

#The function needs to be callable

def fn_entropy(S_i_1, E_i_1, E_i, lnw_i, S_i):
    ''' This is the thing that should be zero '''
    ''' Equation initially is:
            W_i = e^coefficient * (inside) / denominator
            
            introducing ln simplifies to:
                lnW_i = ln(coefficient*inside) - ln(denominator)
                lnW_i = ln(coefficient) + ln(inside) - ln(denominator)
                0 = ln(coefficient) + ln(inside) - ln(denominator) - lnW_i
    '''
    delta_E = E_i_1 - E_i
    assert(delta_E>0)
    delta_S = S_i_1 - S_i
    S_0 = S_i_1
    if abs(delta_S) < 1e-14:
        return np.log(delta_E) + S_0 - lnw_i
    else:
        return np.log(delta_E) + S_0 - np.log((np.exp(delta_S)- 1)/delta_S) - lnw_i

    denominator = np.log((S_i - S_i_1) / (E_i - E_i_1))
    coefficient = (S_i*E_i - S_i_1*E_i_1) / (E_i - E_i_1)
    inside = np.log(
                    np.exp( E_i_1 * ((S_i - S_i_1) / (E_i - E_i_1)) )
                    - np.exp( E_i * ((S_i - S_i_1) / (E_i - E_i_1)) ) 
                    )
    return coefficient + inside - denominator - lnw_i
def optimize_bin_entropy(i, E_bounds, lnw, S_i):
    #i is the bin whose entropy we are calculating
    print('solving E_i_1', E_bounds[i-1], 'E_i', E_bounds[i], 'lnw_i', lnw[i], 'S_i', S_i)
    sol = optimize.root(fn_entropy, [0], args=(E_bounds[i-1], E_bounds[i], lnw[i], S_i))
    # print('   goodness is', fn_entropy(sol.x[0], E_bounds[i-1], E_bounds[i], lnw[i], S_i),
    #       'with entropy', sol.x[0])
    # print('      delta S =', sol.x[0] - S_i, np.exp(sol.x[0] - S_i))
    # print('      delta E =', E_bounds[i-1] - E_bounds[i])
    # print('      S_0 =', S_i)
    # print('      w =', np.exp(lnw[i]))
    # print('      S estimate =', lnw[i] - np.log(E_bounds[i-1] - E_bounds[i]))
    
    # plt.figure('debugging')
    # xxx = np.linspace(-10,100, 10000)
    # answer = np.zeros_like(xxx)
    # for j in range(len(xxx)):
    #     answer[j] = fn_entropy(xxx[j], E_bounds[i-1], E_bounds[i], lnw[i], S_i)
    # plt.plot(xxx, answer)
    # plt.show()
    return sol.x[0]

def compute_entropy_given_Smin(Smin, energy_boundaries, lnw):
    entropy_boundaries = np.zeros_like(energy_boundaries)
    entropy_boundaries[-1] = Smin
    for i in range(len(energy_boundaries)-1, 1, -1):
        print('solving for i =', i)
        entropy_boundaries[i-1] = optimize_bin_entropy(i, energy_boundaries, lnw, entropy_boundaries[i])
    return entropy_boundaries    

def entropy_boundary_badness(Smin, energy_boundaries, lnw):
    entropy_boundaries = compute_entropy_given_Smin(Smin, energy_boundaries, lnw)
    badness = 0
    for i in range(len(energy_boundaries)-2):
        de1 = energy_boundaries[i+1] - energy_boundaries[i]
        de2 = energy_boundaries[i+2] - energy_boundaries[i+1]
        dS1 = entropy_boundaries[i+1] - entropy_boundaries[i]
        dS2 = entropy_boundaries[i+2] - entropy_boundaries[i+1]
        badness += abs((dS2/de2 - dS1/de1)/(de2+de1))
    return badness

def optimize_entropy(energy_boundaries, lnw):
    res = optimize.minimize_scalar(entropy_boundary_badness, args=(energy_boundaries, lnw))
    Smin = res.x
    return compute_entropy_given_Smin(Smin, energy_boundaries, lnw)

def bisect_bin_entropy(i):
    #i is the bin whose entropy we are calculating
    return

def fn_for_beta(x, meanE_over_deltaE):
    if x < 1e-14:
        return 0.5*x - x*meanE_over_deltaE
    return x/(1-np.exp(-x)) - 1 - x*meanE_over_deltaE

def find_beta_deltaE(meanE_over_deltaE):
    # x = np.linspace(-100,100,10000)
    # plt.plot(x, np.vectorize(fn_for_beta)(x, meanE_over_deltaE))
    # plt.show()
    if meanE_over_deltaE == 0:
        x0 = 1e-6
        x1 = -1e-6
    else:
        x0 = 2*meanE_over_deltaE
        x1 = meanE_over_deltaE
    sol = optimize.root_scalar(fn_for_beta, args=(meanE_over_deltaE), x0 = x0, x1 = x1)
    # print(sol)
    return sol.root

def find_entropy_from_beta_and_lnw(beta, lnw, deltaE):
    if abs(beta*deltaE) < 1e-14:
        return lnw - np.log(deltaE)
    return lnw - np.log(deltaE) - np.log((np.exp(beta*deltaE)-1)/(beta*deltaE))


#Read Data
frames = {}
bases = []
for base in args.base:
    if '.cbor' in base or '.yaml' in base:
        base = base[:-5]
    bases.append(base)
    frames[base] = set()
    for f in glob.glob(base+'/*.cbor'):
        f = os.path.splitext(os.path.basename(f))[0]
        frames[base].add(f)

energy_boundaries = {}
mean_energy = {}
lnw = {}
systems = {}
#each file has different path (including extension) so concatenating is easy
for base in bases:
    print(bases)
    energy_boundaries[base] = np.loadtxt(base+'-energy-boundaries.dat')
    mean_energy[base] = np.loadtxt(base+'-mean-energy.dat')
    lnw[base] = np.loadtxt(base+'-lnw.dat')
    with open(base+'-system.dat') as f:
        systems[base] = yaml.safe_load(f)

    if energy_boundaries[base][0] < energy_boundaries[base][-1]:
        energy_boundaries[base] = np.flip(energy_boundaries[base])
        mean_energy[base] = np.flip(mean_energy[base])
        lnw[base] = np.flip(lnw[base])

entropy_boundaries={}
sigma = 1

print('energy_boundaries', energy_boundaries)

which_color = 0
for key in lnw:
    color = colors[which_color]
    which_color += 1

    # sigma.append(np.loadtxt(key+'-sigma.dat'))  #used in D(E) calculation

    #Analysis
    exact_density_of_states = linear_density_of_states
    if 'linear' in key:
        print('\n\n\nusing the linear_density_of_states\n\n\n')
        exact_density_of_states = linear_density_of_states
    elif 'quadratic' in key:
        print('\n\n\nusing the quadratic_density_of_states\n\n\n')
        exact_density_of_states = quadratic_density_of_states
    elif systems[key]['kind'] == 'Gaussian':
        print('\n\n\nusing the gaussian_density_of_states\n\n\n')
        sigma = systems[key]['sigma']
        exact_density_of_states = gaussian_density_of_states
    elif systems[key]['kind'] == 'Erfinv':
        print('\n\n\nusing the erfinv\n\n\n')
        mean_erfinv_energy = systems[key]['mean_energy']
        N = systems[key]['N']
        exact_density_of_states = erfinv_density_of_states
    elif systems[key]['kind'] == 'Pieces':
        print('\n\n\nusing the pieces\n\n\n')
        a = systems[key]['a']
        b = systems[key]['b']
        e1 = systems[key]['e1']
        e2 = systems[key]['e2']
        exact_density_of_states = piecewise_density_of_states
    else:
        print('\n\n\nusing the most bogus density of states\n\n\n', systems[key]['kind'])
        exact_density_of_states = other_density_of_states
