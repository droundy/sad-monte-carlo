#!/usr/bin/python3

import numpy as np
import yaml, cbor, argparse, sys
import scipy.constants as const
import scipy.optimize as optimize
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="fake energies analysis")
parser.add_argument('base', nargs='*', help = 'the yaml or cbor file')

args = parser.parse_args()

def linear_density_of_states(E):
    return np.heaviside(E, 0.5)*np.heaviside(1-E, 0.5)
def quadratic_density_of_states(E):
    return 1.5*np.sqrt(E)*np.heaviside(E, 0)*np.heaviside(1-E, 0)
def gaussian_density_of_states(E):
    sigma = 1
    return np.pi()*(sigma**3)*np.sqrt(32*np.log(E)) / E
def other_density_of_states(E):
    return np.zeros_like(E)

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
    sol = optimize.root_scalar(fn_for_beta, args=(meanE_over_deltaE), x0 = 2*meanE_over_deltaE, x1 = meanE_over_deltaE)
    # print(sol)
    return sol.root

def find_entropy_from_beta_and_lnw(beta, lnw, deltaE):
    if abs(beta*deltaE) < 1e-14:
        return lnw - np.log(deltaE)
    return lnw - np.log(deltaE) - np.log((np.exp(beta*deltaE)-1)/(beta*deltaE))


#Read Data
energy_boundaries = {}
mean_energy = {}
lnw = {}
system = {}
#each file has different path (including extension) so concatenating is easy
for base in args.base:
    print(base)
    energy_boundaries[base] = np.loadtxt(base+'-energy-boundaries.dat')
    mean_energy[base] = np.loadtxt(base+'-mean-energy.dat')
    lnw[base] = np.loadtxt(base+'-lnw.dat')

    with open(base+'.yaml','rb') as stream:
        try:
            system[base] = yaml.full_load(stream)['system']
        except IOError:
            print('An error occurred trying to read the file.')

    if energy_boundaries[base][0] < energy_boundaries[base][-1]:
        energy_boundaries[base] = np.flip(energy_boundaries[base])
        mean_energy[base] = np.flip(mean_energy[base])
        lnw[base] = np.flip(lnw[base])

entropy_boundaries={}

print('energy_boundaries', energy_boundaries)

for key in lnw:
    
    #Analysis
    function_type = list(system[key]['function'].keys())[0]
    exact_density_of_states = linear_density_of_states
    if function_type == 'Linear':
        print('\n\n\nusing the linear_density_of_states\n\n\n')
        exact_density_of_states = linear_density_of_states
    elif function_type == 'Quadratic':
        print('\n\n\nusing the quadratic_density_of_states\n\n\n')
        exact_density_of_states = quadratic_density_of_states
    elif function_type == 'Gaussian':
        print('\n\n\nusing the quadratic_density_of_states\n\n\n')
        exact_density_of_states = gaussian_density_of_states
    else:
        print('\n\n\nusing the most bogus density of states\n\n\n')
        exact_density_of_states = other_density_of_states

    E_lo = energy_boundaries[key][-1]
    dE_lo = energy_boundaries[key][-2] - E_lo
    if np.isnan(mean_energy[key][-1]):
        E = np.linspace(E_lo - 3*dE_lo, energy_boundaries[key].max(), 100000)
    else:
        E = np.linspace(4*mean_energy[key][-1] - 3*E_lo, energy_boundaries[key].max(), 100000)
    S_lo = lnw[key][-1] - np.log(E_lo - mean_energy[key][-1])
    min_T = E_lo - mean_energy[key][-1]
    
    # entropy_boundaries[key] = compute_entropy_given_Smin(S_lo, energy_boundaries[key], lnw[key])
    
    # S = np.interp(E, energy_boundaries[key][::-1], entropy_boundaries[key][::-1])
    # S[E < E_lo] = (S_lo - (E_lo - E)/min_T)[E < E_lo]

    # entropy_boundaries_best = optimize_entropy(energy_boundaries[key], lnw[key])
    # min_T = np.exp(lnw[key][-1] - entropy_boundaries_best[-1])
    # Sbest = np.interp(E, energy_boundaries[key][::-1], entropy_boundaries_best[::-1])
    # Sbest[E < E_lo] = (entropy_boundaries_best[-1] - (E_lo - E) /min_T)[E < E_lo]

    print(len(lnw[key]), len(energy_boundaries[key]))
    middle_entropy = lnw[key][1:-1] - np.log(-np.diff(energy_boundaries[key]))
    Smiddle = np.zeros_like(E)
    Ssloped = np.zeros_like(E)
    for i in range(len(middle_entropy)):
        print(key, 'working on', i, '/', len(middle_entropy))
        here = E<=energy_boundaries[key][i]
        Smiddle[here] = middle_entropy[i]

        deltaE = energy_boundaries[key][i] - energy_boundaries[key][i+1]
        meanE = mean_energy[key][i+1]
        beta = find_beta_deltaE((meanE - energy_boundaries[key][i+1])/deltaE)/deltaE
        S0 = find_entropy_from_beta_and_lnw(beta, lnw[key][i+1], deltaE)
        # print('deltaE', deltaE, 'beta', beta, 'S0', S0, 'dimensionless mean', (meanE - energy_boundaries[key][i+1])/deltaE)
        Ssloped[here] = S0 + beta*(E[here] - energy_boundaries[key][i+1])

    plt.figure('entropy')
    plt.plot(E, Smiddle, '-', label=key + 'simplest')
    plt.plot(E, Ssloped, '--', label=key + 'sloped')
    # plt.plot(E, S, '-', label=key+' optimize_bin_entropy approx.')
    # plt.plot(E, Sbest, '-', label=key + 'smooth')
    # plt.plot(E, np.log(exact_density_of_states(E)), label=key+' exact')
    plt.xlabel('$E$')
    plt.ylabel('$S$')

    plot_dos = False
    if plot_dos:
        plt.figure('density of states')
        # plt.plot(E, np.exp(S), '-', label=key+' optimize_bin_entropy approx.')
        # plt.plot(E, np.exp(Sbest), '-', label=key+' smooth')
        plt.plot(E, exact_density_of_states(E), label=key+' exact')
        plt.xlabel('$E$')
        plt.ylabel('$D(E)$')
        exact = exact_density_of_states(E)
        plt.ylim(0, 1.1*exact[exact == exact].max())


plt.figure('entropy')
plt.tight_layout()
plt.legend(loc='upper right')
    
plt.show()
