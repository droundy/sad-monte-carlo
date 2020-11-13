#!/usr/bin/python3

import numpy as np
import yaml, cbor, argparse, sys, os, glob
import scipy.constants as const
import scipy.optimize as optimize
import matplotlib.pyplot as plt
import colorcet as cc

parser = argparse.ArgumentParser(description="fake energies analysis")
parser.add_argument('base', nargs='*', help = 'the yaml or cbor file')
parser.add_argument('--intensive', action='store_true')

args = parser.parse_args()


prop_cycle = plt.rcParams['axes.prop_cycle']

colors = cc.glasbey_dark

# plt.style.use('dark_background')
# colors = cc.glasbey_light

# prop_cycle = plt.rcParams['axes.prop_cycle']

def read_file(base):
    energy_boundaries = np.loadtxt(base+'-energy-boundaries.dat')
    mean_energy = np.loadtxt(base+'-mean-energy.dat')
    lnw = np.loadtxt(base+'-lnw.dat')
    with open(base+'-system.dat') as f:
        system = yaml.safe_load(f)

    if energy_boundaries[0] < energy_boundaries[-1]:
        energy_boundaries = np.flip(energy_boundaries)
        mean_energy = np.flip(mean_energy)
        lnw = np.flip(lnw)
    lnw -= lnw.max()
    lnw -= np.log(np.sum(np.exp(lnw)))
    return energy_boundaries, mean_energy, lnw, system

def step_entropy(energy_boundaries, mean_energy, lnw):
    step_entropy = []
    step_energy = []
    for i in range(len(energy_boundaries)-1):
        step_energy.append(energy_boundaries[i])
        step_energy.append(energy_boundaries[i+1])
        Shere = lnw[i+1] - np.log(energy_boundaries[i] - energy_boundaries[i+1])
        step_entropy.append(Shere)
        step_entropy.append(Shere)
    step_energy = np.flip(np.array(step_energy))
    step_entropy = np.flip(np.array(step_entropy))
    def entropy(E):
        return np.interp(E, step_energy, step_entropy, left=step_entropy[0], right=step_entropy[-1])
    return entropy, 1*step_energy, 1*step_entropy

def linear_entropy(energy_boundaries, mean_energy, lnw):
    step_entropy = []
    step_energy = []
    for i in range(len(energy_boundaries)-1):
        step_energy.append(energy_boundaries[i])
        step_energy.append(energy_boundaries[i+1])

        deltaE = energy_boundaries[i] - energy_boundaries[i+1]
        meanE = mean_energy[i+1]
        beta = find_beta_deltaE((meanE - energy_boundaries[i+1])/deltaE)/deltaE
        S0 = find_entropy_from_beta_and_lnw(beta, lnw[i+1], deltaE)
        if np.isnan(beta+S0):
            step_entropy.append(lnw[i+1]-np.log(deltaE))
            step_entropy.append(lnw[i+1]-np.log(deltaE))
        else:
            step_entropy.append(S0+beta*deltaE)
            step_entropy.append(S0)
    
    # this is the unbounded low-energy bin, assume exponential DOS
    Tlow = energy_boundaries[-1] - mean_energy[-1]
    Slo = lnw[-1] - np.log(Tlow)
    step_energy.append(energy_boundaries[-1])
    step_energy.append(energy_boundaries[-1] - 20*Tlow)
    step_entropy.append(Slo)
    step_entropy.append(Slo - 10)

    step_energy = np.flip(step_energy)
    step_entropy = np.flip(step_entropy)
    def entropy(E):
        return np.interp(E, step_energy, step_entropy, left=step_entropy[0], right=step_entropy[-1])
    return entropy, 1*step_energy, 1*step_entropy

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
    if E < -e2:
        return a**3 * np.sqrt(E+e1) / 2
    elif E > 0:
        return (b + (b - a)*np.sqrt(E/e2 + 1))**2 / (2*e2*np.sqrt(E/e2+1))
    return a**3 * np.sqrt(E+e1) / 2 + (b - (b - a)*np.sqrt(E/e2 + 1))**2 / (2*e2*np.sqrt(E/e2+1)) + (b + (b - a)*np.sqrt(E/e2 + 1))**2 / (2*e2*np.sqrt(E/e2+1))
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
    if abs(x) < 1e-14:
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

frames = set()
bases = []
print('base', args.base)
#each file has different path (including extension) so concatenating is easy
for base in args.base:
    if '.cbor' in base or '.yaml' in base:
        base = base[:-5]
    print(base)
    bases.append(base)
    for f in glob.glob(base+'/*.cbor'):
        f = os.path.splitext(os.path.basename(f))[0]
        if float(f) > 1e7:
            frames.add(f)

best_energy_boundaries, best_mean_energy, best_lnw, best_system = read_file(args.base[0])
unscaled_best_function, best_energy, best_entropy = linear_entropy(best_energy_boundaries, best_mean_energy, best_lnw)
if args.intensive:
    scale = 1/best_system['N']
    best_energy = scale*best_energy
    best_entropy = scale*best_entropy
    best_function = lambda e: scale*unscaled_best_function(e/scale)
else:
    best_function = unscaled_best_function
energies_to_compare = np.array(best_mean_energy)
energies_reference = best_function(energies_to_compare)
# Only compare energies below the energy with max entropy
# nhi = energies_reference.argmax()
# energies_to_compare = energies_to_compare[nhi:-1]
# energies_reference = energies_reference[nhi:-1]

sigma = 1

frames = sorted(list(frames))

comparison_time = {}
comparison_maxerror = {}
for k in bases:
    comparison_time[k] = []
    comparison_maxerror[k] = []

plt.ion()
for f in frames:
    plt.figure('comparison')
    plt.clf()
    plt.figure('entropy')
    plt.clf()
    main_ax = plt.gca()

    which_color = 0
    for key in bases:
        color = colors[which_color]
        which_color += 1

        final_energy_boundaries, final_mean_energy, final_lnw, final_system = read_file(key)

        frame_filename = f'{key}/{f}' # same as '{}/{}'.format(key, f) or same as '%s/%s' % (key, f)
        print(frame_filename)
        try:
            energy_boundaries, mean_energy, lnw, system = read_file(frame_filename)
        except:
            print('There were no boundaries')
            continue
        
        #Analysis
        exact_density_of_states = linear_density_of_states
        if 'linear' in key:
            print('\n\n\nusing the linear_density_of_states\n\n\n')
            exact_density_of_states = linear_density_of_states
        elif 'quadratic' in key:
            print('\n\n\nusing the quadratic_density_of_states\n\n\n')
            exact_density_of_states = quadratic_density_of_states
        elif system['kind'] == 'Gaussian':
            print('\n\n\nusing the gaussian_density_of_states\n\n\n')
            sigma = system['sigma']
            exact_density_of_states = gaussian_density_of_states
        elif system['kind'] == 'Erfinv':
            print('\n\n\nusing the erfinv\n\n\n')
            mean_erfinv_energy = system['mean_energy']
            N = system['N']
            exact_density_of_states = erfinv_density_of_states
        elif system['kind'] == 'Pieces':
            print('\n\n\nusing the pieces\n\n\n')
            a = system['a']
            b = system['b']
            e1 = system['e1']
            e2 = system['e2']
            exact_density_of_states = piecewise_density_of_states
        else:
            print('\n\n\nusing the most bogus density of states\n\n\n', system['kind'])
            exact_density_of_states = other_density_of_states
        
        s_function, s_energy, s_entropy = step_entropy(energy_boundaries, mean_energy, lnw)
        l_function, l_energy, l_entropy = linear_entropy(energy_boundaries, mean_energy, lnw)

        f_function, f_energy, f_entropy = linear_entropy(final_energy_boundaries, final_mean_energy, final_lnw)

        plt.figure('entropy')
        scale = 1
        if args.intensive:
            scale = 1/system['N']
        print('scale is', scale)
        main_ax.plot(scale*s_energy, scale*s_entropy, '-', label=key + ' step', color=color)
        main_ax.plot(scale*l_energy, scale*l_entropy, '--', label=key + ' linear', color=color)
        main_ax.plot(energies_to_compare, scale*np.log(exact_density_of_states(energies_to_compare/scale)), color='#aaaaaa')
        main_ax.plot(scale*f_energy, scale*f_entropy, '-', linewidth=4, alpha=0.2, color=color, label='final')

        my_s = scale*l_function(energies_to_compare/scale)
        my_s[np.isnan(my_s)] = 1
        #main_ax.plot(energies_to_compare, my_s, 'x', color=color)
        #main_ax.plot(energies_to_compare, energies_reference, '+', color='#999999')

        main_ax.legend(loc='best')
        plt.xlabel('$E$')
        plt.ylabel('$S$')
        plt.title(f.lstrip('0'))
        if args.intensive:
            plt.xlabel('$E/N$')
            plt.ylabel('$S/N$')
        # this is an inset axes over the main axes
        a = plt.axes([.3, .2, .4, .5])
        # a.plot(scale*f_energy, scale*f_entropy, '-', linewidth=4, alpha=0.2, color=color)
        # a.plot(scale*l_energy, scale*l_entropy, '--', color=color)
        # plt.xlim(0, 20)

        a.plot(energies_to_compare[1:-1], (my_s-energies_reference)[1:-1], '-', color=color)
        #a.plot(energies_to_compare, energies_reference, '+', color='#999999')
        a.set_ylabel('error in entropy')

        comparison_time[key].append(float(f))
        me = abs(my_s - energies_reference).max()
        for i in range(len(energies_to_compare)):
            print(energies_to_compare[i], my_s[i])

        comparison_maxerror[key].append(me)

    c = 0
    for key in bases:
        color = colors[c]
        c += 1
        plt.figure('comparison')
        print('max entropy errors', key, comparison_maxerror[key])
        plt.loglog(comparison_time[key], comparison_maxerror[key], '.-', color=color, label=key)
        plt.xlabel('$t$')
        plt.ylabel('max error in $S$')
        plt.legend(loc='best')

    if which_color > 0:
        plt.pause(0.0001)
        plt.draw()

plt.ioff()
plt.show()
