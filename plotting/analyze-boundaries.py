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
def other_density_of_states(E):
    return 2

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

#Read Data
energy_boundaries = {}
mean_energy = {}
lnw = {}
#each file has different path (including extension) so concatenating is easy
for base in args.base:
    print(base)
    energy_boundaries[base] = np.loadtxt(base+'-energy-boundaries.dat')
    mean_energy[base] = np.loadtxt(base+'-mean-energy.dat')
    lnw[base] = np.loadtxt(base+'-lnw.dat')

entropy_boundaries={}

print('energy_boundaries', energy_boundaries)

for key in lnw:
    
    #Analysis
    exact_density_of_states = linear_density_of_states
    if 'linear' in key:
        print('\n\n\nusing the linear_density_of_states\n\n\n')
        exact_density_of_states = linear_density_of_states
    elif 'quadratic' in key:
        print('\n\n\nusing the quadratic_density_of_states\n\n\n')
        exact_density_of_states = quadratic_density_of_states
    else:
        print('\n\n\nusing the most bogus density of states\n\n\n')
        exact_density_of_states = other_density_of_states

    E_lo = energy_boundaries[key][-1]
    E = np.linspace(4*mean_energy[key][-1] - 3*E_lo, energy_boundaries[key].max(), 10000)
    S_lo = lnw[key][-1] - np.log(E_lo - mean_energy[key][-1])
    min_T = E_lo - mean_energy[key][-1]
    
    entropy_boundaries[key] = compute_entropy_given_Smin(S_lo, energy_boundaries[key], lnw[key])
    
    S = np.interp(E, energy_boundaries[key][::-1], entropy_boundaries[key][::-1])
    S[E < E_lo] = (S_lo - (E_lo - E)/min_T)[E < E_lo]

    entropy_boundaries_best = optimize_entropy(energy_boundaries[key], lnw[key])
    min_T = np.exp(lnw[key][-1] - entropy_boundaries_best[-1])
    Sbest = np.interp(E, energy_boundaries[key][::-1], entropy_boundaries_best[::-1])
    Sbest[E < E_lo] = (entropy_boundaries_best[-1] - (E_lo - E) /min_T)[E < E_lo]

    plt.figure('entropy')
    plt.plot(E, S, '-', label=key+' optimize_bin_entropy approx.')
    plt.plot(E, Sbest, '-', label=key + 'smooth')
    plt.plot(E, np.log(exact_density_of_states(E)), label=key+' exact')
    plt.xlabel('$E$')
    plt.ylabel('$S$')

    plt.figure('density of states')
    plt.plot(E, np.exp(S), '-', label=key+' optimize_bin_entropy approx.')
    plt.plot(E, np.exp(Sbest), '-', label=key+' smooth')
    plt.plot(E, exact_density_of_states(E), label=key+' exact')
    plt.xlabel('$E$')
    plt.ylabel('$D(E)$')
    exact = exact_density_of_states(E)
    plt.ylim(0, 1.1*exact[exact == exact].max())


plt.figure('entropy')
plt.tight_layout()
plt.legend(loc='upper right')
    
plt.show()