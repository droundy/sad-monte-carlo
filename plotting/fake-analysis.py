import numpy as np
import yaml, cbor, argparse, sys
import scipy.constants as const
import scipy.optimize as optimize
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="fake energies analysis")
parser.add_argument('fname', help = 'the yaml or cbor file')

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
def bisect_bin_entropy(i):
    #i is the bin whose entropy we are calculating
    return

#Read Data
data_loaded = {}
#each file has different path (including extension) so concatenating is easy
fname = args.fname
print(fname)
with open(fname,'rb') as stream:
    if 'yaml' in fname:
        try:
            data_loaded[fname] = yaml.full_load(stream)
        except IOError:
            print('An error occurred trying to read the file.')
    elif 'cbor' in fname:
        try:
            data_loaded[fname] = cbor.load(stream)
        except IOError:
            print('An error occurred trying to read the file.')
    else:
        print('What kind of file is this?!')
        exit(1)

total_energy={}
histogram={}
rel_bins={}
bin_norm={}
max_energy={}
min_energy={}
lnw={}
for key in data_loaded:
    total_energy[key] = np.array(data_loaded[key]['total_energy'])
    histogram[key] = np.array(data_loaded[key]['histogram'])
    rel_bins[key] = np.array(data_loaded[key]['rel_bins'])
    bin_norm[key] = data_loaded[key]['bin_norm']
    max_energy[key] = data_loaded[key]['max_energy']
    min_energy[key] = data_loaded[key]['min_energy']
    lnw[key] = np.array(data_loaded[key]['lnw'])
    lnw[key] -= lnw[key].max()
    lnw[key] -= np.log(np.sum(np.exp(lnw[key]))) # w = w / sum(w)

#Analysis
exact_density_of_states = linear_density_of_states
if 'linear' in fname:
    print('\n\n\nusing the linear_density_of_states\n\n\n')
    exact_density_of_states = linear_density_of_states
elif 'quadratic' in fname:
    print('\n\n\nusing the quadratic_density_of_states\n\n\n')
    exact_density_of_states = quadratic_density_of_states
else:
    print('\n\n\nusing the most bogus density of states\n\n\n')
    exact_density_of_states = other_density_of_states

energy_boundaries={}
energy_per_rel_bin={}
mean_energy={}
middle_mean_energy={}
energy_width={}
middle_entropies_A={}

for key in data_loaded:

    energy_boundaries[key] = [max_energy[key]]
    energy_per_rel_bin[key] = 1/bin_norm[key]*(max_energy[key]-min_energy[key])
    for b in rel_bins[key]:
        energy_boundaries[key].append( energy_boundaries[key][-1] - b*energy_per_rel_bin[key])
    energy_boundaries[key] = np.array(energy_boundaries[key])
    #print('energy boundarie3 are ', energy_boundaries[key])

    mean_energy[key] = total_energy[key]/histogram[key] #includes unbounded extremes
    #print('mean energies are', mean_energy[key])

    middle_mean_energy[key] = mean_energy[key][1:-1] #excludes unbounded extremes
    #print('middle mean energies are', middle_mean_energy[key])
    
    energy_width[key] = -np.diff(energy_boundaries[key])
    #print('energy widths are', energy_width[key])
    
    middle_entropies_A[key] = lnw[key][1:-1] - np.log(energy_width[key])
    #print('middle entropies are', middle_entropies_A)

entropy_boundaries={}
# #Plotting
# plt.figure('average_energy')
# plt.xlabel('Bin')
# plt.ylabel('Average Energy')
# bins = np.arange(len(middle_mean_energy)) #references bin no. or id eg bin 1, 2...
# plt.plot(bins, middle_mean_energy, label=args.yaml)
# plt.tight_layout()
# plt.legend(loc='upper right')

# plt.figure('bin_sizes')
# plt.xlabel('Bin')
# plt.ylabel('Bin Size')
# bins = np.arange(len(energy_width)) #references bin no. or id eg bin 1, 2...
# plt.plot(bins, energy_width, label=args.yaml)
# plt.tight_layout()
# plt.legend(loc='upper right')

# plt.figure('bin_size v avg_temp')
# plt.xlabel('Bin Size')
# plt.ylabel('Average Energy')
# bins = np.arange(len(middle_mean_energy)) #references bin no. or id eg bin 1, 2
# plt.bar(bins, middle_mean_energy, width=energy_width, label=args.yaml)
# plt.tight_layout()
# plt.legend(loc='upper right')

# plt.figure('entropy')
# plt.xlabel('Average Energy in bin')
# plt.ylabel('Entropy')
# plt.plot(middle_mean_energy, middle_entropies_A, label=args.yaml)

print('energy_boundaries', energy_boundaries)

for key in data_loaded:
    entropy_boundaries[key] = np.zeros_like(energy_boundaries[key])
    E_lo = energy_boundaries[key][-1]
    E = np.linspace(4*mean_energy[key][-1] - 3*E_lo, middle_mean_energy[key].max(), 10000)
    S = np.zeros_like(E)
    S_lo = lnw[key][-1] - np.log(E_lo - mean_energy[key][-1])
    S[E < E_lo] = (S_lo - (E_lo - E) / (E_lo - mean_energy[key][-1]))[E < E_lo]
    entropy_boundaries[key][-1] = S_lo
    
    for i in range(len(energy_boundaries[key])-1, 1, -1):
        print('solving for i =', i)
        entropy_boundaries[key][i-1] = optimize_bin_entropy(i, energy_boundaries[key], lnw[key], entropy_boundaries[key][i])
    
    print(entropy_boundaries)
    
    plt.figure('entropy')
    plt.plot(E, S, label='new approximation')
    plt.plot(energy_boundaries[key], entropy_boundaries[key], '.-',
             label=fname+' optimize_bin_entropy approx.')
    plt.plot(E, np.log(exact_density_of_states(E)), label=fname+' exact')
    plt.xlabel('$E$')
    plt.ylabel('$S$')

    plt.figure('density of states')
    plt.plot(E, np.exp(S), label='new approximation')
    plt.plot(energy_boundaries[key], np.exp(entropy_boundaries[key]), '.-',
             label=fname+' optimize_bin_entropy approx.')
    plt.plot(E, exact_density_of_states(E), label=fname+' exact')
    plt.xlabel('$E$')
    plt.ylabel('$D(E)$')

    exact = exact_density_of_states(E)
    plt.ylim(0, 1.1*exact[exact == exact].max())


plt.figure('entropy')
plt.tight_layout()
plt.legend(loc='upper right')
    
plt.show()