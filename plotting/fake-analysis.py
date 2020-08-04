import numpy as np
import yaml, argparse, sys
import scipy.constants as const
import scipy.optimize as optimize
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="fake energies analysis")
parser.add_argument('yaml', help = 'the yaml file')
args = parser.parse_args()

def linear_density_of_states(E):
    return np.heaviside(E, 0.5)*np.heaviside(1-E, 0.5)
def quadratic_density_of_states(E):
    return 1.5*np.sqrt(E)*np.heaviside(E, 0)*np.heaviside(1-E, 0)
def other_density_of_states(E):
    return 2

#The function needs to be callable

def fn_entropy(S_i_1, E_i_1, E_i, W_i, S_i):
    ''' This is the thing that should be zero '''
    return ((S_i - S_i_1)/(E_i - E_i_1)
            * np.exp( (S_i_1*E_i - S_i*E_i_1) / (E_i - E_i_1) )
            * (
                np.exp( (S_i - S_i_1)/(E_i - E_i_1) * E_i_1 )
                - np.exp( (S_i - S_i_1)/(E_i - E_i_1) * E_i )
            )
            - W_i)
def optimize_bin_entropy(i, E_bounds, W, S_i):
    #i is the bin whose entropy we are calculating
    sol = optimize.root(fn_entropy, [0], args=(E_bounds[i-1], E_bounds[i], W[i-1], S_i))
    return sol.x
def bisect_bin_entropy(i):
    #i is the bin whose entropy we are calculating
    return

#Read From Yaml file
with open(args.yaml,'rb') as stream:
    try:
        data_loaded = yaml.full_load(stream)
    except IOError:
        print('An error occurred trying to read the file.')

total_energy = np.array(data_loaded['total_energy'])
histogram = np.array(data_loaded['histogram'])
rel_bins = np.array(data_loaded['rel_bins'])
bin_norm = data_loaded['bin_norm']
max_energy = data_loaded['max_energy']
min_energy = data_loaded['min_energy']

lnw = np.array(data_loaded['lnw'])
lnw -= lnw.max()
lnw -= np.log(np.sum(np.exp(lnw))) # w = w / sum(w)


#Analysis
exact_density_of_states = linear_density_of_states
if 'linear' in args.yaml:
    exact_density_of_states = linear_density_of_states
elif 'quadratic' in args.yaml:
    exact_density_of_states = quadratic_density_of_states
elif 'other' in args.yaml:
    exact_density_of_states = other_density_of_states
elif 'other' in args.yaml:
    exact_density_of_states = other_density_of_states

energy_boundaries = [max_energy]
energy_per_rel_bin = 1/bin_norm*(max_energy-min_energy)
for b in rel_bins:
    energy_boundaries.append( energy_boundaries[-1] - b*energy_per_rel_bin)
energy_boundaries = np.array(energy_boundaries)
#print('energy boundarie3 are ', energy_boundaries)

mean_energy = total_energy/histogram #includes unbounded extremes
#print('mean energies are', mean_energy)

middle_mean_energy = mean_energy[1:-1] #excludes unbounded extremes
#print('middle mean energies are', middle_mean_energy)

energy_width = -np.diff(energy_boundaries)
#print('energy widths are', energy_width)

middle_entropies_A = lnw[1:-1] - np.log(energy_width)
#print('middle entropies are', middle_entropies_A)

all_entropies = [lnw[0] - (max_energy - energy_boundaries[0])]
all_entropies = np.append(all_entropies, middle_entropies_A)
all_entropies = np.append(all_entropies, lnw[-1] - (energy_boundaries[-1] - min_energy))
#print('all entropies are', all_entropies)


#Plotting
plt.figure('average_energy')
plt.xlabel('Bin')
plt.ylabel('Average Energy')
bins = np.arange(len(middle_mean_energy)) #references bin no. or id eg bin 1, 2...
plt.plot(bins, middle_mean_energy, label=args.yaml)
plt.tight_layout()
plt.legend(loc='upper right')

plt.figure('bin_sizes')
plt.xlabel('Bin')
plt.ylabel('Bin Size')
bins = np.arange(len(energy_width)) #references bin no. or id eg bin 1, 2...
plt.plot(bins, energy_width, label=args.yaml)
plt.tight_layout()
plt.legend(loc='upper right')

plt.figure('bin_size v avg_temp')
plt.xlabel('Bin Size')
plt.ylabel('Average Energy')
bins = np.arange(len(middle_mean_energy)) #references bin no. or id eg bin 1, 2
plt.bar(bins, middle_mean_energy, width=energy_width, label=args.yaml)
plt.tight_layout()
plt.legend(loc='upper right')

plt.figure('entropy')
plt.xlabel('Average Energy in bin')
plt.ylabel('Entropy')
plt.plot(middle_mean_energy, middle_entropies_A, label=args.yaml)

entropy_boundaries = np.zeros_like(energy_boundaries)
E_lo = energy_boundaries[-1]
E = np.linspace(4*mean_energy[-1] - 3*E_lo, middle_mean_energy.max(), 10000)
S = np.zeros_like(E)
S_lo = lnw[-1] - np.log(E_lo - mean_energy[-1])
S[E < E_lo] = (S_lo - (E_lo - E) / (E_lo - mean_energy[-1]))[E < E_lo]
entropy_boundaries[-1] = S_lo


# TODO: Calculate the entropy values in the other bins
# (see top of last page of notes from 7/21)

plt.plot(E, S, label='new approximation')
plt.plot(energy_boundaries, entropy_boundaries, '.-', label='new approximaton in middle')
plt.plot(E, np.log(exact_density_of_states(E)), label='exact')

S_i = 3 #just some random entropy boundary. Its unrealistic
print(optimize_bin_entropy(3, energy_boundaries, np.exp(lnw), S_i))


plt.tight_layout()
plt.legend(loc='upper right')
    
plt.show()