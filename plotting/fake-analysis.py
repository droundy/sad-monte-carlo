import numpy as np
import yaml, argparse, sys
import scipy.constants as scipy
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="fake energies analysis")
parser.add_argument('yaml', help = 'the yaml file')
args = parser.parse_args()

def linear_density_of_states(E):
    return 1
def other_density_of_states(E):
    return 2

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
elif 'other' in args.yaml:
    exact_density_of_states = other_density_of_states

energy_boundaries = [max_energy]
energy_per_rel_bin = 1/bin_norm*(max_energy-min_energy)
for b in rel_bins:
    energy_boundaries.append( energy_boundaries[-1] - b*energy_per_rel_bin)
energy_boundaries = np.array(energy_boundaries)
print('energy boundarie3 are ', energy_boundaries)

mean_energy = total_energy/histogram #includes unbounded extremes
print('mean energies are', mean_energy)

middle_mean_energy = mean_energy[1:-1] #excludes unbounded extremes
print('middle mean energies are', middle_mean_energy)

energy_width = -np.diff(energy_boundaries)
print('energy widths are', energy_width)

middle_entropies_A = lnw[1:-1] - np.log(energy_width)
print('middle entropies are', middle_entropies_A)



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
#plt.plot(middle_mean_energy, np.log(exact_density_of_states), label='Exact')
plt.tight_layout()
plt.legend(loc='upper right')
    
plt.show()