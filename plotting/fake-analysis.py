import numpy as np
import yaml, argparse, sys
import scipy.constants as scipy
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="fake energies analysis")

parser.add_argument('yaml', help = 'the yaml file')
args = parser.parse_args()

with open(args.yaml,'rb') as stream:
    try:
        data_loaded = yaml.full_load(stream)
    except IOError:
        print('An error occurred trying to read the file.')

total_energy = np.array(data_loaded['total_energy'])
histogram = np.array(data_loaded['histogram'])
rel_bins = np.array(data_loaded['rel_bins'])
lnw = np.array(data_loaded['lnw'])
bin_norm = data_loaded['bin_norm']
max_energy = data_loaded['max_energy']
min_energy = data_loaded['min_energy']

energy_boundaries = [max_energy]
energy_per_rel_bin = 1/bin_norm*(max_energy-min_energy)
for b in rel_bins:
    energy_boundaries.append( energy_boundaries[-1] - b*energy_per_rel_bin)
energy_boundaries = np.array(energy_boundaries)
print('energy boundarie3 are ', energy_boundaries)

mean_energy = total_energy/histogram
print('mean energies are', mean_energy)

middle_mean_energy = mean_energy[1:-1]
print('middle mean energies are', middle_mean_energy)

energy_width = -np.diff(energy_boundaries)
print('energy widths are', energy_width)

middle_entropies_A = lnw[1:-1] - np.log(energy_width)
print('middle entropies are', middle_entropies_A)
exit(1)

plt.figure('average_energy')
plt.xlabel('Bin')
plt.ylabel('Average Energy')
bins = np.arange(len(total_energy)) #references bin no. or id eg bin 1, 2...
plt.plot(bins, total_energy/histogram, label=args.yaml)
plt.tight_layout()
plt.legend(loc='upper right')

plt.figure('distribution_of_states')
plt.xlabel('Bin')
plt.ylabel('Proportion in Bin')
bins = np.arange(len(rel_bins)) #references bin no. or id eg bin 1, 2...
plt.plot(bins, rel_bins/bin_norm, label=args.yaml)
plt.tight_layout()
plt.legend(loc='upper right')
    
plt.figure('bin_size v avg_temp')
plt.xlabel('Bin Size')
plt.ylabel('Average Energy')
bins = np.arange(len(total_energy)) #references bin no. or id eg bin 1, 2...
rel_bins = np.append(rel_bins, [0, 0])
plt.bar(bins, total_energy/histogram, width=rel_bins/bin_norm, label=args.yaml)
plt.tight_layout()
plt.legend(loc='upper right')

plt.figure('entropy')
plt.xlabel('Average Energy in bin')
plt.ylabel('Entropy')

dE = rel_bins/bin_norm
E = total_energy/histogram
W = np.array(1 / 2**np.arange(len(dE))[::-1], dtype=np.double)
plt.plot(E, lnw - np.log(dE), label=args.yaml)
plt.plot(E, np.log(W), label='Exact')

print(W)
plt.tight_layout()
plt.legend(loc='upper right')
    
plt.show()