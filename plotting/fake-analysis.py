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
bin_norm = np.array(data_loaded['bin_norm'])

plt.figure('average_energy')
plt.xlabel('Bin')
plt.ylabel('Average Energy')
bins = np.arange(len(total_energy)) #references bin no. or id eg bin 1, 2...
try:
    plt.plot(bins, total_energy/histogram, label=args.yaml)
    plt.tight_layout()
    plt.legend(loc='upper right')
except IndexError:
    plt.plot(bins, total_energy/histogram, label=args.yaml)
    plt.tight_layout()
    plt.legend(loc='upper right')

plt.figure('distribution_of_states')
plt.xlabel('Bin')
plt.ylabel('Proportion in Bin')
bins = np.arange(len(rel_bins)) #references bin no. or id eg bin 1, 2...
try:
    plt.plot(bins, rel_bins/bin_norm, label=args.yaml)
    plt.tight_layout()
    plt.legend(loc='upper right')
except IndexError:
    plt.plot(bins, rel_bins/bin_norm, label=args.yaml)
    plt.tight_layout()
    plt.legend(loc='upper right')
    
plt.figure('bin_size v avg_temp')
plt.xlabel('Bin Size')
plt.ylabel('Average Energy')
bins = np.arange(len(total_energy)) #references bin no. or id eg bin 1, 2...
rel_bins = np.append(rel_bins, [0, 0])
try:
    plt.bar(bins, total_energy/histogram, width=rel_bins/bin_norm, label=args.yaml)
    plt.tight_layout()
    plt.legend(loc='upper right')
except IndexError:
    plt.bar(bins, total_energy/histogram, width=rel_bins/bin_norm, label=args.yaml)
    plt.tight_layout()
    plt.legend(loc='upper right')

plt.figure('entropy')
plt.xlabel('Energy')
plt.ylabel('Entropy')
try:
    dE = rel_bins/bin_norm
    E = total_energy/histogram
    W = np.array(2**np.arange(len(dE))[::-1], dtype = np.double)
    W = np.reciprocal(W)
    plt.plot(E, np.log(W/dE), label=args.yaml)
    plt.plot(E, np.log(W), label='Exact')
    
    print(W)
    plt.tight_layout()
    plt.legend(loc='upper right')
except IndexError:
    plt.plot(dE, rel_bins/bin_norm, label=args.yaml)
    plt.tight_layout()
    plt.legend(loc='upper right')
    
plt.pause(0.6)