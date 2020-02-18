import argparse, sys, yaml
import matplotlib.pyplot as plt
import numpy as np


parser = argparse.ArgumentParser()
parser.add_argument('file')
args = parser.parse_args()
file = args.file


#yaml.warnings({'YAMLLoadWarning': False})
with open(file,'rb') as stream:
    try:
        if 'cbor' in file:
            import cbor
            data_loaded = cbor.load(stream)
        else:
            import yaml
            data_loaded = yaml.full_load(stream)
    except IOError:
        print('An error occurred trying to read the file.')
        
"""
with open('two.yaml','r') as stream:
    try:
        data_loaded = yaml.full_load(stream)
    except yaml.yamlerror as exc:
        print(exc)
"""
time_frame = data_loaded['movies']['time']
entropy_data = data_loaded['movies']['entropy']
hist_data = data_loaded['movies']['histogram']
hist = data_loaded['bins']['histogram']
moves = data_loaded['movies']['time']
N_sites = len(data_loaded['system']['N'])
print('N_sites', N_sites)

number_data = np.array(data_loaded['movies']['number'])
energy_data = np.array(data_loaded['movies']['energy'])

print(len(energy_data))

energy_row = 0
for i in range(len(energy_data)-1):
    if energy_data[i+1] != energy_data[i]:
        energy_row = int(i + 1)
        break

energy_col = int(len(energy_data)/energy_row)

nlist = len(energy_data)

k = 1

energy_data.resize(energy_col, energy_row)
number_data.resize(energy_col, energy_row)

E = np.zeros((energy_col + 1, energy_row + 1))
print(E)
E[:-1,:-1] = energy_data
print(E)
N = np.zeros((energy_col + 1, energy_row + 1))
print(N)
N[:-1,:-1] = number_data
print(N)

dE = abs(E[0,0] - E[1,0])
E -= dE/2
print('dE', dE)

E[-1,:] = E[-2,:] + dE
N[-1,:] = N[-2,:]

E[:,-1] = E[:,-2]
N[:,-1] = N[:,-2] + 1

N -= 0.5

#T_inv = (S(E + dE)-S(E))/dE

S_ideal = number_data*k*(1 + np.log(N_sites/number_data))
S_ideal[number_data==0] = 0

for t in range(len(entropy_data)):
    #print('time', moves[t])
    S_excess = np.array(entropy_data[t])
    S_excess.resize(energy_col, energy_row)
    S0_excess = S_excess[-1,0]
    S_excess = S_excess - S0_excess
    S = S_excess + S_ideal
    hist = np.array(hist_data[t])*1.0 # the multiplication converts it to floating point values
    if hist.max() == 0:
        continue
    hist.resize(energy_col, energy_row)
    S[hist==0] = np.nan
    S_excess[hist==0] = np.nan
    hist[hist==0] = np.nan
    
    plt.figure('entropy')
    plt.clf()
    plt.title(f'{moves[t]} moves')
    plt.pcolor(N,E,S) 
    plt.xlabel('$N$')
    plt.ylabel('$E$')
    plt.colorbar() 
    
    plt.figure('excess entropy')
    plt.clf()
    plt.title(f'{moves[t]} moves')
    plt.pcolor(N,E,S_excess) 
    plt.xlabel('$N$')
    plt.ylabel('$E$')
    plt.colorbar() 
    
    plt.figure('histogram')
    plt.clf()
    plt.title(f'{moves[t]} moves')
    plt.pcolor(N,E,hist)
    plt.colorbar()

    plt.pause(1)

plt.show()
