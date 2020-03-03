import argparse, sys, yaml
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

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

# with open('two.yaml','r') as stream:
#     try:
#         data_loaded = yaml.full_load(stream)
#     except yaml.yamlerror as exc:
#         print(exc)

time_frame = data_loaded['movies']['time']
entropy_data = data_loaded['movies']['entropy']
hist_data = data_loaded['movies']['histogram']
hist = data_loaded['bins']['histogram']
moves = data_loaded['movies']['time']
N_sites = len(data_loaded['system']['N'])
#print('N_sites', N_sites)

number_data = np.array(data_loaded['movies']['number'])
energy_data = np.array(data_loaded['movies']['energy'])

#print(len(energy_data))

row = 0
for i in range(len(energy_data)-1):
    if energy_data[i+1] != energy_data[i]:
        row = int(i + 1)
        break

col = int(len(energy_data)/row)

nlist = len(energy_data)

k = 1

energy_data.resize(col, row)
number_data.resize(col, row)

E = np.zeros((col + 1, row + 1))
#print(E)
E[:-1,:-1] = energy_data
#print(E)
N = np.zeros((col + 1, row + 1))
#print(N)
N[:-1,:-1] = number_data
#print(N)

dE = abs(E[0,0] - E[1,0])
E -= dE/2
#print('dE', dE)

E[-1,:] = E[-2,:] + dE
N[-1,:] = N[-2,:]

E[:,-1] = E[:,-2]
N[:,-1] = N[:,-2] + 1

N -= 0.5

#chemical potential

S_ideal = number_data*k*(1 + np.log(N_sites/number_data))
#is s_ideal correct?
S_ideal[number_data==0] = 0

T = np.zeros((col, row))

for t in range(len(entropy_data)):
    #print('time', moves[t])
    S_excess = np.array(entropy_data[t])
    S_excess.resize(col, row)
    S0_excess = S_excess[-1,0]
    S_excess = S_excess - S0_excess
    S = S_excess + S_ideal
    hist = np.array(hist_data[t])*1.0 # the multiplication converts it to floating point values
    if hist.max() == 0:
        continue
    hist.resize(col, row)
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
    
    plt.figure('temperature')
    plt.clf()
    plt.title(f'{moves[t]} moves')
    plt.pcolor(N,E,S) 
    plt.xlabel('$N$')
    plt.ylabel('$E$')
    plt.colorbar() 
    
    plt.figure('histogram')
    plt.clf()
    plt.title(f'{moves[t]} moves')
    plt.pcolor(N,E,hist)
    plt.colorbar()


    #how to integrate temperature and chemical potential into main loop?
    #temperature
    for i in np.arange(0,col-1,1):
        for j in np.arange(0,row-1,1):
            T[i][j] = dE / (S[i+1][j+1] - S[i][j])

    #T = np.delete(T, col-1, 0)
    #T = np.delete(T, 0, row-1)


    plt.figure('temperature')
    plt.clf()
    plt.pcolor(N,E,T, norm=LogNorm(vmin=0.01, vmax=10)) #, cmap='PuBu_r'vmax=0.5, vmin=0)
    plt.xlabel('$N$')
    plt.ylabel('$E$')
    plt.colorbar() 

    #chemical potential
    chem_pot = np.zeros((col, row))
    for i in np.arange(0,col-2,1):
        for j in np.arange(0,row-2,1):
            chem_pot[i][j] = -T[i][j] * ((S[i+1][j+1] - S[i][j]) / (N[i+1][j+1] - N[i][j]))

    plt.figure('chemical potential')
    plt.clf()
    plt.pcolor(N,E,chem_pot)
    plt.xlabel('$N$')
    plt.ylabel('$E$')
    plt.colorbar()

    plt.pause(1)

plt.show()
